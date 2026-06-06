# =============================================================================
# dyn_abun.R — Dail-Madsen open-population N-mixture family
#
# Latent abundance evolves across primary seasons: N_1 ~ Poisson(lambda); for
# t >= 2, N_t = survivors S_t ~ Binomial(N_{t-1}, omega) plus recruits
# G_t ~ Poisson(gamma); observed via Binomial(N_t, p) over secondary visits. The
# latent abundance sequence is NOT closed-form (unlike the static N-mixture) — it
# is summed out by an exact HMM forward recursion over the abundance states
# (src/dyn_abun_kernel.h), with analytic gradients by forward-mode differentiation
# of the scaled forward algorithm. The fit maximises the exact marginal with that
# analytic gradient (BFGS) and an observed-information covariance; a NUTS path
# samples the same marginal.
#
# Four site-level arms: initial abundance lambda (log; the formula), detection p
# (logit; detection), apparent survival omega (logit; omega_formula), and
# recruitment gamma (log; gamma_formula).
#
#   .tobs_build_dyn_abun()  data binder -> model_type = "dyn_abun"
#   .tobs_fit_dyn_abun()    dispatch to the open N-mixture fit
#   dyn_abun_laplace()      analytic-gradient BFGS over the forward marginal
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a Dail-Madsen open N-mixture. `y` is a 3D array [n_sites x max_visits x
# n_seasons] (or a list of n_sites x max_visits matrices, one per season); a
# missing visit is NA (dropped from its season). All four arms are site-level:
# initial abundance (occ_formula = the tobs formula), detection (det_formula),
# survival (omega_formula), recruitment (gamma_formula).
.tobs_build_dyn_abun <- function(occ_formula, det_formula, data, y,
                                 omega_formula = ~1, gamma_formula = ~1,
                                 mixture = "poisson", K_max = NULL) {
  if (is.list(y) && !is.array(y)) {
    n_seasons <- length(y)
    n_sites <- nrow(y[[1]]); max_visits <- ncol(y[[1]])
    ya <- array(NA_integer_, dim = c(n_sites, max_visits, n_seasons))
    for (t in seq_len(n_seasons)) ya[, , t] <- as.matrix(y[[t]])
    y <- ya
  }
  if (length(dim(y)) != 3) {
    stop("y must be a 3D array [n_sites x max_visits x n_seasons] or a list of ",
         "per-season matrices.", call. = FALSE)
  }
  n_sites <- dim(y)[1]; max_visits <- dim(y)[2]; n_seasons <- dim(y)[3]
  if (n_sites != nrow(data)) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)),
         call. = FALSE)
  }
  if (n_seasons < 2L) {
    stop("dyn_abun() needs >= 2 primary seasons; for a single season use abun().",
         call. = FALSE)
  }
  if (!identical(mixture, "poisson")) {
    stop("dyn_abun() currently supports the Poisson initial abundance only; ",
         "negative-binomial dynamics are not yet wired. (#37)", call. = FALSE)
  }

  bind <- .tobs_bind_formulas(list(lambda = occ_formula, p = det_formula,
                                   omega = omega_formula, gamma = gamma_formula), data)
  X_lambda <- model.matrix(bind$fe$lambda, data)
  X_p      <- model.matrix(bind$fe$p, data)
  X_omega  <- model.matrix(bind$fe$omega, data)
  X_gamma  <- model.matrix(bind$fe$gamma, data)

  # Flat layout site-major then season then visit: j + J*t + J*T*i (0-based),
  # matching compute_dyn_abun_site (dyn_occu's aperm(y, c(2,3,1)) convention).
  y_flat <- as.integer(aperm(y, c(2, 3, 1)))
  y_flat[is.na(y_flat)] <- -1L

  obs <- y[!is.na(y)]
  if (any(obs < 0)) stop("y must contain nonnegative integer counts.", call. = FALSE)
  K <- if (is.null(K_max)) as.integer(max(obs) + 40L) else as.integer(K_max)

  structure(list(
    model_type = "dyn_abun",
    y          = y,
    y_flat     = y_flat,
    X_processes = list(X_lambda, X_p, X_omega, X_gamma),
    formulas   = list(lambda = bind$fe$lambda, p = bind$fe$p,
                      omega = bind$fe$omega, gamma = bind$fe$gamma),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = n_sites, n_seasons = n_seasons, max_visits = max_visits,
    K_max      = K, mixture = mixture,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda), link = "log"),
      list(name = "p",      p = ncol(X_p),      coef_names = colnames(X_p),      link = "logit"),
      list(name = "omega",  p = ncol(X_omega),  coef_names = colnames(X_omega),  link = "logit"),
      list(name = "gamma",  p = ncol(X_gamma),  coef_names = colnames(X_gamma),  link = "log")
    ),
    mean_count = mean(obs), max_count = max(obs)
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter (called from .tobs_fit_model for model_type == "dyn_abun")
# ---------------------------------------------------------------------------

.tobs_fit_dyn_abun <- function(model, max_iter = 300L, tol = 1e-8, verbose = TRUE) {
  raw <- dyn_abun_laplace(
    y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
    J = model$max_visits, K_max = model$K_max,
    X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
    X_omega = model$X_processes[[3]], X_gamma = model$X_processes[[4]],
    max_iter = as.integer(max_iter), tol = as.numeric(tol), verbose = isTRUE(verbose))
  build_dyn_abun_fit(raw, model)
}


# ---------------------------------------------------------------------------
# Laplace fit (analytic-gradient BFGS over the forward marginal)
# ---------------------------------------------------------------------------

#' Maximum-likelihood fit of the Dail-Madsen open N-mixture
#'
#' @description
#' Fits the Dail-Madsen (2011) open-population N-mixture (Poisson initial
#' abundance, binomial survival, Poisson recruitment, binomial detection) by
#' maximising the exact HMM forward marginal likelihood with an analytic gradient
#' (BFGS). The observed-information covariance is the inverse of the negative
#' finite-difference Jacobian of the analytic gradient at the mode. All four arms
#' (initial abundance `lambda`, detection `p`, survival `omega`, recruitment
#' `gamma`) are site-level.
#'
#' @param y_flat Integer vector, flattened detection counts (site-major, then
#'   season, then visit; `-1` marks a missing visit).
#' @param n_sites,T,J Sites, primary seasons, secondary visits per season.
#' @param K_max Abundance-state truncation (states `0..K_max`).
#' @param X_lambda,X_p,X_omega,X_gamma Numeric `[n_sites x p_arm]` design matrices.
#' @param max_iter BFGS iteration budget (default 300).
#' @param tol Convergence tolerance (`optim` `reltol`, default 1e-8).
#' @param verbose Print convergence status.
#'
#' @return A list of class `dyn_abun_fit` with `beta_lambda`, `beta_p`,
#'   `beta_omega`, `beta_gamma`, `log_lik`, `vcov`, `H_obs`, per-site `mean_N1`,
#'   and `converged`.
#'
#' @references Dail, D., Madsen, L. (2011). Models for estimating abundance from
#'   repeated counts of an open metapopulation. *Biometrics* 67, 577-587.
dyn_abun_laplace <- function(y_flat, n_sites, T, J, K_max,
                             X_lambda, X_p, X_omega, X_gamma,
                             max_iter = 300L, tol = 1e-8, verbose = FALSE) {
  y_flat <- as.integer(y_flat)
  K_max <- as.integer(K_max)
  p <- c(ncol(X_lambda), ncol(X_p), ncol(X_omega), ncol(X_gamma))
  off <- cumsum(c(0L, p))
  idx <- list(lambda = off[1] + seq_len(p[1]), p = off[2] + seq_len(p[2]),
              omega = off[3] + seq_len(p[3]), gamma = off[4] + seq_len(p[4]))

  eval_cpp <- function(theta) {
    cpp_dyn_abun_total_log_lik(
      y_flat, n_sites, T, J, K_max,
      as.numeric(X_lambda %*% theta[idx$lambda]),
      as.numeric(X_p      %*% theta[idx$p]),
      as.numeric(X_omega  %*% theta[idx$omega]),
      as.numeric(X_gamma  %*% theta[idx$gamma]))
  }
  grad_design <- function(out) {
    g <- numeric(sum(p))
    g[idx$lambda] <- as.numeric(crossprod(X_lambda, out$grad_eta_lambda))
    g[idx$p]      <- as.numeric(crossprod(X_p,      out$grad_eta_p))
    g[idx$omega]  <- as.numeric(crossprod(X_omega,  out$grad_eta_omega))
    g[idx$gamma]  <- as.numeric(crossprod(X_gamma,  out$grad_eta_gamma))
    g
  }
  neg_ll   <- function(theta) -eval_cpp(theta)$log_lik
  neg_grad <- function(theta) -grad_design(eval_cpp(theta))

  # Warm start: naive initial abundance from the first-season max count, p ~ 0.5,
  # omega ~ 0.6, gamma ~ 0.5.
  J0 <- J; first_counts <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    seg <- y_flat[((i - 1L) * T * J0) + seq_len(J0)]      # season 1 visits
    seg <- seg[seg >= 0]
    first_counts[i] <- if (length(seg)) max(seg) else 0
  }
  theta0 <- numeric(sum(p))
  theta0[idx$lambda[1]] <- log(max(mean(first_counts) / 0.5, 0.5) + 0.5)
  theta0[idx$p[1]]      <- 0
  theta0[idx$omega[1]]  <- stats::qlogis(0.6)
  theta0[idx$gamma[1]]  <- log(0.5)

  # Progress + ETA (gcol33/tulpaObs#43); ON by default. BFGS calls the gradient
  # ~once per quasi-Newton iteration, so ticking there approximates iteration
  # progress (maxit is the ETA denominator); finalised after optim returns.
  .prog <- tulpa:::.tulpa_iter_progress("dyn-abun-laplace", as.integer(max_iter), unit = "iter")
  neg_grad_p <- function(theta) { .prog$tick(); neg_grad(theta) }
  opt <- stats::optim(theta0, neg_ll, neg_grad_p, method = "BFGS",
                      control = list(maxit = as.integer(max_iter), reltol = tol))
  .prog$finish()
  theta <- opt$par
  converged <- opt$convergence == 0L
  out <- eval_cpp(theta)
  Hobs <- -.fp_fd_jacobian(function(th) grad_design(eval_cpp(th)), theta)
  vcov <- tryCatch(solve(Hobs), error = function(e)
    matrix(NA_real_, length(theta), length(theta)))

  nm <- c(paste0("lambda_", colnames(X_lambda)), paste0("p_", colnames(X_p)),
          paste0("omega_", colnames(X_omega)), paste0("gamma_", colnames(X_gamma)))
  dimnames(vcov) <- list(nm, nm); dimnames(Hobs) <- list(nm, nm)
  if (!converged) {
    warning(sprintf("dyn_abun_laplace BFGS did not converge (code %d).", opt$convergence),
            call. = FALSE)
  }

  structure(list(
    beta_lambda = theta[idx$lambda], beta_p = theta[idx$p],
    beta_omega = theta[idx$omega], beta_gamma = theta[idx$gamma],
    means = theta, vcov = vcov, H_obs = Hobs,
    log_lik = out$log_lik, log_lik_site = out$log_lik_site, mean_N1 = out$mean_N1,
    converged = converged, n_iter = opt$counts[[1]], K_max = K_max,
    n_inadmissible = out$n_inadmissible, coef_names = nm),
    class = c("dyn_abun_fit", "list"))
}


# ---------------------------------------------------------------------------
# Fit packer
# ---------------------------------------------------------------------------

build_dyn_abun_fit <- function(raw, model) {
  pi_list <- model$process_info
  nms <- raw$coef_names
  means <- raw$means; names(means) <- nms
  vcov <- as.matrix(raw$vcov); dimnames(vcov) <- list(nms, nms)
  sds <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms
  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov); colnames(draws) <- nms
  ll <- raw$log_lik %||% NA_real_

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    log_prob = rep(ll, n_pseudo), N = sum(model$y >= 0, na.rm = TRUE)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = length(nms), fixed_names = nms,
    process_info = pi_list, model = model, spatial = NULL, method = "laplace",
    log_lik = ll, mean_N1 = raw$mean_N1, K_max = raw$K_max,
    mixture = "poisson",
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "dyn_abun")
# ---------------------------------------------------------------------------

# Per-site arm values and the marginal expected-abundance trajectory
# E[N_t] = omega E[N_{t-1}] + gamma, E[N_1] = lambda.
.tobs_fitted_dyn_abun <- function(object) {
  model <- object$model; means <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  lambda <- exp(as.vector(model$X_processes[[1]] %*% means[off[1] + seq_len(p[1])]))
  pdet   <- stats::plogis(as.vector(model$X_processes[[2]] %*% means[off[2] + seq_len(p[2])]))
  omega  <- stats::plogis(as.vector(model$X_processes[[3]] %*% means[off[3] + seq_len(p[3])]))
  gamma  <- exp(as.vector(model$X_processes[[4]] %*% means[off[4] + seq_len(p[4])]))
  T <- model$n_seasons
  EN <- matrix(0, model$n_sites, T)
  EN[, 1] <- lambda
  for (t in seq_len(T - 1L)) EN[, t + 1L] <- omega * EN[, t] + gamma
  list(lambda = lambda, p = pdet, omega = omega, gamma = gamma, EN = EN)
}

# simulate() for dyn_abun: draw the abundance trajectory and emit counts.
.tobs_simulate_dyn_abun <- function(object, nsim = 1) {
  model <- object$model; draws <- object$draws; n_draws <- nrow(draws)
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  n_sites <- model$n_sites; T <- model$n_seasons; J <- model$max_visits
  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L); th <- draws[di, ]
    lambda <- exp(as.vector(model$X_processes[[1]] %*% th[off[1] + seq_len(p[1])]))
    pdet   <- stats::plogis(as.vector(model$X_processes[[2]] %*% th[off[2] + seq_len(p[2])]))
    omega  <- stats::plogis(as.vector(model$X_processes[[3]] %*% th[off[3] + seq_len(p[3])]))
    gamma  <- exp(as.vector(model$X_processes[[4]] %*% th[off[4] + seq_len(p[4])]))
    ya <- array(0L, dim = c(n_sites, J, T))
    for (i in seq_len(n_sites)) {
      N <- stats::rpois(1L, lambda[i])
      for (t in seq_len(T)) {
        if (t > 1L) N <- stats::rbinom(1L, N, omega[i]) + stats::rpois(1L, gamma[i])
        ya[i, , t] <- stats::rbinom(J, N, pdet[i])
      }
    }
    result[[s]] <- ya
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for dyn_abun, on the marginal expected count mu_itj = E[N_t] * p.
.tobs_residuals_dyn_abun <- function(object, type = c("deviance", "pearson",
                                                    "response")) {
  type  <- match.arg(type)
  fitv  <- .tobs_fitted_dyn_abun(object)
  model <- object$model
  T <- model$n_seasons; J <- model$max_visits
  mu_t <- fitv$EN * fitv$p                       # [n_sites x T]
  y <- model$y                                   # [n_sites x J x T]
  out <- array(NA_real_, dim = dim(y))
  for (t in seq_len(T)) {
    mu <- pmax(matrix(mu_t[, t], model$n_sites, J), 1e-10)
    yt <- y[, , t]
    out[, , t] <- switch(type,
      response = yt - mu,
      pearson  = (yt - mu) / sqrt(mu),
      deviance = {
        d <- 2 * (ifelse(yt > 0, yt * log(yt / mu), 0) - (yt - mu))
        sign(yt - mu) * sqrt(pmax(d, 0))
      })
  }
  out
}

# predict() for dyn_abun: initial abundance lambda at new X (default).
.tobs_predict_dyn_abun <- function(object, X.0 = NULL, type = c("lambda", "gamma")) {
  type  <- match.arg(type)
  model <- object$model
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  k <- if (identical(type, "lambda")) 1L else 4L
  X <- X.0 %||% model$X_processes[[k]]
  exp(as.vector(X %*% object$means[off[k] + seq_len(p[k])]))
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate Dail-Madsen open N-mixture data
#'
#' Latent `N_1 ~ Poisson(lambda)`; for `t >= 2`, `N_t = Binomial(N_{t-1}, omega)
#' + Poisson(gamma)`; observed via `Binomial(N_t, p)` over `J` secondary visits in
#' each of `T` primary seasons. Returns a 3D array `[N x J x T]` suitable for
#' [tobs()] with [dyn_abun()].
#'
#' @param N Number of sites (default 150).
#' @param T Number of primary seasons (default 4).
#' @param J Number of secondary visits per season (default 3).
#' @param n_abund_covs Number of initial-abundance covariates (default 1).
#' @param beta_lambda Initial-abundance coefficients (log). Default
#'   `c(log(5), runif(n_abund_covs, -0.4, 0.4))`.
#' @param p,omega,gamma Detection, apparent-survival, and recruitment parameters
#'   (scalars; defaults 0.5, 0.6, 1.0).
#' @param seed Optional random seed.
#' @return A list with `y` (N x J x T count array), `data` (covariates), and
#'   `truth` (coefficients, per-site `lambda`, scalar `p`/`omega`/`gamma`).
#' @export
simulate_dyn_abun <- function(N = 150, T = 4, J = 3, n_abund_covs = 1,
                              beta_lambda = NULL, p = 0.5, omega = 0.6,
                              gamma = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda)) beta_lambda <- c(log(5), stats::runif(n_abund_covs, -0.4, 0.4))
  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  X_lambda <- stats::model.matrix(~ ., abund_covs)
  lambda <- exp(as.vector(X_lambda %*% beta_lambda))
  y <- array(0L, dim = c(N, J, T))
  Nmat <- matrix(0L, N, T)
  for (i in seq_len(N)) {
    Ni <- stats::rpois(1L, lambda[i])
    for (t in seq_len(T)) {
      if (t > 1L) Ni <- stats::rbinom(1L, Ni, omega) + stats::rpois(1L, gamma)
      Nmat[i, t] <- Ni
      y[i, , t] <- stats::rbinom(J, Ni, p)
    }
  }
  list(y = y, data = abund_covs,
       truth = list(beta_lambda = beta_lambda, lambda = lambda, p = p,
                    omega = omega, gamma = gamma, N = Nmat))
}
