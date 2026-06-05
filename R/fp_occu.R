# =============================================================================
# fp_occu.R — multistate false-positive occupancy family (Miller et al. 2011)
#
# Detections are classified into three states y in {0, 1, 2}: no detection,
# ambiguous detection (a true detection OR a false positive), and certain /
# confirmed detection (only possible when the site is occupied). The latent
# occupancy z marginalises in closed form (two states), so there is no EM: the
# package-internal fit maximises the exact marginal likelihood directly with an
# analytic gradient (BFGS), and the observed-information covariance is the inverse
# of the negative finite-difference Jacobian of that analytic gradient at the mode
# (exact to finite-difference precision). A NUTS path samples the same marginal.
#
# Four site-level logit arms: occupancy psi (formula), true detection p11
# (detection), false-positive rate p10 (fp_formula), and the probability a true
# detection is certain b (b_formula). Per-site math is src/fp_occu_kernel.h.
#
#   .tobs_build_fp_occu()  data binder -> model_type = "fp_occu"
#   .tobs_fit_fp_occu()    dispatch to the false-positive occupancy fit
#   fp_occu_laplace()      analytic-gradient BFGS fit over the exact marginal
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a false-positive occupancy model. `y` is an n_sites x J matrix of
# detection states in {0, 1, 2}; NA visits are dropped (like occu()). All four
# arms are site-level: occupancy (occ_formula), true detection (det_formula),
# false-positive rate (fp_formula), and certain-classification (b_formula).
.tobs_build_fp_occu <- function(occ_formula, det_formula, data, y,
                                fp_formula = ~1, b_formula = ~1) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x J) of detection states in {0, 1, 2}.",
         call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)),
         call. = FALSE)
  }
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  ok <- is.na(y_int) | (y_int %in% 0:2)
  if (!all(ok)) {
    stop("y must contain only the detection states 0 (none), 1 (ambiguous), ",
         "2 (certain), or NA.", call. = FALSE)
  }
  n_sites <- nrow(y_int)

  bind <- .tobs_bind_formulas(list(psi = occ_formula, p11 = det_formula,
                                   p10 = fp_formula, b = b_formula), data)
  X_psi <- model.matrix(bind$fe$psi, data)
  X_p11 <- model.matrix(bind$fe$p11, data)
  X_p10 <- model.matrix(bind$fe$p10, data)
  X_b   <- model.matrix(bind$fe$b,   data)

  # Long form (valid visits only), site-major.
  site_mat <- matrix(seq_len(n_sites), n_sites, ncol(y_int))
  keep <- !is.na(as.vector(t(y_int)))
  y_long   <- as.vector(t(y_int))[keep]
  site_idx <- as.vector(t(site_mat))[keep]

  structure(list(
    model_type = "fp_occu",
    y          = y_int,
    y_long     = as.integer(y_long),
    site_idx   = as.integer(site_idx),
    X_processes = list(X_psi, X_p11, X_p10, X_b),
    formulas   = list(psi = bind$fe$psi, p11 = bind$fe$p11,
                      p10 = bind$fe$p10, b = bind$fe$b),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = n_sites,
    max_visits = ncol(y_int),
    process_info = list(
      list(name = "psi", p = ncol(X_psi), coef_names = colnames(X_psi), link = "logit"),
      list(name = "p11", p = ncol(X_p11), coef_names = colnames(X_p11), link = "logit"),
      list(name = "p10", p = ncol(X_p10), coef_names = colnames(X_p10), link = "logit"),
      list(name = "b",   p = ncol(X_b),   coef_names = colnames(X_b),   link = "logit")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter (called from .tobs_fit_model for model_type == "fp_occu")
# ---------------------------------------------------------------------------

.tobs_fit_fp_occu <- function(model, max_iter = 200L, tol = 1e-8,
                              sigma.beta = NULL, verbose = TRUE) {
  raw <- fp_occu_laplace(
    y        = model$y_long,
    site_idx = model$site_idx,
    X_psi    = model$X_processes[[1]],
    X_p11    = model$X_processes[[2]],
    X_p10    = model$X_processes[[3]],
    X_b      = model$X_processes[[4]],
    max_iter = as.integer(max_iter), tol = as.numeric(tol),
    sigma_beta = sigma.beta, verbose = isTRUE(verbose))
  build_fp_occu_fit(raw, model)
}


# ---------------------------------------------------------------------------
# Laplace fit (analytic-gradient BFGS over the exact marginal)
# ---------------------------------------------------------------------------

# Central-difference Jacobian of a vector-valued function (the analytic gradient),
# used to assemble the observed information at the mode. Exact to FD precision.
.fp_fd_jacobian <- function(fn, x, h = 1e-5) {
  p <- length(x)
  J <- matrix(0, p, p)
  for (i in seq_len(p)) {
    xp <- x; xm <- x; xp[i] <- xp[i] + h; xm[i] <- xm[i] - h
    J[, i] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  0.5 * (J + t(J))      # symmetrise (the marginal Hessian is symmetric)
}

#' Maximum-likelihood fit of the multistate false-positive occupancy model
#'
#' @description
#' Fits the Miller et al. (2011) false-positive occupancy model with confirmed
#' detections (`y` in `{0, 1, 2}`) by maximising the exact two-state marginal
#' likelihood with an analytic gradient (BFGS). The observed-information
#' covariance is the inverse of the negative finite-difference Jacobian of the
#' analytic gradient at the mode. All four arms (occupancy `psi`, true detection
#' `p11`, false-positive `p10`, certain-classification `b`) are site-level logit
#' predictors.
#'
#' @param y Integer vector of detection states (`0`/`1`/`2`), long form (valid
#'   visits only).
#' @param site_idx Integer vector, same length as `y`, 1-based site index.
#' @param X_psi,X_p11,X_p10,X_b Numeric `[n_sites x p_arm]` design matrices for
#'   the four arms.
#' @param sigma_beta Optional Gaussian prior SD on the coefficients (a mild ridge
#'   for stability); `NULL` (default) is the unpenalised MLE.
#' @param max_iter BFGS iteration budget (default 200).
#' @param tol Convergence tolerance (`optim` `reltol`, default 1e-8).
#' @param verbose Print convergence status.
#'
#' @return A list of class `fp_occu_fit` with `beta_psi`, `beta_p11`,
#'   `beta_p10`, `beta_b`, `log_lik`, `vcov`, `H_obs`, per-site posterior
#'   occupancy `w1`, and `converged`.
#'
#' @references
#' Miller, D. A. W., et al. (2011). Improving occupancy estimation when two types
#'   of observational error occur. *Ecology* 92, 1422-1428.
fp_occu_laplace <- function(y, site_idx, X_psi, X_p11, X_p10, X_b,
                            sigma_beta = NULL, max_iter = 200L, tol = 1e-8,
                            verbose = FALSE) {
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  for (nm in c("X_psi", "X_p11", "X_p10", "X_b")) {
    if (!is.matrix(get(nm))) stop(sprintf("`%s` must be a numeric matrix.", nm), call. = FALSE)
  }
  n_sites <- nrow(X_psi)
  p <- c(ncol(X_psi), ncol(X_p11), ncol(X_p10), ncol(X_b))
  off <- cumsum(c(0L, p))
  idx <- list(psi = off[1] + seq_len(p[1]), p11 = off[2] + seq_len(p[2]),
              p10 = off[3] + seq_len(p[3]), b = off[4] + seq_len(p[4]))

  ridge <- if (is.null(sigma_beta)) 0 else 1 / sigma_beta^2
  eval_cpp <- function(theta) {
    eta_psi <- as.numeric(X_psi %*% theta[idx$psi])
    eta_p11 <- as.numeric(X_p11 %*% theta[idx$p11])
    eta_p10 <- as.numeric(X_p10 %*% theta[idx$p10])
    eta_b   <- as.numeric(X_b   %*% theta[idx$b])
    cpp_fp_occu_total_log_lik(y, site_idx, eta_psi, eta_p11, eta_p10, eta_b)
  }
  grad_design <- function(out, theta) {
    g <- numeric(length(theta))
    g[idx$psi] <- as.numeric(crossprod(X_psi, out$grad_eta_psi))
    g[idx$p11] <- as.numeric(crossprod(X_p11, out$grad_eta_p11))
    g[idx$p10] <- as.numeric(crossprod(X_p10, out$grad_eta_p10))
    g[idx$b]   <- as.numeric(crossprod(X_b,   out$grad_eta_b))
    if (ridge > 0) g <- g - ridge * theta
    g
  }
  neg_ll  <- function(theta) {
    val <- -eval_cpp(theta)$log_lik
    if (ridge > 0) val <- val + 0.5 * ridge * sum(theta^2)
    val
  }
  neg_grad <- function(theta) -grad_design(eval_cpp(theta), theta)

  # Warm start: naive occupancy, p11 ~ 0.5, p10 small, b moderate.
  ny <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), function(v) any(v >= 1))
  naive <- mean(ny[!is.na(ny)])
  theta0 <- numeric(sum(p))
  theta0[idx$psi[1]] <- stats::qlogis(min(max(naive, 0.05), 0.95))
  theta0[idx$p11[1]] <- 0                          # p11 ~ 0.5
  theta0[idx$p10[1]] <- stats::qlogis(0.05)        # small false-positive rate
  theta0[idx$b[1]]   <- 0                          # b ~ 0.5

  opt <- stats::optim(theta0, neg_ll, neg_grad, method = "BFGS",
                      control = list(maxit = as.integer(max_iter), reltol = tol))
  theta <- opt$par
  converged <- opt$convergence == 0L
  if (!converged && verbose) {
    warning(sprintf("fp_occu_laplace BFGS did not converge (code %d).", opt$convergence),
            call. = FALSE)
  }

  out  <- eval_cpp(theta)
  Hobs <- -.fp_fd_jacobian(function(th) grad_design(eval_cpp(th), th), theta)
  vcov <- tryCatch(solve(Hobs), error = function(e) matrix(NA_real_, length(theta), length(theta)))

  nm <- c(paste0("psi_", colnames(X_psi)), paste0("p11_", colnames(X_p11)),
          paste0("p10_", colnames(X_p10)), paste0("b_", colnames(X_b)))
  dimnames(vcov) <- list(nm, nm); dimnames(Hobs) <- list(nm, nm)

  structure(list(
    beta_psi = theta[idx$psi], beta_p11 = theta[idx$p11],
    beta_p10 = theta[idx$p10], beta_b = theta[idx$b],
    means = theta, vcov = vcov, H_obs = Hobs,
    log_lik = out$log_lik, log_lik_site = out$log_lik_site, w1 = out$w1,
    converged = converged, n_iter = opt$counts[[1]],
    grad_norm = sqrt(sum(neg_grad(theta)^2)), n_sites = n_sites,
    coef_names = nm), class = c("fp_occu_fit", "list"))
}


# ---------------------------------------------------------------------------
# Fit packer
# ---------------------------------------------------------------------------

build_fp_occu_fit <- function(raw, model) {
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
    log_prob = rep(ll, n_pseudo),
    N = length(model$y_long)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = length(nms), fixed_names = nms,
    process_info = pi_list,
    model = model, spatial = NULL, method = "laplace",
    log_lik = ll, w1 = raw$w1,
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "fp_occu")
# ---------------------------------------------------------------------------

# Per-site fitted occupancy psi, true detection p11, false-positive p10, certain
# prob b, and posterior occupancy z = P(z = 1 | y).
.tobs_fitted_fp_occu <- function(object) {
  model <- object$model; means <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  lp <- function(k) stats::plogis(as.vector(model$X_processes[[k]] %*%
                                            means[off[k] + seq_len(p[k])]))
  list(psi = lp(1), p11 = lp(2), p10 = lp(3), b = lp(4), z = object$w1)
}

# simulate() for fp_occu: draw z, then per visit emit a multistate detection.
.tobs_simulate_fp_occu <- function(object, nsim = 1) {
  model <- object$model; draws <- object$draws; n_draws <- nrow(draws)
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  n_sites <- model$n_sites; J <- model$max_visits
  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L)
    th <- draws[di, ]
    psi <- stats::plogis(as.vector(model$X_processes[[1]] %*% th[off[1] + seq_len(p[1])]))
    p11 <- stats::plogis(as.vector(model$X_processes[[2]] %*% th[off[2] + seq_len(p[2])]))
    p10 <- stats::plogis(as.vector(model$X_processes[[3]] %*% th[off[3] + seq_len(p[3])]))
    b   <- stats::plogis(as.vector(model$X_processes[[4]] %*% th[off[4] + seq_len(p[4])]))
    z <- stats::rbinom(n_sites, 1L, psi)
    y_sim <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) {
      if (z[i] == 1L) {
        det <- stats::rbinom(J, 1L, p11[i])
        cert <- stats::rbinom(J, 1L, b[i])
        y_sim[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
      } else {
        y_sim[i, ] <- stats::rbinom(J, 1L, p10[i])    # false positives -> state 1
      }
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for fp_occu, on the marginal probability of any detection
# (y >= 1): P(y_ij >= 1) = psi * p11 + (1 - psi) * p10.
.tobs_residuals_fp_occu <- function(object, type = c("deviance", "pearson",
                                                   "response")) {
  type  <- match.arg(type)
  fitv  <- .tobs_fitted_fp_occu(object)
  model <- object$model
  mu <- fitv$psi * fitv$p11 + (1 - fitv$psi) * fitv$p10
  y_any <- (model$y >= 1) * 1
  mu_mat <- matrix(mu, model$n_sites, model$max_visits)
  mu_mat <- pmin(pmax(mu_mat, 1e-10), 1 - 1e-10)
  switch(type,
    response = y_any - mu_mat,
    pearson  = (y_any - mu_mat) / sqrt(mu_mat * (1 - mu_mat)),
    deviance = {
      d <- 2 * (ifelse(y_any == 1, log(1 / mu_mat), log(1 / (1 - mu_mat))))
      sign(y_any - mu_mat) * sqrt(pmax(d, 0))
    })
}

# predict() for fp_occu: occupancy psi (default) or true detection p11 at new X.
.tobs_predict_fp_occu <- function(object, X.0 = NULL, type = c("psi", "p11")) {
  type  <- match.arg(type)
  model <- object$model
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  k <- if (identical(type, "psi")) 1L else 2L
  X <- X.0 %||% model$X_processes[[k]]
  stats::plogis(as.vector(X %*% object$means[off[k] + seq_len(p[k])]))
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate multistate false-positive occupancy data
#'
#' Latent occupancy `z_i ~ Bernoulli(psi_i)` with `logit psi = X_psi beta_psi`,
#' observed through `J` replicate visits emitting a detection state `y in {0,1,2}`:
#' at an occupied site a visit detects with probability `p11` and the detection is
#' certain (state 2) with probability `b` else ambiguous (state 1); at an
#' unoccupied site a visit yields a false-positive ambiguous detection (state 1)
#' with probability `p10`. Returns an `N x J` integer matrix in `{0, 1, 2}`
#' suitable for [tobs()] with [fp_occu()].
#'
#' @param N Number of sites (default 300).
#' @param J Number of visits (default 5).
#' @param n_occ_covs Number of occupancy covariates (default 1).
#' @param beta_psi Occupancy coefficients (logit). Default
#'   `c(qlogis(0.5), runif(n_occ_covs, -0.6, 0.6))`.
#' @param p11,p10,b True detection, false-positive, and certain-classification
#'   probabilities (scalars; defaults 0.6, 0.05, 0.5).
#' @param seed Optional random seed.
#' @return A list with `y` (N x J state matrix), `data` (occupancy covariates),
#'   and `truth` (coefficients, per-site `psi`, scalar `p11`/`p10`/`b`, latent `z`).
#' @export
simulate_fp_occu <- function(N = 300, J = 5, n_occ_covs = 1, beta_psi = NULL,
                             p11 = 0.6, p10 = 0.05, b = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_psi)) beta_psi <- c(stats::qlogis(0.5), stats::runif(n_occ_covs, -0.6, 0.6))
  occ_covs <- data.frame(matrix(stats::rnorm(N * n_occ_covs), N, n_occ_covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n_occ_covs))
  X_psi <- stats::model.matrix(~ ., occ_covs)
  psi <- stats::plogis(as.vector(X_psi %*% beta_psi))
  z <- stats::rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- stats::rbinom(J, 1L, p11); cert <- stats::rbinom(J, 1L, b)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else {
      y[i, ] <- stats::rbinom(J, 1L, p10)
    }
  }
  list(y = y, data = occ_covs,
       truth = list(beta_psi = beta_psi, psi = psi, p11 = p11, p10 = p10,
                    b = b, z = z))
}
