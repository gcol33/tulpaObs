# double_observer.R - Double-observer abundance (unmarked multinomPois with the
# independent double-observer pi-function). Two observers survey a site
# independently; each individual is recorded by observer 1 only, observer 2 only,
# or both. With site abundance N_i ~ Poisson(lambda_i) and per-observer detection
# p1_i, p2_i, the observable cell counts are, by Poisson-multinomial thinning,
# INDEPENDENT Poissons -- no latent-N summation is needed:
#
#   n10_i ~ Poisson(lambda_i pi10_i),  pi10 = p1 (1 - p2)     (observer 1 only)
#   n01_i ~ Poisson(lambda_i pi01_i),  pi01 = (1 - p1) p2     (observer 2 only)
#   n11_i ~ Poisson(lambda_i pi11_i),  pi11 = p1 p2           (both observers)
#
# so the marginal is a closed-form product of three Poissons per site. The three
# cell means identify (lambda, p1, p2): pi10/pi11 = (1-p2)/p2, pi01/pi11 =
# (1-p1)/p1, and lambda from pi11. Fit maximises the exact marginal (optim BFGS)
# with an observed-information vcov -- the royle_nichols() recipe, here with a
# closed-form Poisson emission and no latent-state truncation.
#
#   log lambda_i = X_lambda_i . beta_lambda   (abundance)
#   logit p_k_i  = X_det_i    . beta_p_k       (per-observer detection, k = 1, 2)
#
#   .tobs_build_double_observer()   data binder -> model_type = "double_observer"
#   .tobs_fit_double_observer()     optim over the closed-form Poisson marginal
#   .dispatch_double_observer()     tobs() entry (bind + fit + assemble)

# Per-site marginal log-likelihood: three independent Poissons on the observable
# cells. lambda / p1 / p2 length n_sites; n10 / n01 / n11 the per-site cell counts.
.dobs_site_loglik <- function(lambda, p1, p2, n10, n01, n11) {
  p1 <- pmin(pmax(p1, 1e-10), 1 - 1e-10)
  p2 <- pmin(pmax(p2, 1e-10), 1 - 1e-10)
  m10 <- lambda * p1 * (1 - p2)
  m01 <- lambda * (1 - p1) * p2
  m11 <- lambda * p1 * p2
  stats::dpois(n10, m10, log = TRUE) +
    stats::dpois(n01, m01, log = TRUE) +
    stats::dpois(n11, m11, log = TRUE)
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is an n_sites x 3 matrix of cell counts in column order
# (observer-1-only, observer-2-only, both). abund_formula models log lambda;
# det_formula the shared site-level per-observer detection design.
.tobs_build_double_observer <- function(abund_formula, det_formula, data, y) {
  y <- as.matrix(y)
  if (ncol(y) != 3L) {
    stop("double_observer() y must be an N x 3 matrix of cell counts ",
         "(observer-1-only, observer-2-only, both).", call. = FALSE)
  }
  if (any(y < 0 | y != round(y), na.rm = TRUE)) {
    stop("double_observer() y must be non-negative integer counts.",
         call. = FALSE)
  }
  n_sites <- nrow(y)
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  storage.mode(y) <- "integer"

  bind     <- .tobs_bind_formulas(list(lambda = abund_formula, p = det_formula),
                                  data)
  X_lambda <- stats::model.matrix(bind$fe$lambda, data)
  X_det    <- stats::model.matrix(bind$fe$p, data)

  structure(list(
    model_type  = "double_observer",
    y           = y,
    n10         = as.integer(y[, 1L]),
    n01         = as.integer(y[, 2L]),
    n11         = as.integer(y[, 3L]),
    X_processes = list(X_lambda, X_det),
    formulas    = list(lambda = bind$fe$lambda, p = bind$fe$p),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda),
           coef_names = colnames(X_lambda), link = "log"),
      list(name = "p1", p = ncol(X_det), coef_names = colnames(X_det),
           link = "logit"),
      list(name = "p2", p = ncol(X_det), coef_names = colnames(X_det),
           link = "logit")
    )
  ), class = "tobs_model")
}

.dobs_unpack <- function(theta, model) {
  X_lambda <- model$X_processes[[1L]]; X_det <- model$X_processes[[2L]]
  p_lam <- ncol(X_lambda); p_det <- ncol(X_det)
  lambda <- exp(as.vector(X_lambda %*% theta[seq_len(p_lam)]))
  p1 <- stats::plogis(as.vector(X_det %*% theta[p_lam + seq_len(p_det)]))
  p2 <- stats::plogis(as.vector(X_det %*% theta[p_lam + p_det + seq_len(p_det)]))
  list(lambda = lambda, p1 = p1, p2 = p2)
}

# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_double_observer <- function(model, verbose = TRUE, ...) {
  X_lambda <- model$X_processes[[1L]]; X_det <- model$X_processes[[2L]]
  p_lam <- ncol(X_lambda); p_det <- ncol(X_det)
  n10 <- model$n10; n01 <- model$n01; n11 <- model$n11

  nll <- function(theta) {
    up  <- .dobs_unpack(theta, model)
    ll  <- .dobs_site_loglik(up$lambda, up$p1, up$p2, n10, n01, n11)
    val <- -sum(ll)
    if (is.finite(val)) val else 1e10
  }

  # Moment init: total detected / a rough detection guess seed lambda; the
  # both-vs-only ratio seeds the detection probabilities.
  ntot   <- n10 + n01 + n11
  p_hat  <- min(max(2 * mean(n11) / (mean(n10) + mean(n01) + 2 * mean(n11) + 1e-6),
                    0.1), 0.9)
  lam0   <- mean(ntot) / max(1 - (1 - p_hat)^2, 0.1)
  init <- c(log(max(lam0, 1e-2)), rep(0, p_lam - 1L),
            stats::qlogis(p_hat), rep(0, p_det - 1L),
            stats::qlogis(p_hat), rep(0, p_det - 1L))

  opt <- stats::optim(init, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 500L))
  converged <- opt$convergence == 0L

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  means <- opt$par; names(means) <- par_names
  V <- tryCatch(solve(opt$hessian),
                error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = model$n_sites),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = model$process_info,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    convergence  = list(converged = converged, n_iter = opt$counts[[1L]])
  )), class = c("tobs_fit", "tulpa_fit"))
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_double_observer <- function(formula, data, family, detection, y, visits,
                                      engine, priors, control,
                                      approx = "gaussian_laplace",
                                      correction = "none", ...) {
  if (is.null(detection))
    stop("double_observer() requires a `detection` formula (the shared ",
         "site-level per-observer detection model).", call. = FALSE)
  if (is.null(y))
    stop("double_observer() requires `y` (an N x 3 cell-count matrix: ",
         "observer-1-only, observer-2-only, both).", call. = FALSE)
  if (!is.null(visits))
    stop("double_observer() detection is site-level; visit-level detection ",
         "covariates (`visits`) are not yet supported.", call. = FALSE)
  if (!identical(.map_engine(engine, family = "double_observer"), "laplace"))
    stop("double_observer() supports method = \"laplace\" only.", call. = FALSE)
  model <- .tobs_build_double_observer(
    abund_formula = formula, det_formula = detection, data = data, y = y)
  .tobs_fit_double_observer(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# fitted(): per-site abundance lambda, per-observer detection p1 / p2, and the
# expected observable cell means.
.tobs_fitted_double_observer <- function(object) {
  up <- .dobs_unpack(object$means, object$model)
  list(lambda = up$lambda, p1 = up$p1, p2 = up$p2,
       cell10 = up$lambda * up$p1 * (1 - up$p2),
       cell01 = up$lambda * (1 - up$p1) * up$p2,
       cell11 = up$lambda * up$p1 * up$p2)
}

# predict(): abundance (default) or the per-observer detection probabilities.
.tobs_predict_double_observer <- function(object, newdata = NULL,
                                          type = c("abundance", "detection")) {
  type  <- match.arg(type)
  model <- object$model
  Xl <- if (is.null(newdata)) model$X_processes[[1L]]
        else stats::model.matrix(model$formulas$lambda, newdata)
  Xd <- if (is.null(newdata)) model$X_processes[[2L]]
        else stats::model.matrix(model$formulas$p, newdata)
  p_lam <- ncol(model$X_processes[[1L]]); p_det <- ncol(model$X_processes[[2L]])
  b <- object$means
  if (identical(type, "detection")) {
    cbind(p1 = stats::plogis(as.vector(Xd %*% b[p_lam + seq_len(p_det)])),
          p2 = stats::plogis(as.vector(Xd %*% b[p_lam + p_det + seq_len(p_det)])))
  } else {
    exp(as.vector(Xl %*% b[seq_len(p_lam)]))
  }
}

# residuals(): per-site Pearson residual on the total detected count against its
# expected value lambda * (1 - (1 - p1)(1 - p2)).
.tobs_residuals_double_observer <- function(object, type) {
  fv   <- .tobs_fitted_double_observer(object)
  mtot <- fv$cell10 + fv$cell01 + fv$cell11
  obs  <- object$model$n10 + object$model$n01 + object$model$n11
  eps  <- 1e-10
  res <- switch(type,
    response = obs - mtot,
    pearson  = (obs - mtot) / sqrt(mtot + eps),
    deviance = sign(obs - mtot) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pmax(mtot, eps)), 0) - (obs - mtot))))
  list(occ = res, det = NULL)
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws.
.tobs_ploglik_double_observer <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  n10 <- model$n10; n01 <- model$n01; n11 <- model$n11
  t(vapply(seq_len(nrow(draws)), function(d) {
    up <- .dobs_unpack(draws[d, ], model)
    .dobs_site_loglik(up$lambda, up$p1, up$p2, n10, n01, n11)
  }, numeric(model$n_sites)))
}

# Posterior replicate cell counts: draw a coefficient vector, then per site the
# three independent Poisson cell counts.
.tobs_simulate_double_observer <- function(object, nsim = 1) {
  model <- object$model
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    up  <- .dobs_unpack(object$draws[idx, ], model)
    cbind(stats::rpois(model$n_sites, up$lambda * up$p1 * (1 - up$p2)),
          stats::rpois(model$n_sites, up$lambda * (1 - up$p1) * up$p2),
          stats::rpois(model$n_sites, up$lambda * up$p1 * up$p2))
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate a double-observer abundance data set
#'
#' Draws from the [double_observer()] model: site abundance
#' `N ~ Poisson(lambda)` observed by two independent observers with detection
#' `p1` / `p2`, recorded as the three observable cell counts.
#'
#' @param N Number of sites (default 200).
#' @param n_abund_covs,n_det_covs Number of abundance / detection covariates.
#' @param beta_lambda Log-abundance coefficients `c(intercept, slopes...)`.
#'   Default `c(log(8), runif(n_abund_covs, -0.5, 0.5))`.
#' @param beta_p1,beta_p2 Per-observer detection coefficients (logit). Default
#'   moderate detection.
#' @param seed Optional random seed.
#' @return A list with `y` (N x 3 cell counts), `data`, and `truth`.
#' @export
simulate_double_observer <- function(N = 200, n_abund_covs = 1, n_det_covs = 1,
                                     beta_lambda = NULL, beta_p1 = NULL,
                                     beta_p2 = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda))
    beta_lambda <- c(log(8), stats::runif(n_abund_covs, -0.5, 0.5))
  if (is.null(beta_p1))
    beta_p1 <- c(stats::qlogis(0.5), stats::runif(n_det_covs, -0.3, 0.3))
  if (is.null(beta_p2))
    beta_p2 <- c(stats::qlogis(0.45), stats::runif(n_det_covs, -0.3, 0.3))

  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  det_covs <- data.frame(matrix(stats::rnorm(N * n_det_covs), N, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  data <- cbind(abund_covs, det_covs)

  X_lambda <- stats::model.matrix(~ ., abund_covs)
  X_det    <- stats::model.matrix(~ ., det_covs)
  lambda <- exp(as.vector(X_lambda %*% beta_lambda))
  p1 <- plogis(as.vector(X_det %*% beta_p1))
  p2 <- plogis(as.vector(X_det %*% beta_p2))
  N_lat <- stats::rpois(N, lambda)

  y <- matrix(0L, N, 3L)
  for (i in seq_len(N)) {
    if (N_lat[i] == 0L) next
    d1 <- stats::rbinom(N_lat[i], 1L, p1[i])   # seen by observer 1
    d2 <- stats::rbinom(N_lat[i], 1L, p2[i])   # seen by observer 2
    y[i, 1L] <- sum(d1 == 1L & d2 == 0L)
    y[i, 2L] <- sum(d1 == 0L & d2 == 1L)
    y[i, 3L] <- sum(d1 == 1L & d2 == 1L)
  }
  colnames(y) <- c("obs1_only", "obs2_only", "both")
  list(y = y, data = data,
       truth = list(beta_lambda = beta_lambda, beta_p1 = beta_p1,
                    beta_p2 = beta_p2, lambda = lambda, p1 = p1, p2 = p2,
                    N = N_lat))
}
