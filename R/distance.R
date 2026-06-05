# =============================================================================
# distance.R — binned distance-sampling abundance family
#
# Latent abundance N_i ~ Poisson(lambda_i) (or NegBin) in a covered region, each
# individual at distance x with detection g(x; sigma_i[, b]) under a half-normal
# or hazard-rate key. With B distance bins the detected counts are multinomial
# over (bin 1, ..., bin B, undetected) with cell probabilities pi_b = integral of
# g f over the bin and 1 - sum_b pi_b; the latent N is summed out in closed form
# (truncation K_max), so -- as for the static N-mixture and removal families --
# there is no EM: the package-internal `distance_laplace()` fits the marginal
# directly with analytic gradients and observed-Fisher curvature, and a NUTS path
# samples the same marginal. The abundance arm and NB dispersion are the shared
# count-marginal math; the detection arm (a site-level log-scale sigma and an
# optional scalar hazard shape, integrated by quadrature) is distance-specific.
#
#   .tobs_build_distance()  data binder -> model_type = "distance"
#   .tobs_fit_distance()    dispatch to the distance Laplace fit
#   distance_laplace()      R wrapper over cpp_distance_laplace_fixed
#   build_distance_fit()    pack the raw fit into a tobs_fit
# =============================================================================

.dist_key_code      <- function(key) switch(key, halfnorm = 0L, hazard = 1L,
                                            stop("unknown distance key"))
.dist_transect_code <- function(transect) switch(transect, line = 0L, point = 1L,
                                                 stop("unknown transect type"))


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a distance-sampling model. The abundance predictor (X_lambda) and the
# detection-scale predictor (X_sigma) are both site-level (one row per site); the
# response `y` is an n_sites x n_bins matrix of per-bin detected counts. The bin
# cut points and transect geometry define the detection integrals.
.tobs_build_distance <- function(abund_formula, det_formula, data, y,
                                 cutpoints, key = "halfnorm",
                                 transect = "line", mixture = "poisson",
                                 quad_order = 64L) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x n_bins) of integer distance-bin counts.",
         call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)),
         call. = FALSE)
  }
  if (anyNA(y)) {
    stop("distance() needs complete bin counts: y must not contain NA.",
         call. = FALSE)
  }
  if (is.null(cutpoints) || length(cutpoints) != ncol(y) + 1L) {
    stop(sprintf("cutpoints must have length ncol(y) + 1 = %d (the bin edges, ",
                 "0 = c_0 < c_1 < ... < c_B).", ncol(y) + 1L), call. = FALSE)
  }
  cutpoints <- as.numeric(cutpoints)
  if (any(diff(cutpoints) <= 0) || cutpoints[1] < 0) {
    stop("cutpoints must be strictly increasing and start at >= 0.", call. = FALSE)
  }
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  if (any(y_int < 0L)) stop("y must contain nonnegative integer counts.", call. = FALSE)

  bind     <- .tobs_bind_formulas(list(lambda = abund_formula, sigma = det_formula),
                                  data)
  X_lambda <- model.matrix(bind$fe$lambda, data)
  X_sigma  <- model.matrix(bind$fe$sigma, data)

  structure(list(
    model_type = "distance",
    y          = y_int,
    X_processes = list(X_lambda, X_sigma),
    formulas   = list(lambda = bind$fe$lambda, sigma = bind$fe$sigma),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = nrow(y_int),
    n_bins     = ncol(y_int),
    cutpoints  = cutpoints,
    key        = key,
    transect   = transect,
    mixture    = mixture,
    quad_order = as.integer(quad_order),
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda),
           link = "log"),
      list(name = "sigma",  p = ncol(X_sigma),  coef_names = colnames(X_sigma),
           link = "log")
    ),
    mean_count = mean(y_int),
    max_count  = if (length(y_int)) max(rowSums(y_int)) else 0L
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter (called from .tobs_fit_model for model_type == "distance")
# ---------------------------------------------------------------------------

.tobs_fit_distance <- function(model, mixture = "poisson", K_max = NULL,
                               max_iter = 100L, tol = 1e-6, verbose = TRUE) {
  mix_code <- switch(mixture, poisson = "P", negbin = "NB",
                     stop(sprintf("Unknown mixture '%s'.", mixture), call. = FALSE))
  raw <- distance_laplace(
    y          = model$y,
    X_lambda   = model$X_processes[[1]],
    X_sigma    = model$X_processes[[2]],
    cutpoints  = model$cutpoints,
    key        = model$key,
    transect   = model$transect,
    mixture    = mix_code,
    K_max      = K_max,
    quad_order = model$quad_order,
    max_iter   = as.integer(max_iter),
    tol        = as.numeric(tol),
    verbose    = isTRUE(verbose))
  build_distance_fit(raw, model)
}


# ---------------------------------------------------------------------------
# Laplace fit (R wrapper over cpp_distance_laplace_fixed)
# ---------------------------------------------------------------------------

#' Laplace fit of the binned distance-sampling abundance model
#'
#' @description
#' Maximum-likelihood fit (non-spatial, fixed effects only) of a binned
#' distance-sampling abundance model with a Poisson or negative-binomial
#' abundance mixing distribution and a half-normal or hazard-rate detection key.
#' Latent abundance `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, r)`) is summed
#' out in closed form (truncation `K_max`); the detection integrals over distance
#' bins are evaluated by Gauss-Legendre quadrature.
#'
#' @param y Integer matrix `[n_sites x n_bins]` of per-bin detected counts.
#' @param X_lambda Numeric matrix `[n_sites x p_lambda]` of abundance covariates.
#' @param X_sigma Numeric matrix `[n_sites x p_sigma]` of detection-scale
#'   covariates (the `log sigma` linear predictor).
#' @param cutpoints Numeric bin edges, length `n_bins + 1`
#'   (`0 = c_0 < c_1 < ... < c_B`).
#' @param key `"halfnorm"` (default) or `"hazard"` detection key.
#' @param transect `"line"` (default, distances uniform) or `"point"` (radial,
#'   density proportional to distance).
#' @param mixture `"P"` (Poisson, default) or `"NB"` (negative binomial).
#' @param beta_lambda_init,beta_sigma_init,eta_b_init,log_r_init Optional warm
#'   starts.
#' @param r_max Upper bound on the NB size `r` (NB only, default `1e5`).
#' @param K_max Marginal-sum truncation; defaults to `max(site total) + 100`.
#' @param quad_order Gauss-Legendre nodes per bin (default 64).
#' @param max_iter Newton iteration budget (default 100).
#' @param tol Gradient-norm convergence tolerance (default 1e-6).
#' @param verbose Print per-iteration progress.
#'
#' @return A list of class `distance_fit` with `beta_lambda`, `beta_sigma`,
#'   `shape` / `eta_b` (hazard-rate), `log_r` / `r` (NB), `log_lik`, `vcov`,
#'   `H_obs`, per-site `mean_N` / `var_N` / `p_det` / `boundary_weight`.
#'
#' @references
#' Buckland, S. T., et al. (2001). Introduction to Distance Sampling. Oxford.
#' Royle, J. A., Dawson, D. K., Bates, S. (2004). Modeling abundance effects in
#'   distance sampling. *Ecology* 85, 1591-1597.
distance_laplace <- function(y, X_lambda, X_sigma, cutpoints,
                             key = c("halfnorm", "hazard"),
                             transect = c("line", "point"),
                             mixture = c("P", "NB"),
                             beta_lambda_init = NULL, beta_sigma_init = NULL,
                             eta_b_init = NULL, log_r_init = NULL,
                             r_max = 1e5, K_max = NULL, quad_order = 64L,
                             max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  key      <- match.arg(key)
  transect <- match.arg(transect)
  mixture  <- match.arg(mixture)
  nb     <- identical(mixture, "NB")
  hazard <- identical(key, "hazard")
  if (!is.matrix(y))        stop("`y` must be an integer matrix.", call. = FALSE)
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_sigma))  stop("`X_sigma` must be a numeric matrix.", call. = FALSE)
  y <- matrix(as.integer(y), nrow(y), ncol(y))
  n_sites <- nrow(y); n_bins <- ncol(y)
  if (nrow(X_lambda) != n_sites) stop("nrow(X_lambda) must equal nrow(y).", call. = FALSE)
  if (nrow(X_sigma)  != n_sites) stop("nrow(X_sigma) must equal nrow(y).", call. = FALSE)
  if (length(cutpoints) != n_bins + 1L) stop("length(cutpoints) must equal ncol(y)+1.", call. = FALSE)

  p_lambda <- ncol(X_lambda); p_sigma <- ncol(X_sigma)
  site_tot <- rowSums(y)
  R_max    <- max(site_tot)
  if (is.null(beta_lambda_init)) {
    beta_lambda_init <- c(log(max(mean(site_tot), 0.1) + 0.1), rep(0, p_lambda - 1L))
  }
  if (is.null(beta_sigma_init)) {
    # Start sigma near the mid-range of the covered distances.
    beta_sigma_init <- c(log(stats::median(cutpoints[-1])), rep(0, p_sigma - 1L))
  }
  if (length(beta_lambda_init) != p_lambda) stop("length(beta_lambda_init) must equal ncol(X_lambda).", call. = FALSE)
  if (length(beta_sigma_init)  != p_sigma)  stop("length(beta_sigma_init) must equal ncol(X_sigma).", call. = FALSE)
  if (is.null(eta_b_init))  eta_b_init  <- log(2)    # hazard shape ~ 2 (a shoulder)
  if (is.null(log_r_init))  log_r_init  <- 0

  if (is.null(K_max)) {
    # Unlike the N-mixture / removal families, the truncation bounds the latent
    # TOTAL abundance in the covered region, which exceeds the detected total by
    # a factor of 1 / p (the undetected individuals). A "detected total + 100"
    # buffer is therefore too tight when detection is low, so the default uses a
    # multiplicative margin; the boundary-weight warning below still flags any
    # residual truncation so a wider K_max can be set explicitly.
    K_max <- as.integer(3L * R_max + 100L)
  } else {
    K_max <- as.integer(K_max)
    if (K_max < R_max) stop("`K_max` must be >= the largest per-site total.", call. = FALSE)
  }

  fit <- cpp_distance_laplace_fixed(
    y = y, X_lambda_R = X_lambda, X_sigma_R = X_sigma,
    cutpoints = as.numeric(cutpoints),
    transect = .dist_transect_code(transect), key = .dist_key_code(key),
    beta_lambda_init = as.numeric(beta_lambda_init),
    beta_sigma_init  = as.numeric(beta_sigma_init),
    eta_b_init = as.numeric(eta_b_init),
    K_max = K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose), nb = nb,
    log_r_init = as.numeric(log_r_init), theta_max = log(r_max),
    quad_order = as.integer(quad_order))

  nm_lam <- colnames(X_lambda); nm_sig <- colnames(X_sigma)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lambda))
  if (is.null(nm_sig)) nm_sig <- paste0("sig_", seq_len(p_sigma))
  names(fit$beta_lambda) <- nm_lam
  names(fit$beta_sigma)  <- nm_sig
  coef_names <- c(nm_lam, nm_sig, if (hazard) "log_shape", if (nb) "log_r")
  rownames(fit$vcov)  <- colnames(fit$vcov)  <- coef_names
  rownames(fit$H_obs) <- colnames(fit$H_obs) <- coef_names

  fit$mixture <- mixture; fit$key <- key; fit$transect <- transect
  fit$hazard <- hazard; fit$nb <- nb
  fit$cutpoints <- as.numeric(cutpoints); fit$quad_order <- as.integer(quad_order)
  if (!nb) { fit$log_r <- NA_real_; fit$r <- NA_real_ }
  fit$K_max <- K_max; fit$n_sites <- n_sites; fit$n_bins <- n_bins
  fit$call <- match.call()
  if (!fit$converged) {
    warning(sprintf("distance_laplace did not converge in %d iterations (grad_norm = %.2e).",
                    max_iter, fit$grad_norm), call. = FALSE)
  }
  if (nb && isTRUE(fit$dispersion_boundary)) {
    warning(sprintf(paste0("NB dispersion pinned at the boundary (r = r_max = %.3g); ",
            "the data are consistent with Poisson. Consider mixture = \"P\"."), r_max),
            call. = FALSE)
  }
  max_bw <- max(fit$boundary_weight, na.rm = TRUE)
  if (is.finite(max_bw) && max_bw > 1e-4) {
    warning(sprintf("Max posterior weight on N = K_max is %.2e at %d sites; raise K_max.",
                    max_bw, sum(fit$boundary_weight > 1e-4)), call. = FALSE)
  }
  class(fit) <- c("distance_fit", "list")
  fit
}


# ---------------------------------------------------------------------------
# Fit packer
# ---------------------------------------------------------------------------

# Pack a raw distance fit (Laplace `distance_laplace()` or the NUTS summary) into
# a tobs_fit. Coefficient layout: (lambda, sigma[, log_shape][, log_r]); the
# trailing scalars are model coefficients carried with an SE, left on natural
# scale by the per-process unscaler in .tobs_fit_model().
build_distance_fit <- function(raw, model) {
  pi_list <- model$process_info
  p_lam <- pi_list[[1]]$p; p_sig <- pi_list[[2]]$p
  is_nb  <- identical(raw$mixture, "NB") || identical(raw$mixture, "negbin")
  hazard <- isTRUE(raw$hazard) || identical(raw$key, "hazard")

  nms <- c(paste0("lambda_", pi_list[[1]]$coef_names),
           paste0("sigma_",  pi_list[[2]]$coef_names))
  beta_lambda <- raw$beta_lambda; beta_sigma <- raw$beta_sigma
  means <- c(as.numeric(beta_lambda), as.numeric(beta_sigma))
  if (hazard) { nms <- c(nms, "log_shape"); means <- c(means, as.numeric(raw$eta_b)) }
  if (is_nb)  { nms <- c(nms, "log_r");     means <- c(means, as.numeric(raw$log_r)) }
  names(means) <- nms

  p_total <- length(means)
  vcov <- as.matrix(raw$vcov)
  if (any(dim(vcov) != p_total)) vcov <- matrix(NA_real_, p_total, p_total)
  rownames(vcov) <- colnames(vcov) <- nms
  sds <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms

  n_fixed <- length(nms); fixed_names <- nms
  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov)
  colnames(draws) <- nms

  dispersion <- NULL
  if (is_nb && is.finite(raw$log_r %||% NA_real_)) {
    se_logr <- sqrt(pmax(vcov["log_r", "log_r"], 0))
    dispersion <- list(r = as.numeric(raw$r %||% exp(raw$log_r)),
                       log_r = as.numeric(raw$log_r),
                       r_sd = as.numeric(exp(raw$log_r) * se_logr))
  }
  shape <- if (hazard) list(shape = as.numeric(raw$shape %||% exp(raw$eta_b)),
                            eta_b = as.numeric(raw$eta_b)) else NULL

  ll <- raw$log_lik %||% NA_real_
  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    log_prob = rep(ll, n_pseudo),
    N = model$n_sites),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = n_fixed, fixed_names = fixed_names,
    process_info = pi_list,
    model = model, spatial = NULL,
    method = "laplace",
    log_lik = ll,
    K_max = raw$K_max,
    mean_N = raw$mean_N, var_N = raw$var_N, p_det = raw$p_det,
    boundary_weight = raw$boundary_weight,
    key = model$key, transect = model$transect,
    cutpoints = model$cutpoints, quad_order = model$quad_order,
    mixture = if (is_nb) "negbin" else "poisson",
    nmix_dispersion = dispersion,
    distance_shape = shape,
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# Detection-function helpers (R; used by simulate / residuals / fitted)
# ---------------------------------------------------------------------------

# Detection probability g(x) for a key, sigma (and hazard shape b).
.distance_g <- function(x, key, sigma, shape = NULL) {
  if (identical(key, "halfnorm")) return(exp(-x^2 / (2 * sigma^2)))
  # hazard-rate
  out <- 1 - exp(-(x / sigma)^(-shape))
  out[x <= 0] <- 1
  out
}

# Per-bin detection-cell probabilities pi_b = integral_bin g(x) f(x) dx for a
# single sigma, by numerical integration (a diagnostic path; the fit itself uses
# the C++ quadrature). Returns a length-n_bins vector.
.distance_pi <- function(sigma, cutpoints, key, transect, shape = NULL) {
  W <- cutpoints[length(cutpoints)]
  fdens <- if (identical(transect, "point")) function(x) 2 * x / W^2 else function(x) 1 / W
  n_bins <- length(cutpoints) - 1L
  vapply(seq_len(n_bins), function(b) {
    stats::integrate(function(x) .distance_g(x, key, sigma, shape) * fdens(x),
                     lower = cutpoints[b], upper = cutpoints[b + 1],
                     rel.tol = 1e-8)$value
  }, numeric(1))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "distance")
# ---------------------------------------------------------------------------

# Per-site fitted lambda (expected abundance), sigma (detection scale), and p
# (overall detection probability over the covered region).
.tobs_fitted_distance <- function(object) {
  model <- object$model
  means <- object$means
  p_lam <- model$process_info[[1]]$p
  p_sig <- model$process_info[[2]]$p
  beta_lambda <- means[seq_len(p_lam)]
  beta_sigma  <- means[p_lam + seq_len(p_sig)]
  lambda <- exp(as.vector(model$X_processes[[1]] %*% beta_lambda))
  sigma  <- exp(as.vector(model$X_processes[[2]] %*% beta_sigma))
  p <- object$p_det
  if (is.null(p)) {
    shape <- object$distance_shape$shape
    p <- vapply(sigma, function(s)
      sum(.distance_pi(s, model$cutpoints, model$key, model$transect, shape)),
      numeric(1))
  }
  list(lambda = lambda, sigma = sigma, p = p)
}

# simulate() for distance: draw N_i, then assign each individual to a distance
# bin (or "undetected") by the multinomial cell probabilities.
.tobs_simulate_distance <- function(object, nsim = 1) {
  model   <- object$model
  draws   <- object$draws
  n_draws <- nrow(draws)
  p_lam   <- model$process_info[[1]]$p
  p_sig   <- model$process_info[[2]]$p
  n_sites <- model$n_sites
  n_bins  <- model$n_bins
  r_size  <- object$nmix_dispersion$r
  shape   <- object$distance_shape$shape

  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L)
    beta_lambda <- draws[di, seq_len(p_lam)]
    beta_sigma  <- draws[di, p_lam + seq_len(p_sig)]
    lambda <- exp(as.vector(model$X_processes[[1]] %*% beta_lambda))
    sigma  <- exp(as.vector(model$X_processes[[2]] %*% beta_sigma))
    N <- if (!is.null(r_size) && is.finite(r_size))
      stats::rnbinom(n_sites, size = r_size, mu = lambda) else stats::rpois(n_sites, lambda)
    y_sim <- matrix(0L, n_sites, n_bins)
    for (i in seq_len(n_sites)) {
      pi_b <- .distance_pi(sigma[i], model$cutpoints, model$key, model$transect, shape)
      probs <- c(pi_b, max(1 - sum(pi_b), 0))
      if (N[i] > 0) {
        counts <- stats::rmultinom(1L, N[i], probs)
        y_sim[i, ] <- counts[seq_len(n_bins)]
      }
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for distance. The bin-b count is marginally
# y_ib ~ Poisson(lambda_i * pi_ib) under Poisson abundance; score on that mean.
.tobs_residuals_distance <- function(object, type = c("deviance", "pearson",
                                                    "response")) {
  type  <- match.arg(type)
  model <- object$model
  fitv  <- .tobs_fitted_distance(object)
  lambda <- fitv$lambda; sigma <- fitv$sigma
  shape <- object$distance_shape$shape
  n_sites <- model$n_sites; n_bins <- model$n_bins
  pi_mat <- t(vapply(sigma, function(s)
    .distance_pi(s, model$cutpoints, model$key, model$transect, shape),
    numeric(n_bins)))
  mu_mat <- pmax(lambda * pi_mat, 1e-10)
  y <- model$y
  r_mat <- switch(type,
    response = y - mu_mat,
    pearson  = (y - mu_mat) / sqrt(mu_mat),
    deviance = {
      d <- 2 * (ifelse(y > 0, y * log(y / mu_mat), 0) - (y - mu_mat))
      sign(y - mu_mat) * sqrt(pmax(d, 0))
    })
  r_mat
}

# predict() for distance: abundance lambda (density) at new X_lambda, or the
# detection scale sigma at new X_sigma. Mirrors the nmix predictor's lambda mode.
.tobs_predict_distance <- function(object, X.0 = NULL, type = c("lambda", "sigma")) {
  type  <- match.arg(type)
  model <- object$model
  p_lam <- model$process_info[[1]]$p
  p_sig <- model$process_info[[2]]$p
  if (identical(type, "lambda")) {
    X <- X.0 %||% model$X_processes[[1]]
    beta <- object$means[seq_len(p_lam)]
    return(exp(as.vector(X %*% beta)))
  }
  X <- X.0 %||% model$X_processes[[2]]
  beta <- object$means[p_lam + seq_len(p_sig)]
  exp(as.vector(X %*% beta))
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate binned distance-sampling abundance data
#'
#' Latent abundance `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, size)`) with
#' `log lambda_i = X_lambda beta_lambda`, observed through a half-normal or
#' hazard-rate detection function with `log sigma_i = X_sigma beta_sigma`. Each
#' individual lands in a distance bin (or goes undetected) by the multinomial
#' cell probabilities. Returns an `N x n_bins` integer matrix of per-bin counts
#' suitable for [tobs()] with [distance()].
#'
#' @param N Number of sites (default 200).
#' @param cutpoints Distance-bin edges (length `n_bins + 1`). Default
#'   `seq(0, 1, length.out = 6)` (five bins out to 1).
#' @param key `"halfnorm"` (default) or `"hazard"`.
#' @param transect `"line"` (default) or `"point"`.
#' @param n_abund_covs,n_sigma_covs Number of abundance / detection covariates.
#' @param beta_lambda Abundance coefficients (log scale). Default
#'   `c(log(40), runif(n_abund_covs, -0.4, 0.4))`.
#' @param beta_sigma Detection-scale coefficients (log scale). Default
#'   `c(log(0.4), runif(n_sigma_covs, -0.3, 0.3))`.
#' @param shape Hazard-rate shape `b` (`key = "hazard"` only, default 3).
#' @param mixture `"poisson"` (default) or `"negbin"`.
#' @param size Negative-binomial size `r` (`mixture = "negbin"` only, default 5).
#' @param seed Optional random seed.
#' @return A list with `y` (N x n_bins count matrix), `data` (covariates),
#'   `cutpoints`, and `truth` (coefficients, per-site `lambda` / `sigma`, latent
#'   `N`, key/transect/mixture/shape/size).
#' @export
simulate_distance <- function(N = 200, cutpoints = seq(0, 1, length.out = 6),
                              key = c("halfnorm", "hazard"),
                              transect = c("line", "point"),
                              n_abund_covs = 1, n_sigma_covs = 1,
                              beta_lambda = NULL, beta_sigma = NULL,
                              shape = 3, mixture = c("poisson", "negbin"),
                              size = 5, seed = NULL) {
  key      <- match.arg(key)
  transect <- match.arg(transect)
  mixture  <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda)) beta_lambda <- c(log(40), stats::runif(n_abund_covs, -0.4, 0.4))
  if (is.null(beta_sigma))  beta_sigma  <- c(log(0.4), stats::runif(n_sigma_covs, -0.3, 0.3))
  cutpoints <- as.numeric(cutpoints)
  n_bins <- length(cutpoints) - 1L

  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  sigma_covs <- data.frame(matrix(stats::rnorm(N * n_sigma_covs), N, n_sigma_covs))
  names(sigma_covs) <- paste0("sigma_cov", seq_len(n_sigma_covs))
  data <- cbind(abund_covs, sigma_covs)

  X_lambda <- stats::model.matrix(~ ., abund_covs)
  X_sigma  <- stats::model.matrix(~ ., sigma_covs)
  lambda <- exp(as.vector(X_lambda %*% beta_lambda))
  sigma  <- exp(as.vector(X_sigma  %*% beta_sigma))
  Nlat <- if (identical(mixture, "negbin"))
    stats::rnbinom(N, size = size, mu = lambda) else stats::rpois(N, lambda)

  sh <- if (identical(key, "hazard")) shape else NULL
  y <- matrix(0L, N, n_bins)
  for (i in seq_len(N)) {
    pi_b <- .distance_pi(sigma[i], cutpoints, key, transect, sh)
    probs <- c(pi_b, max(1 - sum(pi_b), 0))
    if (Nlat[i] > 0) {
      counts <- stats::rmultinom(1L, Nlat[i], probs)
      y[i, ] <- counts[seq_len(n_bins)]
    }
  }

  list(
    y = y, data = data, cutpoints = cutpoints,
    truth = list(beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                 lambda = lambda, sigma = sigma, N = Nlat,
                 key = key, transect = transect, mixture = mixture,
                 shape = if (identical(key, "hazard")) shape else NA_real_,
                 size = if (identical(mixture, "negbin")) size else NA_real_)
  )
}
