#' Continuous-field (SPDE) Royle (2004) N-mixture model via nested Laplace
#'
#' @description
#' Nested-Laplace fit of the spatial N-mixture model with a continuous Matern
#' (SPDE) field on the abundance arm:
#' \deqn{N_i \sim \mathrm{Poisson}(\lambda_i), \qquad
#'       y_{ij} | N_i \sim \mathrm{Binomial}(N_i, p_{ij}),}
#' \eqn{\log \lambda_i = X_\lambda^{(i)} \beta_\lambda + (A u)_i}, where
#' \eqn{u \sim \mathrm{N}(0, Q(\mathrm{range}, \sigma)^{-1})} is a Gaussian
#' Markov random field on the FEM mesh with the proper Matern precision
#' \eqn{Q}, and \eqn{A} (\eqn{n_{\mathrm{sites}} \times n_{\mathrm{mesh}}})
#' projects mesh nodes onto sites. The detection arm is
#' \eqn{\mathrm{logit}\, p_{ij} = X_p^{(ij)} \beta_p}.
#'
#' The SPDE hyperparameters \eqn{(\mathrm{range}, \sigma)} are integrated out
#' by the SAME outer grid the areal path (`nmix_laplace_icar()` /
#' `nmix_laplace_bym2()` / `nmix_laplace_car_proper()`) uses for
#' \eqn{(\tau[, \rho], \sigma)}: at each grid point the inner Newton finds the
#' joint mode of \eqn{(\beta_\lambda, \beta_p, u)} and the Laplace log-marginal
#' \eqn{\log p(y \mid \mathrm{range}_k, \sigma_k)} is accumulated, then the
#' grid is normalised into posterior weights. The Matern PC priors on
#' \eqn{(\mathrm{range}, \sigma)} (from the `spde()` term) enter each grid
#' point's log-marginal, mirroring `fit_spde()`. The precision \eqn{Q} (and its
#' \eqn{\log|Q|}) is built once per grid point on the R side via the same FEM
#' assembly the occupancy SPDE path uses (`tulpa:::.spde_precision_Q`), so the
#' C++ kernel stays agnostic to the precision parameterisation.
#'
#' Unlike the intrinsic ICAR field, \eqn{Q} is full rank, so the
#' (intercept, field-mean) direction is identified by \eqn{Q} itself: no
#' sum-to-zero centering and no constrained-covariance projection (the proper
#' CAR path is the areal analogue).
#'
#' @param y Integer vector of observed counts (long form, one entry per visit).
#' @param site_idx Integer vector, 1-based site index for each visit.
#' @param X_lambda Numeric matrix `[n_sites x p_lambda]` of abundance covariates.
#' @param X_p Numeric matrix `[n_obs x p_p]` of detection covariates.
#' @param spatial A `tobs_spatial` SPDE term (`type == "spde"`) carrying the
#'   FEM `tulpa_spec` (`A`, `C0_diag`, `G`, `n_mesh`, `nu`, `prior_range`,
#'   `prior_sigma`). The projection `A` must have `n_sites` rows.
#' @param mixture Abundance mixing distribution: `"P"` (Poisson, default) or
#'   `"NB"` (negative binomial). Under `"NB"` the NB size \eqn{r} is integrated
#'   as an additional outer grid dimension; the posterior `r_mean` / `r_sd` are
#'   reported from the grid weights.
#' @param r_grid Optional numeric vector of NB size grid points (NB only).
#'   Defaults to `exp(seq(log(0.5), log(40), length.out = 6))`.
#' @param range_grid,sigma_grid Optional numeric vectors of SPDE range / sigma
#'   grid points. Defaults centre a log-spaced grid on the PC-prior medians.
#' @param beta_lambda_init Optional warm start; default `c(log(mean(y)+0.1), 0, ...)`.
#' @param beta_p_init Optional warm start; default `rep(0, p_p)`.
#' @param u_init Optional warm start for the mesh field; default zeros.
#' @param K_max Truncation for the per-site marginal sum over N. Defaults to
#'   `max(y) + 100`.
#' @param max_iter,tol Inner Newton iteration budget and gradient-norm tol.
#' @param verbose Print per-iteration / per-grid-point progress.
#'
#' @return A list of class `nmix_spatial_fit` with the same shape as the areal
#'   fitters plus `u_mean` (the posterior-mean mesh field) and
#'   `range_mean` / `range_sd` / `sigma_mean` / `sigma_sd`.
#'
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Lindgren, F., Rue, H., Lindstrom, J. (2011). An explicit link between
#'   Gaussian fields and Gaussian Markov random fields: the SPDE approach.
#'   *JRSS-B* 73, 423-498.
#' Fuglstad, G.-A., Simpson, D., Lindgren, F., Rue, H. (2019). Constructing
#'   priors that penalize the complexity of Gaussian random fields. *JASA* 114.
#'
#' @keywords internal
nmix_laplace_spde <- function(y, site_idx, X_lambda, X_p, spatial,
                              mixture = c("P", "NB"), r_grid = NULL,
                              range_grid = NULL, sigma_grid = NULL,
                              beta_lambda_init = NULL, beta_p_init = NULL,
                              u_init = NULL, K_max = NULL,
                              max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  mixture  <- match.arg(mixture)
  y        <- as.integer(y)
  site_idx <- as.integer(site_idx)
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)

  ts <- spatial$tulpa_spec
  if (is.null(ts) || !identical(ts$type, "spde")) {
    stop("nmix_laplace_spde() requires an SPDE tulpa_spec.", call. = FALSE)
  }
  n_sites <- nrow(X_lambda)
  n_obs   <- nrow(X_p)
  p_lam   <- ncol(X_lambda)
  p_p     <- ncol(X_p)
  n_mesh  <- ts$n_mesh
  A_dense <- as.matrix(ts$A)
  if (nrow(A_dense) != n_sites) {
    stop(sprintf("SPDE projection A has %d rows but the model has %d sites.",
                 nrow(A_dense), n_sites), call. = FALSE)
  }
  if (length(y) != n_obs) stop("length(y) must equal nrow(X_p).", call. = FALSE)
  if (length(site_idx) != n_obs) stop("length(site_idx) must equal nrow(X_p).", call. = FALSE)

  if (is.null(beta_lambda_init)) {
    beta_lambda_init <- c(log(mean(y) + 0.1), rep(0, p_lam - 1L))
  }
  if (is.null(beta_p_init)) beta_p_init <- rep(0, p_p)
  if (length(beta_lambda_init) != p_lam) {
    stop("length(beta_lambda_init) must equal ncol(X_lambda).", call. = FALSE)
  }
  if (length(beta_p_init) != p_p) {
    stop("length(beta_p_init) must equal ncol(X_p).", call. = FALSE)
  }
  if (is.null(K_max)) {
    K_max <- as.integer(max(y) + 100L)
  } else {
    K_max <- as.integer(K_max)
    if (K_max < max(y)) stop("K_max must be >= max(y).", call. = FALSE)
  }
  if (!is.null(u_init) && length(u_init) != n_mesh) {
    stop("length(u_init) must equal n_mesh.", call. = FALSE)
  }

  # --- Outer (range, sigma) grid centred on the PC-prior medians ----------
  prior_range <- ts$prior_range
  prior_sigma <- ts$prior_sigma
  if (is.null(range_grid)) {
    # PC prior on range: P(range < U_r) = alpha_r -> rate on r^{-1}; median
    # range ~ U_r at alpha_r = 0.5. Centre a 5-point log grid there.
    r_med <- prior_range[1]
    range_grid <- exp(seq(log(r_med * 0.35), log(r_med * 2.5), length.out = 5L))
  }
  if (is.null(sigma_grid)) {
    s_scale <- prior_sigma[1]
    sigma_grid <- exp(seq(log(s_scale * 0.35), log(s_scale * 2.0), length.out = 5L))
  }
  if (any(range_grid <= 0)) stop("range_grid must be strictly positive.", call. = FALSE)
  if (any(sigma_grid <= 0)) stop("sigma_grid must be strictly positive.", call. = FALSE)
  r_grid_use <- .nmix_resolve_r_grid(mixture, r_grid)

  # --- Per-grid-point precision Q(range, sigma) and log|Q| ----------------
  # kappa / tau_spde from the Matern parameterisation (matches fit_spde).
  build_Q <- function(range_val, sigma_val) {
    kappa    <- sqrt(8 * ts$nu) / range_val
    tau_spde <- 1 / (sqrt(4 * pi) * kappa * sigma_val)
    Q <- tulpa:::.spde_precision_Q(ts, kappa, tau_spde)
    Q <- Matrix::forceSymmetric(Q)
    list(Q = as.matrix(Q), log_det = .spde_logdet_Q(Q))
  }

  # Outer grid: (r [NB], range, sigma). Poisson -> single r = Inf node.
  grid <- expand.grid(range = range_grid, sigma = sigma_grid,
                      r = r_grid_use, KEEP.OUT.ATTRS = FALSE)
  n_grid <- nrow(grid)
  Q_list   <- vector("list", n_grid)
  log_dets <- numeric(n_grid)
  pc_lp    <- numeric(n_grid)
  cache    <- list()
  for (k in seq_len(n_grid)) {
    key <- paste0(grid$range[k], "_", grid$sigma[k])
    if (is.null(cache[[key]])) cache[[key]] <- build_Q(grid$range[k], grid$sigma[k])
    Q_list[[k]]  <- cache[[key]]$Q
    log_dets[k]  <- cache[[key]]$log_det
    pc_lp[k]     <- tulpa:::pc_prior_log_density(grid$range[k], grid$sigma[k],
                                                 prior_range, prior_sigma)
  }
  theta_grid <- as.matrix(grid[, c("range", "sigma", "r"), drop = FALSE])

  fit <- .cpp_nmix_progress(cpp_nested_laplace_nmix_spde,
    y = y, site_idx = site_idx,
    X_lambda_R = X_lambda, X_p_R = X_p, A_R = A_dense,
    Q_list = Q_list, log_det_Q = log_dets,
    theta_grid_R = theta_grid, r_grid = as.numeric(grid$r),
    beta_lambda_init = as.numeric(beta_lambda_init),
    beta_p_init = as.numeric(beta_p_init),
    u_init = if (is.null(u_init)) NULL else as.numeric(u_init),
    K_max = K_max, max_iter = as.integer(max_iter),
    tol = as.numeric(tol), verbose = isTRUE(verbose)
  )

  # Add the PC prior to each grid log-marginal, then normalise.
  lm_post <- fit$log_marginal + pc_lp
  weights <- tulpa:::.nl_normalise_weights_safe(lm_post, "range / sigma grid")

  rng   <- .tobs_weighted_moment(weights, theta_grid[, "range"])
  sigma <- .tobs_weighted_moment(weights, theta_grid[, "sigma"])
  range_mean <- unname(rng["mean"]);   range_sd <- unname(rng["sd"])
  sigma_mean <- unname(sigma["mean"]); sigma_sd <- unname(sigma["sd"])
  disp <- .nmix_dispersion_summary(mixture, theta_grid, weights)

  modes <- fit$modes
  beta_lambda_mean <- as.numeric(crossprod(weights, modes[, seq_len(p_lam), drop = FALSE]))
  beta_p_mean      <- as.numeric(crossprod(
    weights, modes[, p_lam + seq_len(p_p), drop = FALSE]))
  u_mean <- as.numeric(crossprod(
    weights, modes[, p_lam + p_p + seq_len(n_mesh), drop = FALSE]))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lam))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(beta_lambda_mean) <- nm_lam
  names(beta_p_mean)      <- nm_p

  out <- c(fit, list(
    range_grid       = range_grid,
    sigma_grid       = sigma_grid,
    r_grid           = r_grid_use,
    weights          = weights,
    range_mean       = range_mean,
    range_sd         = range_sd,
    sigma_mean       = sigma_mean,
    sigma_sd         = sigma_sd,
    mixture          = mixture,
    r_mean           = disp$r_mean,
    r_sd             = disp$r_sd,
    beta_lambda_mean = beta_lambda_mean,
    beta_p_mean      = beta_p_mean,
    vcov             = .nmix_grid_vcov(fit$cov_blocks, modes, weights,
                                       p_lam, p_p, c(nm_lam, nm_p)),
    u_mean           = u_mean,
    n_sites          = n_sites,
    n_obs            = n_obs,
    n_spatial        = n_mesh,
    prior_type       = "spde",
    call             = match.call()
  ))
  if (any(out$boundary_max > 1e-4, na.rm = TRUE)) {
    warning(sprintf(
      "Max posterior weight on N = K_max is %.2e at one or more grid points; raise K_max.",
      max(out$boundary_max, na.rm = TRUE)), call. = FALSE)
  }
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

# log|Q| for a sparse symmetric PD precision. Matrix::determinant on a
# dsCMatrix returns the log-determinant of the matrix itself (logarithm = TRUE
# by default), computed through a sparse Cholesky.
.spde_logdet_Q <- function(Q) {
  Qs <- methods::as(Matrix::forceSymmetric(Q), "CsparseMatrix")
  as.numeric(Matrix::determinant(Qs, logarithm = TRUE)$modulus)
}
