#' Spatial Royle (2004) N-mixture model via nested Laplace
#'
#' @description
#' Nested-Laplace fit of the spatial N-mixture model
#' \deqn{N_i \sim \mathrm{Poisson}(\lambda_i), \qquad
#'       y_{ij} | N_i \sim \mathrm{Binomial}(N_i, p_{ij}),}
#' with abundance linear predictor
#' \eqn{\log \lambda_i = X_\lambda^{(i)} \beta_\lambda + z_{u(i)}} where
#' \eqn{z \sim \mathrm{ICAR}(\tau)} is an intrinsic conditional autoregressive
#' field on the user-supplied adjacency graph, and detection linear predictor
#' \eqn{\mathrm{logit}\, p_{ij} = X_p^{(ij)} \beta_p}.
#'
#' The hyperparameter \eqn{\tau} (ICAR precision) is integrated out by an
#' outer grid: at each \eqn{\tau_k} the inner Newton finds the joint mode
#' of \eqn{(\beta_\lambda, \beta_p, z)} and the Laplace log-marginal
#' \eqn{\log p(y \mid \tau_k)} is accumulated. Posterior weights over
#' \eqn{\tau} normalise the grid.
#'
#' Inner Newton uses the marginal observed Fisher information matrix for
#' curvature (with the Var\eqn{[N \mid y_i]} rank-1 correction encoding
#' cross-arm coupling) and falls back to the complete-data Fisher block when
#' the observed-info matrix is not PSD. A small diagonal ridge keeps the
#' Cholesky stable through the (intercept, constant-\eqn{z}) structural null
#' direction; \eqn{z} is centered to sum zero after every step.
#'
#' @param y Integer vector of observed counts (long form, one entry per visit).
#' @param site_idx Integer vector, 1-based site index for each visit, same
#'   length as `y`.
#' @param map_site_to_unit Integer vector of length `n_sites`, 1-based spatial
#'   unit index for each site. Sites can share spatial units; the data
#'   contribution to `z[u]` aggregates across all sites that map to `u`.
#' @param X_lambda Numeric matrix `[n_sites x p_lambda]` of abundance covariates.
#' @param X_p Numeric matrix `[n_obs x p_p]` of detection covariates.
#' @param adj_row_ptr,adj_col_idx,n_neighbors CSR adjacency for the ICAR graph
#'   on the `n_spatial` units. `adj_row_ptr` has length `n_spatial + 1`,
#'   `adj_col_idx` lists 0-based neighbours, and `n_neighbors[s]` is the row
#'   degree of unit s.
#' @param n_spatial Number of spatial units.
#' @param tau_grid Optional numeric vector of \eqn{\tau} grid points. Defaults
#'   to `exp(seq(log(0.3), log(30), length.out = 9))`.
#' @param mixture Abundance mixing distribution: `"P"` (Poisson, default) or
#'   `"NB"` (negative binomial). Under `"NB"` the NB size \eqn{r} is integrated
#'   as an additional outer grid dimension (alongside \eqn{\tau}); the posterior
#'   `r_mean` / `r_sd` are reported from the grid weights.
#' @param r_grid Optional numeric vector of NB size grid points (NB only).
#'   Defaults to `exp(seq(log(0.5), log(40), length.out = 6))`. Ignored under
#'   Poisson.
#' @param beta_lambda_init Optional warm start; default `c(log(mean(y)+0.1), 0, ...)`.
#' @param beta_p_init Optional warm start; default `rep(0, p_p)`.
#' @param z_init Optional warm start for the spatial field; default zeros.
#' @param K_max Truncation for the per-site marginal sum over N. Defaults to
#'   `max(y) + 100`. Returned `boundary_max` flags any grid point whose worst
#'   site puts non-trivial mass on `K_max` -- raise `K_max` if it exceeds 1e-4.
#' @param max_iter,tol Inner Newton iteration budget and gradient-norm tol.
#' @param verbose Print per-iteration and per-grid-point progress.
#'
#' @return A list of class `nmix_spatial_fit`:
#'   * `tau_grid` -- input grid
#'   * `log_marginal` -- log marginal at each tau (up to a tau-independent constant)
#'   * `weights` -- normalised grid weights (sum to 1)
#'   * `tau_mean`, `tau_sd` -- posterior moments of tau
#'   * `modes` -- `[n_grid x (p_lambda + p_p + n_spatial)]` matrix of inner modes
#'   * `beta_lambda_mean`, `beta_p_mean` -- weighted-mean coefficient estimates
#'   * `z_mean` -- weighted-mean spatial field
#'   * `n_iter`, `converged`, `grad_norm`, `log_lik`, `boundary_max` -- per-grid diagnostics
#'   * `p_lambda`, `p_p`, `n_spatial`, `K_max` -- echoed dimensions
#'   * `call` -- matched call
#'
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Besag, J., York, J., Mollie, A. (1991). Bayesian image restoration with two
#'   applications in spatial statistics. *Ann. Inst. Statist. Math.* 43, 1-20.
#' Rue, H., Martino, S., Chopin, N. (2009). Approximate Bayesian inference for
#'   latent Gaussian models by using integrated nested Laplace approximations.
#'   *JRSS-B* 71, 319-392.
#'
nmix_laplace_icar <- function(y,
                                    site_idx,
                                    map_site_to_unit,
                                    X_lambda,
                                    X_p,
                                    adj_row_ptr,
                                    adj_col_idx,
                                    n_neighbors,
                                    n_spatial,
                                    tau_grid = NULL,
                                    mixture = c("P", "NB"),
                                    r_grid = NULL,
                                    beta_lambda_init = NULL,
                                    beta_p_init = NULL,
                                    z_init = NULL,
                                    K_max = NULL,
                                    max_iter = 100L,
                                    tol = 1e-6,
                                    verbose = FALSE) {
  mixture <- match.arg(mixture)
  y                <- as.integer(y)
  site_idx         <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .count_spatial_prep(y, site_idx, X_lambda, X_p, mixture,
                            beta_lambda_init, beta_p_init, K_max, r_grid,
                            map_site_to_unit = map_site_to_unit,
                            n_spatial = n_spatial, adj_row_ptr = adj_row_ptr,
                            n_neighbors = n_neighbors,
                            latent = list(z_init = z_init),
                            n_latent = c(n_spatial = n_spatial))
  if (is.null(tau_grid)) tau_grid <- .count_spatial_default_grid("tau_icar")
  tau_grid <- .count_spatial_check_grid(tau_grid, "tau_grid", 0, Inf)

  fit <- .cpp_nmix_progress(cpp_nested_laplace_nmix_icar,
    y                  = y,
    site_idx           = site_idx,
    map_site_to_unit_R = map_site_to_unit,
    X_lambda_R         = X_lambda,
    X_p_R              = X_p,
    adj_row_ptr        = as.integer(adj_row_ptr),
    adj_col_idx        = as.integer(adj_col_idx),
    n_neighbors        = as.integer(n_neighbors),
    n_spatial          = as.integer(n_spatial),
    tau_grid           = as.numeric(tau_grid),
    r_grid             = as.numeric(pp$r_grid),
    beta_lambda_init   = as.numeric(pp$beta_lambda_init),
    beta_p_init        = as.numeric(pp$beta_p_init),
    z_init             = if (is.null(z_init)) NULL else as.numeric(z_init),
    K_max              = pp$K_max,
    max_iter           = as.integer(max_iter),
    tol                = as.numeric(tol),
    verbose            = isTRUE(verbose)
  )

  out <- c(fit, .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                           X_lambda, X_p, mixture),
           list(n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "icar",
                call = match.call()))
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' Proper CAR Royle (2004) N-mixture model via nested Laplace
#'
#' @description
#' Nested-Laplace fit of the spatial N-mixture model with a proper conditional
#' autoregressive prior on the abundance-arm spatial offset:
#' \deqn{N_i \sim \mathrm{Poisson}(\lambda_i), \qquad
#'       y_{ij} | N_i \sim \mathrm{Binomial}(N_i, p_{ij}),}
#' \eqn{\log \lambda_i = X_\lambda^{(i)} \beta_\lambda + z_{u(i)}}, where
#' \eqn{z \mid \tau, \rho \sim \mathrm{N}(0, [\tau (D - \rho W)]^{-1})}.
#' Both hyperparameters are integrated over an outer 2D grid; the inner
#' Newton step shares the kernel with [nmix_laplace_icar()] (ICAR is
#' the \eqn{\rho = 1} limit).
#'
#' Unlike the ICAR fit, no sum-to-zero centering is applied -- \eqn{Q(\rho)}
#' is full rank for \eqn{\rho < 1}. The per-rho \eqn{\log |Q(\rho)|} is
#' precomputed once via a dense Cholesky on the \eqn{n_{\mathrm{spatial}}
#' \times n_{\mathrm{spatial}}} precision matrix.
#'
#' @inheritParams nmix_laplace_icar
#' @param tau_grid Optional numeric vector of \eqn{\tau} grid points
#'   (defaults to `exp(seq(log(0.3), log(30), length.out = 7L))`).
#' @param mixture Abundance mixing distribution `"P"` (default) or `"NB"`; under
#'   `"NB"` the NB size \eqn{r} is integrated as an additional outer grid axis.
#' @param r_grid Optional NB size grid (NB only); see [nmix_laplace_icar()].
#' @param rho_grid Optional numeric vector of \eqn{\rho} grid points in
#'   the valid eigenvalue interval. Defaults to a 5-point grid in
#'   \eqn{(0, 1)} -- callers that want eigenvalue-derived bounds should
#'   compute them via tulpa::spatial_car_proper(adjacency)`$`rho_bounds and
#'   pass an explicit grid in that interval.
#'
#' @return A list of class `nmix_spatial_fit`:
#'   * `theta_grid` -- `[n_grid x 2]` matrix of (tau, rho) per grid point
#'   * `tau_grid`, `rho_grid` -- input axes
#'   * `log_det_Q_rho` -- precomputed log determinants per rho
#'   * `log_marginal`, `weights`, `tau_mean`, `tau_sd`, `rho_mean`, `rho_sd`
#'     -- as in [nmix_laplace_icar()] but with the rho marginal added
#'   * other diagnostic fields and named coefficient means as in ICAR
#'
#' @references
#' Cressie, N. (1993). Statistics for Spatial Data. Wiley.
#' Rue, H., Held, L. (2005). Gaussian Markov Random Fields. CRC.
#'
nmix_laplace_car_proper <- function(y,
                                          site_idx,
                                          map_site_to_unit,
                                          X_lambda,
                                          X_p,
                                          adj_row_ptr,
                                          adj_col_idx,
                                          n_neighbors,
                                          n_spatial,
                                          tau_grid = NULL,
                                          rho_grid = NULL,
                                          mixture = c("P", "NB"),
                                          r_grid = NULL,
                                          beta_lambda_init = NULL,
                                          beta_p_init = NULL,
                                          z_init = NULL,
                                          K_max = NULL,
                                          max_iter = 100L,
                                          tol = 1e-6,
                                          verbose = FALSE) {
  mixture <- match.arg(mixture)
  y                <- as.integer(y)
  site_idx         <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .count_spatial_prep(y, site_idx, X_lambda, X_p, mixture,
                            beta_lambda_init, beta_p_init, K_max, r_grid,
                            map_site_to_unit = map_site_to_unit,
                            n_spatial = n_spatial, adj_row_ptr = adj_row_ptr,
                            n_neighbors = n_neighbors,
                            latent = list(z_init = z_init),
                            n_latent = c(n_spatial = n_spatial))
  if (is.null(tau_grid)) tau_grid <- .count_spatial_default_grid("tau_car")
  if (is.null(rho_grid)) rho_grid <- .count_spatial_default_grid("rho_car")
  tau_grid <- .count_spatial_check_grid(tau_grid, "tau_grid", 0, Inf)
  rho_grid <- .count_spatial_check_grid(
    rho_grid, "rho_grid", 0, 1, open = TRUE,
    hint = "Pass explicit eigenvalue bounds via spatial_car_proper().")

  fit <- .cpp_nmix_progress(cpp_nested_laplace_nmix_car_proper,
    y                  = y,
    site_idx           = site_idx,
    map_site_to_unit_R = map_site_to_unit,
    X_lambda_R         = X_lambda,
    X_p_R              = X_p,
    adj_row_ptr        = as.integer(adj_row_ptr),
    adj_col_idx        = as.integer(adj_col_idx),
    n_neighbors        = as.integer(n_neighbors),
    n_spatial          = as.integer(n_spatial),
    tau_grid           = as.numeric(tau_grid),
    rho_grid           = as.numeric(rho_grid),
    r_grid             = as.numeric(pp$r_grid),
    beta_lambda_init   = as.numeric(pp$beta_lambda_init),
    beta_p_init        = as.numeric(pp$beta_p_init),
    z_init             = if (is.null(z_init)) NULL else as.numeric(z_init),
    K_max              = pp$K_max,
    max_iter           = as.integer(max_iter),
    tol                = as.numeric(tol),
    verbose            = isTRUE(verbose)
  )

  # rho summaries are proper-CAR-specific; the rest of the grid summarisation is
  # shared with the ICAR path.
  common <- .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                       X_lambda, X_p, mixture)
  rho <- .tobs_weighted_moment(common$weights, fit$theta_grid[, "rho"])

  out <- c(fit, common,
           list(rho_mean = unname(rho["mean"]), rho_sd = unname(rho["sd"]),
                n_sites = pp$n_sites, n_obs = pp$n_obs,
                prior_type = "car_proper", call = match.call()))
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' BYM2 Royle (2004) N-mixture model via nested Laplace
#'
#' @description
#' Nested-Laplace fit of the spatial N-mixture model with a BYM2 (Riebler et al.
#' 2016) prior on the abundance-arm spatial offset. The offset decomposes as
#' \deqn{\phi_u = \sigma \left(\sqrt{\rho / s} \, v_u + \sqrt{1 - \rho} \, w_u\right),}
#' with \eqn{v \sim \mathrm{ICAR}} (unscaled, sum-to-zero) and
#' \eqn{w \sim \mathrm{N}(0, I)} iid. \eqn{s} is the Riebler scaling factor,
#' the geometric mean of the marginal variances \eqn{\mathrm{diag}(Q^{+})} of
#' the intrinsic ICAR field; \eqn{\sigma} is then the joint marginal standard
#' deviation of \eqn{\phi}, and \eqn{\rho \in [0, 1]} is the spatial fraction
#' of variance.
#'
#' The inner Newton works in the joint state
#' \eqn{x = (\beta_\lambda, \beta_p, v, w)} (dimension
#' \eqn{p_\lambda + p_p + 2 n_{\mathrm{spatial}}}). At the converged mode the
#' Laplace log-marginal is accumulated; the outer 2D grid integrates over
#' \eqn{(\sigma, \rho)}.
#'
#' @inheritParams nmix_laplace_icar
#' @param mixture Abundance mixing distribution `"P"` (default) or `"NB"`; under
#'   `"NB"` the NB size \eqn{r} is integrated as an additional outer grid axis.
#' @param r_grid Optional NB size grid (NB only); see [nmix_laplace_icar()].
#' @param sigma_grid Numeric vector of \eqn{\sigma} (joint sd) grid points.
#'   Defaults to `exp(seq(log(0.2), log(3), length.out = 5L))`.
#' @param rho_grid Numeric vector of spatial-fraction grid points in
#'   \eqn{[0, 1]}. Defaults to `c(0.05, 0.3, 0.5, 0.7, 0.95)`.
#' @param scale_factor Optional scalar Riebler scaling factor \eqn{s}, the
#'   geometric mean of the marginal variances \eqn{\mathrm{diag}(Q^{+})} of the
#'   intrinsic ICAR field. If `NULL`, it is computed from the adjacency via
#'   dense eigendecomposition.
#' @param v_init,w_init Optional warm starts (each length `n_spatial`) for the
#'   two BYM2 latent components -- the structured (spatial) and unstructured
#'   (iid) effects. `NULL` (default) starts both at zero.
#'
#' @return A list of class `nmix_spatial_fit`:
#'   * `theta_grid` -- `[n_grid x 2]` matrix of (sigma, rho) per grid point
#'   * `sigma_grid`, `rho_grid` -- input axes
#'   * `scale_factor` -- the Riebler scaling factor
#'   * `log_marginal`, `weights`, `sigma_mean`, `sigma_sd`, `rho_mean`, `rho_sd`
#'     -- posterior summaries of the joint hyperparameters
#'   * `modes` -- `[n_grid x (p_lambda + p_p + 2 n_spatial)]` per-grid modes;
#'     columns `(p_lambda + p_p + 1) .. (p_lambda + p_p + n_spatial)` are `v`,
#'     the next `n_spatial` columns are `w`
#'   * `beta_lambda_mean`, `beta_p_mean` -- weighted-mean coefficients
#'   * `v_mean`, `w_mean` -- weighted-mean ICAR and iid components
#'   * `phi_mean` -- weighted-mean total offset \eqn{\phi}
#'   * other diagnostic fields as in [nmix_laplace_icar()]
#'
#' @references
#' Riebler, A., Sorbye, S. H., Simpson, D., Rue, H. (2016). An intuitive
#'   Bayesian spatial model for disease mapping that accounts for scaling.
#'   *Statistical Methods in Medical Research* 25, 1145-1165.
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A.,
#'   DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing
#'   the BYM2 model in Stan. *Spatial and Spatio-temporal Epidemiology* 31.
#'
nmix_laplace_bym2 <- function(y,
                                    site_idx,
                                    map_site_to_unit,
                                    X_lambda,
                                    X_p,
                                    adj_row_ptr,
                                    adj_col_idx,
                                    n_neighbors,
                                    n_spatial,
                                    sigma_grid = NULL,
                                    rho_grid = NULL,
                                    mixture = c("P", "NB"),
                                    r_grid = NULL,
                                    scale_factor = NULL,
                                    beta_lambda_init = NULL,
                                    beta_p_init = NULL,
                                    v_init = NULL,
                                    w_init = NULL,
                                    K_max = NULL,
                                    max_iter = 100L,
                                    tol = 1e-6,
                                    verbose = FALSE) {
  mixture <- match.arg(mixture)
  y                <- as.integer(y)
  site_idx         <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .count_spatial_prep(y, site_idx, X_lambda, X_p, mixture,
                            beta_lambda_init, beta_p_init, K_max, r_grid,
                            map_site_to_unit = map_site_to_unit,
                            n_spatial = n_spatial, adj_row_ptr = adj_row_ptr,
                            n_neighbors = n_neighbors,
                            latent = list(v_init = v_init, w_init = w_init),
                            n_latent = c(n_spatial = n_spatial))
  if (is.null(sigma_grid)) sigma_grid <- .count_spatial_default_grid("sigma_bym2")
  if (is.null(rho_grid))   rho_grid   <- .count_spatial_default_grid("rho_bym2")
  sigma_grid <- .count_spatial_check_grid(sigma_grid, "sigma_grid", 0, Inf)
  rho_grid   <- .count_spatial_check_grid(rho_grid, "rho_grid", 0, 1, open = FALSE)
  scale_factor <- .bym2_resolve_scale(scale_factor, adj_row_ptr, adj_col_idx,
                                      n_spatial)

  fit <- .cpp_nmix_progress(cpp_nested_laplace_nmix_bym2,
    y                  = y,
    site_idx           = site_idx,
    map_site_to_unit_R = map_site_to_unit,
    X_lambda_R         = X_lambda,
    X_p_R              = X_p,
    adj_row_ptr        = as.integer(adj_row_ptr),
    adj_col_idx        = as.integer(adj_col_idx),
    n_neighbors        = as.integer(n_neighbors),
    n_spatial          = as.integer(n_spatial),
    sigma_grid         = as.numeric(sigma_grid),
    rho_grid           = as.numeric(rho_grid),
    r_grid             = as.numeric(pp$r_grid),
    scale_factor       = as.numeric(scale_factor),
    beta_lambda_init   = as.numeric(pp$beta_lambda_init),
    beta_p_init        = as.numeric(pp$beta_p_init),
    v_init             = if (is.null(v_init)) NULL else as.numeric(v_init),
    w_init             = if (is.null(w_init)) NULL else as.numeric(w_init),
    K_max              = pp$K_max,
    max_iter           = as.integer(max_iter),
    tol                = as.numeric(tol),
    verbose            = isTRUE(verbose)
  )

  out <- c(fit, .count_spatial_pack_bym2_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                                X_lambda, X_p, mixture, scale_factor),
           list(n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "bym2",
                call = match.call()))
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' @exportS3Method print nmix_spatial_fit
print.nmix_spatial_fit <- function(x, ...) {
  ptype <- x$prior_type %||% "icar"
  label <- switch(ptype,
                  icar       = "ICAR",
                  car_proper = "CAR(rho)",
                  bym2       = "BYM2",
                  toupper(ptype))
  mix <- x$mixture %||% "P"
  cat(sprintf("tulpa spatial N-mixture (%s, mixture = %s) nested-Laplace fit\n",
              label, mix))
  cat(sprintf("  n_sites = %d   n_obs = %d   n_spatial = %d   K_max = %d\n",
              x$n_sites, x$n_obs, x$n_spatial, x$K_max))
  if (identical(mix, "NB") && is.finite(x$r_mean %||% NA_real_)) {
    cat(sprintf("  NB size r: mean = %.3g   sd = %.3g (grid-integrated)\n",
                x$r_mean, x$r_sd))
  }
  if (ptype == "icar") {
    cat(sprintf("  tau grid = [%.3g, %.3g] over %d points\n",
                min(x$tau_grid), max(x$tau_grid), x$n_grid))
    cat(sprintf("  tau_mean = %.3g   tau_sd = %.3g\n",
                x$tau_mean, x$tau_sd))
  } else if (ptype == "bym2") {
    cat(sprintf("  sigma grid = [%.3g, %.3g] x %d   rho grid = [%.3g, %.3g] x %d\n",
                min(x$sigma_grid), max(x$sigma_grid), length(x$sigma_grid),
                min(x$rho_grid), max(x$rho_grid), length(x$rho_grid)))
    cat(sprintf("  scale_factor = %.4g\n", x$scale_factor))
    cat(sprintf("  sigma_mean = %.3g   sigma_sd = %.3g\n",
                x$sigma_mean, x$sigma_sd))
    cat(sprintf("  rho_mean   = %.3g   rho_sd   = %.3g\n",
                x$rho_mean, x$rho_sd))
  } else {
    cat(sprintf("  tau grid = [%.3g, %.3g] x %d   rho grid = [%.3g, %.3g] x %d\n",
                min(x$tau_grid), max(x$tau_grid), length(x$tau_grid),
                min(x$rho_grid), max(x$rho_grid), length(x$rho_grid)))
    cat(sprintf("  tau_mean = %.3g   tau_sd = %.3g\n",
                x$tau_mean, x$tau_sd))
    if (!is.null(x$rho_mean)) {
      cat(sprintf("  rho_mean = %.3g   rho_sd = %.3g\n",
                  x$rho_mean, x$rho_sd))
    }
  }
  cat("\nabundance (log lambda):\n")
  print(x$beta_lambda_mean)
  cat("\ndetection (logit p):\n")
  print(x$beta_p_mean)
  invisible(x)
}

# Resolve the NB-dispersion grid for the spatial fits. Poisson -> a single
# r = Inf node (the kernel takes its Poisson branch). NB -> a user grid or a
# default log-spaced grid of NB sizes spanning strong overdispersion to nearly
# Poisson. The outer nested-Laplace integration treats r like tau/rho/sigma.
.nmix_resolve_r_grid <- function(mixture, r_grid) {
  if (identical(mixture, "P")) return(Inf)
  if (is.null(r_grid)) {
    return(exp(seq(log(0.5), log(40), length.out = 6L)))
  }
  if (!is.numeric(r_grid) || length(r_grid) < 1L || any(r_grid <= 0) || anyNA(r_grid)) {
    stop("`r_grid` must be a vector of positive NB sizes.", call. = FALSE)
  }
  r_grid
}

# ---------------------------------------------------------------------------
# Shared preamble for the count-marginal spatial fitters
#
# The areal N-mixture wrappers (icar / car_proper / bym2), the SPDE mesh wrapper
# and the removal wrappers fit the same count marginal over the same designs,
# and so validate the same arguments. Everything below is the part they share;
# what differs is passed in, so a default or a check is a property of the axis
# or of the family rather than of which wrapper the caller happened to reach.
# ---------------------------------------------------------------------------

# Default outer-grid axes for the areal count fitters.
#
# The tau axis carries 9 nodes on the 1D ICAR grid and 7 on the 2D proper-CAR
# grid: the proper-CAR grid is the product of both axes, so 7 x 5 = 35 inner
# Newton solves already costs four times the ICAR grid's 9. Both count families
# read these, so the axis is what decides, not the family.
.count_spatial_default_grid <- function(axis) {
  switch(
    axis,
    tau_icar   = exp(seq(log(0.3), log(30), length.out = 9L)),
    tau_car    = exp(seq(log(0.3), log(30), length.out = 7L)),
    rho_car    = c(0.1, 0.3, 0.5, 0.75, 0.95),
    sigma_bym2 = exp(seq(log(0.2), log(3), length.out = 5L)),
    rho_bym2   = c(0.05, 0.3, 0.5, 0.7, 0.95),
    stop(sprintf("Unknown count-spatial grid axis '%s'.", axis), call. = FALSE)
  )
}

# Range check for a resolved outer-grid axis. `open` = the proper-CAR rho
# interval, which excludes both ends (rho = 1 is the intrinsic limit, where
# Q(rho) drops rank and the dense Cholesky behind log|Q(rho)| fails); the closed
# form is the BYM2 mixing weight, where both ends are meaningful (all-iid and
# all-structured). Applied at every door, so a grid valid for one areal count
# family is valid for the other.
.count_spatial_check_grid <- function(x, name, lo, hi, open = TRUE, hint = "") {
  if (!is.numeric(x) || !length(x) || anyNA(x)) {
    stop(sprintf("`%s` must be a non-empty numeric vector.", name), call. = FALSE)
  }
  bad <- if (open) any(x <= lo) || any(x >= hi) else any(x < lo) || any(x > hi)
  if (bad) {
    stop(sprintf("`%s` values must lie %s.%s", name,
                 if (open) sprintf("strictly in (%g, %g)", lo, hi)
                 else sprintf("in [%g, %g]", lo, hi),
                 if (nzchar(hint)) paste0(" ", hint) else ""),
         call. = FALSE)
  }
  as.numeric(x)
}

# Latent-N truncation floor: the largest count at any single site, which the
# N-mixture marginal must be able to represent.
.count_K_floor_max_y <- function(y, site_idx, n_sites) {
  list(value = max(y), label = "max(y)")
}

# Removal's floor: the depleting passes at a site sum, so the latent N clears
# that site's TOTAL removal, not its largest single pass.
.count_K_floor_site_total <- function(y, site_idx, n_sites) {
  tot <- tapply(as.integer(y),
                factor(as.integer(site_idx), levels = seq_len(n_sites)), sum)
  tot[is.na(tot)] <- 0L
  list(value = max(as.integer(tot)),
       label = "the largest per-site removal total")
}

# Resolve the latent-N truncation. An explicit K_max below the family's floor
# cannot represent the data; the default clears the floor by the same 100
# headroom at every door.
.count_spatial_K_max <- function(K_max, floor) {
  if (is.null(K_max)) return(as.integer(floor$value + 100L))
  K_max <- as.integer(K_max)
  if (K_max < floor$value) {
    stop(sprintf("K_max must be >= %s.", floor$label), call. = FALSE)
  }
  K_max
}

# Resolve and validate everything the count-marginal spatial fitters share:
# the two designs against the declared dimensions, the site-to-unit map and the
# CSR graph against `n_spatial`, the coefficient warm starts, the NB-size grid
# and the latent-N truncation, plus any latent-field warm start.
#
# `map_site_to_unit` / `adj_row_ptr` are NULL on the SPDE door (a mesh carries a
# projector, not a graph), and the checks they gate are skipped there.
# `K_floor` names the family's truncation rule; `latent` is the named list of
# field warm starts, each checked against `n_latent` -- whose NAME is what the
# message calls that length (`n_spatial` areal, `n_mesh` on the mesh).
.count_spatial_prep <- function(y, site_idx, X_lambda, X_p, mixture,
                                beta_lambda_init, beta_p_init, K_max, r_grid,
                                K_floor = .count_K_floor_max_y,
                                map_site_to_unit = NULL, n_spatial = NULL,
                                adj_row_ptr = NULL, n_neighbors = NULL,
                                latent = list(), n_latent = NULL) {
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  n_sites <- nrow(X_lambda); n_obs <- nrow(X_p)
  p_lam   <- ncol(X_lambda); p_p   <- ncol(X_p)
  if (length(y) != n_obs) stop("length(y) must equal nrow(X_p).", call. = FALSE)
  if (length(site_idx) != n_obs) {
    stop("length(site_idx) must equal nrow(X_p).", call. = FALSE)
  }
  if (!is.null(map_site_to_unit)) {
    if (length(map_site_to_unit) != n_sites) {
      stop("length(map_site_to_unit) must equal nrow(X_lambda).", call. = FALSE)
    }
    if (any(map_site_to_unit < 1L) || any(map_site_to_unit > n_spatial)) {
      stop("map_site_to_unit values must lie in [1, n_spatial].", call. = FALSE)
    }
  }
  if (!is.null(adj_row_ptr)) {
    if (length(adj_row_ptr) != n_spatial + 1L) {
      stop("length(adj_row_ptr) must equal n_spatial + 1.", call. = FALSE)
    }
    if (length(n_neighbors) != n_spatial) {
      stop("length(n_neighbors) must equal n_spatial.", call. = FALSE)
    }
  }
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
  for (nm in names(latent)) {
    v <- latent[[nm]]
    if (!is.null(v) && length(v) != n_latent) {
      stop(sprintf("length(%s) must equal %s.", nm, names(n_latent)),
           call. = FALSE)
    }
  }
  list(n_sites = n_sites, n_obs = n_obs, p_lam = p_lam, p_p = p_p,
       beta_lambda_init = beta_lambda_init, beta_p_init = beta_p_init,
       K_max = .count_spatial_K_max(K_max, K_floor(y, site_idx, n_sites)),
       r_grid = .nmix_resolve_r_grid(mixture, r_grid))
}

# The boundary diagnostic every count-marginal spatial fitter reports: posterior
# mass sitting on N = K_max at any outer-grid point means the truncation clipped
# the latent-N sum there, so the fit is answering a question the window could
# not hold. One wording, one threshold, one call.
.count_spatial_warn_boundary <- function(out) {
  if (any(out$boundary_max > 1e-4, na.rm = TRUE)) {
    warning(sprintf(paste0("Max posterior weight on N = K_max is %.2e at one or ",
                           "more grid points; raise K_max."),
                    max(out$boundary_max, na.rm = TRUE)), call. = FALSE)
  }
  invisible(out)
}

# Shared grid-summarisation for the single-field areal count fits (ICAR /
# proper-CAR, abundance-arm field on one spatial unit per site). The N-mixture
# and removal nested-Laplace wrappers differ only in the C++ entry point and the
# prior-specific hyperparameter summaries (tau here; rho added by the proper-CAR
# wrapper); the weights, weighted-mean coefficients, weighted field, grid-
# integrated coefficient covariance, and NB-size summary are identical. Single
# source of truth. Returns the common block the wrappers splice into `out`.
.count_spatial_pack_common <- function(fit, p_lam, p_p, n_spatial,
                                       X_lambda, X_p, mixture) {
  weights <- .tobs_grid_weights(fit, "tau_grid / data")
  tau  <- .tobs_weighted_moment(weights, fit$theta_grid[, "tau"])
  disp <- .nmix_dispersion_summary(mixture, fit$theta_grid, weights)

  modes <- fit$modes
  beta_lambda_mean <- as.numeric(crossprod(weights, modes[, seq_len(p_lam), drop = FALSE]))
  beta_p_mean      <- as.numeric(crossprod(
    weights, modes[, p_lam + seq_len(p_p), drop = FALSE]))
  z_mean <- as.numeric(crossprod(
    weights, modes[, p_lam + p_p + seq_len(n_spatial), drop = FALSE]))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lam))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(beta_lambda_mean) <- nm_lam
  names(beta_p_mean)      <- nm_p

  list(
    weights          = weights,
    tau_mean         = unname(tau["mean"]),
    tau_sd           = unname(tau["sd"]),
    mixture          = mixture,
    r_mean           = disp$r_mean,
    r_sd             = disp$r_sd,
    beta_lambda_mean = beta_lambda_mean,
    beta_p_mean      = beta_p_mean,
    vcov             = .nmix_grid_vcov(fit$cov_blocks, modes, weights,
                                       p_lam, p_p, c(nm_lam, nm_p)),
    z_mean           = z_mean
  )
}

# Shared grid-summarisation for the BYM2 areal count fits (two-field v / w on the
# abundance arm). Like .count_spatial_pack_common but with the BYM2-specific
# sigma / rho summaries, the per-component (v, w) weighted means, and the total
# offset phi = sigma (sqrt(rho/scale) v + sqrt(1-rho) w) integrated over the grid.
# Single source of truth for the N-mixture and removal BYM2 wrappers.
.count_spatial_pack_bym2_common <- function(fit, p_lam, p_p, n_spatial,
                                            X_lambda, X_p, mixture, scale_factor) {
  weights <- .tobs_grid_weights(fit)
  sigma_vec <- fit$theta_grid[, "sigma"]; rho_vec <- fit$theta_grid[, "rho"]
  sigma <- .tobs_weighted_moment(weights, sigma_vec)
  rho   <- .tobs_weighted_moment(weights, rho_vec)
  disp  <- .nmix_dispersion_summary(mixture, fit$theta_grid, weights)

  modes <- fit$modes
  beta_lambda_mean <- as.numeric(crossprod(weights, modes[, seq_len(p_lam), drop = FALSE]))
  beta_p_mean      <- as.numeric(crossprod(
    weights, modes[, p_lam + seq_len(p_p), drop = FALSE]))
  v_idx <- p_lam + p_p + seq_len(n_spatial)
  w_idx <- p_lam + p_p + n_spatial + seq_len(n_spatial)
  v_mean <- as.numeric(crossprod(weights, modes[, v_idx, drop = FALSE]))
  w_mean <- as.numeric(crossprod(weights, modes[, w_idx, drop = FALSE]))

  phi_grid <- matrix(0, nrow = nrow(modes), ncol = n_spatial)
  for (k in seq_len(nrow(modes))) {
    sg <- sigma_vec[k]; rg <- rho_vec[k]
    phi_grid[k, ] <- (sg * sqrt(rg / scale_factor)) * modes[k, v_idx] +
                     (sg * sqrt(1 - rg)) * modes[k, w_idx]
  }
  phi_mean <- as.numeric(crossprod(weights, phi_grid))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lam))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(beta_lambda_mean) <- nm_lam; names(beta_p_mean) <- nm_p

  list(
    weights = weights,
    sigma_mean = unname(sigma["mean"]), sigma_sd = unname(sigma["sd"]),
    rho_mean = unname(rho["mean"]), rho_sd = unname(rho["sd"]),
    mixture = mixture,
    r_mean = disp$r_mean, r_sd = disp$r_sd,
    beta_lambda_mean = beta_lambda_mean, beta_p_mean = beta_p_mean,
    vcov = .nmix_grid_vcov(fit$cov_blocks, modes, weights, p_lam, p_p,
                           c(nm_lam, nm_p)),
    v_mean = v_mean, w_mean = w_mean, phi_mean = phi_mean)
}

# Posterior mean / sd of the NB size r from the outer grid weights (the "r"
# column of theta_grid). Poisson -> NA. Weighted moments mirror tau/rho/sigma.
.nmix_dispersion_summary <- function(mixture, theta_grid, weights) {
  if (identical(mixture, "P") || !("r" %in% colnames(theta_grid))) {
    return(list(r_mean = NA_real_, r_sd = NA_real_))
  }
  r <- .tobs_weighted_moment(weights, theta_grid[, "r"])
  list(r_mean = unname(r["mean"]), r_sd = unname(r["sd"]))
}

# Grid-integrated covariance of beta = (beta_lambda, beta_p) over the outer
# hyperparameter grid. `cov_blocks[[k]]` is the within-grid Laplace covariance
# at the k-th mode (the beta block of the joint H^{-1} the C++ kernel returns);
# .tobs_grid_vcov() adds the between-grid mode spread. A grid with no usable
# weight yields an all-NA covariance.
.nmix_grid_vcov <- function(cov_blocks, modes, weights, p_lam, p_p, nm) {
  p_beta     <- p_lam + p_p
  beta_modes <- modes[, seq_len(p_beta), drop = FALSE]
  ok <- is.finite(weights) & weights > 0
  V  <- if (any(ok)) {
    w <- weights; w[!ok] <- 0
    .tobs_grid_vcov(beta_modes, w, cov_blocks)
  } else {
    matrix(NA_real_, p_beta, p_beta)
  }
  dimnames(V) <- list(nm, nm)
  V
}
