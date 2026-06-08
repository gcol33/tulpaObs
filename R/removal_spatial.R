# =============================================================================
# removal_spatial.R - areal-spatial removal-sampling abundance (gcol33/tulpaObs#51)
#
# An ICAR or proper-CAR field on the abundance arm of the removal marginal, fit
# by nested Laplace (outer grid over tau[, rho][, NB size r], inner Newton over
# (beta_lambda, beta_p, z)). The removal per-site marginal shares the count-
# marginal moment interface with the Royle N-mixture, so the C++ inner Newton /
# Laplace driver and the R grid-summarisation (.count_spatial_pack_common) are
# reused verbatim; only the C++ entry point (cpp_nested_laplace_removal_*) and the
# K_max default (the per-site removal TOTAL, since depletion sums passes) differ.
#
#   removal_laplace_icar()         R wrapper over cpp_nested_laplace_removal_icar
#   removal_laplace_car_proper()   R wrapper over cpp_nested_laplace_removal_car_proper
#   .tobs_fit_removal_spatial()    dispatch from .tobs_fit_model (icar / car_proper)
# =============================================================================

# Per-site removal total + buffer: the marginal truncation must clear each site's
# summed removals across passes (the latent N is >= the total removed).
.removal_spatial_K_max <- function(y, site_idx, n_sites, K_max) {
  if (!is.null(K_max)) {
    K_max <- as.integer(K_max)
    site_tot <- tapply(as.integer(y), factor(as.integer(site_idx),
                                             levels = seq_len(n_sites)), sum)
    site_tot[is.na(site_tot)] <- 0L
    if (K_max < max(as.integer(site_tot)))
      stop("K_max must be >= the largest per-site removal total.", call. = FALSE)
    return(K_max)
  }
  site_tot <- tapply(as.integer(y), factor(as.integer(site_idx),
                                           levels = seq_len(n_sites)), sum)
  site_tot[is.na(site_tot)] <- 0L
  as.integer(max(as.integer(site_tot)) + 100L)
}

# Shared input validation for the removal areal wrappers (mirrors the N-mixture
# spatial wrappers; the abundance design is per-site, the detection design and y
# are per-pass long form). Returns the resolved inits / K_max / r_grid.
.removal_spatial_prep <- function(y, site_idx, map_site_to_unit, X_lambda, X_p,
                                  adj_row_ptr, n_neighbors, n_spatial, mixture,
                                  beta_lambda_init, beta_p_init, K_max, r_grid) {
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  n_sites <- nrow(X_lambda); n_obs <- nrow(X_p)
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  if (length(y) != n_obs) stop("length(y) must equal nrow(X_p).", call. = FALSE)
  if (length(site_idx) != n_obs) stop("length(site_idx) must equal nrow(X_p).", call. = FALSE)
  if (length(map_site_to_unit) != n_sites)
    stop("length(map_site_to_unit) must equal nrow(X_lambda).", call. = FALSE)
  if (any(map_site_to_unit < 1L) || any(map_site_to_unit > n_spatial))
    stop("map_site_to_unit values must lie in [1, n_spatial].", call. = FALSE)
  if (length(adj_row_ptr) != n_spatial + 1L)
    stop("length(adj_row_ptr) must equal n_spatial + 1.", call. = FALSE)
  if (length(n_neighbors) != n_spatial)
    stop("length(n_neighbors) must equal n_spatial.", call. = FALSE)
  if (is.null(beta_lambda_init))
    beta_lambda_init <- c(log(mean(y) + 0.1), rep(0, p_lam - 1L))
  if (is.null(beta_p_init)) beta_p_init <- rep(0, p_p)
  if (length(beta_lambda_init) != p_lam)
    stop("length(beta_lambda_init) must equal ncol(X_lambda).", call. = FALSE)
  if (length(beta_p_init) != p_p)
    stop("length(beta_p_init) must equal ncol(X_p).", call. = FALSE)
  list(n_sites = n_sites, n_obs = n_obs, p_lam = p_lam, p_p = p_p,
       beta_lambda_init = beta_lambda_init, beta_p_init = beta_p_init,
       K_max = .removal_spatial_K_max(y, site_idx, n_sites, K_max),
       r_grid = .nmix_resolve_r_grid(mixture, r_grid))
}

#' Areal ICAR removal-sampling abundance via nested Laplace (internal)
#' @keywords internal
#' @noRd
removal_laplace_icar <- function(y, site_idx, map_site_to_unit, X_lambda, X_p,
                                 adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
                                 tau_grid = NULL, mixture = c("P", "NB"),
                                 r_grid = NULL, beta_lambda_init = NULL,
                                 beta_p_init = NULL, z_init = NULL, K_max = NULL,
                                 max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  mixture <- match.arg(mixture)
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .removal_spatial_prep(y, site_idx, map_site_to_unit, X_lambda, X_p,
                              adj_row_ptr, n_neighbors, n_spatial, mixture,
                              beta_lambda_init, beta_p_init, K_max, r_grid)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))

  fit <- .cpp_nmix_progress(cpp_nested_laplace_removal_icar,
    y = y, site_idx = site_idx, map_site_to_unit_R = map_site_to_unit,
    X_lambda_R = X_lambda, X_p_R = X_p,
    adj_row_ptr = as.integer(adj_row_ptr), adj_col_idx = as.integer(adj_col_idx),
    n_neighbors = as.integer(n_neighbors), n_spatial = as.integer(n_spatial),
    tau_grid = as.numeric(tau_grid), r_grid = as.numeric(pp$r_grid),
    beta_lambda_init = as.numeric(pp$beta_lambda_init),
    beta_p_init = as.numeric(pp$beta_p_init),
    z_init = if (is.null(z_init)) NULL else as.numeric(z_init),
    K_max = pp$K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose))

  out <- c(fit, .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                           X_lambda, X_p, mixture),
           list(n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "icar",
                call = match.call()))
  if (any(out$boundary_max > 1e-4, na.rm = TRUE))
    warning(sprintf("Max posterior weight on N = K_max is %.2e at one or more grid points; raise K_max.",
                    max(out$boundary_max, na.rm = TRUE)), call. = FALSE)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' Areal proper-CAR removal-sampling abundance via nested Laplace (internal)
#' @keywords internal
#' @noRd
removal_laplace_car_proper <- function(y, site_idx, map_site_to_unit, X_lambda,
                                       X_p, adj_row_ptr, adj_col_idx, n_neighbors,
                                       n_spatial, tau_grid = NULL, rho_grid = NULL,
                                       mixture = c("P", "NB"), r_grid = NULL,
                                       beta_lambda_init = NULL, beta_p_init = NULL,
                                       z_init = NULL, K_max = NULL,
                                       max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  mixture <- match.arg(mixture)
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .removal_spatial_prep(y, site_idx, map_site_to_unit, X_lambda, X_p,
                              adj_row_ptr, n_neighbors, n_spatial, mixture,
                              beta_lambda_init, beta_p_init, K_max, r_grid)
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))
  if (is.null(rho_grid)) rho_grid <- seq(0.1, 0.95, length.out = 6L)

  fit <- .cpp_nmix_progress(cpp_nested_laplace_removal_car_proper,
    y = y, site_idx = site_idx, map_site_to_unit_R = map_site_to_unit,
    X_lambda_R = X_lambda, X_p_R = X_p,
    adj_row_ptr = as.integer(adj_row_ptr), adj_col_idx = as.integer(adj_col_idx),
    n_neighbors = as.integer(n_neighbors), n_spatial = as.integer(n_spatial),
    tau_grid = as.numeric(tau_grid), rho_grid = as.numeric(rho_grid),
    r_grid = as.numeric(pp$r_grid),
    beta_lambda_init = as.numeric(pp$beta_lambda_init),
    beta_p_init = as.numeric(pp$beta_p_init),
    z_init = if (is.null(z_init)) NULL else as.numeric(z_init),
    K_max = pp$K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose))

  rho_vec <- fit$theta_grid[, "rho"]
  common  <- .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                        X_lambda, X_p, mixture)
  w <- common$weights
  rho_mean <- sum(w * rho_vec, na.rm = TRUE)
  rho_sd   <- sqrt(max(0, sum(w * rho_vec^2, na.rm = TRUE) - rho_mean^2))
  out <- c(fit, common,
           list(rho_mean = rho_mean, rho_sd = rho_sd,
                n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "car_proper",
                call = match.call()))
  if (any(out$boundary_max > 1e-4, na.rm = TRUE))
    warning(sprintf("Max posterior weight on N = K_max is %.2e at one or more grid points; raise K_max.",
                    max(out$boundary_max, na.rm = TRUE)), call. = FALSE)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' Areal BYM2 removal-sampling abundance via nested Laplace (internal)
#' @keywords internal
#' @noRd
removal_laplace_bym2 <- function(y, site_idx, map_site_to_unit, X_lambda, X_p,
                                 adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
                                 sigma_grid = NULL, rho_grid = NULL,
                                 scale_factor = 1, mixture = c("P", "NB"),
                                 r_grid = NULL, beta_lambda_init = NULL,
                                 beta_p_init = NULL, v_init = NULL, w_init = NULL,
                                 K_max = NULL, max_iter = 100L, tol = 1e-6,
                                 verbose = FALSE) {
  mixture <- match.arg(mixture)
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .removal_spatial_prep(y, site_idx, map_site_to_unit, X_lambda, X_p,
                              adj_row_ptr, n_neighbors, n_spatial, mixture,
                              beta_lambda_init, beta_p_init, K_max, r_grid)
  if (is.null(sigma_grid)) sigma_grid <- exp(seq(log(0.2), log(3), length.out = 5L))
  if (is.null(rho_grid))   rho_grid <- c(0.05, 0.3, 0.5, 0.7, 0.95)
  if (scale_factor <= 0) stop("scale_factor must be positive.", call. = FALSE)

  fit <- .cpp_nmix_progress(cpp_nested_laplace_removal_bym2,
    y = y, site_idx = site_idx, map_site_to_unit_R = map_site_to_unit,
    X_lambda_R = X_lambda, X_p_R = X_p,
    adj_row_ptr = as.integer(adj_row_ptr), adj_col_idx = as.integer(adj_col_idx),
    n_neighbors = as.integer(n_neighbors), n_spatial = as.integer(n_spatial),
    sigma_grid = as.numeric(sigma_grid), rho_grid = as.numeric(rho_grid),
    r_grid = as.numeric(pp$r_grid), scale_factor = as.numeric(scale_factor),
    beta_lambda_init = as.numeric(pp$beta_lambda_init),
    beta_p_init = as.numeric(pp$beta_p_init),
    v_init = if (is.null(v_init)) NULL else as.numeric(v_init),
    w_init = if (is.null(w_init)) NULL else as.numeric(w_init),
    K_max = pp$K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose))

  out <- c(fit, .count_spatial_pack_bym2_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                                X_lambda, X_p, mixture, scale_factor),
           list(n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "bym2",
                call = match.call()))
  if (any(out$boundary_max > 1e-4, na.rm = TRUE))
    warning(sprintf("Max posterior weight on N = K_max is %.2e at one or more grid points; raise K_max.",
                    max(out$boundary_max, na.rm = TRUE)), call. = FALSE)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

# Dispatch an areal-spatial removal fit from .tobs_fit_model. icar() / car_proper()
# / bym2() on one spatial unit per site (the field enters log lambda); the
# continuous spde() field is not yet wired for removal, and a weighted (SVC) field
# is rejected. The packed fit reuses build_nmix_fit (the removal coefficient
# layout is the shared (lambda, p) block).
.tobs_fit_removal_spatial <- function(model, spatial, mixture = "P", K_max = NULL,
                                      max_iter = 100L, tol = 1e-6, verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "removal abundance spatial")
  if (spatial$type %in% c("spde", "gp", "multiscale_gp")) {
    stop(sprintf(paste0("removal() areal spatial supports icar() / car_proper() / ",
                        "bym2() under method = \"nested_laplace\"; the '%s' field ",
                        "is not yet wired for removal. (tulpaObs#51)"), spatial$type),
         call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "car_proper", "bym2")) {
    stop(sprintf("removal() areal spatial supports icar() / car_proper() / bym2(); got '%s'.",
                 spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for removal.",
                 spatial$n_units, n_sites), call. = FALSE)
  }
  csr <- if (!is.null(spatial$adj_row_ptr)) {
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors)
  } else {
    adjacency_to_csr(spatial$graph)
  }
  mix_code <- switch(mixture, poisson = "P", negbin = "NB", P = "P", NB = "NB",
                     stop(sprintf("Unknown mixture '%s'.", mixture), call. = FALSE))
  common <- list(
    y = model$y_long, site_idx = model$site_idx,
    map_site_to_unit = seq_len(n_sites),
    X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, n_spatial = spatial$n_units,
    mixture = mix_code, K_max = K_max, max_iter = as.integer(max_iter),
    tol = as.numeric(tol), verbose = isTRUE(verbose))
  raw <- switch(spatial$type,
    icar       = do.call(removal_laplace_icar, common),
    car_proper = do.call(removal_laplace_car_proper, common),
    bym2       = do.call(removal_laplace_bym2, c(common,
                   list(scale_factor = spatial$scale_factor %||%
                          compute_bym2_scale(spatial$graph)))))
  build_nmix_fit(raw, model, spatial = spatial)
}
