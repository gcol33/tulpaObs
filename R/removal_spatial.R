# =============================================================================
# removal_spatial.R
# - areal-spatial removal-sampling abundance
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

# The removal latent-N truncation on its own, for the BFGS driver -- which packs
# its own state and so does not go through the wrapper preamble below.
.removal_spatial_K_max <- function(y, site_idx, n_sites, K_max) {
  .count_spatial_K_max(K_max, .count_K_floor_site_total(y, site_idx, n_sites))
}

# Input validation for the removal areal wrappers. The count-marginal preamble
# is shared with the N-mixture spatial wrappers; removal supplies only its own
# truncation rule, since depleting passes at a site sum and the latent N must
# clear that site's TOTAL removal.
.removal_spatial_prep <- function(y, site_idx, map_site_to_unit, X_lambda, X_p,
                                  adj_row_ptr, n_neighbors, n_spatial, mixture,
                                  beta_lambda_init, beta_p_init, K_max, r_grid,
                                  latent = list()) {
  .count_spatial_prep(y, site_idx, X_lambda, X_p, mixture,
                      beta_lambda_init, beta_p_init, K_max, r_grid,
                      K_floor = .count_K_floor_site_total,
                      map_site_to_unit = map_site_to_unit,
                      n_spatial = n_spatial, adj_row_ptr = adj_row_ptr,
                      n_neighbors = n_neighbors,
                      latent = latent, n_latent = c(n_spatial = n_spatial))
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
                              beta_lambda_init, beta_p_init, K_max, r_grid,
                              latent = list(z_init = z_init))
  if (is.null(tau_grid)) tau_grid <- .count_spatial_default_grid("tau_icar")
  tau_grid <- .count_spatial_check_grid(tau_grid, "tau_grid", 0, Inf)

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
  .count_spatial_warn_boundary(out)
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
                              beta_lambda_init, beta_p_init, K_max, r_grid,
                              latent = list(z_init = z_init))
  if (is.null(tau_grid)) tau_grid <- .count_spatial_default_grid("tau_car")
  if (is.null(rho_grid)) rho_grid <- .count_spatial_default_grid("rho_car")
  tau_grid <- .count_spatial_check_grid(tau_grid, "tau_grid", 0, Inf)
  rho_grid <- .count_spatial_check_grid(
    rho_grid, "rho_grid", 0, 1, open = TRUE,
    hint = "Pass explicit eigenvalue bounds via spatial_car_proper().")

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

  common <- .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                       X_lambda, X_p, mixture)
  rho <- .tobs_weighted_moment(common$weights, fit$theta_grid[, "rho"])
  out <- c(fit, common,
           list(rho_mean = unname(rho["mean"]), rho_sd = unname(rho["sd"]),
                n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = "car_proper",
                call = match.call()))
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' Areal BYM2 removal-sampling abundance via nested Laplace (internal)
#' @keywords internal
#' @noRd
removal_laplace_bym2 <- function(y, site_idx, map_site_to_unit, X_lambda, X_p,
                                 adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
                                 sigma_grid = NULL, rho_grid = NULL,
                                 scale_factor = NULL, mixture = c("P", "NB"),
                                 r_grid = NULL, beta_lambda_init = NULL,
                                 beta_p_init = NULL, v_init = NULL, w_init = NULL,
                                 K_max = NULL, max_iter = 100L, tol = 1e-6,
                                 verbose = FALSE) {
  mixture <- match.arg(mixture)
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  pp <- .removal_spatial_prep(y, site_idx, map_site_to_unit, X_lambda, X_p,
                              adj_row_ptr, n_neighbors, n_spatial, mixture,
                              beta_lambda_init, beta_p_init, K_max, r_grid,
                              latent = list(v_init = v_init, w_init = w_init))
  if (is.null(sigma_grid)) sigma_grid <- .count_spatial_default_grid("sigma_bym2")
  if (is.null(rho_grid))   rho_grid   <- .count_spatial_default_grid("rho_bym2")
  sigma_grid <- .count_spatial_check_grid(sigma_grid, "sigma_grid", 0, Inf)
  rho_grid   <- .count_spatial_check_grid(rho_grid, "rho_grid", 0, 1, open = FALSE)
  scale_factor <- .bym2_resolve_scale(scale_factor, adj_row_ptr, adj_col_idx,
                                      n_spatial)

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
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

# Dispatch an areal-spatial removal fit from .tobs_fit_model. icar() / car_proper()
# / bym2() on one spatial unit per site (the field enters log lambda); the
# continuous spde() field is not yet wired for removal, and a weighted (SVC) field
# is rejected. The packed fit reuses build_nmix_fit (the removal coefficient
# layout is the shared (lambda, p) block).
.tobs_fit_removal_spatial <- function(model, spatial, temporal = NULL, svc = NULL,
                                      mixture = "P", K_max = NULL,
                                      max_iter = 100L, tol = 1e-6, verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "removal abundance spatial")
  # Detection-arm field: a field in the `detection=` formula carries shared =
  # c(abundance, detection) = c(FALSE, TRUE). It loads on the per-pass detection
  # logit (a spatially-varying capture probability) instead of the abundance arm.
  # The C++ nested-Laplace removal kernels carry the field on the abundance arm
  # only, so a detection-arm field routes through the shared areal-BFGS driver
  # (the removal marginal exposes the per-pass detection gradient
  # cpp_removal_total_log_lik$grad_eta_p, summed to a per-site field gradient
  # inside the fitter).
  det_arm <- isTRUE(spatial$shared[2L]) && !isTRUE(spatial$shared[1L])
  if (det_arm)
    return(.tobs_fit_removal_spatial_bfgs(model, spatial, temporal, det_arm = TRUE,
                                          mixture = mixture, K_max = K_max,
                                          max_iter = max_iter, tol = tol,
                                          verbose = verbose))
  if (spatial$type %in% c("spde", "gp", "multiscale_gp")) {
    stop(sprintf(paste0("removal() areal spatial supports icar() / car_proper() / ",
                        "bym2() under method = \"nested_laplace\"; the '%s' field ",
                        "is not yet wired for removal."), spatial$type),
         call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "car_proper", "bym2")) {
    stop(sprintf("removal() areal spatial supports icar() / car_proper() / bym2(); got '%s'.",
                 spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for removal."),
                 spatial$n_units, n_sites), call. = FALSE)
  }
  # Spatial + temporal, or spatial + a continuous varying coefficient: route
  # through the shared areal-BFGS driver (a second latent block alongside the
  # field), reusing the removal marginal's analytic per-site / per-pass gradient
  # (cpp_removal_total_log_lik). The C++ count-spatial driver carries only the
  # field, so the extra blocks live on the BFGS path. Spatial-only fits keep the
  # C++ driver.
  if (!is.null(temporal) || !is.null(svc)) {
    return(.tobs_fit_removal_spatial_bfgs(model, spatial, temporal, svc = svc,
                                          mixture = mixture, K_max = K_max,
                                          max_iter = max_iter, tol = tol,
                                          verbose = verbose))
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
                          .bym2_scale(spatial$graph)))))
  fit <- build_nmix_fit(raw, model, spatial = spatial)
  # The field is part of log lambda, so its per-site contribution is recorded for
  # the post-fit readers (the areal-BFGS branches above record theirs through
  # .tobs_attach_field_results()). Without it fitted() / residuals() / predict()
  # and WAIC / LOO / DIC / CPO score log lambda = X beta. A detection-arm term
  # took the BFGS branch above, so the field is the abundance arm's here.
  fit$spatial_field_arm <- "abundance"
  .tobs_set_field_eta_offset(fit, 1L,
                             .tobs_spatial_field_offset(fit, spatial, model))
}

# Areal field + temporal block on the removal abundance arm via the shared areal-
# BFGS driver. The removal marginal exposes an analytic per-site abundance
# gradient + per-pass detection gradient + the NB dispersion score
# (cpp_removal_total_log_lik), so the driver's eval() contract is met directly:
# it runs BFGS over (beta_lambda, beta_p[, log_r], field_sp, field_tmp)
# + both block priors, integrated over the product of the two blocks'
# hyperparameter grids. The field + temporal block both load onto eta_lambda.
.tobs_fit_removal_spatial_bfgs <- function(model, spatial, temporal,
                                           det_arm = FALSE, svc = NULL,
                                           mixture = "P", K_max = NULL,
                                           max_iter = 200L, tol = 1e-8,
                                           verbose = TRUE) {
  n_sites <- model$n_sites
  map <- seq_len(n_sites)
  temporal_only <- is.null(spatial) && !is.null(temporal)
  X_lam <- model$X_processes[[1]]
  .tobs_check_svc_arm(svc, det_arm, "removal")
  # A detection-arm areal field (det_arm, #114) is a single areal block that loads
  # on eta_p; temporal does not compose with the detection-arm field on this path.
  field <- .tobs_build_field_spec(spatial, temporal, "removal", n_sites, map,
                                  svc = svc, X_svc = X_lam)

  X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lam); p_p <- ncol(X_p)
  y_long   <- as.integer(model$y_long)
  site_idx <- as.integer(model$site_idx)
  is_nb <- mixture %in% c("negbin", "NB")
  K_max <- .removal_spatial_K_max(y_long, site_idx, n_sites, K_max)

  i_lam <- seq_len(p_lam); i_p <- p_lam + seq_len(p_p)
  i_logr <- if (is_nb) p_lam + p_p + 1L else NA_integer_
  n_fixed <- p_lam + p_p + if (is_nb) 1L else 0L

  eval <- function(theta_fix, offset) {
    # `offset` is per-SITE (length n_sites). On the abundance arm it enters eta_lam
    # directly; on the detection arm (det_arm) it is a spatially-varying capture
    # logit shared across a site's passes, so it is expanded to per-pass via
    # site_idx and the per-pass detection gradient is summed back to per-site.
    off_lam <- if (det_arm) numeric(n_sites) else offset
    eta_lam <- as.numeric(X_lam %*% theta_fix[i_lam]) + off_lam
    eta_p   <- as.numeric(X_p %*% theta_fix[i_p]) + (if (det_arm) offset[site_idx] else 0)
    rr <- if (is_nb) exp(theta_fix[i_logr]) else Inf
    out <- cpp_removal_total_log_lik(y_long, site_idx, eta_p, eta_lam, K_max, rr)
    g <- numeric(n_fixed)
    g[i_lam] <- as.numeric(crossprod(X_lam, out$grad_eta_lambda))
    g[i_p]   <- as.numeric(crossprod(X_p,   out$grad_eta_p))
    if (is_nb) g[i_logr] <- out$grad_theta
    grad_eta <- if (det_arm) {           # sum per-pass det gradient into per-site
      gs  <- numeric(n_sites)
      agg <- rowsum(out$grad_eta_p, site_idx, reorder = FALSE)
      gs[as.integer(rownames(agg))] <- agg[, 1L]
      gs
    } else out$grad_eta_lambda
    list(log_lik = out$log_lik, grad_fixed = g, grad_eta = grad_eta)
  }

  warm <- tryCatch(
    removal_laplace(y = y_long, site_idx = site_idx, X_lambda = X_lam, X_p = X_p,
                    mixture = if (is_nb) "NB" else "P", K_max = K_max,
                    max_iter = as.integer(max_iter), tol = as.numeric(tol),
                    verbose = FALSE),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm)) {
    th <- c(unname(warm$beta_lambda), unname(warm$beta_p))
    if (is_nb) th <- c(th, log(if (is.finite(warm$r %||% NA_real_)) warm$r else 2))
    th
  } else {
    th <- c(log(mean(y_long) + 0.1), rep(0, p_lam - 1L), rep(0, p_p))
    if (is_nb) th <- c(th, log(2))
    th
  }
  if (length(theta0_fix) != n_fixed) theta0_fix <- numeric(n_fixed)

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol,
                              label = "removal-spatial", integration = "grid")

  means <- res$beta_mean; V <- res$vcov
  raw <- list(
    mixture = if (is_nb) "NB" else "P",
    beta_lambda = means[i_lam], beta_p = means[i_p],
    log_r = if (is_nb) means[i_logr] else NA_real_,
    r = if (is_nb) exp(means[i_logr]) else NA_real_,
    r_mean = if (is_nb) exp(means[i_logr]) else NA_real_,
    r_sd = if (is_nb && is.finite(V[i_logr, i_logr]))
             exp(means[i_logr]) * sqrt(V[i_logr, i_logr]) else NA_real_,
    vcov = V[c(i_lam, i_p), c(i_lam, i_p), drop = FALSE],
    log_lik = res$log_lik, converged = TRUE, K_max = K_max)
  fit <- build_nmix_fit(raw, model, spatial = spatial)
  # The field loads on the abundance (log lambda) arm by default, or on the
  # per-pass capture logit when the term sits in the detection formula (#114).
  # removal does not surface the outer Pareto-k diagnostic (pareto_k = FALSE).
  .tobs_attach_field_results(fit, res, det_arm, temporal, temporal_only, "abundance",
                             pareto_k = FALSE, svc = svc,
                             has_spatial = !is.null(spatial),
                             X_svc = X_lam, family = "removal")
}

# Areal-spatial removal sampling via NUTS: a FIXED-HYPER non-centered PROPER-CAR
# field on the abundance arm of the removal marginal. The field precision (tau,
# rho) is fixed at the nested-Laplace posterior mean and the whitened raw ~ N(0,
# I) (z = Linv %*% raw) is sampled jointly with the coefficients via the shared
# count-marginal NUTS field block (cpp_removal_nuts over marginal_count_nuts.h /
# nuts_field_block.h, byte-identical to abun's). Same car_proper-only restriction
# (the intrinsic ICAR's flat field-mean needs a sum-to-zero reparameterisation).
# Poisson or NB.
.tobs_fit_removal_nuts_spatial <- function(model, spatial = NULL, temporal = NULL,
                                           mixture = "poisson",
                                           K_max = NULL, sigma.beta = NULL,
                                           sigma.logr = NULL, n.iter = NULL,
                                           n.warmup = NULL, n.chains = NULL, n.thin = NULL,
                                           n.threads = NULL,
                                           max.treedepth = NULL, adapt.delta = NULL,
                                           seed = NULL, verbose = FALSE) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  # This fitter carries the FIXED-HYPER non-centered field on the abundance arm
  # under NUTS from EITHER an areal term (icar/car_proper/bym2; #72/#113) OR a
  # temporal() term (ar1/rw1/rw2/iid; #114). Both reduce to a whitened raw ~ N(0,I)
  # sampled through the shared count-marginal field block; only the loading L, the
  # per-site field map, and the warm-start source differ, so the sampling tail is
  # shared (no duplicate sampler).
  temporal_only <- is.null(spatial) && !is.null(temporal)
  if (!temporal_only) {
    .tobs_reject_weighted_spatial(spatial, "removal NUTS abundance spatial")
    .tobs_reject_det_arm_spatial(spatial, "removal() NUTS", "abundance",
                                 "a spatially-varying capture logit",
                                 "is wired under method = \"nested_laplace\".")
    if (!spatial$type %in% c("icar", "car_proper", "bym2"))
      stop(sprintf(paste0("removal() NUTS + areal spatial supports icar() / ",
                          "car_proper() / bym2() on the abundance arm; got '%s'. ",
                          ""), spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (!temporal_only && spatial$n_units != n_sites)
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for removal NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  is_nb <- mixture %in% c("negbin", "NB")
  mix_code <- if (is_nb) "NB" else "P"
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  y_long <- as.integer(model$y_long); site_idx <- as.integer(model$site_idx)
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  K_max <- .removal_spatial_K_max(y_long, site_idx, n_sites, K_max)

  if (temporal_only) {
    # Warm coefficients + fixed temporal hyper (tau, rho) from the temporal-only
    # nested-Laplace areal-BFGS fit; the whitened temporal-field loading is the
    # eigen-decomposition of tau Q_temporal (rank-deficient for rw1/rw2, so
    # n_raw < n_t, the same sum-to-zero geometry as an intrinsic areal field).
    ti <- as.integer(temporal$time_idx)
    if (length(ti) != n_sites)
      stop(sprintf(paste0("temporal term has %d time indices but the model has %d ",
                          "sites; one time index per site is required for removal ",
                          "NUTS + temporal."), length(ti), n_sites), call. = FALSE)
    n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
           else max(ti, na.rm = TRUE)
    nlfit <- .tobs_fit_removal_spatial_bfgs(model, spatial = NULL, temporal = temporal,
                                            mixture = mix_code, K_max = K_max,
                                            max_iter = 100L, tol = 1e-6, verbose = FALSE)
    hyper <- nlfit$temporal_hyper
    hv <- function(k) suppressWarnings(as.numeric(hyper[[k]]))
    fl <- .tobs_nuts_temporal_loading(temporal$type, n_t, tau = hv("tau"), rho = hv("rho"))
    field_map <- ti
    n_field_units <- n_t
    beta0 <- as.numeric(nlfit$means[seq_len(p_lam + p_p)])
    if (is_nb) beta0 <- c(beta0, log(if (is.finite(nlfit$r %||% NA_real_)) nlfit$r else 2))
  } else {
    adj <- as.matrix(spatial$graph)
    csr <- if (!is.null(spatial$adj_row_ptr))
      list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
           n_neighbors = spatial$n_neighbors) else adjacency_to_csr(spatial$graph)

    # Fixed hyper + warm coefficients from the matching nested-Laplace areal fit,
    # and the whitened-field loading L (square for car_proper, sum-to-zero for
    # icar / bym2/#113).
    common <- list(
      y = y_long, site_idx = site_idx, map_site_to_unit = seq_len(n_sites),
      X_lambda = X_lambda, X_p = X_p, adj_row_ptr = csr$row_ptr,
      adj_col_idx = csr$col_idx, n_neighbors = csr$n_neighbors, n_spatial = n_sites,
      mixture = mix_code, K_max = K_max, max_iter = 100L, tol = 1e-6, verbose = FALSE)
    nl <- switch(spatial$type,
      icar       = do.call(removal_laplace_icar, common),
      car_proper = do.call(removal_laplace_car_proper, common),
      bym2       = do.call(removal_laplace_bym2, c(common,
                     list(scale_factor = spatial$scale_factor %||%
                            .bym2_scale(spatial$graph)))))
    fl <- .tobs_nuts_field_loading(adj, spatial$type, n_sites,
                                   tau = nl$tau_mean, rho = nl$rho_mean,
                                   sigma = nl$sigma_mean,
                                   scale_factor = spatial$scale_factor)
    field_map <- seq_len(n_sites)
    n_field_units <- n_sites
    beta0 <- c(nl$beta_lambda_mean, nl$beta_p_mean)
    if (is_nb && is.finite(nl$r_mean %||% NA_real_)) beta0 <- c(beta0, log(nl$r_mean))
    else if (is_nb) beta0 <- c(beta0, log(2))
  }
  field_load <- fl$field_load; n_raw <- fl$n_raw

  spec <- list(y = y_long, site_idx = site_idx, X_lambda = X_lambda, X_p = X_p,
               n_sites = n_sites, K_max = K_max, is_nb = is_nb,
               n_field_units = n_field_units, field_map = field_map,
               field_load = field_load)

  n_base <- p_lam + p_p + if (is_nb) 1L else 0L
  theta0 <- c(beta0, numeric(n_raw))
  inv_metric <- c(rep(0.1, n_base), rep(1, n_raw))

  run_chain <- function(ch)
    cpp_removal_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                     sigma_logr = sigma.logr, inv_metric = inv_metric,
                     n_iter = as.integer(n.iter + n.warmup),
                     n_warmup = as.integer(n.warmup),
                     max_treedepth = as.integer(max.treedepth),
                     adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                     verbose = isTRUE(verbose))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           if (is_nb) "log_r", paste0("raw_", seq_len(n_raw)))
  run <- .tobs_nuts_field_draws(run_chain, n.chains, nms, n_base, n_raw,
                                field_load, n.thin = n.thin,
                                n.threads = n.threads)
  par <- run$par; cov <- run$cov

  lay  <- .tobs_abun_nuts_layout(p_lam, p_p, is_nb)
  marg <- .tobs_removal_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$p],
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik
  raw_fit <- list(mixture = mix_code,
                  beta_lambda = unname(par[lay$lambda]), beta_p = unname(par[lay$p]),
                  log_r = if (is_nb) unname(par[lay$log_r]) else NA_real_,
                  r = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
                  vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max)
  fit <- build_nmix_fit(raw_fit, model, spatial = if (temporal_only) NULL else spatial)
  .tobs_nuts_field_attach(
    fit, run, ll_mean, n.chains,
    prior_type = if (temporal_only) temporal$type else spatial$type, fl = fl,
    field_map = field_map,
    temporal = if (temporal_only) temporal else NULL)
}
