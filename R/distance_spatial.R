# =============================================================================
# distance_spatial.R - areal-spatial binned distance-sampling abundance (#51)
#
# An ICAR / proper-CAR / BYM2 field on the abundance arm (log lambda) of the
# binned distance marginal, via the shared areal-BFGS nested-Laplace driver
# (R/areal_bfgs.R). The distance kernel exposes the analytic per-site gradient
# (cpp_distance_site_sweep over compute_distance_site); the driver runs BFGS over
# (beta_lambda, beta_sigma[, log_r], field) + the field prior and forms the
# Laplace marginal from an FD-Hessian of that gradient at the mode -- which is the
# distance marginal observed information (the documented diag(info_lam, info_sig)
# - var_N v v' structure, recovered numerically). The field loads onto eta_lambda
# exactly like the abundance intercept (one unit per site). Half-normal key only;
# Poisson or NB (the NB size log_r is jointly estimated, as in the non-spatial
# distance fit).
#
#   .tobs_fit_distance_spatial()   dispatch from .tobs_fit_model
# =============================================================================

.tobs_fit_distance_spatial <- function(model, spatial, temporal = NULL,
                                       svc = NULL, mixture = "poisson",
                                       K_max = NULL, headroom = NULL,
                                       max_iter = 200L, tol = 1e-6,
                                       verbose = TRUE, integration = "grid") {
  temporal_only <- is.null(spatial) && !is.null(temporal)
  # Detection-arm field (gcol33/tulpaObs#114): a field in the `detection=` formula
  # carries shared = c(state, detection) = c(FALSE, TRUE). It loads on the per-site
  # detection scale eta_sigma (a spatially-varying detection scale) instead of the
  # abundance arm; the distance kernel exposes the per-site sigma gradient
  # (sw$grad_sig), so the areal-BFGS driver needs only the arm-routed offset +
  # gradient. Under the hazard-rate key the detection scale sigma still carries the
  # field while the log-shape eta_b stays a single global coordinate (the eval
  # threads both -- the field on eta_sigma, eta_b as a fixed parameter).
  det_arm <- !is.null(spatial) && isTRUE(spatial$shared[2L]) &&
             !isTRUE(spatial$shared[1L])
  if (!is.null(spatial))
    .tobs_reject_weighted_spatial(spatial, "distance abundance spatial")
  hazard <- identical(model$key, "hazard")
  key_code <- .dist_key_code(model$key)
  map <- seq_len(model$n_sites)
  X_lam <- model$X_processes[[1]]
  .tobs_check_svc_arm(svc, det_arm, "distance")
  field <- .tobs_build_field_spec(spatial, temporal, "distance", model$n_sites, map,
                                  svc = svc, X_svc = X_lam)

  X_sig <- model$X_processes[[2]]
  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig)
  y <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  cutpoints <- as.numeric(model$cutpoints)
  transect_code <- .dist_transect_code(model$transect)
  quad_order <- as.integer(model$quad_order %||% 64L)
  # Built ONCE outside `eval()` -- the areal-BFGS driver calls `eval()` every
  # gradient evaluation, so rebuilding the quadrature there would pay the
  # Newton-Raphson root-find on every BFGS step (gcol33/tulpaObs#165).
  quad_xptr <- cpp_distance_build_quad(cutpoints, transect_code, quad_order)
  site_tot <- if (length(y)) rowSums(y) else integer(0)
  trunc    <- .dist_truncation(K_max, site_tot)
  K_max    <- trunc$K_max
  headroom <- if (is.null(headroom)) trunc$headroom else as.integer(headroom)
  is_nb <- mixture %in% c("negbin", "NB")
  # NB's heavier tail is not sized for the per-site window; keep the shared
  # ceiling under NB (mirrors nmix_laplace(), R/nmix_laplace.R).
  if (is_nb) headroom <- -1L

  # Fixed coefficient layout: (beta_lambda, beta_sigma[, eta_b][, log_r]). Under
  # the hazard-rate key the scalar log-shape eta_b is a global coordinate (#79),
  # placed after the detection-scale block and before the NB dispersion.
  i_lam <- seq_len(p_lam); i_sig <- p_lam + seq_len(p_sig)
  off <- p_lam + p_sig
  i_b   <- if (hazard) { off <- off + 1L; off } else NA_integer_
  i_logr <- if (is_nb) off + 1L else NA_integer_
  n_fixed <- off + if (is_nb) 1L else 0L

  eval <- function(theta_fix, offset) {
    eta_lam <- as.numeric(X_lam %*% theta_fix[i_lam]) + (if (det_arm) 0 else offset)
    eta_sig <- as.numeric(X_sig %*% theta_fix[i_sig]) + (if (det_arm) offset else 0)
    eta_b <- if (hazard) theta_fix[i_b] else 0.0
    rr <- if (is_nb) exp(theta_fix[i_logr]) else Inf
    sw <- cpp_distance_site_sweep(y, eta_lam, eta_sig, quad_xptr,
                                  K_max, nb = is_nb, r = rr,
                                  key = key_code, eta_b = eta_b,
                                  headroom = headroom)
    g <- numeric(n_fixed)
    g[i_lam] <- as.numeric(crossprod(X_lam, sw$grad_lam))
    g[i_sig] <- as.numeric(crossprod(X_sig, sw$grad_sig))
    if (hazard) g[i_b] <- sw$grad_b
    if (is_nb) g[i_logr] <- sw$grad_logr
    list(log_lik = sum(sw$log_lik), grad_fixed = g,
         grad_eta = if (det_arm) sw$grad_sig else sw$grad_lam)
  }

  # Warm start from the non-spatial distance Laplace fit (carries the hazard
  # log-shape when the key is hazard-rate).
  warm <- tryCatch(
    .tobs_distance_re_warm(model, mixture = if (is_nb) "NB" else "P",
                           K_max = K_max, max_iter = max_iter, tol = tol),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm)) {
    th <- c(warm$beta_lambda, warm$beta_p)
    if (hazard) th <- c(th, if (is.finite(warm$eta_b %||% NA_real_)) warm$eta_b else 0)
    if (is_nb) th <- c(th, log(if (is.finite(warm$r %||% NA_real_)) warm$r else 2))
    th
  } else {
    th <- c(log(max(mean(rowSums(y)), 0.5) + 0.5), rep(0, p_lam - 1L),
            log(stats::median(cutpoints[-1])), rep(0, p_sig - 1L))
    if (hazard) th <- c(th, 0)
    if (is_nb) th <- c(th, log(2))
    th
  }
  if (length(theta0_fix) != n_fixed) theta0_fix <- numeric(n_fixed)

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol, label = "distance-spatial",
                              integration = integration)

  # Guard the per-site truncation at the BFGS mode (gcol33/tulpaObs#168): the
  # same score-gap comparison distance_laplace() runs, but sandwiched at the
  # field-adjusted eta the areal fit actually converged to (`res$eta_offset`,
  # the SAME per-site offset eval()'s own `offset` argument carries).
  if (headroom >= 0L) {
    gap <- tryCatch({
      means0 <- res$beta_mean
      offset_mode <- res$eta_offset %||% numeric(model$n_sites)
      eta_lam_mode <- as.numeric(X_lam %*% means0[i_lam]) + (if (det_arm) 0 else offset_mode)
      eta_sig_mode <- as.numeric(X_sig %*% means0[i_sig]) + (if (det_arm) offset_mode else 0)
      eta_b_mode <- if (hazard) means0[i_b] else 0
      r_mode <- if (is_nb) exp(means0[i_logr]) else Inf
      sw_h <- cpp_distance_site_sweep(y, eta_lam_mode, eta_sig_mode, quad_xptr, K_max,
                                      nb = is_nb, r = r_mode, key = key_code,
                                      eta_b = eta_b_mode, headroom = headroom)
      sw_u <- cpp_distance_site_sweep(y, eta_lam_mode, eta_sig_mode, quad_xptr, K_max,
                                      nb = is_nb, r = r_mode, key = key_code,
                                      eta_b = eta_b_mode, headroom = -1L)
      .dist_score_gap(sw_h$grad_lam, sw_u$grad_lam, X_lam,
                      sw_h$grad_sig, sw_u$grad_sig, X_sig)
    }, error = function(e) NA_real_)
    if (is.finite(gap) && gap > .NMIX_SCORE_TOL) {
      h_next <- .nmix_widen_headroom(headroom, K_max)
      if (!is.null(h_next)) {
        cl <- match.call()
        cl$headroom <- h_next
        # `eval` is shadowed in this scope by the sweep closure above -- call
        # base::eval explicitly, not the closure.
        return(base::eval(cl, parent.frame()))
      }
    }
  }

  means <- res$beta_mean
  raw <- list(
    mixture = if (is_nb) "negbin" else "poisson",
    beta_lambda = means[i_lam], beta_sigma = means[i_sig],
    eta_b = if (hazard) means[i_b] else NA_real_,
    shape = if (hazard) exp(means[i_b]) else NA_real_,
    log_r = if (is_nb) means[i_logr] else NA_real_,
    r = if (is_nb) exp(means[i_logr]) else NA_real_,
    vcov = res$vcov, log_lik = res$log_lik, converged = TRUE,
    key = model$key, transect = model$transect, hazard = hazard, K_max = K_max)
  fit <- build_distance_fit(raw, model)
  # The field loads on the abundance (log lambda) arm by default, or on the
  # detection scale (log sigma) arm when the term sits in the detection formula.
  .tobs_attach_field_results(fit, res, det_arm, temporal, temporal_only, "abundance",
                             svc = svc, has_spatial = !is.null(spatial),
                             X_svc = X_lam, family = "distance")
}

# Areal-spatial binned distance sampling via NUTS (gcol33/tulpaObs#72): a FIXED-
# HYPER non-centered PROPER-CAR field on the abundance arm (log lambda) of the
# bin-multinomial distance marginal. The field precision (tau, rho) is fixed at the
# nested-Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~
# N(0, I) (z = Linv %*% raw) is sampled jointly with the coefficients via the
# distance NUTS field block (cpp_distance_nuts over nuts_field_block.h). The areal
# Laplace fit supplies warm coefficients + the field hyper; NUTS then samples the
# exact coefficient + whitened-field posterior. Half-normal key only (the spatial
# Laplace path, gcol33/tulpaObs#79); icar / car_proper / bym2 -- the intrinsic
# icar / bym2 fields sample via the #71 sum-to-zero reparameterisation (#113).
# Poisson or NB.
.tobs_fit_distance_nuts_spatial <- function(model, spatial = NULL, temporal = NULL,
                                            mixture = "poisson",
                                            K_max = NULL, sigma.beta = NULL,
                                            sigma.shape = 1.5, sigma.logr = NULL,
                                            n.iter = NULL, n.warmup = NULL,
                                            n.chains = NULL, max.treedepth = NULL,
                                            adapt.delta = NULL, seed = NULL,
                                            verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  # FIXED-HYPER non-centered field on the abundance (log lambda) arm from EITHER an
  # areal term (#72/#113) OR a temporal() term (#114); both share the sampling tail,
  # only the loading / field map / warm source differ. The hazard-rate key adds a
  # single global log-shape coordinate eta_b (#114): it is orthogonal to the field
  # block (the C++ target places it before the whitened field raw and sums its own
  # gradient), so the field path carries it by threading the extra base coordinate.
  temporal_only <- is.null(spatial) && !is.null(temporal)
  hazard <- identical(model$key, "hazard")
  if (!temporal_only) {
    .tobs_reject_weighted_spatial(spatial, "distance NUTS abundance spatial")
    if (!spatial$type %in% c("icar", "car_proper", "bym2"))
      stop(sprintf(paste0("distance() NUTS + areal spatial supports icar() / ",
                          "car_proper() / bym2() on the abundance arm; got '%s'. ",
                          "(tulpaObs#72, #113)"), spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (!temporal_only && spatial$n_units != n_sites)
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for distance NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  is_nb <- mixture %in% c("negbin", "NB")
  mix_code <- if (is_nb) "NB" else "P"
  X_lambda <- model$X_processes[[1]]; X_sigma <- model$X_processes[[2]]
  y <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  p_lam <- ncol(X_lambda); p_sig <- ncol(X_sigma)
  R_max <- if (length(y)) max(rowSums(y)) else 0L
  K_max <- if (is.null(K_max)) as.integer(3L * R_max + 100L) else as.integer(K_max)

  if (temporal_only) {
    ti <- as.integer(temporal$time_idx)
    if (length(ti) != n_sites)
      stop(sprintf(paste0("temporal term has %d time indices but the model has %d ",
                          "sites; one time index per site is required for distance ",
                          "NUTS + temporal."), length(ti), n_sites), call. = FALSE)
    n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
           else max(ti, na.rm = TRUE)
    nl <- .tobs_fit_distance_spatial(model, spatial = NULL, temporal = temporal,
                                     mixture = mixture, K_max = K_max,
                                     max_iter = 200L, tol = 1e-6, verbose = FALSE,
                                     integration = "grid")
    hyper <- nl$temporal_hyper
    hv <- function(k) suppressWarnings(as.numeric(hyper[[k]]))
    fl <- .tobs_nuts_temporal_loading(temporal$type, n_t, tau = hv("tau"), rho = hv("rho"))
    field_map <- ti
    n_field_units <- n_t
  } else {
    adj <- as.matrix(spatial$graph)
    # Warm coefficients + fixed field hyper (tau, rho) from the nested-Laplace fit.
    nl <- .tobs_fit_distance_spatial(model, spatial, mixture = mixture, K_max = K_max,
                                     max_iter = 200L, tol = 1e-6, verbose = FALSE,
                                     integration = "grid")
    hyper <- nl$spatial_hyper
    hv <- function(k) suppressWarnings(as.numeric(hyper[k]))
    fl <- .tobs_nuts_field_loading(adj, spatial$type, n_sites,
                                   tau = hv("tau"), rho = hv("rho"),
                                   sigma = hv("sigma"),
                                   scale_factor = spatial$scale_factor)
    field_map <- seq_len(n_sites)
    n_field_units <- n_sites
  }
  field_load <- fl$field_load; n_raw <- fl$n_raw

  cm <- nl$means
  beta0 <- c(unname(cm[paste0("lambda_", model$process_info[[1]]$coef_names)]),
             unname(cm[paste0("sigma_",  model$process_info[[2]]$coef_names)]))
  # Hazard-rate key: the global log-shape eta_b sits after the detection-scale block
  # and before the NB dispersion, matching the C++ theta layout (#114).
  if (hazard) {
    eta_b0 <- suppressWarnings(as.numeric(cm[["log_shape"]]))
    beta0 <- c(beta0, if (is.finite(eta_b0)) eta_b0 else log(2))
  }
  if (is_nb) beta0 <- c(beta0, log(if (is.finite(nl$r %||% NA_real_)) nl$r else 2))
  n_base <- p_lam + p_sig + (if (hazard) 1L else 0L) + if (is_nb) 1L else 0L
  # Warm-start raw: a full-rank car_proper field starts near the integrated field
  # via the precision Cholesky (z = Linv %*% raw -> raw = L %*% z); the intrinsic
  # icar / bym2 sum-to-zero loadings and the (rank-deficient) temporal loadings are
  # non-square, so their whitened raw starts at 0 (the well-conditioned fixed-hyper
  # geometry adapts in warmup, #71/#113/#114).
  raw0 <- if (!temporal_only && identical(spatial$type, "car_proper")) {
    L <- chol(.areal_Q(as.matrix(spatial$graph), fl$rho) * fl$tau +
              diag(1e-4 * fl$tau, n_sites))
    as.numeric(L %*% (nl$spatial_field %||% numeric(n_sites)))
  } else numeric(n_raw)
  theta0 <- c(beta0, raw0)
  inv_metric <- c(rep(0.1, n_base), rep(1, n_raw))

  spec <- list(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
               quad_xptr = cpp_distance_build_quad(
                 as.numeric(model$cutpoints), .dist_transect_code(model$transect),
                 as.integer(model$quad_order %||% 64L)),
               key = .dist_key_code(model$key), K_max = K_max, is_nb = is_nb,
               n_field_units = n_field_units, field_map = field_map,
               field_load = field_load)

  run_chain <- function(ch)
    cpp_distance_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      sigma_shape = sigma.shape, sigma_logr = sigma.logr,
                      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                      verbose = isTRUE(verbose))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("sigma_",  model$process_info[[2]]$coef_names),
           if (hazard) "log_shape", if (is_nb) "log_r",
           paste0("raw_", seq_len(n_raw)))
  run <- .tobs_nuts_field_draws(run_chain, n.chains, nms, n_base, n_raw, field_load)
  par <- run$par; cov <- run$cov

  lay <- .tobs_distance_nuts_layout(p_lam, p_sig, hazard, is_nb)
  marg <- .tobs_distance_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  eta_b_hat <- if (hazard) as.numeric(par[lay$log_shape]) else 0
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$sigma], eta_b = eta_b_hat,
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik
  raw_fit <- list(
    mixture = mix_code, key = model$key, transect = model$transect,
    hazard = hazard, nb = is_nb,
    beta_lambda = unname(par[lay$lambda]), beta_sigma = unname(par[lay$sigma]),
    eta_b = if (hazard) eta_b_hat else NA_real_,
    shape = if (hazard) exp(eta_b_hat) else NA_real_,
    log_r = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max)
  fit <- build_distance_fit(raw_fit, model)
  .tobs_nuts_field_attach(
    fit, run, ll_mean, n.chains,
    prior_type = if (temporal_only) temporal$type else spatial$type, fl = fl,
    temporal = if (temporal_only) temporal else NULL)
}
