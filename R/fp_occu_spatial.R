# =============================================================================
# fp_occu_spatial.R - areal-spatial multistate false-positive occupancy (#51)
#
# An ICAR or proper-CAR field on the OCCUPANCY (psi) arm of the Miller 2011 two-
# state false-positive model, via the shared areal-BFGS nested-Laplace driver
# (R/areal_bfgs.R): the fp_occu marginal exposes an analytic gradient
# (cpp_fp_occu_total_log_lik) but no analytic per-site Hessian, so the driver runs
# BFGS over (beta_psi, beta_p11, beta_p10, beta_b, z) + the CAR prior and forms the
# Laplace marginal from an FD-Hessian at the mode. The field z loads onto eta_psi
# exactly like the occupancy intercept (one unit per site); the false-positive
# arms (p11, p10, b) carry fixed effects only.
#
#   .tobs_fit_fp_occu_spatial()   dispatch from .tobs_fit_model (icar / car_proper)
# =============================================================================

.tobs_fit_fp_occu_spatial <- function(model, spatial, temporal = NULL,
                                      max_iter = 200L,
                                      tol = 1e-8, verbose = TRUE,
                                      integration = "grid") {
  temporal_only <- is.null(spatial)
  if (!temporal_only)
    .tobs_reject_weighted_spatial(spatial, "fp_occu occupancy spatial")
  # Detection-arm field (gcol33/tulpaObs#114): a field in the `detection=` formula
  # carries shared = c(occupancy, detection) = c(FALSE, TRUE). It loads on the
  # per-visit true-positive detection logit eta_p11 (a spatially-varying detection
  # probability) instead of the psi arm; the marginal exposes the per-visit p11
  # gradient (cpp_fp_occu_total_log_lik$grad_eta_p11), summed to a per-site field
  # gradient. The false-positive arms (p10, b) never carry a structured field.
  det_arm <- !temporal_only && isTRUE(spatial$shared[2L]) && !isTRUE(spatial$shared[1L])
  n_sites <- model$n_sites
  map <- seq_len(model$n_sites)
  # Temporal-only fit (gcol33/tulpaObs#114): the areal-BFGS driver runs the single
  # temporal block; otherwise the areal field is block 1 and the temporal block 2.
  field <- if (temporal_only) {
    list(.tobs_temporal_field_spec(temporal, model$n_sites, "fp_occu"))
  } else {
    field_sp <- .tobs_areal_field_spec(spatial, model$n_sites, "fp_occu", map)
    if (is.null(temporal)) field_sp
    else list(field_sp,
              .tobs_temporal_field_spec(temporal, model$n_sites, "fp_occu"))
  }

  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  p_psi <- ncol(X_psi); p_p11 <- ncol(X_p11); p_p10 <- ncol(X_p10); p_b <- ncol(X_b)
  y_long <- as.integer(model$y_long); site_idx <- as.integer(model$site_idx)

  off <- cumsum(c(0L, p_psi, p_p11, p_p10, p_b))
  i_psi <- off[1] + seq_len(p_psi); i_p11 <- off[2] + seq_len(p_p11)
  i_p10 <- off[3] + seq_len(p_p10); i_b   <- off[4] + seq_len(p_b)
  n_fixed <- off[5]

  eval <- function(theta_fix, offset) {
    # `offset` is per-site (length n_sites). The fp_occu detection design is
    # per-site (one true-positive logit per site, applied across every visit), so
    # a p11-arm field enters eta_p11 directly and its per-site gradient scatters
    # back with no aggregation (as on the psi arm).
    eta_psi <- as.numeric(X_psi %*% theta_fix[i_psi]) + (if (det_arm) 0 else offset)
    eta_p11 <- as.numeric(X_p11 %*% theta_fix[i_p11]) + (if (det_arm) offset else 0)
    out <- cpp_fp_occu_total_log_lik(
      y_long, site_idx, eta_psi, eta_p11,
      as.numeric(X_p10 %*% theta_fix[i_p10]), as.numeric(X_b %*% theta_fix[i_b]))
    g <- numeric(n_fixed)
    g[i_psi] <- as.numeric(crossprod(X_psi, out$grad_eta_psi))
    g[i_p11] <- as.numeric(crossprod(X_p11, out$grad_eta_p11))
    g[i_p10] <- as.numeric(crossprod(X_p10, out$grad_eta_p10))
    g[i_b]   <- as.numeric(crossprod(X_b,   out$grad_eta_b))
    list(log_lik = out$log_lik, grad_fixed = g,
         grad_eta = if (det_arm) out$grad_eta_p11 else out$grad_eta_psi)
  }

  warm <- tryCatch(
    fp_occu_laplace(y = y_long, site_idx = site_idx, X_psi = X_psi, X_p11 = X_p11,
                    X_p10 = X_p10, X_b = X_b, sigma_beta = NULL,
                    max_iter = as.integer(max_iter), tol = as.numeric(tol),
                    verbose = FALSE),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm))
    c(warm$beta_psi, warm$beta_p11, warm$beta_p10, warm$beta_b)
  else c(0, rep(0, p_psi - 1L), rep(0, p_p11),
         stats::qlogis(0.05), rep(0, p_p10 - 1L), rep(0, p_b))

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol, label = "fp-occu-spatial",
                              integration = integration)
  if (!isTRUE(res$ok))
    stop("fp_occu() areal spatial fit produced no usable grid point.", call. = FALSE)

  nm <- c(paste0("psi_", model$process_info[[1]]$coef_names),
          paste0("p11_", model$process_info[[2]]$coef_names),
          paste0("p10_", model$process_info[[3]]$coef_names),
          paste0("b_",   model$process_info[[4]]$coef_names))
  means <- res$beta_mean; names(means) <- nm
  V <- res$vcov; dimnames(V) <- list(nm, nm)
  # Posterior occupancy w1 at the integrated estimate + field (for fitted()). The
  # field enters the psi arm, or the p11 detection arm when det_arm (#114).
  cl <- .tobs_clamp_eta
  fld <- res$field_mean
  eta_psi <- cl(as.numeric(X_psi %*% means[i_psi]) + (if (det_arm) 0 else fld[map]))
  eta_p11 <- as.numeric(X_p11 %*% means[i_p11]) + (if (det_arm) fld[map] else 0)
  ev <- cpp_fp_occu_total_log_lik(
    y_long, site_idx, eta_psi, eta_p11,
    as.numeric(X_p10 %*% means[i_p10]), as.numeric(X_b %*% means[i_b]))
  raw <- list(
    beta_psi = means[i_psi], beta_p11 = means[i_p11],
    beta_p10 = means[i_p10], beta_b = means[i_b],
    means = means, vcov = V, theta_se = sqrt(pmax(diag(V), 0)),
    log_lik = res$log_lik, w1 = ev$w1, converged = TRUE, n_iter = NA_integer_,
    coef_names = nm)
  fit <- build_fp_occu_fit(raw, model)
  fit$method <- "nested_laplace"
  fit$spatial_integration <- res$integration
  fit$spatial_pareto_k <- res$pareto_k
  if (temporal_only) {
    # The single temporal block is reported by the driver as block 1
    # (field_mean / hyper); relabel it as the temporal field (#114).
    fit$temporal <- temporal
    fit$temporal_field <- res$field_mean
    fit$temporal_hyper <- res$hyper
  } else {
    fit$spatial_field <- res$field_mean
    fit$spatial_hyper <- res$hyper
    # The field loads on the occupancy (psi) arm by default, or on the per-visit
    # true-positive detection logit p11 when the term sits in `detection=` (#114).
    fit$spatial_field_arm <- if (det_arm) "detection" else "occupancy"
    if (!is.null(temporal)) {
      fit$temporal <- temporal
      fit$temporal_field <- res$temporal_field
      fit$temporal_hyper <- res$temporal_hyper
    }
  }
  fit
}

# Areal-spatial multistate false-positive occupancy via NUTS (gcol33/tulpaObs#72):
# a FIXED-HYPER non-centered PROPER-CAR field on the occupancy (psi) arm of the
# two-state false-positive marginal. The field precision (tau, rho) is fixed at the
# nested-Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~
# N(0, I) (z = Linv %*% raw) is sampled jointly with the four logit arms' fixed
# effects via the fp_occu NUTS field block (cpp_fp_occu_nuts over
# nuts_field_block.h). The areal Laplace fit supplies warm coefficients + the field
# hyper. The false-positive arms (p11 / p10 / b) carry fixed effects only;
# icar / car_proper / bym2 -- the intrinsic icar / bym2 fields sample via the #71
# sum-to-zero reparameterisation (#113). Occupancy fields are more weakly
# identified than count fields (one binary site per node).
.tobs_fit_fp_occu_nuts_spatial <- function(model, spatial = NULL, temporal = NULL,
                                           sigma.beta = 10,
                                           n.iter = 1000L, n.warmup = 1000L,
                                           n.chains = 1L, max.treedepth = 10L,
                                           adapt.delta = 0.9, seed = 1L,
                                           verbose = FALSE) {
  # FIXED-HYPER non-centered field on the occupancy (psi) arm from EITHER an areal
  # term (#72/#113) OR a temporal() term (#114); both share the sampling tail,
  # differing only in the loading / field map / warm source.
  temporal_only <- is.null(spatial) && !is.null(temporal)
  if (!temporal_only) {
    .tobs_reject_weighted_spatial(spatial, "fp_occu NUTS occupancy spatial")
    if (isTRUE(spatial$shared[2L]) && !isTRUE(spatial$shared[1L]))
      stop(paste0("fp_occu() NUTS carries the areal field on the occupancy (psi) ",
                  "arm; a detection-arm field (a spatially-varying p11 logit) is ",
                  "wired under method = \"nested_laplace\". (tulpaObs#114)"),
           call. = FALSE)
    if (!spatial$type %in% c("icar", "car_proper", "bym2"))
      stop(sprintf(paste0("fp_occu() NUTS + areal spatial supports icar() / ",
                          "car_proper() / bym2() on the psi arm; got '%s'. ",
                          "(tulpaObs#72, #113)"), spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (!temporal_only && spatial$n_units != n_sites)
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for fp_occu NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]

  if (temporal_only) {
    ti <- as.integer(temporal$time_idx)
    if (length(ti) != n_sites)
      stop(sprintf(paste0("temporal term has %d time indices but the model has %d ",
                          "sites; one time index per site is required for fp_occu ",
                          "NUTS + temporal."), length(ti), n_sites), call. = FALSE)
    n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
           else max(ti, na.rm = TRUE)
    nl <- .tobs_fit_fp_occu_spatial(model, spatial = NULL, temporal = temporal,
                                    max_iter = 200L, tol = 1e-8, verbose = FALSE,
                                    integration = "grid")
    hyper <- nl$temporal_hyper
    hv <- function(k) suppressWarnings(as.numeric(hyper[[k]]))
    fl <- .tobs_nuts_temporal_loading(temporal$type, n_t, tau = hv("tau"), rho = hv("rho"))
    field_map <- ti
    n_field_units <- n_t
  } else {
    adj <- as.matrix(spatial$graph)
    # Warm coefficients + fixed field hyper (tau, rho) from the nested-Laplace fit.
    nl <- .tobs_fit_fp_occu_spatial(model, spatial, max_iter = 200L, tol = 1e-8,
                                    verbose = FALSE, integration = "grid")
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

  cm <- as.numeric(nl$means)
  n_base <- length(cm)
  # car_proper warm-starts raw near the integrated field; the intrinsic icar / bym2
  # sum-to-zero loadings and the (rank-deficient) temporal loadings are non-square,
  # so their whitened raw starts at 0 (#71/#113/#114).
  raw0 <- if (!temporal_only && identical(spatial$type, "car_proper")) {
    L <- chol(.areal_Q(as.matrix(spatial$graph), fl$rho) * fl$tau +
              diag(1e-4 * fl$tau, n_sites))
    as.numeric(L %*% (nl$spatial_field %||% numeric(n_sites)))
  } else numeric(n_raw)
  theta0 <- c(cm, raw0)
  inv_metric <- c(rep(0.5, n_base), rep(1, n_raw))

  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
               n_sites = n_sites, n_field_units = n_field_units,
               field_map = field_map, field_load = field_load)

  run_chain <- function(ch)
    cpp_fp_occu_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                     inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                     n_warmup = as.integer(n.warmup),
                     max_treedepth = as.integer(max.treedepth),
                     adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                     verbose = isTRUE(verbose))
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("psi_", model$process_info[[1]]$coef_names),
           paste0("p11_", model$process_info[[2]]$coef_names),
           paste0("p10_", model$process_info[[3]]$coef_names),
           paste0("b_",   model$process_info[[4]]$coef_names),
           paste0("raw_", seq_len(n_raw)))
  colnames(draws) <- nms
  b_idx <- seq_len(n_base)
  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws[, b_idx, drop = FALSE])
  raw_idx <- n_base + seq_len(n_raw)
  z_mean <- as.numeric(field_load %*% colMeans(draws[, raw_idx, drop = FALSE]))

  lay <- .tobs_fp_occu_nuts_layout(ncol(X_psi), ncol(X_p11), ncol(X_p10), ncol(X_b))
  marg <- .tobs_fp_occu_nuts_marginal(model)
  ev_mean <- marg$eval_beta(par[lay$psi], par[lay$p11], par[lay$p10], par[lay$b])
  raw_fit <- list(means = unname(par[b_idx]), vcov = cov, coef_names = nms[b_idx],
                  log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
                  w1 = ev_mean$w1, converged = TRUE)
  fit <- build_fp_occu_fit(raw_fit, model)
  fit$draws <- draws[, b_idx, drop = FALSE]
  fit$means <- par[b_idx]; fit$sds <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms[b_idx]
  fit$vcov <- cov
  fit$n_samples <- nrow(draws); fit$log_prob <- rep(ev_mean$log_lik, nrow(draws))
  accept <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  fit$accept_prob <- accept; fit$divergent <- divergent
  fit$method <- "nuts"
  if (temporal_only) { fit$temporal <- temporal; fit$temporal_field <- z_mean }
  else               fit$spatial_field <- z_mean
  fit$nuts <- list(accept_prob = accept, divergent = divergent,
                   treedepth = as.integer(unlist(lapply(chains, `[[`, "treedepth"))),
                   epsilon = chains[[1L]]$epsilon, n_chains = as.integer(n.chains),
                   divergent_total = sum(divergent), tau = fl$tau, rho = fl$rho,
                   prior_type = if (temporal_only) temporal$type else spatial$type,
                   fixed_hyper = TRUE)
  fit
}
