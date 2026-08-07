# =============================================================================
# dyn_abun_spatial.R - areal-spatial Dail-Madsen open N-mixture (tulpaObs#51)
#
# An ICAR / proper-CAR / BYM2 field on the INITIAL-abundance arm (log lambda_1)
# via the shared areal-BFGS nested-Laplace driver (R/areal_bfgs.R): the forward-
# HMM marginal exposes an analytic gradient (cpp_dyn_abun_total_log_lik) but no
# analytic per-site Hessian, so the driver runs BFGS over (betas, field) + the
# field prior and forms the Laplace marginal from an FD-Hessian at the mode. This
# file supplies only the family eval (log-lik + per-arm gradient + the per-site
# initial-abundance eta-gradient the field scatters). The field loads onto
# eta_lambda exactly like the initial-abundance intercept (one unit per site).
#
#   .tobs_fit_dyn_abun_spatial()   dispatch from .tobs_fit_model
# =============================================================================

.tobs_fit_dyn_abun_spatial <- function(model, spatial, temporal = NULL,
                                       svc = NULL, mixture = "poisson",
                                       K_max = NULL, max_iter = 300L, tol = 1e-8,
                                       verbose = TRUE, integration = "grid") {
  temporal_only <- is.null(spatial) && !is.null(temporal)
  if (!is.null(spatial))
    .tobs_reject_weighted_spatial(spatial, "dyn_abun abundance spatial")
  # Detection-arm field (gcol33/tulpaObs#114): a field in the `detection=` formula
  # carries shared = c(abundance, detection) = c(FALSE, TRUE). The dyn_abun
  # detection design is per-site ([n_sites x p]), so a site-level detection field
  # (a spatially-varying detection probability applied across every season's obs
  # pmf) loads on eta_p directly and the marginal's per-site grad_eta_p scatters
  # back with no aggregation. omega / gamma never carry a structured field.
  det_arm <- !is.null(spatial) && isTRUE(spatial$shared[2L]) &&
             !isTRUE(spatial$shared[1L])
  map <- seq_len(model$n_sites)
  X_lam <- model$X_processes[[1]]
  .tobs_check_svc_arm(svc, det_arm, "dyn_abun")
  field <- .tobs_build_field_spec(spatial, temporal, "dyn_abun", model$n_sites, map,
                                  svc = svc, X_svc = X_lam)

  X_p <- model$X_processes[[2]]
  X_om  <- model$X_processes[[3]]; X_gm <- model$X_processes[[4]]
  p_lam <- ncol(X_lam); p_p <- ncol(X_p); p_om <- ncol(X_om); p_gm <- ncol(X_gm)
  y_flat <- as.integer(model$y_flat); N <- model$n_sites
  T <- model$n_seasons; J <- model$max_visits; K <- model$K_max
  is_nb <- mixture %in% c("negbin", "NB")

  off <- cumsum(c(0L, p_lam, p_p, p_om, p_gm))
  i_lam <- off[1] + seq_len(p_lam); i_p <- off[2] + seq_len(p_p)
  i_om  <- off[3] + seq_len(p_om);  i_gm <- off[4] + seq_len(p_gm)
  i_logr <- if (is_nb) off[5] + 1L else NA_integer_
  n_fixed <- off[5] + if (is_nb) 1L else 0L

  eval <- function(theta_fix, offset) {
    # `offset` is per-site (length N). On the abundance arm it enters eta_lambda;
    # on the detection arm (det_arm) it enters the per-site eta_p directly (both
    # arms are per-site here, so no per-observation aggregation is needed).
    eta_lam <- as.numeric(X_lam %*% theta_fix[i_lam]) + (if (det_arm) 0 else offset)
    eta_p   <- as.numeric(X_p %*% theta_fix[i_p]) + (if (det_arm) offset else 0)
    out <- cpp_dyn_abun_total_log_lik(
      y_flat, N, T, J, K, eta_lam, eta_p,
      as.numeric(X_om %*% theta_fix[i_om]),
      as.numeric(X_gm %*% theta_fix[i_gm]),
      use_nb = is_nb, eta_logr = if (is_nb) theta_fix[i_logr] else 0.0)
    g <- numeric(n_fixed)
    g[i_lam] <- as.numeric(crossprod(X_lam, out$grad_eta_lambda))
    g[i_p]   <- as.numeric(crossprod(X_p,   out$grad_eta_p))
    g[i_om]  <- as.numeric(crossprod(X_om,  out$grad_eta_omega))
    g[i_gm]  <- as.numeric(crossprod(X_gm,  out$grad_eta_gamma))
    if (is_nb) g[i_logr] <- as.numeric(out$grad_eta_logr)
    list(log_lik = out$log_lik, grad_fixed = g,
         grad_eta = if (det_arm) out$grad_eta_p else out$grad_eta_lambda)
  }

  warm <- tryCatch(dyn_abun_laplace(
    y_flat = y_flat, n_sites = N, T = T, J = J, K_max = K,
    X_lambda = X_lam, X_p = X_p, X_omega = X_om, X_gamma = X_gm,
    mixture = if (is_nb) "negbin" else "poisson", verbose = FALSE),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm) && length(warm$means) == n_fixed) warm$means
                else numeric(n_fixed)

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol, label = "dyn-abun-spatial",
                              integration = integration)

  nm <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
          paste0("p_",      model$process_info[[2]]$coef_names),
          paste0("omega_",  model$process_info[[3]]$coef_names),
          paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (is_nb) nm <- c(nm, "log_r")
  means <- res$beta_mean; names(means) <- nm
  V <- res$vcov; dimnames(V) <- list(nm, nm)
  raw <- list(
    beta_lambda = means[i_lam], beta_p = means[i_p],
    beta_omega = means[i_om], beta_gamma = means[i_gm],
    log_r = if (is_nb) means[i_logr] else NA_real_,
    r = if (is_nb) exp(means[i_logr]) else NA_real_,
    mixture = if (is_nb) "negbin" else "poisson",
    means = means, vcov = V, log_lik = res$log_lik, mean_N1 = NULL,
    K_max = K, converged = TRUE, n_iter = NA_integer_, coef_names = nm)
  fit <- build_dyn_abun_fit(raw, model)
  # The field loads on the initial-abundance (log lambda_1) arm by default, or on
  # the per-site detection logit eta_p when the term sits in `detection=` (#114).
  .tobs_attach_field_results(fit, res, det_arm, temporal, temporal_only, "abundance",
                             svc = svc, has_spatial = !is.null(spatial),
                             X_svc = X_lam, family = "dyn_abun")
}

# Areal-spatial Dail-Madsen open N-mixture via NUTS (gcol33/tulpaObs#72): a FIXED-
# HYPER non-centered PROPER-CAR field on the initial-abundance (log lambda_1) arm of
# the forward-HMM marginal. The field precision (tau, rho) is fixed at the nested-
# Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~ N(0, I)
# (z = Linv %*% raw) is sampled jointly with the four arms' coefficients via the
# dyn_abun NUTS field block (cpp_dyn_abun_nuts over nuts_field_block.h). The areal
# Laplace fit supplies warm coefficients + the field hyper. icar / car_proper /
# bym2 -- the intrinsic icar / bym2 fields sample via the #71 sum-to-zero
# reparameterisation (#113); Poisson or NB initial abundance.
.tobs_fit_dyn_abun_nuts_spatial <- function(model, spatial, mixture = "poisson",
                                            K_max = NULL, sigma.beta = NULL,
                                            n.iter = NULL, n.warmup = NULL,
                                            n.chains = NULL, max.treedepth = NULL,
                                            adapt.delta = NULL, seed = NULL,
                                            verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  .tobs_reject_weighted_spatial(spatial, "dyn_abun NUTS abundance spatial")
  if (isTRUE(spatial$shared[2L]) && !isTRUE(spatial$shared[1L]))
    stop(paste0("dyn_abun() NUTS carries the areal field on the initial-abundance ",
                "arm; a detection-arm field (a spatially-varying detection logit) ",
                "is wired under method = \"nested_laplace\". (tulpaObs#114)"),
         call. = FALSE)
  if (!spatial$type %in% c("icar", "car_proper", "bym2"))
    stop(sprintf(paste0("dyn_abun() NUTS + areal spatial supports icar() / ",
                        "car_proper() / bym2() on the initial-abundance arm; got ",
                        "'%s'. (tulpaObs#72, #113)"), spatial$type), call. = FALSE)
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites)
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for dyn_abun NUTS."),
                 spatial$n_units, n_sites), call. = FALSE)
  use_nb <- identical(model$mixture %||% mixture, "negbin") ||
            mixture %in% c("negbin", "NB")
  X_lam <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_om  <- model$X_processes[[3]]; X_gm <- model$X_processes[[4]]
  adj <- as.matrix(spatial$graph)

  # Warm coefficients + fixed field hyper (tau, rho) from the nested-Laplace fit.
  nl <- .tobs_fit_dyn_abun_spatial(model, spatial,
                                   mixture = if (use_nb) "negbin" else "poisson",
                                   K_max = K_max, max_iter = 300L, tol = 1e-8,
                                   verbose = FALSE, integration = "grid")
  hyper <- nl$spatial_hyper
  hv <- function(k) suppressWarnings(as.numeric(hyper[k]))
  fl <- .tobs_nuts_field_loading(adj, spatial$type, n_sites,
                                 tau = hv("tau"), rho = hv("rho"),
                                 sigma = hv("sigma"),
                                 scale_factor = spatial$scale_factor)
  field_load <- fl$field_load; n_raw <- fl$n_raw

  cm <- as.numeric(nl$means)
  n_base <- length(cm)
  # car_proper warm-starts raw near the integrated field; icar / bym2 (non-square
  # sum-to-zero loadings) start raw at 0 (#71/#113).
  raw0 <- if (identical(spatial$type, "car_proper")) {
    L <- chol(.areal_Q(adj, fl$rho) * fl$tau + diag(1e-4 * fl$tau, n_sites))
    as.numeric(L %*% (nl$spatial_field %||% numeric(n_sites)))
  } else numeric(n_raw)
  theta0 <- c(cm, raw0)
  inv_metric <- c(rep(0.2, n_base), rep(1, n_raw))

  spec <- list(y = as.integer(model$y_flat), n_sites = n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = X_lam, X_p = X_p, X_omega = X_om, X_gamma = X_gm,
               use_nb = use_nb, n_field_units = n_sites,
               field_map = seq_len(n_sites), field_load = field_load)

  run_chain <- function(ch)
    cpp_dyn_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                      verbose = isTRUE(verbose))
  nms <- .tobs_dyn_abun_nuts_names(model, use_nb, n_raw)
  run <- .tobs_nuts_field_draws(run_chain, n.chains, nms, n_base, n_raw, field_load)

  ev  <- .tobs_dyn_abun_nuts_eval(model, run$par, X_lam, X_p, X_om, X_gm, use_nb)
  fit <- build_dyn_abun_fit(
    .tobs_dyn_abun_nuts_raw(run, nms, ev, model, use_nb), model)
  .tobs_nuts_field_attach(fit, run, ev$log_lik, n.chains,
                          prior_type = spatial$type, fl = fl)
}

# Coefficient + whitened-field column names for a field dyn_abun NUTS run.
.tobs_dyn_abun_nuts_names <- function(model, use_nb, n_raw) {
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           paste0("omega_",  model$process_info[[3]]$coef_names),
           paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (use_nb) nms <- c(nms, "log_r")
  c(nms, paste0("raw_", seq_len(n_raw)))
}

# Forward-HMM marginal at the posterior-mean coefficients.
.tobs_dyn_abun_nuts_eval <- function(model, par, X_lam, X_p, X_om, X_gm, use_nb) {
  lay <- .tobs_dyn_abun_nuts_layout(ncol(X_lam), ncol(X_p), ncol(X_om), ncol(X_gm),
                                    use_nb = use_nb)
  log_r <- if (use_nb) as.numeric(par[lay$logr]) else NA_real_
  ev <- .tobs_dyn_abun_nuts_marginal(model)$eval_beta(
    par[lay$lambda], par[lay$p], par[lay$omega], par[lay$gamma],
    if (use_nb) log_r else 0)
  c(ev, list(log_r = log_r))
}

# The raw fit list build_dyn_abun_fit() consumes, from the posterior mean.
.tobs_dyn_abun_nuts_raw <- function(run, nms, ev, model, use_nb) {
  list(means = unname(run$par[run$b_idx]), vcov = run$cov,
       coef_names = nms[run$b_idx],
       log_lik = ev$log_lik, log_lik_site = ev$log_lik_site,
       mean_N1 = ev$mean_N1, K_max = model$K_max, converged = TRUE,
       mixture = if (use_nb) "negbin" else "poisson",
       log_r = ev$log_r, r = if (use_nb) exp(ev$log_r) else NA_real_)
}

# Temporal-field Dail-Madsen open N-mixture via NUTS (gcol33/tulpaObs#114): a
# FIXED-HYPER non-centered ar1 / rw1 / rw2 / iid field on the initial-abundance
# (log lambda_1) arm of the forward-HMM marginal. Structurally identical to the
# areal NUTS path -- the temporal field is a GMRF whose whitened loading is fixed
# at the nested-Laplace temporal-only posterior mean and whose per-site field_map
# is the period index (many sites share a period), so it rides the SAME dyn_abun
# NUTS field block (cpp_dyn_abun_nuts over nuts_field_block.h) with no engine
# change. Poisson or NB initial abundance; temporal-only (no simultaneous areal
# field on this path).
.tobs_fit_dyn_abun_nuts_temporal <- function(model, temporal, mixture = "poisson",
                                             K_max = NULL, sigma.beta = NULL,
                                             n.iter = NULL, n.warmup = NULL,
                                             n.chains = NULL, max.treedepth = NULL,
                                             adapt.delta = NULL, seed = NULL,
                                             verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  n_sites <- model$n_sites
  ti <- as.integer(temporal$time_idx)
  if (length(ti) != n_sites)
    stop(sprintf(paste0("temporal term has %d time indices but the model has %d ",
                        "sites; one time index per site is required for dyn_abun ",
                        "NUTS + temporal."), length(ti), n_sites), call. = FALSE)
  n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
         else max(ti, na.rm = TRUE)
  use_nb <- identical(model$mixture %||% mixture, "negbin") ||
            mixture %in% c("negbin", "NB")
  X_lam <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_om  <- model$X_processes[[3]]; X_gm <- model$X_processes[[4]]

  # Warm coefficients + fixed temporal field hyper (tau[, rho]) from the
  # nested-Laplace temporal-only fit (spatial = NULL).
  nl <- .tobs_fit_dyn_abun_spatial(model, spatial = NULL, temporal = temporal,
                                   mixture = if (use_nb) "negbin" else "poisson",
                                   K_max = K_max, max_iter = 300L, tol = 1e-8,
                                   verbose = FALSE, integration = "grid")
  hyper <- nl$temporal_hyper
  hv <- function(k) suppressWarnings(as.numeric(hyper[k]))
  fl <- .tobs_nuts_temporal_loading(temporal$type, n_t,
                                    tau = hv("tau"), rho = hv("rho"))
  field_load <- fl$field_load; n_raw <- fl$n_raw

  cm <- as.numeric(nl$means)
  n_base <- length(cm)
  theta0 <- c(cm, numeric(n_raw))                 # raw starts flat (sum-to-zero safe)
  inv_metric <- c(rep(0.2, n_base), rep(1, n_raw))

  spec <- list(y = as.integer(model$y_flat), n_sites = n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = X_lam, X_p = X_p, X_omega = X_om, X_gamma = X_gm,
               use_nb = use_nb, n_field_units = n_t,
               field_map = ti, field_load = field_load)

  run_chain <- function(ch)
    cpp_dyn_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                      verbose = isTRUE(verbose))
  nms <- .tobs_dyn_abun_nuts_names(model, use_nb, n_raw)
  run <- .tobs_nuts_field_draws(run_chain, n.chains, nms, n_base, n_raw, field_load)

  ev  <- .tobs_dyn_abun_nuts_eval(model, run$par, X_lam, X_p, X_om, X_gm, use_nb)
  fit <- build_dyn_abun_fit(
    .tobs_dyn_abun_nuts_raw(run, nms, ev, model, use_nb), model)
  .tobs_nuts_field_attach(fit, run, ev$log_lik, n.chains,
                          prior_type = temporal$type, fl = fl,
                          temporal = temporal)
}
