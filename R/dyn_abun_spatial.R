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
                                       mixture = "poisson",
                                       K_max = NULL, max_iter = 300L, tol = 1e-8,
                                       verbose = TRUE, integration = "grid") {
  .tobs_reject_weighted_spatial(spatial, "dyn_abun abundance spatial")
  map <- seq_len(model$n_sites)
  field_sp <- .tobs_areal_field_spec(spatial, model$n_sites, "dyn_abun", map)
  field <- if (is.null(temporal)) field_sp else {
    list(field_sp,
         .tobs_temporal_field_spec(temporal, model$n_sites, "dyn_abun"))
  }

  X_lam <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
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
    eta_lam <- as.numeric(X_lam %*% theta_fix[i_lam]) + offset
    out <- cpp_dyn_abun_total_log_lik(
      y_flat, N, T, J, K, eta_lam,
      as.numeric(X_p %*% theta_fix[i_p]), as.numeric(X_om %*% theta_fix[i_om]),
      as.numeric(X_gm %*% theta_fix[i_gm]),
      use_nb = is_nb, eta_logr = if (is_nb) theta_fix[i_logr] else 0.0)
    g <- numeric(n_fixed)
    g[i_lam] <- as.numeric(crossprod(X_lam, out$grad_eta_lambda))
    g[i_p]   <- as.numeric(crossprod(X_p,   out$grad_eta_p))
    g[i_om]  <- as.numeric(crossprod(X_om,  out$grad_eta_omega))
    g[i_gm]  <- as.numeric(crossprod(X_gm,  out$grad_eta_gamma))
    if (is_nb) g[i_logr] <- as.numeric(out$grad_eta_logr)
    list(log_lik = out$log_lik, grad_fixed = g, grad_eta = out$grad_eta_lambda)
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
  if (!isTRUE(res$ok))
    stop("dyn_abun() areal spatial fit produced no usable grid point.", call. = FALSE)

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
  fit$method <- "nested_laplace"
  fit$spatial_field <- res$field_mean
  fit$spatial_hyper <- res$hyper
  fit$spatial_integration <- res$integration
  fit$spatial_pareto_k <- res$pareto_k
  if (!is.null(temporal)) {
    fit$temporal <- temporal
    fit$temporal_field <- res$temporal_field
    fit$temporal_hyper <- res$temporal_hyper
  }
  fit
}

# Areal-spatial Dail-Madsen open N-mixture via NUTS (gcol33/tulpaObs#72): a FIXED-
# HYPER non-centered PROPER-CAR field on the initial-abundance (log lambda_1) arm of
# the forward-HMM marginal. The field precision (tau, rho) is fixed at the nested-
# Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~ N(0, I)
# (z = Linv %*% raw) is sampled jointly with the four arms' coefficients via the
# dyn_abun NUTS field block (cpp_dyn_abun_nuts over nuts_field_block.h). The areal
# Laplace fit supplies warm coefficients + the field hyper. car_proper only
# (intrinsic icar = #71); Poisson or NB initial abundance.
.tobs_fit_dyn_abun_nuts_spatial <- function(model, spatial, mixture = "poisson",
                                            K_max = NULL, sigma.beta = 10,
                                            n.iter = 1000L, n.warmup = 1000L,
                                            n.chains = 1L, max.treedepth = 10L,
                                            adapt.delta = 0.9, seed = 1L,
                                            verbose = FALSE) {
  .tobs_reject_weighted_spatial(spatial, "dyn_abun NUTS abundance spatial")
  if (!identical(spatial$type, "car_proper"))
    stop(sprintf(paste0("dyn_abun() NUTS + areal spatial supports the proper-CAR ",
                        "field car_proper(); the intrinsic '%s' field needs a ",
                        "sum-to-zero reparameterisation for NUTS -- use method = ",
                        "\"nested_laplace\" for the icar()/bym2() areal fit. ",
                        "(tulpaObs#72)"), spatial$type), call. = FALSE)
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites)
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for dyn_abun NUTS.",
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
  tau <- max(unname(hyper[["tau"]]), 1e-3)
  rho <- min(max(unname(hyper[["rho"]]), 0.01), 0.99)
  Linv <- .tobs_field_linv(adj, tau, rho, n_sites)

  cm <- as.numeric(nl$means)
  n_base <- length(cm)
  L <- chol(.areal_Q(adj, rho) * tau + diag(1e-4 * tau, n_sites))
  raw0 <- as.numeric(L %*% (nl$spatial_field %||% numeric(n_sites)))
  theta0 <- c(cm, raw0)
  inv_metric <- c(rep(0.2, n_base), rep(1, n_sites))

  spec <- list(y = as.integer(model$y_flat), n_sites = n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = X_lam, X_p = X_p, X_omega = X_om, X_gamma = X_gm,
               use_nb = use_nb, n_field_units = n_sites,
               field_map = seq_len(n_sites), field_Linv = Linv)

  run_chain <- function(ch)
    cpp_dyn_abun_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                      verbose = isTRUE(verbose))
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           paste0("omega_",  model$process_info[[3]]$coef_names),
           paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (use_nb) nms <- c(nms, "log_r")
  nms <- c(nms, paste0("raw_", seq_len(n_sites)))
  colnames(draws) <- nms
  b_idx <- seq_len(n_base)
  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws[, b_idx, drop = FALSE])
  raw_idx <- n_base + seq_len(n_sites)
  z_mean <- as.numeric(Linv %*% colMeans(draws[, raw_idx, drop = FALSE]))

  lay <- .tobs_dyn_abun_nuts_layout(ncol(X_lam), ncol(X_p), ncol(X_om), ncol(X_gm),
                                    use_nb = use_nb)
  marg <- .tobs_dyn_abun_nuts_marginal(model)
  log_r <- if (use_nb) as.numeric(par[lay$logr]) else NA_real_
  ev_mean <- marg$eval_beta(par[lay$lambda], par[lay$p], par[lay$omega],
                            par[lay$gamma], if (use_nb) log_r else 0)
  raw_fit <- list(means = unname(par[b_idx]), vcov = cov, coef_names = nms[b_idx],
                  log_lik = ev_mean$log_lik, log_lik_site = ev_mean$log_lik_site,
                  mean_N1 = ev_mean$mean_N1, K_max = model$K_max, converged = TRUE,
                  mixture = if (use_nb) "negbin" else "poisson",
                  log_r = log_r, r = if (use_nb) exp(log_r) else NA_real_)
  fit <- build_dyn_abun_fit(raw_fit, model)
  fit$draws <- draws[, b_idx, drop = FALSE]
  fit$means <- par[b_idx]; fit$sds <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms[b_idx]
  fit$vcov <- cov
  fit$n_samples <- nrow(draws); fit$log_prob <- rep(ev_mean$log_lik, nrow(draws))
  accept <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(chains, `[[`, "divergent"))
  fit$accept_prob <- accept; fit$divergent <- divergent
  fit$method <- "nuts"; fit$spatial_field <- z_mean
  fit$nuts <- list(accept_prob = accept, divergent = divergent,
                   treedepth = as.integer(unlist(lapply(chains, `[[`, "treedepth"))),
                   epsilon = chains[[1L]]$epsilon, n_chains = as.integer(n.chains),
                   divergent_total = sum(divergent), tau = tau, rho = rho,
                   prior_type = spatial$type, fixed_hyper = TRUE)
  fit
}
