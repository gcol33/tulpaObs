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
  .tobs_reject_weighted_spatial(spatial, "fp_occu occupancy spatial")
  map <- seq_len(model$n_sites)
  field_sp <- .tobs_areal_field_spec(spatial, model$n_sites, "fp_occu", map)
  field <- if (is.null(temporal)) field_sp else {
    list(field_sp,
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
    eta_psi <- as.numeric(X_psi %*% theta_fix[i_psi]) + offset
    out <- cpp_fp_occu_total_log_lik(
      y_long, site_idx, eta_psi,
      as.numeric(X_p11 %*% theta_fix[i_p11]), as.numeric(X_p10 %*% theta_fix[i_p10]),
      as.numeric(X_b %*% theta_fix[i_b]))
    g <- numeric(n_fixed)
    g[i_psi] <- as.numeric(crossprod(X_psi, out$grad_eta_psi))
    g[i_p11] <- as.numeric(crossprod(X_p11, out$grad_eta_p11))
    g[i_p10] <- as.numeric(crossprod(X_p10, out$grad_eta_p10))
    g[i_b]   <- as.numeric(crossprod(X_b,   out$grad_eta_b))
    list(log_lik = out$log_lik, grad_fixed = g, grad_eta = out$grad_eta_psi)
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
  # Posterior occupancy w1 at the integrated estimate + field (for fitted()).
  cl <- .tobs_clamp_eta
  eta_psi <- cl(as.numeric(X_psi %*% means[i_psi]) + res$field_mean[map])
  ev <- cpp_fp_occu_total_log_lik(
    y_long, site_idx, eta_psi, as.numeric(X_p11 %*% means[i_p11]),
    as.numeric(X_p10 %*% means[i_p10]), as.numeric(X_b %*% means[i_b]))
  raw <- list(
    beta_psi = means[i_psi], beta_p11 = means[i_p11],
    beta_p10 = means[i_p10], beta_b = means[i_b],
    means = means, vcov = V, theta_se = sqrt(pmax(diag(V), 0)),
    log_lik = res$log_lik, w1 = ev$w1, converged = TRUE, n_iter = NA_integer_,
    coef_names = nm)
  fit <- build_fp_occu_fit(raw, model)
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

# Areal-spatial multistate false-positive occupancy via NUTS (gcol33/tulpaObs#72):
# a FIXED-HYPER non-centered PROPER-CAR field on the occupancy (psi) arm of the
# two-state false-positive marginal. The field precision (tau, rho) is fixed at the
# nested-Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~
# N(0, I) (z = Linv %*% raw) is sampled jointly with the four logit arms' fixed
# effects via the fp_occu NUTS field block (cpp_fp_occu_nuts over
# nuts_field_block.h). The areal Laplace fit supplies warm coefficients + the field
# hyper. The false-positive arms (p11 / p10 / b) carry fixed effects only;
# car_proper only (intrinsic icar = #71). Occupancy fields are more weakly
# identified than count fields (one binary site per node).
.tobs_fit_fp_occu_nuts_spatial <- function(model, spatial, sigma.beta = 10,
                                           n.iter = 1000L, n.warmup = 1000L,
                                           n.chains = 1L, max.treedepth = 10L,
                                           adapt.delta = 0.9, seed = 1L,
                                           verbose = FALSE) {
  .tobs_reject_weighted_spatial(spatial, "fp_occu NUTS occupancy spatial")
  if (!identical(spatial$type, "car_proper"))
    stop(sprintf(paste0("fp_occu() NUTS + areal spatial supports the proper-CAR ",
                        "field car_proper(); the intrinsic '%s' field needs a ",
                        "sum-to-zero reparameterisation for NUTS -- use method = ",
                        "\"nested_laplace\" for the icar()/bym2() areal fit. ",
                        "(tulpaObs#72)"), spatial$type), call. = FALSE)
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites)
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for fp_occu NUTS.",
                 spatial$n_units, n_sites), call. = FALSE)
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  adj <- as.matrix(spatial$graph)

  # Warm coefficients + fixed field hyper (tau, rho) from the nested-Laplace fit.
  nl <- .tobs_fit_fp_occu_spatial(model, spatial, max_iter = 200L, tol = 1e-8,
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
  inv_metric <- c(rep(0.5, n_base), rep(1, n_sites))

  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = X_psi, X_p11 = X_p11, X_p10 = X_p10, X_b = X_b,
               n_sites = n_sites, n_field_units = n_sites,
               field_map = seq_len(n_sites), field_Linv = Linv)

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
           paste0("raw_", seq_len(n_sites)))
  colnames(draws) <- nms
  b_idx <- seq_len(n_base)
  par <- colMeans(draws); names(par) <- nms
  cov <- stats::cov(draws[, b_idx, drop = FALSE])
  raw_idx <- n_base + seq_len(n_sites)
  z_mean <- as.numeric(Linv %*% colMeans(draws[, raw_idx, drop = FALSE]))

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
  fit$method <- "nuts"; fit$spatial_field <- z_mean
  fit$nuts <- list(accept_prob = accept, divergent = divergent,
                   treedepth = as.integer(unlist(lapply(chains, `[[`, "treedepth"))),
                   epsilon = chains[[1L]]$epsilon, n_chains = as.integer(n.chains),
                   divergent_total = sum(divergent), tau = tau, rho = rho,
                   prior_type = spatial$type, fixed_hyper = TRUE)
  fit
}
