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

.tobs_fit_distance_spatial <- function(model, spatial, mixture = "poisson",
                                       K_max = NULL, max_iter = 200L, tol = 1e-6,
                                       verbose = TRUE, integration = "grid") {
  .tobs_reject_weighted_spatial(spatial, "distance abundance spatial")
  if (!identical(model$key, "halfnorm")) {
    stop("distance() areal spatial supports the half-normal key only; the ",
         "hazard-rate shape is a global coordinate not yet wired into the ",
         "spatial path. (tulpaObs#51)", call. = FALSE)
  }
  map <- seq_len(model$n_sites)
  field <- .tobs_areal_field_spec(spatial, model$n_sites, "distance", map)

  X_lam <- model$X_processes[[1]]; X_sig <- model$X_processes[[2]]
  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig)
  y <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  cutpoints <- as.numeric(model$cutpoints)
  transect_code <- .dist_transect_code(model$transect)
  quad_order <- as.integer(model$quad_order %||% 64L)
  R_max <- if (length(y)) max(rowSums(y)) else 0L
  K_max <- if (is.null(K_max)) as.integer(3L * R_max + 100L) else as.integer(K_max)
  is_nb <- mixture %in% c("negbin", "NB")

  i_lam <- seq_len(p_lam); i_sig <- p_lam + seq_len(p_sig)
  i_logr <- if (is_nb) p_lam + p_sig + 1L else NA_integer_
  n_fixed <- p_lam + p_sig + if (is_nb) 1L else 0L

  eval <- function(theta_fix, offset) {
    eta_lam <- as.numeric(X_lam %*% theta_fix[i_lam]) + offset
    eta_sig <- as.numeric(X_sig %*% theta_fix[i_sig])
    rr <- if (is_nb) exp(theta_fix[i_logr]) else Inf
    sw <- cpp_distance_site_sweep(y, eta_lam, eta_sig, cutpoints, transect_code,
                                  quad_order, K_max, nb = is_nb, r = rr)
    g <- numeric(n_fixed)
    g[i_lam] <- as.numeric(crossprod(X_lam, sw$grad_lam))
    g[i_sig] <- as.numeric(crossprod(X_sig, sw$grad_sig))
    if (is_nb) g[i_logr] <- sw$grad_logr
    list(log_lik = sum(sw$log_lik), grad_fixed = g, grad_eta = sw$grad_lam)
  }

  # Warm start from the non-spatial distance Laplace fit.
  warm <- tryCatch(
    .tobs_distance_re_warm(model, mixture = if (is_nb) "NB" else "P",
                           K_max = K_max, max_iter = max_iter, tol = tol),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm)) {
    th <- c(warm$beta_lambda, warm$beta_p)
    if (is_nb) th <- c(th, log(if (is.finite(warm$r %||% NA_real_)) warm$r else 2))
    th
  } else {
    th <- c(log(max(mean(rowSums(y)), 0.5) + 0.5), rep(0, p_lam - 1L),
            log(stats::median(cutpoints[-1])), rep(0, p_sig - 1L))
    if (is_nb) th <- c(th, log(2))
    th
  }
  if (length(theta0_fix) != n_fixed) theta0_fix <- numeric(n_fixed)

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol, label = "distance-spatial",
                              integration = integration)
  if (!isTRUE(res$ok))
    stop("distance() areal spatial fit produced no usable grid point.", call. = FALSE)

  means <- res$beta_mean
  raw <- list(
    mixture = if (is_nb) "negbin" else "poisson",
    beta_lambda = means[i_lam], beta_sigma = means[i_sig],
    log_r = if (is_nb) means[i_logr] else NA_real_,
    r = if (is_nb) exp(means[i_logr]) else NA_real_,
    vcov = res$vcov, log_lik = res$log_lik, converged = TRUE,
    key = model$key, transect = model$transect, hazard = FALSE, K_max = K_max)
  fit <- build_distance_fit(raw, model)
  fit$method <- "nested_laplace"
  fit$spatial_field <- res$field_mean
  fit$spatial_hyper <- res$hyper
  fit$spatial_integration <- res$integration
  fit$spatial_pareto_k <- res$pareto_k
  fit
}

# Areal-spatial binned distance sampling via NUTS (gcol33/tulpaObs#72): a FIXED-
# HYPER non-centered PROPER-CAR field on the abundance arm (log lambda) of the
# bin-multinomial distance marginal. The field precision (tau, rho) is fixed at the
# nested-Laplace areal posterior mean (fit$spatial_hyper) and the whitened raw ~
# N(0, I) (z = Linv %*% raw) is sampled jointly with the coefficients via the
# distance NUTS field block (cpp_distance_nuts over nuts_field_block.h). The areal
# Laplace fit supplies warm coefficients + the field hyper; NUTS then samples the
# exact coefficient + whitened-field posterior. Half-normal key only (the spatial
# Laplace path, gcol33/tulpaObs#79); car_proper only (intrinsic icar = #71).
# Poisson or NB.
.tobs_fit_distance_nuts_spatial <- function(model, spatial, mixture = "poisson",
                                            K_max = NULL, sigma.beta = 10,
                                            sigma.shape = 1.5, sigma.logr = 1.5,
                                            n.iter = 1000L, n.warmup = 1000L,
                                            n.chains = 1L, max.treedepth = 10L,
                                            adapt.delta = 0.9, seed = 1L,
                                            verbose = FALSE) {
  .tobs_reject_weighted_spatial(spatial, "distance NUTS abundance spatial")
  if (!identical(model$key, "halfnorm"))
    stop("distance() NUTS + areal spatial supports the half-normal key only ",
         "(the hazard-rate shape is not wired into the spatial path). ",
         "(tulpaObs#72/#79)", call. = FALSE)
  if (!identical(spatial$type, "car_proper"))
    stop(sprintf(paste0("distance() NUTS + areal spatial supports the proper-CAR ",
                        "field car_proper(); the intrinsic '%s' field needs a ",
                        "sum-to-zero reparameterisation for NUTS -- use method = ",
                        "\"nested_laplace\" for the icar()/bym2() areal fit. ",
                        "(tulpaObs#72)"), spatial$type), call. = FALSE)
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites)
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for distance NUTS.",
                 spatial$n_units, n_sites), call. = FALSE)
  is_nb <- mixture %in% c("negbin", "NB")
  mix_code <- if (is_nb) "NB" else "P"
  X_lambda <- model$X_processes[[1]]; X_sigma <- model$X_processes[[2]]
  y <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  p_lam <- ncol(X_lambda); p_sig <- ncol(X_sigma)
  R_max <- if (length(y)) max(rowSums(y)) else 0L
  K_max <- if (is.null(K_max)) as.integer(3L * R_max + 100L) else as.integer(K_max)
  adj <- as.matrix(spatial$graph)

  # Warm coefficients + fixed field hyper (tau, rho) from the nested-Laplace fit.
  nl <- .tobs_fit_distance_spatial(model, spatial, mixture = mixture, K_max = K_max,
                                   max_iter = 200L, tol = 1e-6, verbose = FALSE,
                                   integration = "grid")
  hyper <- nl$spatial_hyper
  tau <- max(unname(hyper[["tau"]]), 1e-3)
  rho <- min(max(unname(hyper[["rho"]]), 0.01), 0.99)
  Linv <- .tobs_field_linv(adj, tau, rho, n_sites)

  cm <- nl$means
  beta0 <- c(unname(cm[paste0("lambda_", model$process_info[[1]]$coef_names)]),
             unname(cm[paste0("sigma_",  model$process_info[[2]]$coef_names)]))
  if (is_nb) beta0 <- c(beta0, log(if (is.finite(nl$r %||% NA_real_)) nl$r else 2))
  n_base <- p_lam + p_sig + if (is_nb) 1L else 0L
  # Warm-start raw from the Laplace field mean via the precision Cholesky:
  # z = Linv %*% raw -> raw = L %*% z (raw0 = forward-solve), so the sampler starts
  # near the integrated field.
  L <- chol(.areal_Q(adj, rho) * tau + diag(1e-4 * tau, n_sites))   # upper L'L = Qr
  raw0 <- as.numeric(L %*% (nl$spatial_field %||% numeric(n_sites)))
  theta0 <- c(beta0, raw0)
  inv_metric <- c(rep(0.1, n_base), rep(1, n_sites))

  spec <- list(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
               cutpoints = as.numeric(model$cutpoints),
               transect = .dist_transect_code(model$transect),
               key = .dist_key_code(model$key), K_max = K_max, is_nb = is_nb,
               quad_order = as.integer(model$quad_order %||% 64L),
               n_field_units = n_sites, field_map = seq_len(n_sites),
               field_Linv = Linv)

  run_chain <- function(ch)
    cpp_distance_nuts(spec, theta0 = theta0, sigma_beta = sigma.beta,
                      sigma_shape = sigma.shape, sigma_logr = sigma.logr,
                      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
                      verbose = isTRUE(verbose))
  chains <- lapply(seq_len(as.integer(n.chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("sigma_",  model$process_info[[2]]$coef_names),
           if (is_nb) "log_r", paste0("raw_", seq_len(n_sites)))
  colnames(draws) <- nms
  b_idx <- seq_len(n_base)
  par <- colMeans(draws); cov <- stats::cov(draws[, b_idx, drop = FALSE])
  raw_idx <- n_base + seq_len(n_sites)
  z_mean <- as.numeric(Linv %*% colMeans(draws[, raw_idx, drop = FALSE]))

  lay <- .tobs_distance_nuts_layout(p_lam, p_sig, FALSE, is_nb)
  marg <- .tobs_distance_nuts_marginal(model, mixture = mix_code, K_max = K_max)
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$sigma], eta_b = 0,
                            r = if (is_nb) exp(par[lay$log_r]) else Inf)$log_lik
  raw_fit <- list(
    mixture = mix_code, key = model$key, transect = model$transect,
    hazard = FALSE, nb = is_nb,
    beta_lambda = unname(par[lay$lambda]), beta_sigma = unname(par[lay$sigma]),
    eta_b = NA_real_, shape = NA_real_,
    log_r = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max)
  fit <- build_distance_fit(raw_fit, model)
  fit$draws <- draws[, b_idx, drop = FALSE]
  fit$means <- par[b_idx]; fit$sds <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms[b_idx]
  fit$vcov <- cov
  fit$n_samples <- nrow(draws); fit$log_prob <- rep(ll_mean, nrow(draws))
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
