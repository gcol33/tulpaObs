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

.tobs_fit_dyn_abun_spatial <- function(model, spatial, mixture = "poisson",
                                       K_max = NULL, max_iter = 300L, tol = 1e-8,
                                       verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "dyn_abun abundance spatial")
  map <- seq_len(model$n_sites)
  field <- .tobs_areal_field_spec(spatial, model$n_sites, "dyn_abun", map)

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
                              max_iter = max_iter, tol = tol, label = "dyn-abun-spatial")
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
  fit
}
