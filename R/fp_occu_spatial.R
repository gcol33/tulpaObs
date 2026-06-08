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

.tobs_fit_fp_occu_spatial <- function(model, spatial, max_iter = 200L,
                                      tol = 1e-8, verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "fp_occu occupancy spatial")
  map <- seq_len(model$n_sites)
  field <- .tobs_areal_field_spec(spatial, model$n_sites, "fp_occu", map)

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
                              max_iter = max_iter, tol = tol, label = "fp-occu-spatial")
  if (!isTRUE(res$ok))
    stop("fp_occu() areal spatial fit produced no usable grid point.", call. = FALSE)

  nm <- c(paste0("psi_", model$process_info[[1]]$coef_names),
          paste0("p11_", model$process_info[[2]]$coef_names),
          paste0("p10_", model$process_info[[3]]$coef_names),
          paste0("b_",   model$process_info[[4]]$coef_names))
  means <- res$beta_mean; names(means) <- nm
  V <- res$vcov; dimnames(V) <- list(nm, nm)
  # Posterior occupancy w1 at the integrated estimate + field (for fitted()).
  cl <- function(e) pmin(pmax(e, -30), 30)
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
  fit
}
