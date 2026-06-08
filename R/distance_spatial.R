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
                                       verbose = TRUE) {
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
                              max_iter = max_iter, tol = tol, label = "distance-spatial")
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
  fit
}
