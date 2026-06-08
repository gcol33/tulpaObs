# =============================================================================
# dyn_abun_spatial.R - areal-spatial Dail-Madsen open N-mixture (tulpaObs#51)
#
# An ICAR or proper-CAR field on the INITIAL-abundance arm (log lambda_1) of the
# open N-mixture, fit by nested Laplace: outer grid over the field precision tau
# (+ rho for proper CAR), inner optimisation of (beta_lambda, beta_p, beta_omega,
# beta_gamma, z[, log_r]) by BFGS with the exact analytic gradient the forward-HMM
# kernel already returns (cpp_dyn_abun_total_log_lik), plus the CAR log-prior. The
# Laplace marginal at each grid point uses a finite-difference Hessian of that
# analytic gradient at the mode -- the same observed-information route the non-
# spatial dyn_abun fit uses (R/dyn_abun.R). The field z loads onto eta_lambda
# exactly like the initial-abundance intercept (one unit per site).
#
#   .tobs_fit_dyn_abun_spatial()   dispatch from .tobs_fit_model (icar / car_proper)
# =============================================================================

.dyn_abun_spatial_Q <- function(adj, rho) {
  deg <- rowSums(adj != 0)
  diag(deg) - rho * adj
}

.tobs_fit_dyn_abun_spatial <- function(model, spatial, mixture = "poisson",
                                       K_max = NULL, max_iter = 300L, tol = 1e-8,
                                       verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "dyn_abun abundance spatial")
  if (spatial$type %in% c("bym2", "spde", "gp", "multiscale_gp")) {
    stop(sprintf(paste0("dyn_abun() areal spatial supports icar() / car_proper() ",
                        "under method = \"nested_laplace\"; the '%s' field is not ",
                        "yet wired for dyn_abun. (tulpaObs#51)"), spatial$type),
         call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "car_proper")) {
    stop(sprintf("dyn_abun() areal spatial supports icar() / car_proper(); got '%s'.",
                 spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for dyn_abun.",
                 spatial$n_units, n_sites), call. = FALSE)
  }
  adj <- if (!is.null(spatial$graph)) as.matrix(spatial$graph) else
    stop("dyn_abun() spatial term must carry an adjacency graph.", call. = FALSE)
  n_sp <- spatial$n_units; map <- seq_len(n_sites)

  X_lam <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_om  <- model$X_processes[[3]]; X_gm <- model$X_processes[[4]]
  p_lam <- ncol(X_lam); p_p <- ncol(X_p); p_om <- ncol(X_om); p_gm <- ncol(X_gm)
  y_flat <- as.integer(model$y_flat); N <- model$n_sites
  T <- model$n_seasons; J <- model$max_visits; K <- model$K_max
  is_nb <- mixture %in% c("negbin", "NB")

  off <- cumsum(c(0L, p_lam, p_p, p_om, p_gm))
  i_lam <- off[1] + seq_len(p_lam); i_p <- off[2] + seq_len(p_p)
  i_om  <- off[3] + seq_len(p_om);  i_gm <- off[4] + seq_len(p_gm)
  p_fix <- off[5]; i_logr <- if (is_nb) p_fix + 1L else NA_integer_
  p_fixed <- p_fix + if (is_nb) 1L else 0L
  z_idx <- p_fixed + seq_len(n_sp)
  n_x <- p_fixed + n_sp

  ev <- function(theta) {
    eta_lam <- as.numeric(X_lam %*% theta[i_lam]) + theta[z_idx][map]
    cpp_dyn_abun_total_log_lik(
      y_flat, N, T, J, K, eta_lam,
      as.numeric(X_p %*% theta[i_p]), as.numeric(X_om %*% theta[i_om]),
      as.numeric(X_gm %*% theta[i_gm]),
      use_nb = is_nb, eta_logr = if (is_nb) theta[i_logr] else 0.0)
  }
  grad_design <- function(out, theta, tau, Q) {
    g <- numeric(n_x)
    g[i_lam] <- as.numeric(crossprod(X_lam, out$grad_eta_lambda))
    g[i_p]   <- as.numeric(crossprod(X_p,   out$grad_eta_p))
    g[i_om]  <- as.numeric(crossprod(X_om,  out$grad_eta_omega))
    g[i_gm]  <- as.numeric(crossprod(X_gm,  out$grad_eta_gamma))
    if (is_nb) g[i_logr] <- as.numeric(out$grad_eta_logr)
    gz <- numeric(n_sp)
    for (s in seq_len(N)) gz[map[s]] <- gz[map[s]] + out$grad_eta_lambda[s]
    g[z_idx] <- gz - tau * as.numeric(Q %*% theta[z_idx])
    g
  }
  car_logprior <- function(z, tau, rho, kind, log_det_Q, Q) {
    quad <- as.numeric(t(z) %*% Q %*% z)
    if (kind == "icar") -0.5 * tau * quad + 0.5 * (n_sp - 1) * log(tau)
    else 0.5 * log_det_Q + 0.5 * n_sp * log(tau) - 0.5 * tau * quad
  }

  warm <- tryCatch(dyn_abun_laplace(
    y_flat = y_flat, n_sites = N, T = T, J = J, K_max = K,
    X_lambda = X_lam, X_p = X_p, X_omega = X_om, X_gamma = X_gm,
    mixture = if (is_nb) "negbin" else "poisson", verbose = FALSE),
    error = function(e) NULL)
  theta0_fix <- if (!is.null(warm)) warm$means else numeric(p_fixed)
  if (length(theta0_fix) != p_fixed) theta0_fix <- numeric(p_fixed)

  tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))
  kind <- if (identical(spatial$type, "icar")) "icar" else "car_proper"
  rho_grid <- if (kind == "car_proper") seq(0.1, 0.95, length.out = 6L) else 1.0
  log_det_Q <- vapply(rho_grid, function(rho) {
    if (kind == "icar") 0 else {
      ch <- tryCatch(chol(.dyn_abun_spatial_Q(adj, rho)), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
  }, numeric(1)); names(log_det_Q) <- as.character(rho_grid)

  grid <- expand.grid(tau = tau_grid, rho = rho_grid)
  n_grid <- nrow(grid)
  modes <- matrix(NA_real_, n_grid, p_fixed); fields <- matrix(NA_real_, n_grid, n_sp)
  cov_blocks <- vector("list", n_grid); logm <- rep(-Inf, n_grid)

  .prog <- tulpa:::.tulpa_iter_progress("dyn-abun-spatial", n_grid, unit = "cells")
  for (k in seq_len(n_grid)) {
    tau <- grid$tau[k]; rho <- grid$rho[k]; ldQ <- log_det_Q[[as.character(rho)]]
    if (!is.finite(ldQ)) { .prog$tick(); next }
    Q <- .dyn_abun_spatial_Q(adj, rho)
    neg_lp <- function(theta) {
      o <- ev(theta)
      -(o$log_lik + car_logprior(theta[z_idx], tau, rho, kind, ldQ, Q))
    }
    neg_gr <- function(theta) -grad_design(ev(theta), theta, tau, Q)
    th0 <- c(theta0_fix, numeric(n_sp))
    opt <- tryCatch(stats::optim(th0, neg_lp, neg_gr, method = "BFGS",
                   control = list(maxit = as.integer(max_iter), reltol = tol)),
                   error = function(e) NULL)
    if (is.null(opt)) { .prog$tick(); next }
    th <- opt$par
    if (kind == "icar") th[z_idx] <- th[z_idx] - mean(th[z_idx])  # sum-to-zero
    out <- ev(th)
    # Observed-info Hessian: FD-Jacobian of the (loglik + prior) gradient at mode.
    H <- tryCatch(-.fp_fd_jacobian(function(t) grad_design(ev(t), t, tau, Q), th),
                  error = function(e) NULL)
    if (is.null(H)) { .prog$tick(); next }
    H <- 0.5 * (H + t(H))
    ridge <- max(1e-8 * mean(abs(diag(H))), 1e-10); diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(ch)) { .prog$tick(); next }
    log_det_H <- 2 * sum(log(diag(ch)))
    Hc <- H
    if (kind == "icar") {
      pen <- 1e6 * mean(abs(diag(H)))
      Hc[z_idx, z_idx] <- Hc[z_idx, z_idx] + pen
    }
    cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
    if (is.null(cov_full)) { .prog$tick(); next }
    modes[k, ] <- th[seq_len(p_fixed)]; fields[k, ] <- th[z_idx]
    cov_blocks[[k]] <- cov_full[seq_len(p_fixed), seq_len(p_fixed), drop = FALSE]
    logm[k] <- out$log_lik + car_logprior(th[z_idx], tau, rho, kind, ldQ, Q) -
               0.5 * log_det_H
    .prog$tick()
  }
  .prog$finish()
  ok <- is.finite(logm)
  if (!any(ok)) stop("dyn_abun() areal spatial fit produced no usable grid point.", call. = FALSE)
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data"); w[!ok] <- 0; w <- w / sum(w)

  beta_mean <- as.numeric(crossprod(w, ifelse(is.finite(modes), modes, 0)))
  field_mean <- as.numeric(crossprod(w, ifelse(is.finite(fields), fields, 0)))
  nm <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
          paste0("p_",      model$process_info[[2]]$coef_names),
          paste0("omega_",  model$process_info[[3]]$coef_names),
          paste0("gamma_",  model$process_info[[4]]$coef_names))
  if (is_nb) nm <- c(nm, "log_r")
  names(beta_mean) <- nm
  V <- matrix(0, p_fixed, p_fixed)
  for (k in which(ok)) {
    Vk <- cov_blocks[[k]]; if (is.null(Vk)) next
    dk <- modes[k, ] - beta_mean; V <- V + w[k] * (Vk + outer(dk, dk))
  }
  dimnames(V) <- list(nm, nm)

  raw <- list(
    beta_lambda = beta_mean[i_lam], beta_p = beta_mean[i_p],
    beta_omega = beta_mean[i_om], beta_gamma = beta_mean[i_gm],
    log_r = if (is_nb) beta_mean[i_logr] else NA_real_,
    r = if (is_nb) exp(beta_mean[i_logr]) else NA_real_,
    mixture = if (is_nb) "negbin" else "poisson",
    means = beta_mean, vcov = V, log_lik = sum(w * ifelse(ok, logm, 0)),
    mean_N1 = NULL, K_max = K, converged = TRUE, n_iter = NA_integer_,
    coef_names = nm)
  fit <- build_dyn_abun_fit(raw, model)
  fit$method <- "nested_laplace"
  fit$spatial_field <- field_mean
  fit
}
