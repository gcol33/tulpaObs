# =============================================================================
# areal_bfgs.R - shared nested-Laplace driver for areal fields on a gradient-only
# family marginal (gcol33/tulpaObs#51).
#
# Families whose marginal exposes an analytic gradient but NO analytic per-site
# Hessian (the open N-mixture dyn_abun, the false-positive occupancy fp_occu)
# carry an areal ICAR / proper-CAR field on one arm by: an outer grid over the
# field precision tau (+ rho for proper CAR); per cell, BFGS over (fixed params,
# field z) with the family's analytic gradient + the CAR log-prior; the Laplace
# marginal from a finite-difference Hessian of that gradient at the mode (the same
# observed-information route the families' non-spatial fits use). The driver owns
# the grid walk, the prior, the FD-Hessian + constrained covariance, and the grid
# integration; the family supplies `eval(theta) -> list(log_lik, grad)` over
# theta = (fixed[1:n_fixed], z[n_fixed + 1:n_sp]) with the field scatter folded in
# and NO prior. Single source of truth for the two areal-BFGS families.
# =============================================================================

.areal_Q <- function(adj, rho) { deg <- rowSums(adj != 0); diag(deg) - rho * adj }

.tobs_areal_bfgs_fit <- function(eval, n_fixed, n_sp, adj, kind,
                                 tau_grid, rho_grid, theta0_fix,
                                 max_iter = 300L, tol = 1e-8, label = "areal-bfgs") {
  z_idx <- n_fixed + seq_len(n_sp)
  car_lp <- function(z, tau, rho, ldQ, Q) {
    quad <- as.numeric(t(z) %*% Q %*% z)
    if (kind == "icar") -0.5 * tau * quad + 0.5 * (n_sp - 1) * log(tau)
    else 0.5 * ldQ + 0.5 * n_sp * log(tau) - 0.5 * tau * quad
  }
  log_det_Q <- vapply(rho_grid, function(rho) {
    if (kind == "icar") 0 else {
      ch <- tryCatch(chol(.areal_Q(adj, rho)), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
  }, numeric(1)); names(log_det_Q) <- as.character(rho_grid)

  grid <- expand.grid(tau = tau_grid, rho = rho_grid); n_grid <- nrow(grid)
  modes <- matrix(NA_real_, n_grid, n_fixed); fields <- matrix(NA_real_, n_grid, n_sp)
  cov_blocks <- vector("list", n_grid); logm <- rep(-Inf, n_grid)

  prog <- tulpa:::.tulpa_iter_progress(label, n_grid, unit = "cells")
  for (k in seq_len(n_grid)) {
    tau <- grid$tau[k]; rho <- grid$rho[k]; ldQ <- log_det_Q[[as.character(rho)]]
    if (!is.finite(ldQ)) { prog$tick(); next }
    Q <- .areal_Q(adj, rho)
    neg_lp <- function(theta) {
      o <- eval(theta); -(o$log_lik + car_lp(theta[z_idx], tau, rho, ldQ, Q))
    }
    neg_gr <- function(theta) {
      g <- -eval(theta)$grad
      g[z_idx] <- g[z_idx] + tau * as.numeric(Q %*% theta[z_idx]); g
    }
    th0 <- c(theta0_fix, numeric(n_sp))
    opt <- tryCatch(stats::optim(th0, neg_lp, neg_gr, method = "BFGS",
                   control = list(maxit = as.integer(max_iter), reltol = tol)),
                   error = function(e) NULL)
    if (is.null(opt)) { prog$tick(); next }
    th <- opt$par
    if (kind == "icar") th[z_idx] <- th[z_idx] - mean(th[z_idx])
    o <- eval(th)
    grad_wp <- function(t) {
      g <- eval(t)$grad; g[z_idx] <- g[z_idx] - tau * as.numeric(Q %*% t[z_idx]); g
    }
    H <- tryCatch(-.fp_fd_jacobian(grad_wp, th), error = function(e) NULL)
    if (is.null(H)) { prog$tick(); next }
    H <- 0.5 * (H + t(H)); ridge <- max(1e-8 * mean(abs(diag(H))), 1e-10)
    diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(ch)) { prog$tick(); next }
    ldH <- 2 * sum(log(diag(ch)))
    Hc <- H
    if (kind == "icar") { pen <- 1e6 * mean(abs(diag(H))); Hc[z_idx, z_idx] <- Hc[z_idx, z_idx] + pen }
    cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
    if (is.null(cov_full)) { prog$tick(); next }
    modes[k, ] <- th[seq_len(n_fixed)]; fields[k, ] <- th[z_idx]
    cov_blocks[[k]] <- cov_full[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE]
    logm[k] <- o$log_lik + car_lp(th[z_idx], tau, rho, ldQ, Q) - 0.5 * ldH
    prog$tick()
  }
  prog$finish()
  ok <- is.finite(logm)
  if (!any(ok)) return(list(ok = FALSE))
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data"); w[!ok] <- 0; w <- w / sum(w)
  beta_mean <- as.numeric(crossprod(w, ifelse(is.finite(modes), modes, 0)))
  field_mean <- as.numeric(crossprod(w, ifelse(is.finite(fields), fields, 0)))
  V <- matrix(0, n_fixed, n_fixed)
  for (k in which(ok)) {
    Vk <- cov_blocks[[k]]; if (is.null(Vk)) next
    dk <- modes[k, ] - beta_mean; V <- V + w[k] * (Vk + outer(dk, dk))
  }
  list(ok = TRUE, beta_mean = beta_mean, field_mean = field_mean, vcov = V,
       log_lik = sum(w * ifelse(ok, logm, 0)))
}

# Resolve the standard areal grid axes + adjacency for the areal-BFGS families.
.tobs_areal_setup <- function(spatial, n_sites, family) {
  if (spatial$type %in% c("bym2", "spde", "gp", "multiscale_gp"))
    stop(sprintf(paste0("%s() areal spatial supports icar() / car_proper() under ",
                        "method = \"nested_laplace\"; the '%s' field is not yet ",
                        "wired for %s. (tulpaObs#51)"), family, spatial$type, family),
         call. = FALSE)
  if (!spatial$type %in% c("icar", "car_proper"))
    stop(sprintf("%s() areal spatial supports icar() / car_proper(); got '%s'.",
                 family, spatial$type), call. = FALSE)
  if (spatial$n_units != n_sites)
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for %s.",
                 spatial$n_units, n_sites, family), call. = FALSE)
  adj <- if (!is.null(spatial$graph)) as.matrix(spatial$graph) else
    stop(sprintf("%s() spatial term must carry an adjacency graph.", family), call. = FALSE)
  kind <- if (identical(spatial$type, "icar")) "icar" else "car_proper"
  list(adj = adj, n_sp = spatial$n_units, kind = kind,
       tau_grid = exp(seq(log(0.3), log(30), length.out = 9L)),
       rho_grid = if (kind == "car_proper") seq(0.1, 0.95, length.out = 6L) else 1.0)
}
