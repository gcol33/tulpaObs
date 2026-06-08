# =============================================================================
# areal_bfgs.R - shared nested-Laplace driver for areal fields on a gradient-only
# family marginal (gcol33/tulpaObs#51).
#
# Families whose marginal exposes an analytic gradient but NO analytic per-site
# Hessian (the open N-mixture dyn_abun, the false-positive occupancy fp_occu)
# carry an areal field on one arm by: an outer grid over the field hyperparameters;
# per cell, BFGS over (fixed params, field params) with the family's analytic
# gradient + the field log-prior; the Laplace marginal from a finite-difference
# Hessian of that gradient at the mode (the observed-information route the
# families' non-spatial fits use).
#
# The driver owns ALL field handling through a `field` spec, so a new field
# structure (ICAR / proper-CAR single block, or BYM2 two-block v/w) is a new spec
# rather than a copied driver. The family supplies only
#   eval(theta_fixed, offset) -> list(log_lik, grad_fixed, grad_eta)
# where `offset` is the per-site eta offset on the field arm and `grad_eta` the
# per-site d log L / d eta_field. The field spec maps its parameters to that
# per-site offset and scatters grad_eta back, owns the prior, the sum-to-zero
# constraint, and the reported field. Single source of truth for the areal-BFGS
# families.
# =============================================================================

.areal_Q <- function(adj, rho) { deg <- rowSums(adj != 0); diag(deg) - rho * adj }

# ICAR / proper-CAR single-block field (parameter z, length n_sp; eta += z[map]).
.areal_field_car <- function(adj, kind, map, n_sp) {
  tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))
  rho_grid <- if (kind == "car_proper") seq(0.1, 0.95, length.out = 6L) else 1.0
  cells <- list()
  for (rho in rho_grid) {
    ldQ <- if (kind == "icar") 0 else {
      ch <- tryCatch(chol(.areal_Q(adj, rho)), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
    Q <- .areal_Q(adj, rho)
    for (tau in tau_grid) cells[[length(cells) + 1L]] <- list(tau = tau, rho = rho, ldQ = ldQ, Q = Q)
  }
  list(
    n_field = n_sp, n_sp = n_sp, cells = cells,
    valid = function(cell) is.finite(cell$ldQ),
    offset = function(fp, cell) fp[map],
    scatter = function(grad_eta) {
      g <- numeric(n_sp); for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]; g
    },
    prior_logp = function(fp, cell) {
      quad <- as.numeric(t(fp) %*% cell$Q %*% fp)
      if (kind == "icar") -0.5 * cell$tau * quad + 0.5 * (n_sp - 1) * log(cell$tau)
      else 0.5 * cell$ldQ + 0.5 * n_sp * log(cell$tau) - 0.5 * cell$tau * quad
    },
    prior_grad = function(fp, cell) cell$tau * as.numeric(cell$Q %*% fp),  # d(-logp)/dfp
    center = function(fp) if (kind == "icar") fp - mean(fp) else fp,
    constrain = if (kind == "icar") rep(TRUE, n_sp) else rep(FALSE, n_sp),
    to_phi = function(fp, cell) fp
  )
}

# BYM2 two-block field (v = ICAR, w = iid; eta += a v[map] + b w[map],
# a = sigma sqrt(rho/scale), b = sigma sqrt(1-rho)). The (v, w) priors are
# independent of (sigma, rho) (Riebler 2016), so the prior is constant across
# cells; the hyperparameters enter only the likelihood through (a, b).
.areal_field_bym2 <- function(adj, scale_factor, map, n_sp) {
  sigma_grid <- exp(seq(log(0.2), log(3), length.out = 5L))
  rho_grid   <- c(0.05, 0.3, 0.5, 0.7, 0.95)
  Q <- .areal_Q(adj, 1.0)                       # ICAR precision for v
  cells <- list()
  for (sg in sigma_grid) for (rho in rho_grid) {
    a <- sg * sqrt(rho / scale_factor); b <- sg * sqrt(1 - rho)
    cells[[length(cells) + 1L]] <- list(sigma = sg, rho = rho, a = a, b = b)
  }
  scat1 <- function(grad_eta) {
    g <- numeric(n_sp); for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]; g
  }
  list(
    n_field = 2L * n_sp, n_sp = n_sp, cells = cells,
    valid = function(cell) cell$sigma > 0 && cell$rho >= 0 && cell$rho <= 1,
    offset = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      cell$a * v[map] + cell$b * w[map]
    },
    scatter = function(grad_eta, cell) {
      s <- scat1(grad_eta); c(cell$a * s, cell$b * s)
    },
    prior_logp = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      -0.5 * as.numeric(t(v) %*% Q %*% v) - 0.5 * sum(w^2)
    },
    prior_grad = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      c(as.numeric(Q %*% v), w)
    },
    center = function(fp) { fp[seq_len(n_sp)] <- fp[seq_len(n_sp)] - mean(fp[seq_len(n_sp)]); fp },
    constrain = c(rep(TRUE, n_sp), rep(FALSE, n_sp)),  # sum-to-zero on v only
    to_phi = function(fp, cell) cell$a * fp[seq_len(n_sp)] + cell$b * fp[n_sp + seq_len(n_sp)],
    bym2 = TRUE
  )
}

.tobs_areal_bfgs_fit <- function(eval, n_fixed, field, theta0_fix,
                                 max_iter = 300L, tol = 1e-8, label = "areal-bfgs") {
  nf <- field$n_field; fi <- n_fixed + seq_len(nf); n_x <- n_fixed + nf
  is_bym2 <- isTRUE(field$bym2)
  cells <- field$cells; n_grid <- length(cells)
  modes <- matrix(NA_real_, n_grid, n_fixed); phis <- matrix(NA_real_, n_grid, field$n_sp)
  cov_blocks <- vector("list", n_grid); logm <- rep(-Inf, n_grid)

  prog <- tulpa:::.tulpa_iter_progress(label, n_grid, unit = "cells")
  for (k in seq_len(n_grid)) {
    cell <- cells[[k]]
    if (!field$valid(cell)) { prog$tick(); next }
    # gradient of (log_lik + log_prior) over theta = (fixed, field).
    grad_wp <- function(theta) {
      fp <- theta[fi]; e <- eval(theta[seq_len(n_fixed)], field$offset(fp, cell))
      gf <- if (is_bym2) field$scatter(e$grad_eta, cell) else field$scatter(e$grad_eta)
      c(e$grad_fixed, gf - field$prior_grad(fp, cell))
    }
    ll_fn <- function(theta) {
      fp <- theta[fi]; e <- eval(theta[seq_len(n_fixed)], field$offset(fp, cell))
      e$log_lik + field$prior_logp(fp, cell)
    }
    th0 <- c(theta0_fix, numeric(nf))
    opt <- tryCatch(stats::optim(th0, function(t) -ll_fn(t), function(t) -grad_wp(t),
                   method = "BFGS", control = list(maxit = as.integer(max_iter), reltol = tol)),
                   error = function(e) NULL)
    if (is.null(opt)) { prog$tick(); next }
    th <- opt$par; th[fi] <- field$center(th[fi])
    H <- tryCatch(-.fp_fd_jacobian(grad_wp, th), error = function(e) NULL)
    if (is.null(H)) { prog$tick(); next }
    H <- 0.5 * (H + t(H)); ridge <- max(1e-8 * mean(abs(diag(H))), 1e-10); diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL); if (is.null(ch)) { prog$tick(); next }
    ldH <- 2 * sum(log(diag(ch)))
    Hc <- H
    cc <- which(field$constrain)
    if (length(cc)) {                     # sum-to-zero penalty on the constrained field block
      pen <- 1e6 * mean(abs(diag(H))); idx <- n_fixed + cc
      Hc[idx, idx] <- Hc[idx, idx] + pen
    }
    cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
    if (is.null(cov_full)) { prog$tick(); next }
    modes[k, ] <- th[seq_len(n_fixed)]; phis[k, ] <- field$to_phi(th[fi], cell)
    cov_blocks[[k]] <- cov_full[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE]
    logm[k] <- ll_fn(th) - 0.5 * ldH
    prog$tick()
  }
  prog$finish()
  ok <- is.finite(logm)
  if (!any(ok)) return(list(ok = FALSE))
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data"); w[!ok] <- 0; w <- w / sum(w)
  beta_mean <- as.numeric(crossprod(w, ifelse(is.finite(modes), modes, 0)))
  field_mean <- as.numeric(crossprod(w, ifelse(is.finite(phis), phis, 0)))
  V <- matrix(0, n_fixed, n_fixed)
  for (k in which(ok)) {
    Vk <- cov_blocks[[k]]; if (is.null(Vk)) next
    dk <- modes[k, ] - beta_mean; V <- V + w[k] * (Vk + outer(dk, dk))
  }
  list(ok = TRUE, beta_mean = beta_mean, field_mean = field_mean, vcov = V,
       log_lik = sum(w * ifelse(ok, logm, 0)))
}

# Resolve the field spec for an areal-BFGS family from the spatial term.
.tobs_areal_field_spec <- function(spatial, n_sites, family, map) {
  if (spatial$type %in% c("spde", "gp", "multiscale_gp"))
    stop(sprintf(paste0("%s() areal spatial supports icar() / car_proper() / bym2() ",
                        "under method = \"nested_laplace\"; the '%s' field is not yet ",
                        "wired for %s. (tulpaObs#51)"), family, spatial$type, family),
         call. = FALSE)
  if (!spatial$type %in% c("icar", "car_proper", "bym2"))
    stop(sprintf("%s() areal spatial supports icar() / car_proper() / bym2(); got '%s'.",
                 family, spatial$type), call. = FALSE)
  if (spatial$n_units != n_sites)
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for %s.",
                 spatial$n_units, n_sites, family), call. = FALSE)
  adj <- if (!is.null(spatial$graph)) as.matrix(spatial$graph) else
    stop(sprintf("%s() spatial term must carry an adjacency graph.", family), call. = FALSE)
  if (identical(spatial$type, "bym2")) {
    sf <- spatial$scale_factor %||% compute_bym2_scale(spatial$graph)
    .areal_field_bym2(adj, sf, map, spatial$n_units)
  } else {
    .areal_field_car(adj, if (identical(spatial$type, "icar")) "icar" else "car_proper",
                     map, spatial$n_units)
  }
}
