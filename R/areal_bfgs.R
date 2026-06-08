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
  # One cell from physical hyperparameters: ICAR -> (tau); proper CAR -> (tau, rho).
  make_cell <- function(theta) {
    tau <- theta[1L]
    rho <- if (kind == "car_proper") theta[2L] else 1.0
    ldQ <- if (kind == "icar") 0 else {
      ch <- tryCatch(chol(.areal_Q(adj, rho)), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
    list(tau = tau, rho = rho, ldQ = ldQ, Q = .areal_Q(adj, rho))
  }
  cells <- list()
  for (rho in rho_grid) for (tau in tau_grid)
    cells[[length(cells) + 1L]] <- make_cell(if (kind == "car_proper") c(tau, rho) else tau)
  axes <- c(list(.tobs_ccd_axis("tau", "log", lower = 0.3, upper = 30, start = 3)),
            if (kind == "car_proper")
              list(.tobs_ccd_axis("rho", "identity", lower = 0.1, upper = 0.95, start = 0.5)))
  list(
    n_field = n_sp, n_sp = n_sp, cells = cells, axes = axes, make_cell = make_cell,
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
  make_cell <- function(theta) {
    sg <- theta[1L]; rho <- theta[2L]
    list(sigma = sg, rho = rho,
         a = sg * sqrt(rho / scale_factor), b = sg * sqrt(1 - rho))
  }
  cells <- list()
  for (sg in sigma_grid) for (rho in rho_grid)
    cells[[length(cells) + 1L]] <- make_cell(c(sg, rho))
  axes <- list(
    .tobs_ccd_axis("sigma", "log",      lower = 0.2,  upper = 3,    start = 0.77),
    .tobs_ccd_axis("rho",   "identity", lower = 0.05, upper = 0.95, start = 0.5))
  scat1 <- function(grad_eta) {
    g <- numeric(n_sp); for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]; g
  }
  list(
    n_field = 2L * n_sp, n_sp = n_sp, cells = cells, axes = axes, make_cell = make_cell,
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
                                 max_iter = 300L, tol = 1e-8, label = "areal-bfgs",
                                 integration = c("auto", "ccd", "grid")) {
  integration <- match.arg(integration)
  nf <- field$n_field; fi <- n_fixed + seq_len(nf)
  is_bym2 <- isTRUE(field$bym2)

  # One Laplace fit at a fixed field-hyperparameter cell. Returns the joint
  # marginal of (fixed params, field) plus the fixed-parameter mode / cov block /
  # reported field, or NULL on an invalid or numerically failed cell.
  solve_cell <- function(cell) {
    if (!field$valid(cell)) return(NULL)
    grad_wp <- function(theta) {              # grad of (log_lik + log_prior) over (fixed, field)
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
    if (is.null(opt)) return(NULL)
    th <- opt$par; th[fi] <- field$center(th[fi])
    H <- tryCatch(-.fp_fd_jacobian(grad_wp, th), error = function(e) NULL)
    if (is.null(H)) return(NULL)
    H <- 0.5 * (H + t(H)); ridge <- max(1e-8 * mean(abs(diag(H))), 1e-10); diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL); if (is.null(ch)) return(NULL)
    ldH <- 2 * sum(log(diag(ch)))
    Hc <- H
    cc <- which(field$constrain)
    if (length(cc)) {                     # sum-to-zero penalty on the constrained field block
      pen <- 1e6 * mean(abs(diag(H))); idx <- n_fixed + cc
      Hc[idx, idx] <- Hc[idx, idx] + pen
    }
    cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
    if (is.null(cov_full)) return(NULL)
    list(logm = ll_fn(th) - 0.5 * ldH,
         mode = th[seq_len(n_fixed)],
         phi  = field$to_phi(th[fi], cell),
         cov  = cov_full[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE])
  }

  # Weighted (marginalised) summaries from a set of evaluated cells + weights.
  summarise <- function(res, w, method, pareto_k = NA_real_) {
    ok  <- vapply(res, Negate(is.null), TRUE) & is.finite(w) & w > 0
    if (!any(ok)) return(list(ok = FALSE))
    w[!ok] <- 0; w <- w / sum(w)
    modes <- t(vapply(which(ok), function(k) res[[k]]$mode, numeric(n_fixed)))
    phis  <- t(vapply(which(ok), function(k) res[[k]]$phi,  numeric(field$n_sp)))
    wk    <- w[ok]
    beta_mean  <- as.numeric(crossprod(wk, modes))
    field_mean <- as.numeric(crossprod(wk, phis))
    V <- matrix(0, n_fixed, n_fixed)
    for (j in seq_along(wk)) {
      dk <- modes[j, ] - beta_mean
      V <- V + wk[j] * (res[[which(ok)[j]]]$cov + outer(dk, dk))
    }
    logm <- vapply(which(ok), function(k) res[[k]]$logm, numeric(1))
    list(ok = TRUE, beta_mean = beta_mean, field_mean = field_mean, vcov = V,
         log_lik = sum(wk * logm), integration = method, pareto_k = pareto_k)
  }

  # ---- outer integration: opt-in mode-centred CCD (gcol33/tulpaObs#60), silently
  # declining to the fixed tensor grid when the outer curvature is ill-conditioned.
  if (identical(integration, "ccd") && !is.null(field$axes)) {
    eval_logm <- function(theta_phys) {
      r <- solve_cell(field$make_cell(theta_phys))
      if (is.null(r) || !is.finite(r$logm)) NA_real_ else r$logm
    }
    cc <- tryCatch(.tobs_ccd_outer_grid(eval_logm, field$axes),
                   error = function(e) NULL)
    if (!is.null(cc)) {
      nn  <- nrow(cc$nodes)
      prog <- tulpa:::.tulpa_iter_progress(label, nn, unit = "cells")
      res <- vector("list", nn); logm <- rep(-Inf, nn)
      for (k in seq_len(nn)) {
        r <- solve_cell(field$make_cell(cc$nodes[k, ]))
        if (!is.null(r) && is.finite(r$logm)) { res[[k]] <- r; logm[k] <- r$logm }
        prog$tick()
      }
      prog$finish()
      if (any(is.finite(logm))) {
        w <- cc$dnode * exp(logm - max(logm[is.finite(logm)]))
        out <- summarise(res, w, "ccd", cc$pareto_k)
        if (isTRUE(out$ok)) return(out)
      }
    }
  }

  cells <- field$cells; n_grid <- length(cells)
  prog <- tulpa:::.tulpa_iter_progress(label, n_grid, unit = "cells")
  res <- vector("list", n_grid); logm <- rep(-Inf, n_grid)
  for (k in seq_len(n_grid)) {
    r <- solve_cell(cells[[k]])
    if (!is.null(r) && is.finite(r$logm)) { res[[k]] <- r; logm[k] <- r$logm }
    prog$tick()
  }
  prog$finish()
  if (!any(is.finite(logm))) return(list(ok = FALSE))
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data")
  summarise(res, w, "grid")
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
