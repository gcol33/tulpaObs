# =============================================================================
# sla_cover_hurdle_joint.R -- Simplified-Laplace correction on the joint
# nested-Laplace path for the cover hurdle family.
#
# The joint engine `tulpa_nested_laplace_joint(..., store_Q = TRUE)` returns
# per-grid joint modes `modes[k, ]` over `x = (beta_occ, beta_pos, w_field)`
# and per-grid sparse joint precision `Q_k`. We need:
#
#   * per-grid inner skewness gamma_kj via FD on the joint inner log-lik;
#   * mixture third-moment combiner across the outer grid (variance-of-means
#     contribution from per-grid mode shifts + skewness contribution from
#     gamma_kj scaled by sigma_kj^3).
#
# Per-grid inner variance is reused from `.joint_inner_var()` (already
# implements the sum-to-zero constraint correction on phi / theta blocks);
# the constraint-corrected Sigma-column solve here applies the *same*
# correction to the off-diagonal column we step along.
# =============================================================================


# ---------------------------------------------------------------------------
# Per-unit unit-variance field from the joint mode at grid point k.
#
# Under the (sigma_occ, sigma_pos) reparam (gcol33/tulpa#18) the stored
# modes are unit-variance latent z. Each arm scales it by its own
# `sigma_arm_k`. For BYM2 the z-field is
#     z[s] = sqrt(rho) * scale_factor * phi[s] + sqrt(1 - rho) * theta[s];
# for ICAR / CAR_proper it is just `phi[s]`.
#
# Inputs:
#   x_lat -- length-`n_x` latent vector with the same block layout as
#            `fit$modes[k, ]`. Allows callers to substitute displaced
#            betas while keeping the field block at the mode.
#   k     -- grid index.
#   fit   -- joint fit (from `tulpa_nested_laplace_joint`).
#
# Returns: numeric length `n_s` (number of spatial units).
# ---------------------------------------------------------------------------
.joint_field_unit_z <- function(x_lat, k, fit) {
  layout <- fit$arm_layout
  theta_grid <- fit$theta_grid
  phi_start   <- layout$phi_start
  theta_start <- layout$theta_start

  # n_s: spatial-block length, inferred from layout (phi block always present
  # under any of bym2 / icar / car_proper). For BYM2 we have phi + theta of
  # equal length; for ICAR / CAR_proper just phi.
  n_s_phi <- if (!is.null(theta_start)) {
    as.integer(theta_start - phi_start)
  } else {
    as.integer(layout$n_x - phi_start)
  }
  if (is.null(phi_start) || n_s_phi <= 0L) {
    return(numeric(0))
  }
  phi_k <- x_lat[phi_start + seq_len(n_s_phi)]

  has_rho <- "rho" %in% colnames(theta_grid)
  if (!is.null(theta_start) && has_rho) {
    n_s_theta    <- as.integer(layout$n_x - theta_start)
    theta_k      <- x_lat[theta_start + seq_len(n_s_theta)]
    rho_k        <- theta_grid[k, "rho"]
    scale_factor <- as.numeric(attr(fit, "scale_factor") %||% 1.0)
    sqrt(max(rho_k, 0) + 1e-10) * scale_factor * phi_k +
      sqrt(max(1 - rho_k, 0) + 1e-10) * theta_k
  } else {
    phi_k
  }
}


# ---------------------------------------------------------------------------
# Joint inner log-lik at fixed grid point k (no field prior; cancels under
# FD since the prior is Gaussian -> third derivative is zero).
#
# x is a length-n_x latent vector (same layout as fit$modes[k, ]). The
# log-lik sums the occurrence-arm Bernoulli + the positive-arm Beta or
# Gaussian (lognormal on log y -- the Jacobian is beta-independent and
# also cancels under FD along the beta direction).
# ---------------------------------------------------------------------------
.loglik_cover_joint_at_grid <- function(x, k, fit, enc, positive) {
  L <- fit$arm_layout
  p_occ <- L$p[1]
  p_pos <- L$p[2]
  beta_occ <- x[L$beta_start[1] + seq_len(p_occ)]
  beta_pos <- x[L$beta_start[2] + seq_len(p_pos)]
  field_z  <- .joint_field_unit_z(x, k, fit)

  spi_full <- enc$..spi_full
  spi_pos  <- enc$..spi_pos
  s_occ_k  <- fit$theta_grid[k, "sigma_occ"]
  s_pos_k  <- fit$theta_grid[k, "sigma_pos"]

  eta_occ <- as.numeric(enc$occ_data$X %*% beta_occ)
  if (length(field_z) > 0L && length(spi_full) > 0L) {
    eta_occ <- eta_occ + s_occ_k * field_z[spi_full]
  }
  eta_pos <- as.numeric(enc$pos_data$X %*% beta_pos)
  if (length(field_z) > 0L && length(spi_pos) > 0L) {
    eta_pos <- eta_pos + s_pos_k * field_z[spi_pos]
  }

  # Occurrence Bernoulli (logit link), log1p-stable.
  y_occ <- enc$occ_data$y
  ll_occ <- sum(y_occ * eta_occ - log1p(exp(eta_occ)))

  # Positive arm. Beta on natural scale, or Gaussian on log y (Jacobian
  # `-log y` is beta-independent; we drop it).
  ll_pos <- if (identical(positive, "beta")) {
    y_pos <- enc$pos_data$y
    if (length(y_pos) == 0L) {
      0
    } else {
      phi_k <- if ("phi_pos" %in% colnames(fit$theta_grid)) {
        fit$theta_grid[k, "phi_pos"]
      } else {
        # Profiled beta path (legacy): phi held fixed at the prefit value.
        as.numeric(attr(fit, "phi_fixed") %||% 1.0)
      }
      mu_pos <- plogis(eta_pos)
      a <- mu_pos * phi_k
      b <- (1 - mu_pos) * phi_k
      sum(lgamma(phi_k) - lgamma(a) - lgamma(b) +
            (a - 1) * log(y_pos) + (b - 1) * log1p(-y_pos))
    }
  } else {
    # Lognormal: enc$pos_data$y is already log(cover). Noise SD lives on
    # the per-grid `phi_pos` axis (gaussian arm's phi is the residual SD);
    # legacy `sigma_pos_noise` name is still recognized for forward compat,
    # with a final fallback to the posterior-mean `sigma_pos_fixed` attr.
    z <- enc$pos_data$y
    if (length(z) == 0L) {
      0
    } else {
      sig_k <- if ("phi_pos" %in% colnames(fit$theta_grid)) {
        fit$theta_grid[k, "phi_pos"]
      } else if ("sigma_pos_noise" %in% colnames(fit$theta_grid)) {
        fit$theta_grid[k, "sigma_pos_noise"]
      } else {
        as.numeric(attr(fit, "sigma_pos_fixed") %||% 1.0)
      }
      sum(stats::dnorm(z, eta_pos, sig_k, log = TRUE))
    }
  }
  ll_occ + ll_pos
}


# ---------------------------------------------------------------------------
# Constraint-corrected Sigma columns at beta-block indices.
#
# Returns Sigma_k[, beta_idx_arm] (n_x x p_arm) with the same sum-to-zero
# constraint correction as `.joint_inner_var()`:
#   Sigma_c = Q^{-1} - Q^{-1} A^T (A Q^{-1} A^T)^{-1} A Q^{-1}
# where A picks the per-block sums of phi (and theta, for BYM2).
#
# Reuses the constraint-row construction from `.joint_inner_var` (cf.
# R/family_cover_hurdle.R:1163-1234) so the FD steps along the *same*
# directions as the variance estimator.
# ---------------------------------------------------------------------------
.joint_constrained_solve_columns <- function(Qk, layout, beta_idx_arm) {
  n_x <- as.integer(layout$n_x)
  p_arm <- length(beta_idx_arm)
  if (p_arm == 0L) return(matrix(0, n_x, 0))

  A_cols <- list()
  if (!is.null(layout$phi_start)) {
    n_s_phi <- (layout$theta_start %||% layout$n_x) - layout$phi_start
    A_cols[[length(A_cols) + 1L]] <- layout$phi_start + seq_len(n_s_phi)
  }
  if (!is.null(layout$theta_start)) {
    n_s_theta <- layout$n_x - layout$theta_start
    A_cols[[length(A_cols) + 1L]] <- layout$theta_start + seq_len(n_s_theta)
  }
  k_constr <- length(A_cols)

  E <- Matrix::sparseMatrix(
    i = beta_idx_arm, j = seq_len(p_arm), x = 1,
    dims = c(n_x, p_arm)
  )
  V <- tryCatch(as.matrix(Matrix::solve(Qk, E)), error = function(e) NULL)
  if (is.null(V)) return(NULL)

  if (k_constr == 0L) return(V)

  ii <- unlist(A_cols)
  jj <- rep(seq_len(k_constr), vapply(A_cols, length, integer(1)))
  A_t <- Matrix::sparseMatrix(i = ii, j = jj, x = 1,
                              dims = c(n_x, k_constr))

  W <- tryCatch(as.matrix(Matrix::solve(Qk, A_t)), error = function(e) NULL)
  if (is.null(W)) return(V)
  AV <- as.matrix(Matrix::crossprod(A_t, V))           # k_constr x p_arm
  M  <- as.matrix(Matrix::crossprod(A_t, W))           # k_constr x k_constr
  M_inv_AV <- tryCatch(solve(M, AV), error = function(e) NULL)
  if (is.null(M_inv_AV)) return(V)
  V - W %*% M_inv_AV
}


# ---------------------------------------------------------------------------
# Per-grid SLA gamma via 5-point central FD along the constrained Sigma
# columns.
#
# Step size derivation. Parameterize the FD path by the standardised
# variable `s = h * sigma_j`, so one unit of `s` is one marginal SD along
# v_j (`d^3 L/ds^3 |_{s=0}` is exactly `gamma_j`). The 5-point central
# rule's optimal step in `s` is `s* = eps^(1/5)`, giving
#     h_j = eps^(1/5) / sigma_j.
#
# Equivalently: pin the beta_j-component displacement at the natural
# scale of the j-th marginal posterior. Along direction `v_j = Sigma[, j]`
# we have `v_j[beta_j] = sigma_j^2`, so the beta_j-coordinate moves by
# `h * v_j[beta_j] = (eps^(1/5)/sigma_j) * sigma_j^2 = eps^(1/5) * sigma_j`
# -- a fraction of one marginal SD, well-conditioned for FD.
#
# Why not `h = eps^(1/5) * sigma_j / ||v_j||` (the standalone formula).
# When `v_j` is a column of a *beta-only* Sigma, `||v_j|| ~ sigma_j^2` and
# the two formulas agree. When `v_j` is a column of the *full* joint
# Sigma, `||v_j||` is inflated by field-block components (proportional
# to `sigma_pos_k * field_z_scale`), making `h` orders of magnitude
# smaller than optimal. Roundoff in the Beta `lgamma` evaluations then
# dominates the FD signal -- see dev_notes/probe_sla_fd_step.R for an
# h-sweep on the failing grid point.
# ---------------------------------------------------------------------------
.sla_inner_gamma_joint <- function(k, beta_idx_arm, fit, enc, positive) {
  layout <- fit$arm_layout
  n_x <- as.integer(layout$n_x)
  p_arm <- length(beta_idx_arm)
  if (p_arm == 0L) return(numeric(0))

  Qp <- fit$Q_csc_p_per_grid
  Qi <- fit$Q_csc_i_per_grid
  Qx <- fit$Q_csc_x_per_grid
  if (is.null(Qp) || is.null(Qi) || is.null(Qx)) {
    return(rep(NA_real_, p_arm))
  }
  if (is.null(Qp[[k]]) || length(Qx[[k]]) == 0L) {
    return(rep(NA_real_, p_arm))
  }

  Qk_lt <- Matrix::sparseMatrix(
    i = as.integer(Qi[[k]]) + 1L,
    p = as.integer(Qp[[k]]),
    x = as.numeric(Qx[[k]]),
    dims = c(n_x, n_x),
    symmetric = FALSE,
    index1 = TRUE
  )
  Qk <- Matrix::forceSymmetric(Qk_lt, uplo = "L")

  V <- .joint_constrained_solve_columns(Qk, layout, beta_idx_arm)
  if (is.null(V)) return(rep(NA_real_, p_arm))

  beta_hat <- fit$modes[k, ]
  eps_h <- .Machine$double.eps^(1 / 5)

  gamma <- numeric(p_arm)
  for (j in seq_len(p_arm)) {
    v_j <- V[, j]
    sigma2_j <- max(v_j[beta_idx_arm[j]], 0)
    sigma_j  <- sqrt(sigma2_j)
    if (sigma_j <= 0) { gamma[j] <- 0; next }
    if (sum(v_j^2) <= 0) { gamma[j] <- 0; next }
    # Standardised step: `s = h * sigma_j` is one marginal SD; pick
    # `s = eps^(1/5)` for the 5-point central rule => `h = eps^(1/5) / sigma_j`.
    # See header comment for the field-coupling rationale.
    h <- eps_h / sigma_j
    L_p2 <- .loglik_cover_joint_at_grid(beta_hat + 2 * h * v_j, k, fit, enc, positive)
    L_p1 <- .loglik_cover_joint_at_grid(beta_hat +     h * v_j, k, fit, enc, positive)
    L_m1 <- .loglik_cover_joint_at_grid(beta_hat -     h * v_j, k, fit, enc, positive)
    L_m2 <- .loglik_cover_joint_at_grid(beta_hat - 2 * h * v_j, k, fit, enc, positive)
    d3 <- (L_p2 - 2 * L_p1 + 2 * L_m1 - L_m2) / (2 * h^3)
    gamma[j] <- d3 / sigma_j^3
  }
  gamma
}


# ---------------------------------------------------------------------------
# Mixture third-moment combiner.
#
#   M3_j = Sum_k w_k * [ (mu_kj - mu_j)^3
#                      + 3 * (mu_kj - mu_j) * sigma_kj^2
#                      + gamma_kj * sigma_kj^3 ]
#   gamma_j = M3_j / sigma_j^3
#
# Inputs:
#   modes_j     -- length n_grid, per-grid posterior mean of beta_j.
#   weights     -- length n_grid, outer-grid posterior weights (sum to 1).
#   sigma2_kj   -- length n_grid, per-grid inner variance (from
#                  `.joint_inner_var()`).
#   gamma_kj    -- length n_grid, per-grid standardised inner skewness.
#   mu_j        -- scalar marginal mean (already computed).
#   sigma2_j    -- scalar marginal variance (already computed).
# ---------------------------------------------------------------------------
.sla_combine_grid_skewness <- function(modes_j, weights, sigma2_kj,
                                       gamma_kj, mu_j, sigma2_j) {
  dmu <- modes_j - mu_j
  sigma_kj <- sqrt(pmax(sigma2_kj, 0))
  # Replace non-finite gammas (e.g. failed FD) with zero so they neither
  # contribute nor propagate NaN through the mixture sum.
  g <- gamma_kj
  g[!is.finite(g)] <- 0
  M3 <- sum(weights * (dmu^3 + 3 * dmu * sigma2_kj +
                       g * sigma_kj^3))
  sigma_j <- sqrt(max(sigma2_j, 0))
  if (sigma_j <= 0) return(0)
  M3 / sigma_j^3
}


# ---------------------------------------------------------------------------
# Orchestrator: per-coefficient marginal SLA gamma for both arms.
#
# Returns `list(gamma_occ, gamma_pos, valid, reason)` shaped identically
# to `.sla_compute_cover_hurdle()` (standalone-Laplace path), so the
# existing `.sla_build_cover_hurdle_draws()` consumer can stay unchanged.
# ---------------------------------------------------------------------------
.sla_compute_cover_hurdle_joint <- function(fit, enc, positive) {
  layout <- fit$arm_layout
  if (is.null(layout) || is.null(fit$Q_csc_p_per_grid)) {
    return(list(gamma_occ = NULL, gamma_pos = NULL, valid = FALSE,
                reason = "joint fit missing arm_layout / store_Q"))
  }
  p_occ <- layout$p[1]
  p_pos <- layout$p[2]
  bocc_idx <- layout$beta_start[1] + seq_len(p_occ)
  bpos_idx <- layout$beta_start[2] + seq_len(p_pos)

  weights <- fit$weights
  n_grid  <- length(weights)
  if (n_grid == 0L) {
    return(list(gamma_occ = NULL, gamma_pos = NULL, valid = FALSE,
                reason = "joint fit has no grid points"))
  }

  # Per-grid inner variance (constraint-corrected). Returns
  # n_grid x (p_occ + p_pos); first `p_occ` cols are occ, rest pos.
  beta_idx_all <- c(bocc_idx, bpos_idx)
  inner_var <- .joint_inner_var(fit, beta_idx_all)
  if (is.null(inner_var)) {
    return(list(gamma_occ = NULL, gamma_pos = NULL, valid = FALSE,
                reason = "could not invert joint Q (no inner variance)"))
  }

  gamma_grid_occ <- matrix(NA_real_, n_grid, p_occ)
  gamma_grid_pos <- matrix(NA_real_, n_grid, p_pos)
  for (k in seq_len(n_grid)) {
    g_occ <- tryCatch(
      .sla_inner_gamma_joint(k, bocc_idx, fit, enc, positive),
      error = function(e) rep(NA_real_, p_occ)
    )
    g_pos <- tryCatch(
      .sla_inner_gamma_joint(k, bpos_idx, fit, enc, positive),
      error = function(e) rep(NA_real_, p_pos)
    )
    if (length(g_occ) == p_occ) gamma_grid_occ[k, ] <- g_occ
    if (length(g_pos) == p_pos) gamma_grid_pos[k, ] <- g_pos
  }

  # Posterior-weighted marginal means / variances (var-of-means +
  # mean-of-var). Reuses the same arithmetic as the SE path in
  # fit_cover_hurdle_joint_nested.
  modes_occ <- fit$modes[, bocc_idx, drop = FALSE]
  modes_pos <- fit$modes[, bpos_idx, drop = FALSE]
  mu_occ <- as.numeric(crossprod(weights, modes_occ))
  mu_pos <- as.numeric(crossprod(weights, modes_pos))
  vom_occ <- as.numeric(crossprod(weights, modes_occ^2)) - mu_occ^2
  vom_pos <- as.numeric(crossprod(weights, modes_pos^2)) - mu_pos^2
  mov_occ <- as.numeric(crossprod(weights, inner_var[, seq_len(p_occ), drop = FALSE]))
  mov_pos <- as.numeric(crossprod(weights,
    inner_var[, p_occ + seq_len(p_pos), drop = FALSE]))
  sigma2_occ <- pmax(vom_occ + mov_occ, 0)
  sigma2_pos <- pmax(vom_pos + mov_pos, 0)

  gamma_occ <- vapply(seq_len(p_occ), function(j) {
    .sla_combine_grid_skewness(
      modes_j   = modes_occ[, j],
      weights   = weights,
      sigma2_kj = inner_var[, j],
      gamma_kj  = gamma_grid_occ[, j],
      mu_j      = mu_occ[j],
      sigma2_j  = sigma2_occ[j]
    )
  }, numeric(1))
  gamma_pos <- vapply(seq_len(p_pos), function(j) {
    .sla_combine_grid_skewness(
      modes_j   = modes_pos[, j],
      weights   = weights,
      sigma2_kj = inner_var[, p_occ + j],
      gamma_kj  = gamma_grid_pos[, j],
      mu_j      = mu_pos[j],
      sigma2_j  = sigma2_pos[j]
    )
  }, numeric(1))

  names(gamma_occ) <- colnames(enc$occ_data$X)
  names(gamma_pos) <- colnames(enc$pos_data$X)

  if (!all(is.finite(gamma_occ)) || !all(is.finite(gamma_pos))) {
    return(list(gamma_occ = gamma_occ, gamma_pos = gamma_pos, valid = FALSE,
                reason = "non-finite marginal gamma on at least one coefficient"))
  }
  list(gamma_occ = gamma_occ, gamma_pos = gamma_pos, valid = TRUE,
       reason = "ok")
}
