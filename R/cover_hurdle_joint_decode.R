
# Pre-fit the lognormal residual SD on the positive subset before handing
# control to the joint engine. The joint integrand reads `phi` as the noise
# SD; without a sensible pre-fit it sees scale 1 regardless of truth and the
# log-marginal across the alpha grid becomes near-flat (issue #4).
#
# Strategy: non-spatial Gaussian fit on the positive subset, residual SD as
# the point estimate. This is an upper bound on the true noise SD (it
# includes the alpha-mediated field variance), but it sits inside the same
# order of magnitude as the truth, which is enough to restore the joint
# engine's discrimination across the alpha grid. The post-hoc sigma_pos in
# `fit_cover_hurdle_joint_nested` then refines this by subtracting the
# alpha-scaled posterior field.
.prefit_lognormal_sigma <- function(enc, control) {
  y <- enc$pos_data$y
  X <- enc$pos_data$X
  n <- length(y); p <- ncol(X)
  if (n <= p) return(1.0)
  beta_init  <- tryCatch(qr.solve(X, y), error = function(e) NULL)
  if (is.null(beta_init)) return(1.0)
  resid_init <- as.numeric(y - X %*% beta_init)
  sigma_init <- sqrt(sum(resid_init^2) / max(n - p, 1L))
  if (!is.finite(sigma_init) || sigma_init <= 0) return(1.0)
  sigma_init
}

# Positive-arm family + dispersion (phi) grid for the joint nested-Laplace cover
# paths -- single source of truth for the shared-field, coupled-trend, MCAR, and
# arm-specific fitters. `phi` is the gaussian residual SD (lognormal), the latent
# log-cover SD (ordinal interval-censored Gaussian, integrated on the same phi
# axis as the lognormal sigma -- its discrete-class sibling), or the beta
# precision. The ordinal grid runs a touch wider on the low side because the
# midpoint-residual prefit over-disperses the censored latent.
.cover_pos_family_grid <- function(positive, enc, control) {
  if (positive == "lognormal") {
    sigma_hat <- .prefit_lognormal_sigma(enc, control)
    list(pos_family = "gaussian", phi_hat = sigma_hat,
         phi_grid_pos = control$phi.grid %||%
           exp(seq(log(sigma_hat / 3), log(sigma_hat * 3), length.out = 7)))
  } else if (positive == "lognormal_trunc") {
    # Upper-truncated Gaussian on log-cover (cover <= 1). The midpoint prefit SD
    # ignores truncation and so over-disperses a touch; the grid runs a little
    # wider on the low side, as for the ordinal censored latent.
    sigma_hat <- .prefit_lognormal_sigma(enc, control)
    list(pos_family = "truncated_gaussian", phi_hat = sigma_hat,
         phi_grid_pos = control$phi.grid %||%
           exp(seq(log(sigma_hat / 4), log(sigma_hat * 3), length.out = 7)))
  } else if (positive == "ordinal") {
    sigma_hat <- .prefit_lognormal_sigma(enc, control)
    list(pos_family = "interval_gaussian", phi_hat = sigma_hat,
         phi_grid_pos = control$phi.grid %||%
           exp(seq(log(sigma_hat / 4), log(sigma_hat * 3), length.out = 7)))
  } else {
    list(pos_family = "beta", phi_hat = 1.0,
         phi_grid_pos = control$phi.grid %||%
           exp(seq(log(2), log(300), length.out = 7)))
  }
}

# Append the per-plot positive-arm bound(s) for families that carry them: the
# interval (lower, upper] for ordinal, the truncation ceiling for lognormal_trunc.
# A no-op for beta / lognormal. The bounds in `enc$pos_data` already run over the
# positive subset, matching the arm rows.
.cover_arm_pos_bounds <- function(arm_pos, enc, positive) {
  if (identical(positive, "ordinal")) {
    arm_pos$lower <- as.numeric(enc$pos_data$lower)
    arm_pos$upper <- as.numeric(enc$pos_data$upper)
  } else if (identical(positive, "lognormal_trunc")) {
    arm_pos$trunc_upper <- as.numeric(enc$pos_data$trunc_upper)
  }
  arm_pos
}

#' Decode a joint-nested-Laplace cover-hurdle fit into a `cover_fit`.
#'
#' Lighter-weight than the single-Laplace decode: the joint engine has
#' already produced posterior moments for beta and the spatial hyperparameters,
#' so we just shape them into the existing `cover_fit` structure.
#'
#' Under an SLA method (`method = "nested_laplace_sla"`), the per-arm marginal
#' skewness is
#' computed via `.sla_compute_cover_hurdle_joint()` (mixture third-moment
#' over the outer grid; per-grid FD of the joint inner log-lik along the
#' constraint-corrected Sigma columns), and per-arm pseudo-draws are
#' resampled from moment-matched skew-normals via
#' [`.sla_build_cover_hurdle_draws()`].
#'
#' @keywords internal
decode_cover_hurdle_joint <- function(fits, enc, family,
                                      approx = "gaussian_laplace") {
  beta_occ <- fits$beta_occ
  beta_pos <- fits$beta_pos
  names(beta_occ) <- colnames(enc$occ_data$X)
  names(beta_pos) <- colnames(enc$pos_data$X)

  se_occ <- fits$se_occ
  se_pos <- fits$se_pos
  if (length(se_occ)) names(se_occ) <- names(beta_occ)
  if (length(se_pos)) names(se_pos) <- names(beta_pos)

  hyperpar <- list(
    spatial = fits$joint$theta_mean,
    engine  = "nested_laplace"
  )
  if (fits$positive %in% c("lognormal", "lognormal_trunc", "ordinal")) {
    hyperpar$sigma_pos    <- fits$sigma_pos
    hyperpar$sigma_pos_sd <- fits$sigma_pos_sd
  } else {
    hyperpar$phi_pos    <- fits$phi_pos
    hyperpar$phi_pos_sd <- fits$phi_pos_sd
  }

  # Simplified-Laplace marginal correction on the joint path. The SLA
  # orchestrator wants `enc$..spi_full` / `enc$..spi_pos` for per-arm field
  # gather inside the inner log-lik evaluator; we attach them here from
  # the fits list (computed once inside `fit_cover_hurdle_joint_nested`).
  skew_occ <- NULL
  skew_pos <- NULL
  draws_occ <- NULL
  draws_pos <- NULL
  sla_status <- "off"
  if (identical(approx, "simplified_laplace") && isTRUE(fits$mcar)) {
    # The simplified-Laplace marginal skew correction over a correlated MCAR
    # field is not wired (the per-arm field gather assumes a single-field copy).
    # Record the no-op status rather than mis-applying the single-field path
    # (gcol33/tulpaObs#64); the Gaussian-Laplace MCAR fit stands on its own.
    sla_status <- "mcar_unsupported"
  } else if (identical(approx, "simplified_laplace") && isTRUE(fits$armspecific)) {
    # The simplified-Laplace marginal skew correction over arm-specific separate
    # latents is not wired (the per-arm field gather assumes a single shared
    # copied field). Record the no-op; the Gaussian-Laplace fit stands on its own
    # (gcol33/tulpaObs#65).
    sla_status <- "armspecific_unsupported"
  } else if (identical(approx, "simplified_laplace") &&
             identical(fits$positive, "ordinal")) {
    # The simplified-Laplace marginal skew correction is not wired for the
    # interval-censored Gaussian arm (its per-arm FD assumes a continuous
    # point-response density). Record the no-op; the Gaussian-Laplace ordinal
    # fit stands on its own.
    sla_status <- "ordinal_unsupported"
  } else if (identical(approx, "simplified_laplace")) {
    enc_sla <- enc
    enc_sla$..spi_full <- as.integer(fits$spi_full %||% integer(0))
    enc_sla$..spi_pos  <- as.integer(fits$spi_pos  %||% integer(0))
    sla_res <- .sla_compute_cover_hurdle_joint(fits$joint, enc_sla,
                                               fits$positive)
    sla_draws <- .sla_build_cover_hurdle_draws(
      beta_occ, se_occ, beta_pos, se_pos, sla_res
    )
    draws_occ <- sla_draws$draws_occ
    draws_pos <- sla_draws$draws_pos
    sla_status <- sla_draws$sla_status
    if (isTRUE(sla_res$valid)) {
      skew_occ <- sla_res$gamma_occ
      skew_pos <- sla_res$gamma_pos
    } else {
      # The orchestrator may still return numeric (possibly non-finite)
      # gamma vectors alongside `valid = FALSE`; surface them only when
      # they are finite so downstream consumers can inspect them.
      if (!is.null(sla_res$gamma_occ) && all(is.finite(sla_res$gamma_occ))) {
        skew_occ <- sla_res$gamma_occ
      }
      if (!is.null(sla_res$gamma_pos) && all(is.finite(sla_res$gamma_pos))) {
        skew_pos <- sla_res$gamma_pos
      }
    }
  }

  out <- structure(
    list(
      occ          = fits$m_occ,
      pos          = fits$m_pos,
      beta_occ     = beta_occ,
      beta_pos     = beta_pos,
      se_occ       = se_occ,
      se_pos       = se_pos,
      positive     = fits$positive,
      sigma_pos    = fits$sigma_pos,
      sigma_pos_sd = fits$sigma_pos_sd,
      phi_pos      = fits$phi_pos,
      phi_pos_sd   = fits$phi_pos_sd,
      pi_one       = enc$oi$pi_one    %||% NA_real_,
      pi_one_sd    = enc$oi$pi_one_sd %||% NA_real_,
      n_ceiling    = enc$oi$n_ceiling %||% NA_integer_,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = enc$oi$n_positive %||% length(enc$idx_pos),
      converged    = TRUE,
      # Unified convergence record (gcol33/tulpaObs#88); see the non-spatial
      # assembly above. The joint nested-Laplace outer grid has no iteration
      # count, so n_iter is NA, matching the other joint_coupled paths.
      convergence  = list(converged = TRUE, n_iter = NA_integer_,
                          sla_status = sla_status),
      log_marginal = c(joint = max(fits$joint$log_marginal)),
      joint        = fits$joint,
      spi_full     = fits$spi_full,
      spi_pos      = fits$spi_pos,
      n_cells      = fits$n_cells,
      n_fields     = fits$n_fields,
      trend_weight  = fits$trend_weight,
      trend_weights = fits$trend_weights,
      trend_w_occ   = fits$trend_w_occ,
      trend_w_pos   = fits$trend_w_pos,
      sigma_trend   = fits$sigma_trend,
      alpha_trend   = fits$alpha_trend,
      mcar             = isTRUE(fits$mcar),
      mcar_field_names = fits$mcar_field_names,
      sigma_mcar       = fits$sigma_mcar,
      rho_mcar         = fits$rho_mcar,
      alpha_mcar       = fits$alpha_mcar,
      alpha_mcar_sd    = fits$alpha_mcar_sd,
      armspecific       = isTRUE(fits$armspecific),
      armspec_blocks    = fits$armspec_blocks,
      sigma_armspecific = fits$sigma_armspecific,
      skew_occ     = skew_occ,
      skew_pos     = skew_pos,
      draws_occ    = draws_occ,
      draws_pos    = draws_pos,
      sla_status   = sla_status
    ),
    class = c("cover_fit", "tobs_multiarm_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


# ---------------------------------------------------------------------------
# Multi-block prior assembly (Phase J-D)
# ---------------------------------------------------------------------------

# Build a multi-block joint prior (spatial + optional temporal + optional
# IID RE blocks) for tulpa::tulpa_nested_laplace_joint() under
# cover-hurdle copy semantics (copy on the spatial block).
#
# Non-spatial blocks are shared identically across the two arms — no
# per-arm scaling (INLA convention). This matches the typical cover
# hurdle use case: a year RE that influences both occurrence and cover
# magnitude in the same way, an observer RE that introduces a shared
# offset on both arms.
.cover_build_multi_prior <- function(prior_spatial, spi_full, spi_pos,
                                     idx_pos, temporal, re,
                                     control, sigma_pos_grid) {
  # Spatial block — fill missing grids with defaults and attach per-arm
  # spatial_idx vectors. (Single-block path stores spi inside the arms;
  # multi-block puts it in the block.)
  sp <- prior_spatial
  if (is.null(sp$sigma_grid)) {
    sp$sigma_grid <- exp(seq(log(0.1), log(3), length.out = 5))
  }
  if (tolower(sp$type) == "bym2" && is.null(sp$rho_grid)) {
    sp$rho_grid <- c(0.25, 0.5, 0.75)
  }
  sp$spatial_idx <- list(as.integer(spi_full), as.integer(spi_pos))

  blocks <- list(sp)

  if (!is.null(temporal)) {
    blocks[[length(blocks) + 1L]] <- .cover_temporal_block(
      temporal, idx_pos, control
    )
  }

  if (!is.null(re) && length(re) > 0L) {
    for (re_i in re) {
      blocks[[length(blocks) + 1L]] <- .cover_re_block(
        re_i, idx_pos, control
      )
    }
  }

  list(
    prior = blocks,
    copy  = list(block = 1L, arm = "pos",
                 sigma_pos_grid = as.numeric(sigma_pos_grid))
  )
}

.cover_temporal_block <- function(temporal, idx_pos, control) {
  if (!inherits(temporal, "tobs_temporal")) {
    stop("`temporal` must be a tobs_temporal() object.", call. = FALSE)
  }
  # tobs_temporal()'s `shared = c(TRUE, FALSE)` default was designed for
  # occupancy + detection (state vs. observation). cover() has two
  # likelihood arms (occurrence + cover magnitude) and the temporal term
  # enters both identically -- the `shared` field is ignored here.
  # Index codes resolved when the temporal() term was constructed (against
  # the same NA-dropped observations these arms use).
  t_full <- as.integer(temporal$time_idx)
  t_pos  <- t_full[idx_pos]
  n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times) else max(t_full)
  type <- temporal$type

  if (type == "ar1") {
    list(
      type         = "ar1",
      n_times      = as.integer(n_t),
      tau_grid     = as.numeric(control$tau.temporal.grid %||%
                                 c(1, 4, 16)),
      rho_grid     = as.numeric(control$rho.temporal.grid %||%
                                 c(0.3, 0.7)),
      temporal_idx = list(as.integer(t_full), as.integer(t_pos))
    )
  } else if (type == "iid") {
    list(
      type       = "iid",
      n_units    = as.integer(n_t),
      sigma_grid = as.numeric(control$sigma.temporal.grid %||%
                               exp(seq(log(0.1), log(1), length.out = 3))),
      obs_idx    = list(as.integer(t_full), as.integer(t_pos))
    )
  } else if (type %in% c("rw1", "rw2")) {
    list(
      type         = type,
      n_times      = as.integer(n_t),
      tau_grid     = as.numeric(control$tau.temporal.grid %||%
                                 c(1, 4, 16)),
      temporal_idx = list(as.integer(t_full), as.integer(t_pos))
    )
  } else {
    stop(sprintf("Unsupported tobs_temporal$type: '%s'", type),
         call. = FALSE)
  }
}

.cover_re_block <- function(re_i, idx_pos, control) {
  if (!inherits(re_i, "tobs_re")) {
    stop("`re` elements must be tobs_re() objects.", call. = FALSE)
  }
  if (!identical(re_i$type, "intercept") && !identical(re_i$type, "iid")) {
    stop("cover() multi-block: tobs_re(type = 'intercept' | 'iid') is the ",
         "only supported config. Random slopes land in a later phase.",
         call. = FALSE)
  }
  if (!identical(re_i$model, "iid")) {
    stop("cover() multi-block: tobs_re(model = 'iid') is the only ",
         "supported temporal structure on RE blocks. AR1/RW1/RW2 on RE ",
         "land in a later phase.", call. = FALSE)
  }
  # Same as in .cover_temporal_block: `shared` is ignored in cover-hurdle
  # context. The RE term enters both arms identically.
  # Group codes resolved when the re() term was constructed.
  g_full <- as.integer(re_i$group_idx)
  g_pos  <- g_full[idx_pos]
  n_g <- if (!is.null(re_i$n_groups)) as.integer(re_i$n_groups) else max(g_full)
  list(
    type       = "iid",
    n_units    = as.integer(n_g),
    sigma_grid = as.numeric(control$sigma.re.grid %||%
                             exp(seq(log(0.1), log(1.5), length.out = 3))),
    obs_idx    = list(as.integer(g_full), as.integer(g_pos))
  )
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Per-grid inner posterior variance for selected latent coordinates,
# applying sum-to-zero constraints on the BYM2/ICAR/CAR_proper spatial
# blocks (phi, theta) so the fixed-effect intercept is data-identified
# rather than prior-bounded.
#
# The joint-Laplace precision Q_k is near-singular along the
# (intercept, mean(phi)) direction whenever the prior on the spatial
# block has a sum-to-zero soft-null (ICAR is rank-deficient, BYM2's
# phi block likewise). The unconstrained inverse maps that direction
# onto the weak beta prior (1e-4 in the kernel = sd 100), producing
# meaningless intercept SEs. The fix is the standard INLA constraint
# correction:
#
#   Sigma_c = Q^{-1} - Q^{-1} A^T (A Q^{-1} A^T)^{-1} A Q^{-1}
#
# where A picks the per-block sums of phi (and theta, for BYM2). With
# A = 0 (no spatial block) this reduces to the unconstrained inverse,
# which is the right behaviour for the SPDE case (Q non-singular).
#
# IMPORTANT: the same constraint is applied to BOTH BYM2 sub-blocks
# (phi: rank-deficient ICAR; theta: proper IID). This is a modelling
# choice rather than a mathematical necessity for theta — the IID
# prior identifies mean(theta) at precision n_s — but matching INLA's
# `f(..., model = "bym2", constr = TRUE)` default keeps the reported
# intercept comparable. Any simulator generating BYM2 data for the
# joint engine must demean both phi_f and theta_f before scaling,
# otherwise the constrained-intercept estimator targets
# `beta_pos_0_truth + alpha * mean(w_s_sim)` rather than the
# population truth and coverage of the population truth collapses
# with alpha (see `simulate_cover_joint()` for a demeaned simulator;
# diagnosed in INLAabun `example/validation/SUMMARY.md` Demo 3).
#
# Returns an `n_grid x length(beta_idx)` matrix of constrained
# Var(beta_j | data, theta_k), or NULL when no Q matrices were stored.
# Build the sum-to-zero constraint columns for the joint-Laplace field
# block(s). Returns a list of 1-based column-index vectors, one all-ones
# constraint per structured spatial field block. Constraint columns index into
# the joint latent vector.
#
# Two layouts are handled:
#   * Multi-block layout (`field_starts` reported): one constraint per ICAR /
#     CAR_proper / BYM2 field, derived from the engine's per-block offsets and
#     sizes (BYM2 contributes the phi sub-block constraint plus the theta IID
#     constraint, matching INLA's bym2 constr = TRUE default).
#   * Single-block layout (`phi_start` / `theta_start`): the original
#     one-(ICAR)-or-two-(BYM2) constraint behaviour, kept byte-identical.
.joint_field_constraint_cols <- function(layout) {
  A_cols <- list()
  if (!is.null(layout$field_starts)) {
    starts <- layout$field_starts
    types  <- layout$field_block_types %||% rep("icar", length(starts))
    bstart <- layout$block_start
    bsize  <- layout$block_size
    for (i in seq_along(starts)) {
      s0   <- starts[i]
      type <- tolower(types[i])
      b    <- match(s0, bstart)
      sz   <- if (is.na(b)) NA_integer_ else bsize[b]
      if (type == "bym2") {
        n_units <- as.integer(sz / 2L)
        A_cols[[length(A_cols) + 1L]] <- s0 + seq_len(n_units)
        A_cols[[length(A_cols) + 1L]] <- s0 + n_units + seq_len(n_units)
      } else {
        n_units <- as.integer(sz)
        A_cols[[length(A_cols) + 1L]] <- s0 + seq_len(n_units)
      }
    }
    return(A_cols)
  }
  if (!is.null(layout$phi_start)) {
    n_s_phi <- (layout$theta_start %||% layout$n_x) - layout$phi_start
    A_cols[[length(A_cols) + 1L]] <- layout$phi_start + seq_len(n_s_phi)
  }
  if (!is.null(layout$theta_start)) {
    n_s_theta <- layout$n_x - layout$theta_start
    A_cols[[length(A_cols) + 1L]] <- layout$theta_start + seq_len(n_s_theta)
  }
  A_cols
}

.joint_inner_var <- function(fit, beta_idx) {
  Qp <- fit$Q_csc_p_per_grid
  Qi <- fit$Q_csc_i_per_grid
  Qx <- fit$Q_csc_x_per_grid
  n_x <- fit$Q_csc_n
  if (is.null(Qp) || is.null(Qi) || is.null(Qx) || is.null(n_x)) return(NULL)

  layout <- fit$arm_layout
  # Build constraint matrix A (k x n_x): one row of all-ones per structured
  # spatial field block. layout offsets are 0-based.
  A_cols <- .joint_field_constraint_cols(layout)
  k_constr <- length(A_cols)

  n_grid <- length(Qp)
  p <- length(beta_idx)
  out <- matrix(NA_real_, n_grid, p)

  E <- Matrix::sparseMatrix(
    i = beta_idx, j = seq_len(p), x = 1,
    dims = c(n_x, p)
  )
  A_t <- if (k_constr > 0L) {
    ii <- unlist(A_cols)
    jj <- rep(seq_len(k_constr), vapply(A_cols, length, integer(1)))
    Matrix::sparseMatrix(i = ii, j = jj, x = 1,
                         dims = c(n_x, k_constr))
  } else NULL

  for (k in seq_len(n_grid)) {
    if (is.null(Qp[[k]]) || length(Qx[[k]]) == 0L) next
    Qk_lt <- Matrix::sparseMatrix(
      i = as.integer(Qi[[k]]) + 1L,
      p = as.integer(Qp[[k]]),
      x = as.numeric(Qx[[k]]),
      dims = c(n_x, n_x),
      symmetric = FALSE,
      index1 = TRUE
    )
    Qk <- Matrix::forceSymmetric(Qk_lt, uplo = "L")
    V <- tryCatch(Matrix::solve(Qk, E), error = function(e) NULL)
    if (is.null(V)) next
    var_uncon <- vapply(seq_len(p),
      function(j) as.numeric(V[beta_idx[j], j]), numeric(1))

    if (k_constr > 0L) {
      W <- tryCatch(Matrix::solve(Qk, A_t), error = function(e) NULL)
      if (!is.null(W)) {
        AV <- as.matrix(Matrix::crossprod(A_t, V))     # k_constr x p
        M  <- as.matrix(Matrix::crossprod(A_t, W))     # k_constr x k_constr
        corr <- vapply(seq_len(p), function(j) {
          v <- AV[, j]
          as.numeric(crossprod(v, solve(M, v)))
        }, numeric(1))
        out[k, ] <- pmax(var_uncon - corr, 0)
        next
      }
    }
    out[k, ] <- pmax(var_uncon, 0)
  }
  out
}

.se_from_hessian <- function(H, scale = 1) {
  if (is.null(H)) return(numeric(0))
  cov <- tryCatch(scale * solve(H), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, nrow(H)))
  sqrt(pmax(diag(cov), 0))
}

# Per-grid constrained covariance block for the selected latent
# coordinates. Same conditioning-by-kriging constraint correction as
# `.joint_inner_var()` on the BYM2/ICAR/CAR_proper spatial blocks, returning
# the full `length(beta_idx) x length(beta_idx)` sub-block per grid cell. Used
# by the joint-engine autoscale unscaling so the intercept SE carries the
# correct cross-covariance contribution from the centered+scaled slopes
# (gcol33/tulpaObs#9).
#
# `beta_idx` stacks `n_dense` leading fixed-effect (betas) coordinates followed
# by the latent field coordinates. With a field present (`n_dense <
# length(beta_idx)`) the per-cell extraction takes the cheap selected-inversion
# recipe (gcol33/tulpa#113): the dense betas block and the betas x field cross
# are exact, the field marginal variances come from one Takahashi pass, and the
# field x field off-diagonal -- never read by the SD summary (it consumes the
# betas block + the diagonal) nor by predict (which draws each cell directly
# from `Q_k`) -- is left at zero. The cells run concurrently in the engine over
# `n_threads` (gcol33/tulpaObs#93). When `beta_idx` is betas-only
# (`n_dense == length(beta_idx)`, the cover()-only callers) the full block is
# formed. The whole loop is the single C++ source
# `tulpa:::cpp_joint_inner_vcov_blocks`, replacing the former serial R
# `solve(Qk, E)` over ~`length(beta_idx)` right-hand sides per cell.
#
# Returns a list of length n_grid, each element either NULL (when the per-cell
# sparse Cholesky failed or the cell stored no Q) or a dense `p x p` matrix.
# Returns NULL when no Q matrices were stored at all.
.joint_inner_vcov_block <- function(fit, beta_idx, n_dense = length(beta_idx),
                                    n_threads = 1L) {
  Qp <- fit$Q_csc_p_per_grid
  Qi <- fit$Q_csc_i_per_grid
  Qx <- fit$Q_csc_x_per_grid
  n_x <- fit$Q_csc_n
  if (is.null(Qp) || is.null(Qi) || is.null(Qx) || is.null(n_x)) return(NULL)

  A_cols <- .joint_field_constraint_cols(fit$arm_layout)
  field_marginal <- n_dense < length(beta_idx)
  nthr <- max(1L, as.integer(n_threads %||% 1L))

  tulpa:::cpp_joint_inner_vcov_blocks(
    Q_p_per_grid = Qp, Q_i_per_grid = Qi, Q_x_per_grid = Qx,
    n_x          = as.integer(n_x),
    idx          = as.integer(beta_idx),
    n_dense      = as.integer(n_dense),
    A_cols_list  = lapply(A_cols, as.integer),
    field_marginal = field_marginal,
    n_threads    = nthr
  )
}

.coef_table <- function(beta, se) {
  if (length(se) != length(beta)) se <- rep(NA_real_, length(beta))
  z <- beta / se
  data.frame(
    estimate = beta,
    std.err  = se,
    z.value  = z,
    row.names = names(beta)
  )
}

.extract_spatial_hyperpar <- function(fit, spec) {
  if (is.null(spec)) return(NULL)
  out <- list()
  for (nm in c("range", "sigma", "sigma2_gp", "phi_gp", "tau_spatial",
               "sigma_spatial", "rho")) {
    if (!is.null(fit[[nm]])) out[[nm]] <- fit[[nm]]
  }
  if (length(out) == 0) NULL else out
}
