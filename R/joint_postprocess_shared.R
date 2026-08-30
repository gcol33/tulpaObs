# =============================================================================
# joint_postprocess_shared.R -- the arm-count-agnostic spine of the joint
# direct-grid postprocess.
#
# `.occu_cover_jc_postprocess()` (3 arms: psi, p, pos) and
# `.occu_jc_postprocess()` (2 arms: psi, p) shape the same object out of the
# same engine return. Everything that does not mention the cover arm is the same
# computation in both: dropping non-converged outer-grid cells, the law-of-
# total-covariance (betas + field) block, the parameter-surface vcov, and the
# per-field split into the intercept field plus the trend fields. The helpers
# below carry that spine, so a correction to the reported joint posterior lands
# on both routes.
#
# Every helper takes the arms as data -- an index vector or a list of per-arm
# (idx, mean) pairs -- so the 2-arm and 3-arm callers differ only in what they
# pass.
# =============================================================================


# Outer-grid weights for a driver result that carries a `theta_grid` but no
# reconciled `weights` of its own (or whose `log_marginal` a caller has since
# adjusted). Each cell carries the prior mass of the node it sits at, the same
# rule the engine weights its own grids with, so a post-processor never falls
# back to weighing every node equally: on an evenly spaced grid the two agree,
# and on an uneven one the spacing is what differs.
.tobs_grid_weights <- function(fit, what = "outer grid", log_marginal = NULL) {
  lm <- if (is.null(log_marginal)) fit$log_marginal else log_marginal
  lq <- tulpa:::.nl_grid_log_quad(fit$theta_grid)
  if (!is.null(lq) && length(lq) != length(lm)) lq <- NULL
  tulpa:::.nl_normalise_weights_safe(lm, what, log_quad = lq)
}


# Outer-grid cells whose inner Newton converged, their softmax weights, and the
# reconciled `fit$weights`.
#
# Cells with a non-finite `log_marginal` hold NaN modes that would poison every
# weighted sum downstream; zero-mass cells stay represented in `fit$weights`, the
# moments just route around them. When the engine left NO usable weight (an
# unguarded upstream normalization collapses `fit$weights` to all NaN once any
# cell is non-finite) the sampler `tulpa_posterior_draws()` uses for predict /
# WAIC would find no positive-weight cell, so fall back to the same pure softmax
# the reported moments use and keep the two consistent. Untouched when the engine
# weights are already usable (every finite-grid fit).
#
# `label` names the route in the error / warning text.
.tobs_joint_ok_cells <- function(fit, label) {
  ok_cells <- which(is.finite(fit$log_marginal))
  if (length(ok_cells) == 0L) {
    stop(sprintf("%s: inner Newton failed at every grid cell. ", label),
         "Bump control$max.iter or tighten control$tol.", call. = FALSE)
  }
  if (length(ok_cells) < length(fit$log_marginal)) {
    n_bad <- length(fit$log_marginal) - length(ok_cells)
    warning(sprintf("%s: dropping %d / %d outer-grid cell(s) ",
                    label, n_bad, length(fit$log_marginal)),
            "whose inner Newton did not converge.", call. = FALSE)
  }
  # The engine's cell weights carry the outer-grid quadrature: node spacing and
  # the declared hyperparameter prior, not the likelihood alone. Posterior
  # summaries are built by renormalising THOSE over the converged cells. A
  # softmax of `log_marginal` would weigh every node equally and so report a
  # posterior against a measure the fit never integrated.
  ew <- fit$weights
  ew_ok <- !is.null(ew) && length(ew) == length(fit$log_marginal) &&
           any(is.finite(ew[ok_cells]) & ew[ok_cells] > 0)
  if (ew_ok) {
    w <- ew[ok_cells]
    w[!is.finite(w) | w < 0] <- 0
    w <- w / sum(w)
  } else {
    w_raw <- exp(fit$log_marginal[ok_cells] - max(fit$log_marginal[ok_cells]))
    w     <- w_raw / sum(w_raw)
  }

  if (!any(is.finite(fit$weights) & fit$weights > 0)) {
    w_full <- numeric(length(fit$log_marginal))
    w_full[ok_cells] <- w
    fit$weights <- w_full
  }
  list(ok_cells = ok_cells, w = w, fit = fit)
}


# Joint (betas + field) posterior covariance by the law of total covariance over
# the outer hyperparameter grid:
#
#   Cov(x | y) = sum_k w_k [ Cov(x | y, theta_k) + (m_k - mbar)(m_k - mbar)' ]
#
# where x = (betas, field) stacked, Cov(x | y, theta_k) is the inner-Laplace
# covariance at grid cell k (the ICAR sum-to-zero constrained sub-block of
# Q_k^-1, so the intercept covariance is data-identified rather than collapsing
# along the (intercept, mean(phi)) near-null direction of the improper field
# prior), and (m_k - mbar) is the between-grid mode deviation. Carrying the field
# block (not just the betas) means downstream derived quantities can marginalize
# the joint betas+field posterior instead of a marginal-only diagonal.
#
# `arms` is a list of one `list(idx =, mean =)` per response arm, in latent-
# vector order; `field_idx` the stacked field columns. Older tulpa without stored
# per-grid Q returns no dense block, and the fallback is the marginal-only
# diagonal (var-of-means plus the diagonal inner-Laplace variance) with no
# betas+field cross-covariance -- `Vj` is NULL there, which every caller branches
# on.
.tobs_joint_beta_field_vcov <- function(fit, modes, w, ok_cells, arms, field_idx,
                                        n_cells, n_fields, n_threads = 1L) {
  beta_idx  <- as.integer(unlist(lapply(arms, `[[`, "idx")))
  p_beta    <- length(beta_idx)
  idx_joint <- c(beta_idx, field_idx)
  blocks    <- .joint_inner_vcov_block(fit, idx_joint, n_dense = p_beta,
                                       n_threads = n_threads)

  if (is.null(blocks)) {
    inner_var <- .joint_inner_var(fit, beta_idx)
    total_var <- function(modes_block, mean_vec, iv_block) {
      vom <- as.numeric(crossprod(w, modes_block^2)) - mean_vec^2
      mov <- if (is.null(iv_block)) {
        0
      } else {
        iv_k <- iv_block[ok_cells, , drop = FALSE]
        iv_k[!is.finite(iv_k)] <- 0  # rank-deficient Q -> var-of-means only
        as.numeric(crossprod(w, iv_k))
      }
      pmax(vom + mov, 0)
    }
    sds_beta <- numeric(0)
    off <- 0L
    for (a in arms) {
      p_a <- length(a$idx)
      iv  <- if (is.null(inner_var)) NULL
             else inner_var[, off + seq_len(p_a), drop = FALSE]
      sds_beta <- c(sds_beta,
                    sqrt(total_var(modes[, a$idx, drop = FALSE], a$mean, iv)))
      off <- off + p_a
    }
    beta_block <- diag(sds_beta^2, nrow = p_beta)

    field_modes    <- modes[, field_idx, drop = FALSE]
    field_at_cell  <- as.numeric(crossprod(w, field_modes))
    field_var      <- as.numeric(crossprod(w, field_modes^2)) - field_at_cell^2
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
    Vj <- NULL  # no joint covariance available
  } else {
    modes_joint <- modes[, idx_joint, drop = FALSE]
    mbar_joint  <- as.numeric(crossprod(w, modes_joint))
    # symmetrize off floating-point constraint residuals
    Vj <- .tobs_grid_vcov(modes_joint, w, blocks[ok_cells],
                          center = mbar_joint, on_missing = "zero",
                          symmetrize = TRUE)
    diag_Vj    <- diag(Vj)
    sds_beta   <- sqrt(pmax(diag_Vj[seq_len(p_beta)], 0))
    beta_block <- Vj[seq_len(p_beta), seq_len(p_beta), drop = FALSE]

    # Field summary uses the full (within + between) variance, demeaned to the
    # sum-to-zero convention the field-block covariance already sits under. One
    # block of n_cells columns per coupled field, in block order.
    n_field_cols   <- n_fields * n_cells
    field_at_cell  <- mbar_joint[p_beta + seq_len(n_field_cols)]
    field_var      <- diag_Vj[p_beta + seq_len(n_field_cols)]
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
  }

  list(beta_idx = beta_idx, p_beta = p_beta,
       sds_beta = sds_beta, beta_block = beta_block,
       field_demeaned = field_demeaned, field_sd = sqrt(pmax(field_var, 0)),
       Vj = Vj)
}


# Parameter-surface vcov by the law of total covariance over the outer grid.
#
# The betas carry within + between (`beta_block`); the hyperparameters ARE the
# grid coordinates, so within a cell their variance -- and their covariance with
# the betas -- is exactly zero, leaving only the between (variance- and
# covariance-of-modes) term. That between term is computed jointly over
# (betas, hyperparameters) so the cross-covariance and the hyper-hyper covariance
# are both retained rather than dropped to a diagonal block.
#
# On the older-tulpa fallback (`Vj = NULL`) the betas are a marginal-only
# diagonal, so a beta-hyper cross block could break PSD; keep the hyper block
# diagonal there, matching the betas' treatment.
.tobs_joint_param_vcov <- function(modes, w, beta_idx, beta_block, p_beta,
                                   hyper_names, hyper_vals, hyper_means,
                                   means, sds, par_names, Vj) {
  n_par <- length(means)
  V <- matrix(0, n_par, n_par)
  V[seq_len(p_beta), seq_len(p_beta)] <- beta_block
  if (length(hyper_names) > 0L) {
    hyper_idx <- p_beta + seq_along(hyper_names)
    H    <- do.call(cbind, hyper_vals[hyper_names])          # n_ok x n_hyper
    H_dm <- sweep(H, 2L, unlist(hyper_means)[hyper_names], "-")
    if (!is.null(Vj)) {
      beta_modes <- modes[, beta_idx, drop = FALSE]
      B_dm   <- sweep(beta_modes, 2L, means[seq_len(p_beta)], "-")
      cross  <- crossprod(B_dm * w, H_dm)                    # p_beta x n_hyper
      hyhy   <- crossprod(H_dm * w, H_dm)                    # n_hyper x n_hyper
      hyhy   <- (hyhy + t(hyhy)) / 2
      V[seq_len(p_beta), hyper_idx] <- cross
      V[hyper_idx, seq_len(p_beta)] <- t(cross)
      V[hyper_idx, hyper_idx]       <- hyhy
    } else {
      diag(V)[hyper_idx] <- sds[hyper_idx]^2
    }
  }
  dimnames(V) <- list(par_names, par_names)
  V
}


# Split the stacked per-field summaries into one block of n_cells per coupled
# field. Field 1 is the intercept field (the back-compat `spatial_field`); the
# next `n_own_fields - 1` are the spatially-varying trend fields, in block order,
# labelled by their weight column. Any blocks beyond `n_own_fields` belong to a
# different arm (the arm-specific cover fields) and are returned unsliced in
# `blocks` for the caller to name.
.tobs_joint_field_split <- function(field_demeaned, field_sd, n_cells, n_fields,
                                    n_own_fields, coupled_trends) {
  field_block <- function(b) {
    idx <- (b - 1L) * n_cells + seq_len(n_cells)
    list(mean = field_demeaned[idx], sd = field_sd[idx])
  }
  fblocks <- lapply(seq_len(n_fields), field_block)

  trend_blocks <- if (n_own_fields >= 2L) fblocks[2:n_own_fields] else list()
  trend_means  <- lapply(trend_blocks, function(b) b$mean)
  trend_tables <- lapply(trend_blocks, .tobs_joint_field_z_table)
  trend_labels <- vapply(coupled_trends, function(tf) tf$weight_label,
                         character(1))
  if (length(trend_labels) == length(trend_means)) {
    names(trend_means)  <- trend_labels
    names(trend_tables) <- trend_labels
  }
  list(blocks       = fblocks,
       field_z_table = .tobs_joint_field_z_table,
       intercept    = fblocks[[1L]]$mean,
       field_table  = .tobs_joint_field_z_table(fblocks[[1L]]),
       trend_means  = trend_means,
       trend_tables = trend_tables,
       trend_labels = trend_labels,
       # Back-compat single-trend accessors (the first trend field).
       field_trend       = if (length(trend_means))  trend_means[[1L]]  else NULL,
       trend_field_table = if (length(trend_tables)) trend_tables[[1L]] else NULL)
}

# Per-cell posterior summary of one field block.
.tobs_joint_field_z_table <- function(blk) {
  data.frame(cell = seq_len(length(blk$mean)), z_mean = blk$mean, z_sd = blk$sd,
             z_lower = blk$mean - 1.96 * blk$sd,
             z_upper = blk$mean + 1.96 * blk$sd)
}


# Names for the stacked field coordinates of the joint betas+field posterior, in
# block order: the intercept field, then the trend fields, then any arm-specific
# cover fields. A lone field of a kind drops the numeric suffix.
.tobs_joint_field_par_names <- function(n_cells, n_trend_fields,
                                        n_pos_fields = 0L) {
  run <- function(prefix, n) lapply(seq_len(n), function(j) {
    suffix <- if (n == 1L) "" else as.character(j)
    paste0(prefix, suffix, "_", seq_len(n_cells))
  })
  unlist(c(list(paste0("field_", seq_len(n_cells))),
           run("trend_field", n_trend_fields),
           run("pos_field",  n_pos_fields)))
}
