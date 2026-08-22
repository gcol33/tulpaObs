# ---------------------------------------------------------------------------
# Joint nested-Laplace fit (shared field across the occurrence and cover arms)
# ---------------------------------------------------------------------------

# Collapse the occurrence (binomial) arm to its exact sufficient statistic.
# Observations agreeing on the occurrence design row AND every per-observation
# component of the linear predictor -- the spatial cell, any temporal / RE block
# index, any SVC weight -- are exchangeable Bernoulli trials, so replacing them
# with one Binomial row (n = count, y = successes) leaves the gradient and
# Hessian pointwise unchanged, and with them the mode, the SEs and every
# posterior weight, while cutting the row count.
#
# The log-likelihood shifts by a parameter-free constant: the engine's binomial
# kernel DOES carry lchoose(n, y), which is 0 on each Bernoulli row and
# lchoose(n, y) on the group that replaces them. `lconst` returns that shift so
# the caller can subtract it and keep the reported marginal on the scale of the
# observed per-plot sequence. Each input row is one Bernoulli trial here (one
# plot, y in {0, 1}), so the pre-aggregation constant is 0.
#
# `keys` is the named list of those per-observation vectors; the X-row plus every
# key forms the grouping key. Returns the aggregated `y`, `n` (n_trials), `X`,
# the per-group representative of each key, `rep_i` (the representative row per
# group) and `lconst`.
.cover_aggregate_occ <- function(y, X, keys) {
  parts <- c(as.data.frame(X, stringsAsFactors = FALSE), keys)
  gid   <- as.integer(factor(do.call(paste, c(parts, list(sep = "\r")))))
  ord   <- order(gid)
  rep_i <- ord[!duplicated(gid[ord])]               # first row per group, group order
  y_agg <- as.numeric(rowsum(as.numeric(y), gid))   # rowsum orders by sorted gid
  n_agg <- as.integer(tabulate(gid, nbins = max(gid)))
  list(
    y      = y_agg,
    n      = n_agg,
    X      = X[rep_i, , drop = FALSE],
    keys   = lapply(keys, function(v) v[rep_i]),
    rep_i  = rep_i,
    lconst = sum(lchoose(n_agg, y_agg))
  )
}

# Aggregate the occurrence arm in place against `keys` and reset its trial /
# RE bookkeeping. Returns the updated arm, the per-group key representatives for
# the caller to scatter back onto the latent blocks, and the binomial
# combinatorial constant the reduction introduces.
.cover_apply_occ_agg <- function(arm_occ, keys) {
  og <- .cover_aggregate_occ(arm_occ$y, arm_occ$X, keys)
  arm_occ$y        <- og$y
  arm_occ$n_trials <- og$n
  arm_occ$X        <- og$X
  arm_occ$re_idx   <- rep(0, length(og$y))
  list(arm_occ = arm_occ, keys = og$keys, lconst = og$lconst)
}

# Collapse the positive (beta) arm to its exact grouped sufficient statistics.
# Beta is not a count family, so there is no single-row collapse: a group of
# plots sharing the positive design row AND every per-observation component of
# its linear predictor (cell, trend weight, RE/time index) are exchangeable
# Beta(mu*phi, (1-mu)*phi) draws, and the beta log-density is linear in log(y)
# and log(1-y). One row carrying (n = count, slog_y = sum log(y), slog_1my = sum
# log(1-y)) therefore leaves the log-likelihood, gradient and (Fisher) Hessian
# pointwise unchanged -- the per-arm `n_trials`, `slog_y`, `slog_1my` are read by
# tulpa's built-in beta spec. `keys` is the named list of per-observation latent
# components; the X-row plus every key forms the grouping key.
.cover_aggregate_pos <- function(y, X, keys) {
  parts <- c(as.data.frame(X, stringsAsFactors = FALSE), keys)
  gid   <- as.integer(factor(do.call(paste, c(parts, list(sep = "\r")))))
  ord   <- order(gid)
  rep_i <- ord[!duplicated(gid[ord])]               # first row per group, group order
  ly    <- log(y)
  l1my  <- log1p(-y)
  list(
    y        = as.numeric(y[rep_i]),                  # representative (length only;
    n        = as.integer(tabulate(gid, nbins = max(gid))),  # the grouped beta spec
    slog_y   = as.numeric(rowsum(ly,   gid)),         # reads n/slog_y/slog_1my, not y)
    slog_1my = as.numeric(rowsum(l1my, gid)),
    X        = X[rep_i, , drop = FALSE],
    keys     = lapply(keys, function(v) v[rep_i]),
    rep_i    = rep_i
  )
}

# Aggregate the positive arm in place against `keys`, attaching the grouped beta
# sufficient statistics and resetting trial / RE bookkeeping. Returns the updated
# arm plus the per-group key representatives for the caller to scatter back onto
# the latent blocks' positive slot.
.cover_apply_pos_agg <- function(arm_pos, keys) {
  ag <- .cover_aggregate_pos(arm_pos$y, arm_pos$X, keys)
  arm_pos$y        <- ag$y
  arm_pos$n_trials <- ag$n
  arm_pos$slog_y   <- ag$slog_y
  arm_pos$slog_1my <- ag$slog_1my
  arm_pos$X        <- ag$X
  arm_pos$re_idx   <- rep(0, length(ag$y))
  list(arm_pos = arm_pos, keys = ag$keys)
}

# The per-arm index field a joint prior block keys on: structured spatial blocks
# carry `spatial_idx`, AR1/RW temporal blocks `temporal_idx`, IID temporal / RE
# blocks `obs_idx`. Each is a list(occ_idx, pos_idx).
.cover_block_idx_field <- function(block) {
  for (f in c("spatial_idx", "temporal_idx", "obs_idx")) {
    if (!is.null(block[[f]])) return(f)
  }
  NULL
}

# Gather every component of one arm's linear predictor carried by the prior
# blocks -- each block's per-arm index and, for a weighted (SVC) block, its
# per-arm weight -- as the exchangeability key (paired with that arm's design
# row). `slot` selects the arm within each block's per-arm lists: 1 = occurrence,
# 2 = positive.
.cover_arm_keys_from_blocks <- function(blocks, slot) {
  keys <- list()
  for (b in seq_along(blocks)) {
    f <- .cover_block_idx_field(blocks[[b]])
    if (!is.null(f)) keys[[sprintf("idx%d", b)]] <- blocks[[b]][[f]][[slot]]
    if (!is.null(blocks[[b]]$svc_weight)) {
      keys[[sprintf("w%d", b)]] <- blocks[[b]]$svc_weight[[slot]]
    }
  }
  keys
}

.cover_occ_keys_from_blocks <- function(blocks) {
  .cover_arm_keys_from_blocks(blocks, 1L)
}
.cover_pos_keys_from_blocks <- function(blocks) {
  .cover_arm_keys_from_blocks(blocks, 2L)
}

# Write the per-group key representatives back onto one arm's slot of every
# block, leaving the other arm untouched. `slot`: 1 = occurrence, 2 = positive.
.cover_scatter_arm_keys <- function(blocks, keys, slot) {
  for (b in seq_along(blocks)) {
    f <- .cover_block_idx_field(blocks[[b]])
    if (!is.null(f)) {
      blocks[[b]][[f]][[slot]] <- as.integer(keys[[sprintf("idx%d", b)]])
    }
    if (!is.null(blocks[[b]]$svc_weight)) {
      blocks[[b]]$svc_weight[[slot]] <- as.numeric(keys[[sprintf("w%d", b)]])
    }
  }
  blocks
}

.cover_scatter_occ_keys <- function(blocks, keys) {
  .cover_scatter_arm_keys(blocks, keys, 1L)
}
.cover_scatter_pos_keys <- function(blocks, keys) {
  .cover_scatter_arm_keys(blocks, keys, 2L)
}

# Reconstruct the lower-Cholesky factor L (Sigma = L L') from one row of
# log-Cholesky coordinates (column-major lower triangle: diagonal stored as
# log L_ii, strict-lower as raw L_ij), matching the engine's
# mcar_sigma_inv_from_logchol layout. Small, self-contained (no tulpa internals).
.cover_mcar_logchol_to_L <- function(theta, p) {
  L <- matrix(0, p, p)
  idx <- 1L
  for (j in seq_len(p)) for (i in j:p) {
    L[i, j] <- if (i == j) exp(theta[idx]) else theta[idx]
    idx <- idx + 1L
  }
  L
}

# ---------------------------------------------------------------------------
# Shared spine of the three joint cover fitters (arm-specific `||`, correlated
# MCAR `|`, and the shared-field nested path). They differ only in how they
# build their latent prior blocks and what derived quantities they report; the
# arm build, the engine control list, and the per-arm beta post-processing are
# the same computation in all three. Extracted here so a correction to the
# reported se_occ / se_pos cannot land on one route and miss the other two.
# ---------------------------------------------------------------------------

# Both response arms for the joint engine, plus the positive-arm family and its
# dispersion grid. `spi_occ` / `spi_pos` add the per-observation `spatial_idx`
# the single-block shared-field path carries; the multi-block routes pass NULL
# and declare each block's own per-arm index instead. Opt-in fixed-effect priors
# (cover_priors()) are attached here so every route penalises identically: the
# natural-scale numbers apply at face value to the (autoscaled) design, where
# every predictor is O(1), and precisions are floored at the engine's own weak
# default so an Inf-sd bucket reproduces the weak ridge rather than dropping the
# diagonal. NULL / FALSE / "none" leaves both arms unpenalised.
.cover_joint_arms <- function(enc, positive, control, priors,
                              spi_occ = NULL, spi_pos = NULL) {
  N     <- enc$N
  N_pos <- length(enc$pos_data$y)
  pfg <- .cover_pos_family_grid(positive, enc, control)

  arm_occ <- c(
    list(y = as.numeric(enc$occ_data$y), n_trials = enc$occ_data$n_trials,
         X = enc$occ_data$X),
    if (!is.null(spi_occ)) list(spatial_idx = as.integer(spi_occ)),
    list(re_idx = rep(0, N), n_re_groups = 0L, sigma_re = 1.0,
         family = "binomial", phi = 1.0))
  arm_pos <- c(
    list(y = as.numeric(enc$pos_data$y), n_trials = rep(1L, N_pos),
         X = enc$pos_data$X),
    if (!is.null(spi_pos)) list(spatial_idx = as.integer(spi_pos)),
    list(re_idx = rep(0, N_pos), n_re_groups = 0L, sigma_re = 1.0,
         family = pfg$pos_family, phi = pfg$phi_hat))
  arm_pos <- .cover_arm_pos_bounds(arm_pos, enc, positive)

  cprior <- .resolve_cover_priors(priors)
  if (!is.null(cprior)) {
    to_prec <- function(pr) {
      if (is.null(pr) || length(pr$sd) == 0L) return(NULL)
      list(mean = as.numeric(pr$mean), prec = pmax(1 / pr$sd^2, 1e-4))
    }
    occ_ap <- to_prec(.cover_arm_prior(cprior, "occ", colnames(arm_occ$X)))
    pos_ap <- to_prec(.cover_arm_prior(cprior, "pos", colnames(arm_pos$X)))
    if (!is.null(occ_ap)) {
      arm_occ$beta_prior_mean <- occ_ap$mean
      arm_occ$beta_prior_prec <- occ_ap$prec
    }
    if (!is.null(pos_ap)) {
      arm_pos$beta_prior_mean <- pos_ap$mean
      arm_pos$beta_prior_prec <- pos_ap$prec
    }
  }

  list(arm_occ = arm_occ, arm_pos = arm_pos,
       pos_family = pfg$pos_family, phi_hat = pfg$phi_hat,
       phi_grid_pos = pfg$phi_grid_pos)
}

# The joint engine's control list. `integration` is the route's default outer-
# grid layout, used when the caller does not set `control$integration`: the
# multi-block routes carry >= 3 latent axes and default to the mode-centred CCD
# (the dense tensor would blow up; the engine declines back to the tensor grid
# on a ridge), while the single-block path passes NULL and lets the engine
# choose. `prune` adds the opt-in cheap-pass screen axes, which only the
# shared-field path exposes.
#
# Notes on the individual entries, which every route inherits:
#   * hessian -- the beta positive arm's observed mixture Hessian is indefinite
#     away from the mode, so observed-curvature Newton steps stall and the inner
#     Newton hits max.iter in every grid cell. Expected/Fisher curvature is PSD
#     by construction and converges in ~12 steps; the final mode-pass always
#     re-factorizes with the observed Hessian, so the reported SEs, log_det and
#     grid weights are unchanged. The lognormal arm is exactly quadratic (one
#     inner step), so observed curvature is already optimal -> "lm".
#   * n_threads_outer -- outer-grid parallelism, one replicated cell-solve
#     state per thread. Preferred over inner per-obs threads on many-core
#     hardware, where the mode-region cells dominate.
#   * adaptive_grid -- brackets the mode with FULL inner solves and densifies
#     near it, so it never approximates the marginal and cannot drop the true
#     mode.
#   * progress / progress.file -- two independent channels, both ON by default.
#     `progress` gates the console bar (NOT tied to `verbose`); `progress.file`
#     is emitted whenever non-empty and is the only channel that survives a
#     detached Start-Process stdout buffer. `[[` (exact) not `$`:
#     `control$progress` prefix-matches `progress.file`.
#   * checkpoint -- grid-cell checkpoint/resume, forwarded verbatim so a
#     killed run resumes instead of restarting.
.cover_joint_control <- function(control, positive, integration = NULL,
                                 prune = FALSE) {
  head <- list(
    max_iter  = control$max.iter  %||% 50L,
    tol       = control$tol       %||% 1e-6,
    n_threads = control$n.threads %||% 1L,
    n_threads_outer = control$n.threads.outer %||% 1L,
    store_Q   = TRUE,
    hessian   = control$hessian   %||% (if (positive == "beta") "fisher" else "lm"))
  screen <- if (prune) {
    list(prune     = control$prune     %||% FALSE,
         prune_tol = control$prune.tol %||% 1e-4)
  } else list()
  tail <- list(
    adaptive_grid             = control$adaptive.grid             %||% TRUE,
    adaptive_grid_edge_thresh = control$adaptive.grid.edge.thresh %||% 0.02,
    adaptive_grid_max_passes  = control$adaptive.grid.max.passes  %||% 1L,
    progress          = control[["progress"]]     %||% TRUE,
    progress.every    = control$progress.every    %||% 0L,
    progress.throttle = control$progress.throttle %||% 2,
    progress.file     = control$progress.file     %||% "",
    checkpoint        = control$checkpoint,
    integration       = control$integration %||% integration)
  c(head, screen, tail)
}

# Per-arm natural-scale beta posterior moments from a joint fit, by the law of
# total covariance over the outer grid.
#
# The engine returns per-cell modes and per-cell precision blocks in the SCALED
# design's parameterization (`encode_cover_hurdle()`), so each cell is
# transformed to the natural scale first and aggregated after. Doing it
# cell-by-cell on the full constrained vcov block preserves the intercept's
# cross-covariance contribution; a diag-only approach would underestimate the
# intercept SE. Returns the beta means, their SEs, and the layout pieces (`p_occ`
# / `p_pos`, the per-arm index vectors, the scale transforms) the callers report
# or reuse.
.cover_joint_beta_moments <- function(fit, enc) {
  layout <- fit$arm_layout
  p_occ  <- layout$p[1]; p_pos <- layout$p[2]
  bocc_idx <- layout$beta_start[1] + seq_len(p_occ)
  bpos_idx <- layout$beta_start[2] + seq_len(p_pos)

  scale_occ <- enc$scale_occ %||% .scale_meta(enc$occ_data$X)
  scale_pos <- enc$scale_pos %||% .scale_meta(enc$pos_data$X)
  T_occ <- .scale_transform(scale_occ); T_pos <- .scale_transform(scale_pos)
  modes_occ <- fit$modes[, bocc_idx, drop = FALSE] %*% t(T_occ)
  modes_pos <- fit$modes[, bpos_idx, drop = FALSE] %*% t(T_pos)
  beta_occ <- as.numeric(crossprod(fit$weights, modes_occ))
  beta_pos <- as.numeric(crossprod(fit$weights, modes_pos))
  var_of_means_occ <- as.numeric(crossprod(fit$weights, modes_occ^2)) - beta_occ^2
  var_of_means_pos <- as.numeric(crossprod(fit$weights, modes_pos^2)) - beta_pos^2

  inner_blocks <- .joint_inner_vcov_block(fit, c(bocc_idx, bpos_idx))
  if (is.null(inner_blocks)) {
    mean_of_var_occ <- rep(0, p_occ); mean_of_var_pos <- rep(0, p_pos)
  } else {
    occ_rows <- seq_along(bocc_idx)
    pos_rows <- length(bocc_idx) + seq_along(bpos_idx)
    n_grid_eff <- length(inner_blocks)
    diag_occ <- matrix(0, n_grid_eff, p_occ)
    diag_pos <- matrix(0, n_grid_eff, p_pos)
    for (k in seq_len(n_grid_eff)) {
      V_block <- inner_blocks[[k]]
      if (is.null(V_block)) next
      V_occ_nat <- T_occ %*% V_block[occ_rows, occ_rows, drop = FALSE] %*% t(T_occ)
      V_pos_nat <- T_pos %*% V_block[pos_rows, pos_rows, drop = FALSE] %*% t(T_pos)
      diag_occ[k, ] <- pmax(diag(V_occ_nat), 0)
      diag_pos[k, ] <- pmax(diag(V_pos_nat), 0)
    }
    w_eff <- fit$weights[seq_len(n_grid_eff)]
    mean_of_var_occ <- as.numeric(crossprod(w_eff, diag_occ))
    mean_of_var_pos <- as.numeric(crossprod(w_eff, diag_pos))
  }

  list(p_occ = p_occ, p_pos = p_pos, bocc_idx = bocc_idx, bpos_idx = bpos_idx,
       T_occ = T_occ, T_pos = T_pos,
       beta_occ = beta_occ, beta_pos = beta_pos,
       se_occ = sqrt(pmax(0, var_of_means_occ + mean_of_var_occ)),
       se_pos = sqrt(pmax(0, var_of_means_pos + mean_of_var_pos)))
}

# Build one areal latent block restricted to a SINGLE arm with NO cross-arm copy.
# The block carries the field's own precision axis (tau for icar/car/car_proper,
# sigma + rho for bym2) integrated on the outer grid; the OTHER arm's per-arm
# spatial_idx is the all-zero sentinel, which the joint multi-block scatter reads
# as "this arm's rows do not see this block" (the `l_b > 0` guard in
# nested_laplace_joint_multi.h), so the field contributes to exactly one arm.
# This is the no-copy assembler; the shared/copied intercept + trend path uses
# the copy spec instead (one assembler, copy off here, copy on there). `slot` is
# the active arm (1 = occ, 2 = pos); `n_occ` / `n_pos` size the sentinel;
# `idx_active` is the per-obs node code on the active arm. `svc_weight` (NULL for
# the intercept field) is the per-obs design column, placed on the active arm
# with a zero sentinel weight on the inactive arm.
.cover_armspecific_block <- function(type, graph, slot, idx_active,
                                     n_occ, n_pos, svc_weight, control,
                                     block_label) {
  csr <- adjacency_to_csr(graph)
  n_nodes <- nrow(graph)
  if (anyNA(idx_active) || min(idx_active) < 1L || max(idx_active) > n_nodes) {
    stop(sprintf(paste0(
      "cover() arm-specific field (%s): node index must be a 1..%d code into ",
      "the graph."), block_label, n_nodes), call. = FALSE)
  }
  n_active <- if (slot == 1L) n_occ else n_pos
  if (length(idx_active) != n_active) {
    stop(sprintf(paste0(
      "internal: arm-specific field (%s) active-arm idx length %d != arm n %d."),
      block_label, length(idx_active), n_active), call. = FALSE)
  }
  # Per-arm spatial_idx: active arm gets the node codes, the other arm all-zero.
  idx_occ <- if (slot == 1L) as.integer(idx_active) else integer(n_occ)
  idx_pos <- if (slot == 2L) as.integer(idx_active) else integer(n_pos)

  # Field precision grid (own amplitude per field). icar/car/car_proper use the
  # single-arm tau parameterization (sigma = 1/sqrt(tau)); bym2 adds rho. The
  # sigma grid default mirrors the shared path's amplitude range. Override via
  # control$sigma.grid (translated to tau for the intrinsic backends). Whether
  # this axis is ours to move: the marker goes on the vector finally written
  # into the block, so track the provenance here and apply it after the tau
  # translation / bym2 pairing below.
  sigma_auto <- is.null(control$sigma.grid)
  sigma_grid <- as.numeric(control$sigma.grid %||%
                           exp(seq(log(0.2), log(2.5), length.out = 7)))

  block <- list(
    type            = if (tolower(type) == "car") "icar" else tolower(type),
    n_spatial_units = as.integer(n_nodes),
    adj_row_ptr     = as.integer(csr$row_ptr),
    adj_col_idx     = as.integer(csr$col_idx),
    n_neighbors     = as.integer(csr$n_neighbors),
    spatial_idx     = list(idx_occ, idx_pos)
  )
  if (block$type %in% c("icar", "car_proper")) {
    block$tau_grid <- .tobs_mark_auto(sort(1.0 / sigma_grid^2), sigma_auto)
    if (block$type == "car_proper") {
      block$rho_car_grid <- .tobs_num_auto(
        control$rho.car.grid %||% .tobs_default_rho_car_grid())
    }
  } else if (block$type == "bym2") {
    # BYM2 fits as a non-copied length-2 block (structured phi ICAR + iid theta
    # on n_nodes), parameterized by (sigma, rho) on the outer grid. The
    # arm-specific draw projection reconstructs the rho-mixed unit field z =
    # sqrt(rho) * scale_factor * phi + sqrt(1 - rho) * theta from the two
    # sub-blocks, so predict / WAIC see the full mix. The block's bym2 grid is
    # the PAIRED (cartesian-expanded) (sigma, rho) vectors the registry
    # consumes, not two separate axes. The mixing-weight axis is READ from the
    # engine rather than restated here: this block is meant to reach the
    # registry as any other bym2 block does, and the engine's own `bym2_rho`
    # default is what every other route gets. A copy of its nodes goes stale
    # the moment the engine moves them, which upstream has already done once.
    rho_auto <- is.null(control$rho.grid)
    rho_vals <- as.numeric(control$rho.grid %||% tulpa:::.nl_grid_axis("bym2_rho"))
    gr <- expand.grid(sigma = sort(sigma_grid), rho = rho_vals,
                      KEEP.OUT.ATTRS = FALSE)
    block$sigma_grid   <- .tobs_mark_auto(gr$sigma, sigma_auto)
    block$rho_grid     <- .tobs_mark_auto(gr$rho, rho_auto)
    # The joint nested-Laplace engine is tulpa's: its `scale_factor` MULTIPLIES
    # the structured block, so it takes 1 / sqrt(s). The draw substrate and the
    # SLA reconstruction read this same stored value back and mirror the
    # engine's spelling, so they need no conversion.
    block$scale_factor <- .bym2_engine_scale(.bym2_scale(graph))
  } else {
    stop(sprintf(paste0(
      "cover() arm-specific field (%s): areal type '%s' is not supported. ",
      "Arm-specific separate fields (single-arm `to`) use icar / car / ",
      "car_proper / bym2."),
      block_label, type), call. = FALSE)
  }

  # Per-arm per-row design weight (SVC column on a non-intercept field): the
  # active arm carries the covariate column, the inactive arm a zero placeholder
  # (its rows skip the block anyway via the 0 spatial_idx, so the weight is never
  # read; the engine still validates per-arm length).
  if (!is.null(svc_weight)) {
    w_occ <- if (slot == 1L) as.numeric(svc_weight) else numeric(n_occ)
    w_pos <- if (slot == 2L) as.numeric(svc_weight) else numeric(n_pos)
    block$svc_weight <- list(w_occ, w_pos)
  }
  block
}

# Fit the cover hurdle with arm-specific separate spatial latent field(s): one
# or more single-arm `||` bars, each a per-arm areal field with its own
# precision and NO cross-arm copy. Each bar's design columns (intercept +
# covariates) become independent non-copied areal blocks placed on that bar's
# single arm via the 0-sentinel spatial_idx on the other arm; each field's
# precision is integrated on the outer nested-Laplace grid. Separate single-arm
# bars (one to = "presence", one to = "positive") are independent per-arm
# fields with no coupling between them. Same output contract as
# fit_cover_hurdle_joint_nested so decode_cover_hurdle_joint consumes it
# unchanged; the arm-specific block layout is recorded on the fit
# (`armspec_blocks`) so the joint-draw projection scatters each block onto its
# own arm only.
.fit_cover_hurdle_joint_armspecific <- function(enc, data, positive = enc$positive,
                                                control = list(), priors = NULL) {
  arms <- enc$armspec
  N     <- enc$N
  N_pos <- length(enc$pos_data$y)
  idx_pos <- enc$idx_pos

  # Both arms + the positive-arm family / dispersion grid (same regime as the
  # single-field path), via the shared spine.
  .arms        <- .cover_joint_arms(enc, positive, control, priors)
  arm_occ      <- .arms$arm_occ
  arm_pos      <- .arms$arm_pos
  phi_grid_pos <- .arms$phi_grid_pos

  # Build one NON-copied block per field column, restricted to its arm. The
  # positive arm's node codes are the occ codes subset by idx_pos. `armspec_meta`
  # records each block's active arm and field role so the draw projection
  # scatters the block onto its arm only.
  blocks <- list(); armspec_meta <- list()
  for (a in arms) {
    slot <- if (a$arm == "presence") 1L else 2L
    idx_active <- if (slot == 1L) a$idx_obs else a$idx_obs[idx_pos]
    for (fi in seq_along(a$fields)) {
      f <- a$fields[[fi]]
      svc <- if (isTRUE(f$is_intercept)) NULL else {
        if (slot == 1L) f$weight else f$weight[idx_pos]
      }
      label <- paste(a$arm, f$column_name, sep = ".")
      blk <- .cover_armspecific_block(
        type = a$type, graph = a$graph, slot = slot, idx_active = idx_active,
        n_occ = N, n_pos = N_pos, svc_weight = svc, control = control,
        block_label = label)
      blocks[[length(blocks) + 1L]] <- blk
      armspec_meta[[length(armspec_meta) + 1L]] <- list(
        arm = a$arm, slot = slot, column_name = f$column_name,
        is_intercept = isTRUE(f$is_intercept),
        weight_occ = if (slot == 1L && !isTRUE(f$is_intercept)) as.numeric(f$weight) else NULL,
        weight_pos = if (slot == 2L && !isTRUE(f$is_intercept)) as.numeric(f$weight[idx_pos]) else NULL,
        idx_active = as.integer(idx_active),
        n_nodes = as.integer(nrow(a$graph)),
        type = blk$type,
        scale_factor = blk$scale_factor %||% 1.0,
        label = label)
    }
  }

  # Each field carries 1 (icar/car) or 2 (bym2/car_proper) latent axes; with the
  # pos-arm phi axis the dense outer tensor grows fast, so the mode-centred CCD
  # is this route's default.
  joint_control <- .cover_joint_control(control, positive, integration = "ccd")

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = blocks,         # list-of-blocks -> multi-block path
    copy      = NULL,           # NO copy: each field on one arm only
    phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
    prior_sigma = control$prior.sigma,
    prior_phi = control$prior.phi,
    control = joint_control
  )

  # Shared per-arm beta post-processing (identical to the single-field path).
  bm <- .cover_joint_beta_moments(fit, enc)
  p_pos <- bm$p_pos
  beta_occ <- bm$beta_occ; beta_pos <- bm$beta_pos
  se_occ   <- bm$se_occ;   se_pos   <- bm$se_pos

  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive %in% c("lognormal", "lognormal_trunc", "gaussian")) {
    sigma_pos <- phi_mu; sigma_pos_sd <- phi_sd
    phi_pos <- NA_real_; phi_pos_sd <- NA_real_
  } else {
    sigma_pos <- NA_real_; sigma_pos_sd <- NA_real_
    phi_pos <- phi_mu; phi_pos_sd <- phi_sd
  }

  # Per-field amplitude (sigma) posterior, marginalized over the outer grid. Each
  # block b carries its own axis b<b>.tau (icar/car_proper) or b<b>.sigma (bym2);
  # sigma = 1/sqrt(tau) for the intrinsic backends. Report the grid-weighted mean
  # per field (the derived-quantity rule: marginalize, do not plug in a MAP).
  tg <- fit$theta_grid; w <- fit$weights
  fin <- is.finite(w) & w > 0
  tg <- tg[fin, , drop = FALSE]; w <- w[fin]; w <- w / sum(w)
  sigma_fields <- numeric(length(blocks))
  field_names  <- character(length(blocks))
  for (b in seq_along(blocks)) {
    field_names[b] <- armspec_meta[[b]]$label
    tau_col <- sprintf("b%d.tau", b)
    sig_col <- sprintf("b%d.sigma", b)
    if (tau_col %in% colnames(tg)) {
      sigma_fields[b] <- sum(w * (1.0 / sqrt(as.numeric(tg[, tau_col]))))
    } else if (sig_col %in% colnames(tg)) {
      sigma_fields[b] <- sum(w * as.numeric(tg[, sig_col]))
    } else {
      sigma_fields[b] <- NA_real_
    }
  }
  names(sigma_fields) <- field_names

  m_occ <- list(mode = beta_occ, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  m_pos <- list(mode = beta_pos, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  attr(fit, "scale_factor") <- 1.0

  list(
    m_occ = m_occ, m_pos = m_pos, positive = positive,
    sigma_pos = sigma_pos, sigma_pos_sd = sigma_pos_sd,
    phi_pos = phi_pos, phi_pos_sd = phi_pos_sd,
    pos_fit_n = N_pos, pos_fit_p = p_pos,
    beta_occ = beta_occ, beta_pos = beta_pos, se_occ = se_occ, se_pos = se_pos,
    n_fields = length(blocks),
    armspecific = TRUE,
    armspec_blocks = armspec_meta,
    sigma_armspecific = sigma_fields,
    joint = fit
  )
}

# Fit the cover hurdle with a correlated separable-MCAR coefficient field shared
# onto the positive arm. The p design columns of the bar become p coupled areal
# fields with a free cross-covariance Sigma (x) Q^-1 (within-arm covariance
# among the fields, integrated over the outer CCD in log-Cholesky coordinates);
# the whole correlated field is copied onto the positive arm with one estimated
# amplitude alpha (the cross-arm transfer). Same output contract as
# fit_cover_hurdle_joint_nested so decode_cover_hurdle_joint consumes it
# unchanged.
.fit_cover_hurdle_joint_mcar <- function(enc, data, positive = enc$positive,
                                         control = list(), priors = NULL) {
  mc    <- enc$mcar
  p     <- mc$n_fields
  N     <- enc$N
  N_pos <- length(enc$pos_data$y)
  idx_pos <- enc$idx_pos

  # Both arms + the positive-arm family / dispersion grid (same regime as the
  # single-field path), via the shared spine.
  .arms        <- .cover_joint_arms(enc, positive, control, priors)
  arm_occ      <- .arms$arm_occ
  arm_pos      <- .arms$arm_pos
  phi_grid_pos <- .arms$phi_grid_pos

  # MCAR block: per-arm cell index (occ, pos) and per-field per-arm design
  # weight (occ, pos). The positive arm slices the occ weights / cell index by
  # `idx_pos`, exactly as the `||` trend weight is subset. With both cover arms
  # on `to` the field is anchored on occ and COPIED to pos with one estimated
  # amplitude alpha (#64); with a single arm it sits on that arm alone via the
  # 0-sentinel cell index on the other arm (no cross-arm copy, #109).
  to_arms  <- mc$to %||% .tobs_cover_arms
  both_arm <- setequal(to_arms, .tobs_cover_arms)
  on_occ   <- "presence" %in% to_arms
  on_pos   <- "positive" %in% to_arms
  idx_occ      <- mc$idx_occ
  idx_pos_cell <- idx_occ[idx_pos]
  mcar_spi_occ <- if (on_occ) as.integer(idx_occ)      else integer(N)
  mcar_spi_pos <- if (on_pos) as.integer(idx_pos_cell) else integer(N_pos)
  field_weight <- lapply(mc$field_weight_occ, function(w_occ) list(
    if (on_occ) as.numeric(w_occ)          else numeric(N),
    if (on_pos) as.numeric(w_occ[idx_pos]) else numeric(N_pos)))

  mcar_block <- list(
    type            = "mcar",
    n_spatial_units = mc$n_spatial_units,
    n_fields        = as.integer(p),
    adj_row_ptr     = mc$adj_row_ptr,
    adj_col_idx     = mc$adj_col_idx,
    n_neighbors     = mc$n_neighbors,
    spatial_idx     = list(mcar_spi_occ, mcar_spi_pos),
    field_weight    = field_weight
  )

  copy_spec <- if (both_arm) {
    alpha_grid <- control$alpha.grid %||%
      .tobs_default_alpha_grid()
    list(arm = "pos", block = 1L,
         alpha_grid = .tobs_num_auto(alpha_grid))
  } else NULL

  # The MCAR block carries p(p+1)/2 + 1 latent axes (log-Cholesky + alpha), so
  # the outer grid uses the mode-centred CCD by default; the dense tensor would
  # blow up.
  joint_control <- .cover_joint_control(control, positive, integration = "ccd")

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = list(mcar_block),
    copy      = copy_spec,
    phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
    prior_sigma = control$prior.sigma,
    prior_phi = control$prior.phi,
    prior_alpha = control$prior.alpha,
    control = joint_control
  )

  # Shared per-arm beta post-processing (identical to the single-field path).
  bm <- .cover_joint_beta_moments(fit, enc)
  p_pos <- bm$p_pos
  beta_occ <- bm$beta_occ; beta_pos <- bm$beta_pos
  se_occ   <- bm$se_occ;   se_pos   <- bm$se_pos

  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive %in% c("lognormal", "lognormal_trunc", "gaussian")) {
    sigma_pos <- phi_mu; sigma_pos_sd <- phi_sd
    phi_pos <- NA_real_; phi_pos_sd <- NA_real_
  } else {
    sigma_pos <- NA_real_; sigma_pos_sd <- NA_real_
    phi_pos <- phi_mu; phi_pos_sd <- phi_sd
  }

  # Cross-covariance Sigma derived quantities, marginalized over the outer grid
  # (the marginalize-derived-quantities rule): reconstruct Sigma per grid cell
  # from the log-Cholesky axes b1.L<i><j>, derive (sigma_a, rho_ab), then take
  # grid-weighted moments. The alpha copy amplitude is its own posterior axis.
  m <- p * (p + 1L) / 2L
  axis_nm <- character(m); tt <- 1L
  for (j in seq_len(p)) for (i in j:p) {
    axis_nm[tt] <- sprintf("b1.L%d%d", i, j); tt <- tt + 1L
  }
  tg <- fit$theta_grid; w <- fit$weights
  fin <- is.finite(w) & w > 0
  tg <- tg[fin, , drop = FALSE]; w <- w[fin]; w <- w / sum(w)
  sd_mat  <- matrix(NA_real_, nrow(tg), p)
  rho_mat <- matrix(NA_real_, nrow(tg), p * (p - 1L) / 2L)
  for (k in seq_len(nrow(tg))) {
    L <- .cover_mcar_logchol_to_L(as.numeric(tg[k, axis_nm]), p)
    Sig <- L %*% t(L)
    sds <- sqrt(pmax(diag(Sig), 0))
    sd_mat[k, ] <- sds
    cc <- 1L
    for (a in seq_len(p - 1L)) for (b in (a + 1L):p) {
      rho_mat[k, cc] <- Sig[a, b] / max(sds[a] * sds[b], 1e-12)
      cc <- cc + 1L
    }
  }
  sigma_mcar <- as.numeric(crossprod(w, sd_mat))
  rho_mcar   <- as.numeric(crossprod(w, rho_mat))
  rho_names  <- character(0)
  for (a in seq_len(p - 1L)) for (b in (a + 1L):p)
    rho_names <- c(rho_names, sprintf("rho_%d%d", a, b))
  names(sigma_mcar) <- mc$field_names
  names(rho_mcar)   <- rho_names
  # Multi-block axis names carry the block prefix: the copy amplitude is b1.alpha.
  # A single-arm correlated field has no copy, so no alpha axis -- report NA
  # rather than indexing a missing name (which errors).
  has_alpha <- "b1.alpha" %in% names(fit$theta_mean)
  alpha_mu <- if (has_alpha) as.numeric(fit$theta_mean[["b1.alpha"]]) else NA_real_
  alpha_sd <- if (has_alpha) as.numeric(fit$theta_sd[["b1.alpha"]])   else NA_real_

  m_occ <- list(mode = beta_occ, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  m_pos <- list(mode = beta_pos, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  attr(fit, "scale_factor") <- 1.0

  list(
    m_occ = m_occ, m_pos = m_pos, positive = positive,
    sigma_pos = sigma_pos, sigma_pos_sd = sigma_pos_sd,
    phi_pos = phi_pos, phi_pos_sd = phi_pos_sd,
    pos_fit_n = N_pos, pos_fit_p = p_pos,
    beta_occ = beta_occ, beta_pos = beta_pos, se_occ = se_occ, se_pos = se_pos,
    spi_full = as.integer(idx_occ), spi_pos = as.integer(idx_pos_cell),
    n_cells = as.integer(mc$n_spatial_units), n_fields = as.integer(p),
    mcar = TRUE, mcar_field_names = mc$field_names,
    sigma_mcar = sigma_mcar, rho_mcar = rho_mcar,
    alpha_mcar = alpha_mu, alpha_mcar_sd = alpha_sd,
    joint = fit
  )
}

#' Fit cover_hurdle as a joint binomial+(gaussian|beta) model with shared
#' spatial field via [tulpa::tulpa_nested_laplace_joint()].
#'
#' For both positive parts the dispersion scalar is integrated on the outer
#' joint hyperparameter grid (per-arm `phi_pos` axis):
#'
#' * `positive = "lognormal"` (gaussian arm): the residual SD is the per-grid
#'   phi. The default 7-point log-spaced grid is centred on the non-spatial
#'   prefit from `.prefit_lognormal_sigma()` and spans `[sigma_hat / 3,
#'   sigma_hat * 3]`. The posterior mean and SD across that axis are
#'   surfaced as `sigma_pos` / `sigma_pos_sd` on the returned `cover_fit`.
#' * `positive = "beta"`: the beta precision is the per-grid phi. The
#'   default 7-point log-spaced grid spans `[2, 300]`; posterior mean and
#'   SD are surfaced as `phi_pos` / `phi_pos_sd`.
#'
#' Override the per-arm phi grid via `control$phi.grid`.
#'
#' @param enc Output of [encode_cover_hurdle()].
#' @param data The original (un-subsetted) data frame — required to resolve
#'   the spatial spec (group_var lookup, n_spatial_units check).
#' @param positive `"lognormal"` or `"beta"`.
#' @param control List with optional `max_iter`, `tol`, `n_threads`,
#'   `sigma_grid`, `rho_grid`, `rho_car_grid`, `alpha_grid`,
#'   `phi_init`, `phi_bounds` (the last two are forwarded to the beta
#'   pre-fit when `positive = "beta"`). For ICAR / CAR_proper backends
#'   `tau_grid` is also accepted and translated to `sigma_grid` as
#'   `sigma = 1 / sqrt(tau)`. The cover arm sees the shared field at
#'   amplitude `alpha * sigma`, so `alpha_grid` is the cover-arm
#'   amplitude axis and `sigma_grid` the donor's. Regularizing hyperpriors on the joint
#'   (sigma, alpha) axes can be set via `prior_sigma` (donor amplitude)
#'   and `prior_alpha` (copy coefficient) — each a length-2 list
#'   `list(family, params)` matching tulpa's `prior_sigma` / `prior_alpha`
#'   args. The prior on alpha directly regularizes the copy scalar at small
#'   `n_pos`, replacing the per-arm `prior_sigma_pos` of the pre-reparam
#'   API. `prior.phi` puts the same kind of regularizing hyperprior on the
#'   cover-arm dispersion grid (the beta precision under `positive = "beta"`,
#'   the log-scale SD under `lognormal`), re-weighting the `phi.grid` axis by
#'   the chosen density instead of an implicit flat prior; same
#'   `list(family, params)` form, forwarded to tulpa's `prior_phi`.
#' @param temporal,re Structured `temporal()` / `re()` blocks from the formula,
#'   stacked onto the shared spatial block via the multi-block joint engine.
#' @param priors Optional [cover_priors()] object (or coercible list). When
#'   supplied, the per-arm fixed-effect prior reaches the joint engine as a
#'   `beta_prior_mean` / `beta_prior_prec` on the occurrence and positive
#'   responses, mirroring the separate-Laplace path. `NULL` / `FALSE` /
#'   `"none"` leave both arms unpenalised.
#' @return List shaped like the single-Laplace fit output but with extra
#'   `joint` field carrying the raw `tulpa_nested_laplace_joint` result.
#' @keywords internal
fit_cover_hurdle_joint_nested <- function(enc, data, positive = enc$positive,
                                          control = list(),
                                          temporal = NULL, re = NULL,
                                          priors = NULL) {
  if (!positive %in% c("lognormal", "lognormal_trunc", "beta", "beta_oi", "ordinal", "gaussian")) {
    stop("method = 'nested_laplace'/'nested_laplace_sla' for cover() supports positive = ",
         "'lognormal', 'lognormal_trunc', 'beta', 'beta_oi', or 'ordinal'. Got '",
         positive, "'.", call. = FALSE)
  }
  has_mcar <- !is.null(enc$mcar)
  if (has_mcar) {
    if (!is.null(temporal) || (!is.null(re) && length(re) > 0L)) {
      stop("cover(): a correlated spatial bar (single `|`) cannot yet be ",
           "combined with temporal()/re() blocks in the same fit.",
           call. = FALSE)
    }
    return(.fit_cover_hurdle_joint_mcar(enc, data, positive, control,
                                        priors = priors))
  }
  if (!is.null(enc$armspec)) {
    if (!is.null(temporal) || (!is.null(re) && length(re) > 0L)) {
      stop("cover(): arm-specific spatial fields (single-arm `to`) cannot be ",
           "combined with temporal()/re() blocks in the same fit.",
           call. = FALSE)
    }
    return(.fit_cover_hurdle_joint_armspecific(enc, data, positive, control,
                                               priors = priors))
  }
  if (is.null(enc$spatial_spec)) {
    stop("method = 'nested_laplace'/'nested_laplace_sla' for cover() requires ",
         "a spatial term in the formula. Add one of `bym2(graph = adj)`, ",
         "`icar(graph = adj)`, `car(graph = adj)`, or `car_proper(graph = adj)` ",
         "to the latent-presence formula.", call. = FALSE)
  }
  spec <- enc$spatial_spec
  if (!inherits(spec, "tulpa_spatial")) {
    stop("method = 'nested_laplace'/'nested_laplace_sla' for cover(): `spatial` must be a ",
         "tulpa_spatial spec.", call. = FALSE)
  }
  spec_type <- tolower(spec$type)
  # tulpa::spatial_car() returns type = "car" but prior_from_spec maps it
  # to backend = "icar"; treat the two as equivalent at dispatch time.
  supported <- c("bym2", "icar", "car", "car_proper")
  if (!spec_type %in% supported) {
    stop("method = 'nested_laplace'/'nested_laplace_sla' for cover() supports spatial types: ",
         paste(shQuote(supported), collapse = ", "),
         ". Got type = '", spec$type, "'.", call. = FALSE)
  }

  # Resolve obs -> spatial unit via tulpa's prior_from_spec. The dropped-NA
  # rows in encode_cover_hurdle (obs_keep) shrink the obs set; subset the
  # spatial_idx vector accordingly.
  data_obs <- data[enc$obs_keep, , drop = FALSE]
  # Replicated CAR: a `by` factor on the shared bar replicates the field across the
  # factor's levels. Build I_L (x) Q and offset each observation's node into its
  # level's copy (tulpa::tulpa_bar_field_replicate), then resolve the field over
  # the replicated graph so its precision is the block-diagonal Kronecker and the
  # field hyperparameters are shared across levels (one sigma[, rho_car]). The
  # coupled trend block copies the same replicated structure and the copy onto the
  # positive arm carries the whole replicated field at the one estimated alpha, so
  # `by` composes with the shared field, the trend, and the cross-arm copy
  # unchanged. No `by` is the identity. Injecting the offset index as the field's
  # node column lets the SAME prior_from_spec build the precision and resolve the
  # index over the replicated graph in one pass.
  by_replicated <- FALSE
  if (!is.null(spec$by_var)) {
    if (is.null(data_obs[[spec$by_var]])) {
      stop(sprintf(paste0(
        "spatial(<bar>, by = \"%s\"): the replication-factor column was not ",
        "found in the data."), spec$by_var), call. = FALSE)
    }
    base_idx <- tulpa::prior_from_spec(spec, data_obs)$spatial_idx
    rep_info <- tulpa::tulpa_bar_field_replicate(spec$adjacency, base_idx,
                                                 data_obs[[spec$by_var]])
    spec$adjacency              <- rep_info$adjacency
    # n_spatial is the spec's cached node count; validate_spatial() checks the
    # replicated node index against it (not against nrow(adjacency)), so it must
    # grow to L * n_nodes alongside the replicated graph or a 1..(L*n) index
    # reads as out of range on the base graph.
    spec$n_spatial              <- rep_info$n_levels * rep_info$n_nodes
    data_obs[[".tobs_by_node"]] <- rep_info$index
    spec$group_var              <- ".tobs_by_node"
    by_replicated               <- rep_info$n_levels > 1L
  }
  # The replicated graph I_L (x) Q is L disjoint copies by construction, so the
  # generic "graph not fully connected" identifiability warning is a false alarm
  # for a `by` field (each level is its own connected component, sharing one
  # precision). Muffle only that message, only when replication actually expanded
  # the graph; any other warning passes through.
  prior <- if (by_replicated) {
    withCallingHandlers(
      tulpa::prior_from_spec(spec, data_obs),
      warning = function(w) {
        if (grepl("not fully connected", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })
  } else {
    tulpa::prior_from_spec(spec, data_obs)
  }
  spi_full <- prior$spatial_idx                 # length N (post-NA-drop)
  spi_pos  <- spi_full[enc$idx_pos]             # length N_pos

  N     <- enc$N
  N_pos <- length(enc$pos_data$y)

  has_multi <- !is.null(temporal) || (!is.null(re) && length(re) > 0L)

  # Coupled spatially-varying trend (SVC) field on the cover hurdle. The
  # unweighted areal formula term is the shared intercept field; a SECOND,
  # weighted areal term (`icar(graph = adj, weight = col, group_var = ...)`)
  # adds a shared areal field on the same graph, weighted per observation by
  # `col` and copied onto the positive arm with its own alpha axis. This is
  # the analogue of the INLA joint model's `f(cell.slope, time, model =
  # "besag") + f(cell.slope.ab, time, copy =)` spatially-varying time trend
  # (two coupled besag fields, two copy coefficients). The per-observation
  # weight enters each arm's field contribution via the engine's per-block
  # `svc_weight`. `enc$trend` carries the per-observation weight over
  # `data_obs`; subset it to the positive arm.
  trend_spec <- if (!is.null(enc$trend)) {
    list(w_occ = enc$trend$w_occ,
         w_pos = enc$trend$w_occ[enc$idx_pos],
         label = enc$trend$label)
  } else {
    NULL
  }
  has_trend  <- !is.null(trend_spec)
  if (has_trend && has_multi) {
    stop("cover(): a coupled trend field (a weighted areal term) cannot yet ",
         "be combined with temporal()/re() blocks in the same fit.", call. = FALSE)
  }

  # Both arms + the positive-arm family / dispersion, via the shared spine. Both
  # dispersion regimes integrate phi on the outer joint hyperparameter grid;
  # `arm_pos$phi` is a placeholder overridden per grid point by the joint engine.
  # The family / phi grid (lognormal noise SD, ordinal latent log-cover SD, or
  # beta precision; 7 log-spaced points, densified by the engine's mode-tracked
  # refinement near the peak) comes from `.cover_pos_family_grid()`; override it
  # via `control$phi.grid`. This is the single-block path, so each arm carries
  # its own per-observation `spatial_idx` (the multi-block branches below replace
  # it with per-block indices). Attaching the fixed-effect priors here lets them
  # ride through the aggregation / scatter steps, which mutate rows rather than
  # the design columns the prior keys on.
  .arms        <- .cover_joint_arms(enc, positive, control, priors,
                                    spi_occ = spi_full, spi_pos = spi_pos)
  arm_occ      <- .arms$arm_occ
  arm_pos      <- .arms$arm_pos
  phi_grid_pos <- .arms$phi_grid_pos

  # Strip the per-obs spatial_idx (tulpa_nested_laplace_joint takes it per
  # arm) and the legacy rho_bounds field (joint car_proper uses rho_car_grid).
  # Forward control-grid overrides per backend.
  #
  # The joint engine parameterizes the copy as (sigma, alpha): the donor
  # (occurrence) arm sees the field at amplitude `sigma`, the cover arm at
  # `alpha * sigma`. For ICAR / CAR_proper the spatial field is unit-precision,
  # so the donor amplitude is `prior$sigma_grid` (`control$sigma.grid`, or
  # `control$tau.grid` as `sigma = 1/sqrt(tau)`), and the coupling axis is the
  # copy's own `alpha_grid` (`control$alpha.grid`).
  prior_for_joint <- prior
  prior_for_joint$spatial_idx <- NULL
  prior_for_joint$rho_bounds  <- NULL
  if (!is.null(control$sigma.grid))   prior_for_joint$sigma_grid   <- control$sigma.grid
  if (!is.null(control$rho.grid))     prior_for_joint$rho_grid     <- control$rho.grid
  if (!is.null(control$tau.grid)) {
    prior_for_joint$sigma_grid <- 1.0 / sqrt(as.numeric(control$tau.grid))
  }
  if (!is.null(control$rho.car.grid)) prior_for_joint$rho_car_grid <- control$rho.car.grid

  # Direct (sigma, alpha) copy axis: the cover arm sees the shared field at
  # amplitude alpha * sigma_donor. The single-block path declares alpha on the
  # pos arm via field_coef (the engine takes no single-block `copy`); the
  # multi-block branches carry it on their copy spec(s).
  alpha_grid <- control$alpha.grid %||%
    .tobs_default_alpha_grid()

  # Outer joint-grid integration controls, shared by the multi-block and
  # single-block dispatch. The dense outer tensor (sigma x [rho] x alpha x
  # phi_pos) concentrates almost all posterior mass on a handful of cells, but
  # the inner latent mode moves substantially across the grid, so the
  # cheap-pass prune is OFF by default: the full-grid solve is the correct
  # default. The rank-safe speed path is the adaptive grid (`adaptive_grid =
  # TRUE`): it brackets the mode with FULL inner solves and densifies near it,
  # so it never approximates the marginal and cannot drop the true mode. The
  # cheap-pass prune is available opt-in via control$prune = TRUE; it is now
  # rank-faithful (a neighbour- warm-started lattice sweep) and gated (a
  # safety check falls back to the full grid if the screen's ranking looks
  # unreliable), but the correct full grid remains the default. Override via
  # control$prune / control$prune.tol. This is the only route that exposes the
  # opt-in cheap-pass screen (`control$prune`), and the only one with no
  # outer-grid layout default: the coupled-trend multi-block path (>= 3 latent
  # axes: intercept + trend sigma/alpha) can request "ccd" via
  # control$integration, and NULL falls through to the engine default.
  joint_control <- .cover_joint_control(control, positive, prune = TRUE)

  # Exact sufficient-statistic reduction of the occurrence (binomial) arm,
  # default ON. The collapse is pointwise exact -- observations sharing the
  # occurrence design row AND every per-observation latent component (cell, trend
  # weight, RE/time index) are exchangeable Bernoulli trials, so one Binomial row
  # (n = count, y = successes) leaves the gradient and Hessian unchanged, and
  # shifts the log-likelihood by the parameter-free binomial constant the
  # reduction introduces (`agg_lconst` below, removed from the reported
  # marginal). Multi-seed parameter recovery on the aggregated path holds against
  # simulated truth (test-cover-hurdle-aggregate-recovery.R), so the reduction is
  # the default; set control$aggregate.occ = FALSE for the full per-plot
  # occurrence arm. `[[` (exact), never `$` (prefix-matching).
  do_agg_occ <- !isFALSE(control[["aggregate.occ"]])

  # sum lchoose(n_g, y_g) over the groups the occurrence reduction forms, in
  # whichever branch below runs it. Subtracted from the engine's log-marginal so
  # the reported value stays the marginal of the observed per-plot sequence and
  # does not depend on an internal performance switch. Every grid cell shifts by
  # the same amount, so posterior weights, the mode and the SEs are untouched.
  agg_lconst <- 0

  # Exact grouped sufficient-statistic reduction of the positive (beta) arm,
  # default ON for the beta arm. Beta has no single-row collapse, so plots
  # sharing the positive design row AND every per-observation latent component
  # are collapsed to one row carrying (n, sum log y, sum log(1-y)); tulpa's
  # built-in beta spec reads those sufficient statistics. Byte-identical to the
  # full per-plot beta arm (test-cover-hurdle-aggregate-pos.R), with the
  # both-arms-aggregated default behind a multi-seed parameter-recovery suite
  # (test-cover-hurdle-aggregate-recovery.R), so the reduction is the default;
  # set control$aggregate.pos = FALSE for the full per-plot positive arm. The
  # collapse is beta-only -- a lognormal positive arm would need its own (n,
  # sum, sum-of-squares) statistics, so an EXPLICIT aggregate.pos = TRUE errors
  # there rather than silently no-op, while the default leaves a non-beta arm
  # untouched. `[[` (exact), never `$` (prefix-matching).
  agg_pos_req <- control[["aggregate.pos"]]
  do_agg_pos  <- if (positive == "beta") !isFALSE(agg_pos_req) else isTRUE(agg_pos_req)
  if (isTRUE(agg_pos_req) && positive != "beta") {
    stop("control$aggregate.pos = TRUE is implemented for positive = \"beta\" ",
         "only (grouped beta sufficient statistics). Got positive = '",
         positive, "'.", call. = FALSE)
  }

  # ---- Multi-block path (Phase J-D) -----------------------------------
  # When `temporal` or `re` components are supplied, stack the spatial
  # block with AR1/RW/IID blocks and dispatch through the multi-block
  # joint engine. Copy semantics remain on the spatial block, carried by its
  # own `alpha` axis (`control$alpha.grid`); other blocks are shared
  # identically across the two arms (no per-arm scale).
  if (has_trend) {
    # Coupled trend path: two shared areal blocks on the same graph -- block 1
    # the unweighted intercept field, block 2 the per-observation-weighted SVC
    # field -- each copied onto the positive arm with its own alpha axis. The
    # engine's multi-block driver (list-valued prior + list-valued copy) carries
    # the per-block svc_weight and per-block alpha; the returned theta_grid axes
    # are b<k>.sigma / b<k>.alpha ( on the cover hurdle).
    base_block <- prior_for_joint
    if (is.null(base_block$sigma_grid)) {
      base_block$sigma_grid <- .tobs_default_sigma_grid()
    }
    if (tolower(base_block$type) == "bym2" && is.null(base_block$rho_grid)) {
      base_block$rho_grid <- .tobs_default_bym2_rho_grid()
    }
    base_block$spatial_idx <- list(as.integer(spi_full), as.integer(spi_pos))

    trend_block <- base_block
    trend_block$svc_weight <- list(as.numeric(trend_spec$w_occ),
                                   as.numeric(trend_spec$w_pos))

    # Exact sufficient-statistic reduction of the occurrence arm. The positive
    # arm is untouched; the two arms couple only through the shared cell field
    # (intercept + trend), and the grouping keys on the cell index and the trend
    # weight, so both fields' per-occ contributions are preserved.
    if (do_agg_occ) {
      blocks <- list(base_block, trend_block)
      ag      <- .cover_apply_occ_agg(arm_occ, .cover_occ_keys_from_blocks(blocks))
      arm_occ <- ag$arm_occ
      agg_lconst <- ag$lconst
      blocks  <- .cover_scatter_occ_keys(blocks, ag$keys)
      base_block  <- blocks[[1L]]
      trend_block <- blocks[[2L]]
    }
    if (do_agg_pos) {
      blocks  <- list(base_block, trend_block)
      agp     <- .cover_apply_pos_agg(arm_pos, .cover_pos_keys_from_blocks(blocks))
      arm_pos <- agp$arm_pos
      blocks  <- .cover_scatter_pos_keys(blocks, agp$keys)
      base_block  <- blocks[[1L]]
      trend_block <- blocks[[2L]]
    }

    alpha_grid_base  <- control$alpha.grid %||%
      .tobs_default_alpha_grid()
    alpha_grid_trend <- control$alpha.grid.trend %||% alpha_grid_base

    prior_coupled <- list(base_block, trend_block)
    copy_coupled  <- list(
      list(arm = "pos", block = 1L,
           alpha_grid = .tobs_num_auto(alpha_grid_base)),
      list(arm = "pos", block = 2L,
           alpha_grid = .tobs_num_auto(alpha_grid_trend))
    )
    arm_occ$spatial_idx <- NULL
    arm_pos$spatial_idx <- NULL
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = prior_coupled,
      copy      = copy_coupled,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma = control$prior.sigma,
      prior_phi = control$prior.phi,
      prior_alpha = control$prior.alpha,
      control = joint_control
    )
  } else if (has_multi) {
    multi <- .cover_build_multi_prior(
      prior_spatial = prior_for_joint,
      spi_full      = spi_full,
      spi_pos       = spi_pos,
      idx_pos       = enc$idx_pos,
      temporal      = temporal,
      re            = re,
      control       = control,
      alpha_grid    = alpha_grid
    )
    # Strip spatial_idx from the arms — it lives inside the spatial
    # block's per-arm spatial_idx list in the multi-block prior.
    arm_occ$spatial_idx <- NULL
    arm_pos$spatial_idx <- NULL
    # Exact sufficient-statistic reduction of the occurrence arm. Every block's
    # per-occ index (spatial cell, AR1/RW/IID time, RE group) enters the
    # grouping key, so only observations sharing the FULL linear predictor merge;
    # the representatives are scattered back onto each block's occ arm.
    if (do_agg_occ) {
      ag          <- .cover_apply_occ_agg(arm_occ, .cover_occ_keys_from_blocks(multi$prior))
      arm_occ     <- ag$arm_occ
      agg_lconst  <- ag$lconst
      multi$prior <- .cover_scatter_occ_keys(multi$prior, ag$keys)
    }
    if (do_agg_pos) {
      agp         <- .cover_apply_pos_agg(arm_pos, .cover_pos_keys_from_blocks(multi$prior))
      arm_pos     <- agp$arm_pos
      multi$prior <- .cover_scatter_pos_keys(multi$prior, agp$keys)
    }
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = multi$prior,
      copy      = multi$copy,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma = control$prior.sigma,
      prior_phi = control$prior.phi,
      prior_alpha = control$prior.alpha,
      control = joint_control
    )
  } else {
    # Exact sufficient-statistic reduction of the occurrence arm. The single
    # spatial cell index (carried on the arm here, not in the block) is the only
    # per-occ latent component, so it is the sole grouping key beyond the design.
    if (do_agg_occ) {
      ag      <- .cover_apply_occ_agg(arm_occ, list(idx1 = arm_occ$spatial_idx))
      arm_occ <- ag$arm_occ
      agg_lconst <- ag$lconst
      arm_occ$spatial_idx <- as.integer(ag$keys$idx1)
    }
    if (do_agg_pos) {
      agp     <- .cover_apply_pos_agg(arm_pos, list(idx1 = arm_pos$spatial_idx))
      arm_pos <- agp$arm_pos
      arm_pos$spatial_idx <- as.integer(agp$keys$idx1)
    }
    # Adaptive grid forwarding. Defaults match the joint engine's defaults
    # (`adaptive_grid = TRUE`, threshold 0.02, one pass) and triggered the
    # under-coverage fix in INLAabun D3. Pass `control$adaptive.grid
    # = FALSE` to recover the legacy fixed-grid behaviour for
    # reproducibility checks.
    arm_pos$field_coef <- list(
      name = "alpha",
      grid = .tobs_num_auto(alpha_grid))
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = prior_for_joint,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma = control$prior.sigma,
      prior_phi = control$prior.phi,
      prior_alpha = control$prior.alpha,
      control = joint_control
    )
  }

  # Put the reported marginal back on the observed per-plot scale (see
  # `agg_lconst` above). A constant shared by every grid cell, so the integrated
  # posterior is untouched -- this only stops the reported value from moving when
  # control$aggregate.occ flips.
  if (agg_lconst != 0 && !is.null(fit$log_marginal)) {
    fit$log_marginal <- fit$log_marginal - agg_lconst
  }

  # Posterior-weighted mean / SE for the per-arm beta blocks.
  bm <- .cover_joint_beta_moments(fit, enc)
  p_pos <- bm$p_pos
  beta_occ <- bm$beta_occ; beta_pos <- bm$beta_pos
  se_occ   <- bm$se_occ;   se_pos   <- bm$se_pos

  # Dispersion summary on the positive arm. Both regimes integrate the
  # dispersion scalar on the outer joint hyperparameter grid; read the
  # posterior mean and SD from the engine's `theta_mean` / `theta_sd`. Those
  # are computed against the phi-axis marginal (foreign-axis slice cells
  # filtered out by `.joint_recalibrate_axis_moments`) with Laplace-at-mode
  # SD at the modal cell, so they are grid-spacing- independent.
  # Hand-rolling `sum(weights * theta_grid^2) - mean^2` against
  # `theta_grid[, "phi_pos"]` underestimates SD on sharply peaked axes and
  # additionally collapses on slice cells that pin phi at the modal value
  # while varying other axes -- that's the legacy pattern were added to
  # replace.
  #
  # The phi axis carries the gaussian residual SD for lognormal and the
  # beta precision for beta; surface under the respective slot names.
  # The phi axis carries the gaussian residual SD (lognormal) or the latent
  # log-cover SD (ordinal interval-censored Gaussian) -- both surfaced as
  # sigma_pos -- and the beta precision otherwise (phi_pos).
  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive %in% c("lognormal", "lognormal_trunc", "ordinal", "gaussian")) {
    sigma_pos    <- phi_mu
    sigma_pos_sd <- phi_sd
    phi_pos      <- NA_real_
    phi_pos_sd   <- NA_real_
  } else {
    sigma_pos    <- NA_real_
    sigma_pos_sd <- NA_real_
    phi_pos      <- phi_mu
    phi_pos_sd   <- phi_sd
  }

  m_occ <- list(mode = beta_occ, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)
  m_pos <- list(mode = beta_pos, H_beta = NULL, converged = TRUE,
                log_marginal = NA_real_)

  # Stash the field-decomposition scale_factor (BYM2 Riebler scaling) on the
  # joint fit so the SLA path can reconstruct per-grid field amplitude
  # without re-deriving it. Dispersion is always integrated on `phi_pos`
  # (both lognormal and beta regimes), so the SLA path reads it directly
  # from `fit$theta_grid[k, "phi_pos"]` and needs no attr fallback.
  if (has_trend) {
    sf_attr <- as.numeric(base_block$scale_factor %||% 1.0)
  } else if (has_multi) {
    sf_attr <- as.numeric(multi$prior[[1L]]$scale_factor %||% 1.0)
  } else {
    sf_attr <- as.numeric(prior_for_joint$scale_factor %||% 1.0)
  }
  attr(fit, "scale_factor") <- sf_attr

  # Trend-field hyperparameter summaries (block 2: sigma_trend, alpha_trend),
  # read off the multi-block (sigma, alpha) axes of the integrated posterior.
  sigma_trend <- if (has_trend) as.numeric(fit$theta_mean[["b2.sigma"]] %||% NA) else NULL
  alpha_trend <- if (has_trend) as.numeric(fit$theta_mean[["b2.alpha"]] %||% NA) else NULL

  list(
    m_occ        = m_occ,
    m_pos        = m_pos,
    positive     = positive,
    sigma_pos    = sigma_pos,
    sigma_pos_sd = sigma_pos_sd,
    phi_pos      = phi_pos,
    phi_pos_sd   = phi_pos_sd,
    pos_fit_n    = N_pos,
    pos_fit_p    = p_pos,
    beta_occ     = beta_occ,
    beta_pos     = beta_pos,
    se_occ       = se_occ,
    se_pos       = se_pos,
    spi_full     = as.integer(spi_full),
    spi_pos      = as.integer(spi_pos),
    n_cells      = as.integer(prior_for_joint$n_spatial_units %||% NA),
    n_fields     = if (has_trend) 2L else 1L,
    trend_weight  = if (has_trend) trend_spec$label else NULL,
    trend_weights = if (has_trend) list(trend_spec$label) else NULL,
    trend_w_occ   = if (has_trend) trend_spec$w_occ else NULL,
    trend_w_pos   = if (has_trend) trend_spec$w_pos else NULL,
    sigma_trend   = sigma_trend,
    alpha_trend   = alpha_trend,
    joint        = fit
  )
}
