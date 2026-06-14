# =============================================================================
# joint_substrate.R -- shared access layer for the joint nested-Laplace fits
# behind cover() (2-arm hurdle) and occu_cover() (3-arm occupancy-cover).
#
# Both families fit a `tulpa_nested_laplace_joint` object carrying the same
# substrate: per-grid weights / modes / sparse precision, an `arm_layout`, a
# shared latent field, and a working `tulpa::tulpa_posterior_draws()`. They
# differ only in (a) which slot holds the joint object, (b) the arm roster
# (occ/pos vs psi/p/pos), and (c) the field-amplitude convention -- cover()'s
# independent (sigma_occ, sigma_pos) reparam versus occu_cover()'s donor sigma
# scaled onto the cover arm by the copy coefficient alpha. This file normalizes
# those three differences ONCE into a family-agnostic draw bundle so predict()
# and the pointwise log-likelihood share a single draws -> linear-predictor
# path (gcol33/tulpaObs#24).
# =============================================================================


# The joint nested-Laplace object regardless of family slot: occu_cover() stores
# it at `$joint_fit`, cover() at `$joint`. NULL when neither is present (a
# non-spatial / separate-Laplace fit that carries no joint object).
.tobs_joint_fit <- function(object) {
  object$joint_fit %||% object$joint
}

# Resolve a grid amplitude axis to a per-draw vector. A multi-block fit prefixes
# the axis with its block (`b<k>.sigma`); a single-block fit uses the bare name.
# `cells` is the outer-grid cell each draw came from, so the returned length-n
# vector carries the amplitude active for that draw. A missing axis returns a
# constant `default`.
.tobs_joint_amp <- function(theta_grid, cells, block, name, default = 1) {
  cn <- colnames(theta_grid)
  j  <- match(paste0("b", block, ".", name), cn)
  if (is.na(j)) j <- match(name, cn)
  if (is.na(j)) return(rep(default, length(cells)))
  as.numeric(theta_grid[cells, j])
}

# Draw the grid-integrated joint posterior and normalize it into a family-
# agnostic bundle. Returns:
#   $n         number of draws
#   $positive  "lognormal" / "beta"
#   $cells     length-n outer-grid cell each draw came from
#   $disp      length-n per-draw positive-arm dispersion (residual SD for
#              lognormal, precision for beta), read off the `phi_pos` grid axis
#   $b         list(occ = [n x p_occ], det = [n x p_det] | NULL,
#                   pos = [n x p_pos]) of per-arm coefficient draws (in the
#              fitted, scaled design space)
#   $blocks    list of shared-field blocks, each:
#                $z       [n x n_cells] unit-variance field draws
#                $amp_occ [n] occupancy-arm field amplitude per draw
#                $amp_pos [n] positive-arm field amplitude per draw
#                $weight  NULL (intercept field) or the per-cell covariate name
#                         weighting this (trend / SVC) field
#   $n_cells   spatial-unit count
.tobs_joint_draws <- function(object, n = 1000L) {
  jf <- .tobs_joint_fit(object)
  if (is.null(jf)) {
    stop("This fit carries no joint nested-Laplace object to sample ",
         "(`$joint_fit` / `$joint`).", call. = FALSE)
  }
  layout <- jf$arm_layout
  n_arms <- layout$n_arms %||% length(layout$p)
  if (isTRUE(object$occu_only_joint)) {
    .tobs_joint_draws_occu(object, jf, layout, n)
  } else if (n_arms == 3L) {
    .tobs_joint_draws_occu_cover(object, jf, layout, n)
  } else if (isTRUE(object$armspecific)) {
    .tobs_joint_draws_cover_armspecific(object, jf, layout, n)
  } else {
    .tobs_joint_draws_cover(object, jf, layout, n)
  }
}

# occu single-arm (2-arm psi/p): one or more independent ICAR fields on the
# OCCUPANCY (psi) arm only -- no cover arm, no copy. Each block is a unit-variance
# latent z scaled on the occupancy arm by its own amplitude (b<b>.sigma); the
# detection (p) arm carries no field (amp_pos = 0). Block 1 is the unweighted
# intercept field; blocks 2.. carry the per-cell trend weight. The bundle reuses
# the cover roster slot names ("occ" for psi, "pos" empty) so the shared
# `.tobs_joint_arm_eta` accumulator and the predict path read it unchanged: the
# occupancy psi arm is "occ", and there is no positive arm. (gcol33/tulpaObs#81)
.tobs_joint_draws_occu <- function(object, jf, layout, n) {
  tg      <- jf$theta_grid
  n_cells <- object$model$n_cells %||% object$model$n_sites
  p       <- layout$p
  bstart  <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_det <- bstart[2L] + seq_len(p[2L])
  starts  <- layout$field_starts %||% layout$phi_start
  n_field <- length(starts)
  field_idx <- lapply(starts, function(s0) s0 + seq_len(n_cells))

  idx   <- c(idx_occ, idx_det, unlist(field_idx))
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
  b_occ <- take(p[1L]); b_det <- take(p[2L])

  cn <- colnames(tg)
  trend_cols <- object$trend_weights %||% object$trend_weight
  blocks <- lapply(seq_len(n_field), function(b) {
    z <- take(n_cells)
    # Each ICAR block grids on precision tau (axis b<b>.tau); the occupancy-arm
    # field amplitude is the SD sigma = 1 / sqrt(tau) per draw cell. The
    # detection arm carries no field (amp_pos = 0, the 0-sentinel node index at
    # fit time already excluded it).
    tau_col <- sprintf("b%d.tau", b)
    amp <- if (tau_col %in% cn) 1.0 / sqrt(as.numeric(tg[cells, tau_col]))
           else rep(1.0, length(cells))
    list(z = z, amp_occ = amp, amp_pos = rep(0, length(cells)),
         weight = if (b == 1L) NULL else trend_cols[[b - 1L]])
  })

  list(n = n, positive = NA_character_, cells = cells,
       disp = rep(1, length(cells)),
       b = list(occ = b_occ, det = b_det, pos = NULL),
       blocks = blocks, n_cells = n_cells)
}

# cover (2-arm occ/pos) with arm-specific separate latents (gcol33/tulpaObs#65):
# one or more NON-copied areal blocks, each placed on exactly ONE arm. Block b is
# stored as a unit-precision latent z; its amplitude on the active arm is sigma_b
# (= 1/sqrt(tau_b) for icar/car_proper, the bym2 mixed amplitude otherwise), and
# its amplitude on the OTHER arm is 0 (no cross-arm copy). The bundle's per-block
# amp_occ / amp_pos encode this directly -- one is the field amplitude, the other
# is 0 -- so the shared `.tobs_joint_arm_eta` accumulator scatters each block onto
# its own arm only. A non-intercept (slope) field also carries its per-arm weight.
.tobs_joint_draws_cover_armspecific <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  cn       <- colnames(tg)
  positive <- object$positive %||% "lognormal"
  p        <- layout$p
  bstart   <- layout$beta_start
  meta     <- object$armspec_blocks

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_pos <- bstart[2L] + seq_len(p[2L])

  field_starts <- layout$field_starts
  n_field <- length(field_starts %||% integer(0))
  if (n_field != length(meta)) {
    stop("internal: arm-specific cover fit has ", length(meta), " field block(s) ",
         "but the joint layout reports ", n_field, ".", call. = FALSE)
  }
  # Every block spans n_nodes latent entries (the graph it sits on). Each block
  # records its own n_nodes; on a single shared graph these are equal, but read
  # per-block so a future multi-graph fit stays correct.
  n_nodes <- vapply(meta, function(m) as.integer(m$n_nodes), integer(1))
  field_idx <- lapply(seq_len(n_field), function(b)
    field_starts[b] + seq_len(n_nodes[b]))

  idx <- c(idx_occ, idx_pos, unlist(field_idx))
  D   <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")
  off <- 0L
  take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
  b_occ <- take(p[1L]); b_pos <- take(p[2L])

  blocks <- lapply(seq_len(n_field), function(b) {
    z <- take(n_nodes[b])
    # Field amplitude on the outer grid: b<b>.sigma (bym2) or 1/sqrt(b<b>.tau)
    # (icar / car_proper), per draw cell. The active arm scales z by this; the
    # inactive arm by 0 (no copy).
    sig_col <- sprintf("b%d.sigma", b)
    tau_col <- sprintf("b%d.tau", b)
    if (sig_col %in% cn) {
      amp <- as.numeric(tg[cells, sig_col])
    } else if (tau_col %in% cn) {
      amp <- 1.0 / sqrt(as.numeric(tg[cells, tau_col]))
    } else {
      amp <- rep(1.0, length(cells))
    }
    slot <- meta[[b]]$slot
    amp_occ <- if (slot == 1L) amp else rep(0, length(cells))
    amp_pos <- if (slot == 2L) amp else rep(0, length(cells))
    # The block sits on ONE arm, encoded by the zero amplitude on the other arm.
    # Its node map and per-cell weight are NOT stored here: the consumer supplies
    # them via `.tobs_joint_arm_eta`'s `units` / `wfun` -- predict() the newdata
    # cell map and column, the pointwise-loglik consumer the per-observation map
    # and weight built from the fit's armspec_blocks (gcol33/tulpaObs#95). A
    # non-intercept (slope) field carries its covariate column name to look up.
    wt <- if (isTRUE(meta[[b]]$is_intercept)) NULL else meta[[b]]$column_name
    list(z = z, amp_occ = amp_occ, amp_pos = amp_pos, weight = wt)
  })

  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
       b = list(occ = b_occ, det = NULL, pos = b_pos),
       blocks = blocks, n_cells = n_nodes[1L])
}

# occu_cover (3-arm psi/p/pos): one or more ICAR fields stored as unit-variance
# z. The occupancy arm scales each block by `sigma`, the positive arm by
# `alpha * sigma`, the detection arm not at all (the field's `field_coef = 0` on
# p). A trend fit adds weighted blocks, each with its own per-cell weight column.
.tobs_joint_draws_occu_cover <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  positive <- object$positive %||% "lognormal"
  # Field nodes (cells), not occupancy units (sites): under group_var the shared
  # field is sized by the graph, so n_sites can exceed it. Each block holds
  # n_cells field entries; consumers map sites -> cells via model$site_cell.
  n_cells  <- object$model$n_cells %||% object$model$n_sites
  p        <- layout$p
  bstart   <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_det <- bstart[2L] + seq_len(p[2L])
  idx_pos <- bstart[3L] + seq_len(p[3L])
  starts  <- layout$field_starts %||% layout$phi_start
  n_field <- length(starts)
  field_idx <- lapply(starts, function(s0) s0 + seq_len(n_cells))

  idx   <- c(idx_occ, idx_det, idx_pos, unlist(field_idx))
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) {
    v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v
  }
  b_occ <- take(p[1L]); b_det <- take(p[2L]); b_pos <- take(p[3L])

  trend_cols <- object$trend_weights %||% object$trend_weight
  blocks <- lapply(seq_len(n_field), function(b) {
    z     <- take(n_cells)
    sigma <- .tobs_joint_amp(tg, cells, b, "sigma")
    alpha <- .tobs_joint_amp(tg, cells, b, "alpha")
    list(z = z, amp_occ = sigma, amp_pos = alpha * sigma,
         weight = if (b == 1L) NULL else trend_cols[[b - 1L]])
  })

  # Pos-arm dispersion: the `phi_pos` axis when it is integrated on the outer
  # grid (control$phi.grid.pos, or the latent path's sigma_u); otherwise the
  # dispersion the fit held FIXED in the cell-coupling spec. Falling back to a
  # bare 1 would score every spatial occu_cover fit at unit dispersion regardless
  # of the value the spec used (gcol33/tulpaObs#34).
  fixed_disp <- object$model$cover_pos_disp %||% 1
  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = fixed_disp),
       b = list(occ = b_occ, det = b_det, pos = b_pos),
       blocks = blocks, n_cells = n_cells)
}

# cover (2-arm occ/pos): a single shared field. Under the (sigma_occ, sigma_pos)
# reparam the occupancy arm scales the unit-variance field z by `sigma_occ`, the
# positive arm by `sigma_pos`. The field is stored as phi (ICAR / proper CAR) or
# phi + theta (BYM2); for BYM2 the unit-variance z is
#   z = sqrt(rho) * scale_factor * phi + sqrt(1 - rho) * theta
# reconstructed per draw from the draw's grid rho.
.tobs_joint_draws_cover <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  cn       <- colnames(tg)
  positive <- object$positive %||% "lognormal"
  p        <- layout$p
  bstart   <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_pos <- bstart[2L] + seq_len(p[2L])

  # Coupled multi-block path (intercept field + one or more SVC trend fields):
  # ICAR blocks under the (sigma, alpha) per-block copy convention, axes named
  # b<k>.sigma / b<k>.alpha. The occupancy arm scales block k by b<k>.sigma, the
  # positive arm by b<k>.alpha * b<k>.sigma; block 1 is the unweighted intercept
  # field, blocks 2.. carry the per-observation trend weight. (gcol33/tulpaObs#15)
  field_starts <- layout$field_starts
  n_field <- length(field_starts %||% integer(0))
  if (n_field > 1L) {
    n_cells <- object$n_cells %||% as.integer(field_starts[2L] - field_starts[1L])
    field_idx <- lapply(field_starts, function(s0) s0 + seq_len(n_cells))
    idx <- c(idx_occ, idx_pos, unlist(field_idx))
    D   <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
    cells <- attr(D, "cells")
    off <- 0L
    take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
    b_occ <- take(p[1L]); b_pos <- take(p[2L])
    trend_cols <- object$trend_weights %||% list(object$trend_weight)
    blocks <- lapply(seq_len(n_field), function(b) {
      z     <- take(n_cells)
      sigma <- .tobs_joint_amp(tg, cells, b, "sigma")
      alpha <- .tobs_joint_amp(tg, cells, b, "alpha")
      list(z = z, amp_occ = sigma, amp_pos = alpha * sigma,
           weight = if (b == 1L) NULL else trend_cols[[b - 1L]])
    })
    return(list(n = n, positive = positive, cells = cells,
                disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
                b = list(occ = b_occ, det = NULL, pos = b_pos),
                blocks = blocks, n_cells = n_cells))
  }

  phi_start   <- layout$phi_start
  theta_start <- layout$theta_start
  n_phi <- if (!is.null(theta_start)) as.integer(theta_start - phi_start)
           else as.integer(layout$n_x - phi_start)
  phi_idx   <- phi_start + seq_len(n_phi)
  has_theta <- !is.null(theta_start)
  theta_idx <- if (has_theta) theta_start + seq_len(as.integer(layout$n_x - theta_start))
               else integer(0)

  idx   <- c(idx_occ, idx_pos, phi_idx, theta_idx)
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) {
    v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v
  }
  b_occ <- take(p[1L]); b_pos <- take(p[2L])
  phi   <- take(n_phi)

  # Unit-variance field z. BYM2 mixes the structured (phi) and unstructured
  # (theta) components by the draw's grid rho; ICAR / CAR use phi directly.
  has_rho <- "rho" %in% cn
  if (has_theta && has_rho) {
    theta <- take(length(theta_idx))
    rho   <- as.numeric(tg[cells, "rho"])
    sf    <- as.numeric(attr(jf, "scale_factor") %||% 1.0)
    z <- sweep(phi, 1L, sqrt(pmax(rho, 0) + 1e-10) * sf, "*") +
         sweep(theta, 1L, sqrt(pmax(1 - rho, 0) + 1e-10), "*")
  } else {
    z <- phi
  }

  amp_occ <- .tobs_joint_amp(tg, cells, 1L, "sigma_occ")
  amp_pos <- .tobs_joint_amp(tg, cells, 1L, "sigma_pos")
  block <- list(z = z, amp_occ = amp_occ, amp_pos = amp_pos, weight = NULL)

  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
       b = list(occ = b_occ, det = NULL, pos = b_pos),
       blocks = list(block), n_cells = n_phi)
}

# Pointwise-log-likelihood consumer helpers for arm-specific cover() fits
# (gcol33/tulpaObs#65): the per-arm fields store no node map / weight on their
# bundle blocks (gcol33/tulpaObs#95), so the loglik consumer -- which runs over
# the fit's observations -- rebuilds them from `object$armspec_blocks` and hands
# them to `.tobs_joint_arm_eta` as `units` / `wfun`. The predict consumer instead
# supplies the newdata cell map and a newdata-column lookup. Every block on an
# arm shares one per-observation node map (they index the same graph), so the
# first block on the arm carries it.
#
# Per-arm per-observation node map; `seq_len(n_rows)` when the arm has no field
# (the amplitude check in `.tobs_joint_arm_eta` then never indexes it).
.tobs_armspec_obs_units <- function(object, slot, n_rows) {
  for (m in object$armspec_blocks) {
    if (isTRUE(m$slot == slot)) return(as.integer(m$idx_active))
  }
  seq_len(n_rows)
}

# Per-arm weight lookup (column name -> per-observation weight) for the trend /
# SVC blocks on this arm; NULL when the arm carries no weighted block.
.tobs_armspec_obs_wfun <- function(object, slot) {
  lut <- list()
  for (m in object$armspec_blocks) {
    if (!isTRUE(m$slot == slot) || isTRUE(m$is_intercept)) next
    w <- if (slot == 1L) m$weight_occ else m$weight_pos
    if (!is.null(w)) lut[[m$column_name]] <- as.numeric(w)
  }
  if (length(lut) == 0L) return(NULL)
  function(nm) {
    v <- lut[[nm]]
    if (is.null(v)) {
      stop("internal: arm-specific weight column '", nm,
           "' is absent from the fit's armspec_blocks.", call. = FALSE)
    }
    v
  }
}

# Per-arm linear predictor for every draw: [nrow(X) x n]. `X` is the arm's
# design (scaled, when the family autoscales), `arm` is "occ" / "det" / "pos",
# `units` maps each design row to its spatial unit (1..n_cells), and `wfun`
# resolves a weighted field block's per-cell covariate to a length-nrow(X)
# vector (only consulted for weighted blocks). The detection arm sees no field.
# This is the single source of truth for the field accumulation across both
# families and both consumers (prediction and pointwise log-likelihood).
#
# A block's arm membership is carried entirely by its per-arm amplitude
# (`amp_occ` / `amp_pos`): a shared field scales both arms, an arm-specific field
# (gcol33/tulpaObs#65) has zero amplitude on the arm it does not sit on (no
# cross-arm copy). Every block on an arm then reads the SAME consumer-supplied
# `units` / `wfun` -- predict() passes the newdata cell map and a newdata-column
# lookup; the pointwise-log-likelihood consumer passes the per-observation node
# map and weight. Blocks store neither map nor weight themselves, so an
# arm-specific field maps correctly at predict time (gcol33/tulpaObs#95).
.tobs_joint_arm_eta <- function(bundle, X, arm, units, wfun = NULL) {
  B <- bundle$b[[arm]]
  if (is.null(B)) {
    stop("Arm '", arm, "' is absent from this joint fit.", call. = FALSE)
  }
  eta <- tcrossprod(X, B)                       # [nrow x n]
  if (identical(arm, "det")) return(eta)        # detection arm carries no field
  is_occ <- identical(arm, "occ")
  for (blk in bundle$blocks) {
    amp <- if (is_occ) blk$amp_occ else blk$amp_pos
    # Zero amplitude on this arm = block not on this arm; skip before any node
    # indexing, so a block sized to its own graph is never over-indexed.
    if (all(amp == 0)) next
    z_unit <- blk$z[, units, drop = FALSE]      # [n x nrow]; `units` maps each
    contr  <- t(z_unit * amp)                   # design row to its field cell
    if (!is.null(blk$weight)) {
      # Weighted (trend / SVC) field: scale each row by its per-cell covariate,
      # resolved by the consumer (predict: from newdata; loglik: per-observation).
      if (is.null(wfun)) {
        stop("Weighted field block '", blk$weight,
             "' needs a per-cell weight lookup.", call. = FALSE)
      }
      contr <- contr * wfun(blk$weight)         # row c scaled by weight[c]
    }
    eta <- eta + contr
  }
  eta
}
