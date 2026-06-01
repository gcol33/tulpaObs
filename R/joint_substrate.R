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
  if (n_arms == 3L) {
    .tobs_joint_draws_occu_cover(object, jf, layout, n)
  } else {
    .tobs_joint_draws_cover(object, jf, layout, n)
  }
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

  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
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

# Per-arm linear predictor for every draw: [nrow(X) x n]. `X` is the arm's
# design (scaled, when the family autoscales), `arm` is "occ" / "det" / "pos",
# `units` maps each design row to its spatial unit (1..n_cells), and `wfun`
# resolves a weighted field block's per-cell covariate to a length-nrow(X)
# vector (only consulted for weighted blocks). The detection arm sees no field.
# This is the single source of truth for the field accumulation across both
# families and both consumers (prediction and pointwise log-likelihood).
.tobs_joint_arm_eta <- function(bundle, X, arm, units, wfun = NULL) {
  B <- bundle$b[[arm]]
  if (is.null(B)) {
    stop("Arm '", arm, "' is absent from this joint fit.", call. = FALSE)
  }
  eta <- tcrossprod(X, B)                       # [nrow x n]
  if (identical(arm, "det")) return(eta)        # detection arm carries no field
  for (blk in bundle$blocks) {
    amp    <- if (identical(arm, "occ")) blk$amp_occ else blk$amp_pos
    z_unit <- blk$z[, units, drop = FALSE]      # [n x nrow]
    contr  <- t(z_unit * amp)                   # [nrow x n], draw d scaled by amp[d]
    if (!is.null(blk$weight)) {
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
