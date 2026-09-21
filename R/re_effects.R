# =============================================================================
# Shared random-effect parameter layout, naming, and BLUP reconstruction.
#
# The NUTS sampler lays the random-effect parameters out TYPE-BLOCKED (see
# tulpa hmc_param_layout.cpp): the engine's block right after the fixed-effect
# columns, and a likelihood's own random effects (occu() visit-level terms)
# among its extra parameters at the end, in the same order:
#   [ all terms' log_sigma(n_coefs) ]
#   [ all correlated terms' chol_raw(n_coefs*(n_coefs-1)/2) ]
#   [ all terms' z(n_groups * n_coefs, group-major: index (g-1)*n_coefs + c) ]
# with the non-centred group effects recovered as
#   b_{g,c} = sigma_c * (L %*% z_g)_c,   L = build_L_from_raw(chol_raw)
# (L = I for an uncorrelated block). This module names those columns and
# reconstructs the per-group BLUP table `re_effects` consumed by
# ranef.tobs_fit(). The deterministic Laplace path builds the same
# `re_effects` shape directly in .tobs_re_param_block().
# =============================================================================

# Lower-triangular correlation Cholesky factor from its strictly-lower raw
# parameters (row-major), through the engine's partial-correlation map
# (tulpa::build_L_from_raw), so the rebuilt effects use the L the sampler drew.
.tobs_re_chol_factor <- function(raw, k) {
  cpp_re_chol_factor(as.numeric(raw), as.integer(k))
}

# The random effects over visit rows as `cpp_occu_fit` reads them: per term the
# 0-based group of each unit-major visit row (-1 where the visit is absent), the
# visit-row design, the group count, the correlation flag and the SD prior
# scale. `design` is .tobs_re_design() of the visit-level terms.
.tobs_visit_re_spec <- function(design, re_list) {
  lapply(seq_along(design), function(t) {
    d <- design[[t]]
    group <- d$idx - 1L
    group[is.na(group)] <- -1L
    Z <- d$Z
    Z[is.na(Z)] <- 0
    list(group = as.integer(group), Z = Z, n_groups = d$n_groups,
         correlated = isTRUE(d$correlated),
         sigma_scale = as.numeric(re_list[[t]]$sigma_scale %||% 1))
  })
}

# Section sizes per term, used by both the namer and the reconstructor so the
# type-blocked offsets are computed in one place.
.tobs_re_section_sizes <- function(design) {
  ncs   <- vapply(design, function(d) as.integer(d$n_coefs), integer(1))
  ngs   <- vapply(design, function(d) as.integer(d$n_groups), integer(1))
  nchol <- vapply(design, function(d)
    if (isTRUE(d$correlated)) as.integer(d$n_coefs * (d$n_coefs - 1L) / 2L) else 0L,
    integer(1))
  list(ncs = ncs, ngs = ngs, nchol = nchol)
}

# NUTS RE column names in the engine's TYPE-BLOCKED order. `design` is the
# list from .tobs_re_design().
.tobs_re_nuts_param_names <- function(design) {
  sig <- character(0); chl <- character(0); zz <- character(0)
  for (d in design) {
    nc <- d$n_coefs; ng <- d$n_groups; g <- d$group_label
    sig <- c(sig, sprintf("log_sigma_%s_%s", g, d$coef_names))
    if (isTRUE(d$correlated)) {
      n_chol <- nc * (nc - 1L) / 2L
      if (n_chol > 0L) chl <- c(chl, sprintf("chol_%s_%d", g, seq_len(n_chol)))
    }
    z_nms <- character(ng * nc)
    for (gi in seq_len(ng)) for (c in seq_len(nc)) {
      z_nms[(gi - 1L) * nc + c] <- sprintf("z_%s_%s[%d]", g, d$coef_names[c], gi)
    }
    zz <- c(zz, z_nms)
  }
  c(sig, chl, zz)
}

# Reconstruct the per-group BLUP table from NUTS draws (marginalising over the
# posterior: b is reconstructed per draw, then summarised). `n_lead` is the
# number of columns before the RE block. Returns
# a named list of per-term data.frames with group/level/term/estimate/std.error.
.tobs_re_nuts_effects <- function(draws, design, n_lead) {
  sz <- .tobs_re_section_sizes(design)
  sig_base  <- n_lead
  chol_base <- n_lead + sum(sz$ncs)
  z_base    <- n_lead + sum(sz$ncs) + sum(sz$nchol)
  sc <- 0L; cc <- 0L; zc <- 0L
  re_effects <- list()
  for (t in seq_along(design)) {
    d <- design[[t]]; nc <- sz$ncs[t]; ng <- sz$ngs[t]; g <- d$group_label
    sig_cols  <- sig_base + sc + seq_len(nc)
    chol_cols <- if (sz$nchol[t] > 0L) chol_base + cc + seq_len(sz$nchol[t]) else integer(0)
    z_cols    <- z_base + zc + seq_len(ng * nc)

    sig_draws  <- exp(draws[, sig_cols, drop = FALSE])      # n_draws x nc
    z_draws    <- draws[, z_cols, drop = FALSE]             # n_draws x ng*nc
    chol_draws <- if (length(chol_cols)) draws[, chol_cols, drop = FALSE] else NULL
    nd <- nrow(draws)

    B <- array(0, dim = c(nd, ng, nc))
    for (s in seq_len(nd)) {
      L <- if (!is.null(chol_draws)) .tobs_re_chol_factor(chol_draws[s, ], nc) else diag(nc)
      sig_s <- sig_draws[s, ]
      for (gi in seq_len(ng)) {
        zg <- z_draws[s, (gi - 1L) * nc + seq_len(nc)]
        B[s, gi, ] <- sig_s * as.numeric(L %*% zg)
      }
    }
    est <- apply(B, c(2, 3), mean)            # ng x nc
    se  <- apply(B, c(2, 3), stats::sd)
    re_effects[[g]] <- data.frame(
      group = g,
      level = rep(seq_len(ng), times = nc),
      term  = rep(d$coef_names, each = ng),
      estimate = as.numeric(est),
      std.error = as.numeric(se),
      stringsAsFactors = FALSE)
    sc <- sc + nc; cc <- cc + sz$nchol[t]; zc <- zc + ng * nc
  }
  re_effects
}
