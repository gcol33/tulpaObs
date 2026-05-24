# =============================================================================
# Shared random-effect parameter layout, naming, and BLUP reconstruction.
#
# The NUTS sampler lays the random-effect parameters out TYPE-BLOCKED (see
# tulpa hmc_param_layout.cpp), after the fixed-effect + visit-level columns:
#   [ all terms' log_sigma(n_coefs) ]
#   [ all correlated terms' chol_raw(n_coefs*(n_coefs-1)/2) ]
#   [ all terms' z(n_groups * n_coefs, group-major: index (g-1)*n_coefs + c) ]
# with the non-centred group effects recovered as
#   b_{g,c} = sigma_c * (L %*% z_g)_c,   L = tanh-Cholesky(chol_raw)
# (L = I for an uncorrelated block). This module names those columns and
# reconstructs the per-group BLUP table `re_effects` consumed by
# ranef.tobs_fit(). The deterministic Laplace path builds the same
# `re_effects` shape directly in .tobs_re_param_block().
# =============================================================================

# Lower-triangular tanh-Cholesky factor from strictly-lower raw parameters
# (row-major), mirroring tulpa_priors_re.h. The diagonal is derived from the
# unit-norm row constraint.
.tobs_tanh_chol <- function(raw, k) {
  L <- matrix(0, k, k)
  idx <- 1L
  for (r in seq_len(k)) {
    s2 <- 0
    for (cc in seq_len(r - 1L)) {
      L[r, cc] <- tanh(raw[idx]); s2 <- s2 + L[r, cc]^2; idx <- idx + 1L
    }
    L[r, r] <- sqrt(max(1 - s2, 1e-10))
  }
  L
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
# number of columns before the RE block (fixed effects + visit-level). Returns
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
      L <- if (!is.null(chol_draws)) .tobs_tanh_chol(chol_draws[s, ], nc) else diag(nc)
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
