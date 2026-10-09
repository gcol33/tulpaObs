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

# Per-draw natural-scale random effects from the sampler's coordinates. For
# each term: the SDs (n_draws x nc), the correlations of a correlated block
# (n_draws x nc(nc-1)/2, pairs (ci, cj) with ci < cj, ci-major) read off the
# correlation Cholesky factor, and the group effects B (n_draws x ng x nc).
# `n_lead` is the number of columns before the RE block.
.tobs_re_nuts_blocks <- function(draws, design, n_lead) {
  sz <- .tobs_re_section_sizes(design)
  sig_base  <- n_lead
  chol_base <- n_lead + sum(sz$ncs)
  z_base    <- n_lead + sum(sz$ncs) + sum(sz$nchol)
  sc <- 0L; cc <- 0L; zc <- 0L
  nd <- nrow(draws)
  lapply(seq_along(design), function(t) {
    nc <- sz$ncs[t]; ng <- sz$ngs[t]
    sig_cols  <- sig_base + sc + seq_len(nc)
    chol_cols <- if (sz$nchol[t] > 0L) chol_base + cc + seq_len(sz$nchol[t]) else integer(0)
    z_cols    <- z_base + zc + seq_len(ng * nc)
    sc <<- sc + nc; cc <<- cc + sz$nchol[t]; zc <<- zc + ng * nc

    sig_draws <- exp(draws[, sig_cols, drop = FALSE])
    z_draws   <- draws[, z_cols, drop = FALSE]
    pairs <- if (nc > 1L && length(chol_cols))
      do.call(rbind, lapply(seq_len(nc - 1L), function(ci)
        cbind(ci, seq(ci + 1L, nc)))) else matrix(integer(0), 0L, 2L)
    cor_draws <- matrix(NA_real_, nd, nrow(pairs))
    B <- array(0, dim = c(nd, ng, nc))
    for (s in seq_len(nd)) {
      L <- if (length(chol_cols))
        .tobs_re_chol_factor(draws[s, chol_cols], nc) else diag(nc)
      if (nrow(pairs)) cor_draws[s, ] <- tcrossprod(L)[pairs]
      sig_s <- sig_draws[s, ]
      for (gi in seq_len(ng)) {
        zg <- z_draws[s, (gi - 1L) * nc + seq_len(nc)]
        B[s, gi, ] <- sig_s * as.numeric(L %*% zg)
      }
    }
    list(sigma = sig_draws, cor = cor_draws, pairs = pairs, B = B,
         cols = c(sig_cols, chol_cols, z_cols))
  })
}

# Reconstruct the per-group BLUP table from NUTS draws (marginalising over the
# posterior: b is reconstructed per draw, then summarised). Returns a named
# list of per-term data.frames with group/level/term/estimate/std.error.
.tobs_re_nuts_effects <- function(draws, design, n_lead,
                                  blocks = .tobs_re_nuts_blocks(draws, design, n_lead)) {
  re_effects <- list()
  for (t in seq_along(design)) {
    d <- design[[t]]; B <- blocks[[t]]$B
    ng <- dim(B)[2L]; nc <- dim(B)[3L]; g <- d$group_label
    est <- apply(B, c(2, 3), mean)            # ng x nc
    se  <- apply(B, c(2, 3), stats::sd)
    re_effects[[g]] <- .tobs_re_effects_table(d, as.numeric(est),
                                              as.numeric(se))
  }
  re_effects
}

# One term's per-group effect table (group / level / term / estimate /
# std.error; coefficient-major, groups in code order), the shape ranef()
# stacks. It carries the term's grouping expression and process as attributes,
# which predict(newdata = ) reads to add a row's group effect.
.tobs_re_effects_table <- function(d, estimate, std.error) {
  ng <- d$n_groups; nc <- length(d$coef_names)
  tab <- data.frame(
    group = d$group_label,
    level = rep(d$levels, times = nc),
    term  = rep(d$coef_names, each = ng),
    estimate = estimate, std.error = std.error,
    stringsAsFactors = FALSE)
  attr(tab, "group_expr") <- d$group_expr
  attr(tab, "process")    <- d$process
  tab
}

# The RE block of a NUTS fit on the natural scale, in the layout and names the
# Laplace path reports (.tobs_re_param_block()): per term `sigma_<g>_<coef>`,
# `cor_<g>_<ci>_<cj>` on a correlated block, then `re_<g>_<coef>[k]` per
# coefficient. It has as many columns as the sampler's
# log_sigma / chol / z block, so it replaces that block in place.
.tobs_re_nuts_natural <- function(draws, design, n_lead,
                                  blocks = .tobs_re_nuts_blocks(draws, design, n_lead)) {
  out <- list(); nms <- character(0)
  for (t in seq_along(design)) {
    d <- design[[t]]; bl <- blocks[[t]]; g <- d$group_label; cn <- d$coef_names
    ng <- dim(bl$B)[2L]; nc <- dim(bl$B)[3L]
    out <- c(out, list(bl$sigma))
    nms <- c(nms, sprintf("sigma_%s_%s", g, cn))
    if (nrow(bl$pairs)) {
      out <- c(out, list(bl$cor))
      nms <- c(nms, sprintf("cor_%s_%s_%s", g, cn[bl$pairs[, 1L]], cn[bl$pairs[, 2L]]))
    }
    for (c in seq_len(nc)) {
      out <- c(out, list(matrix(bl$B[, , c], nrow = nrow(draws))))
      nms <- c(nms, sprintf("re_%s_%s[%d]", g, cn[c], seq_len(ng)))
    }
  }
  M <- do.call(cbind, out)
  colnames(M) <- nms
  M
}


# =============================================================================
# One ranef() layout for every fit.
#
# Every branch of ranef.tobs_fit() -- a family handler, the per-term BLUP blocks
# a joint fit stores on `fit$re`, the `re_effects` tables of the single-species
# paths, and the engine's own table -- emits through .tobs_ranef_table(), so a
# consumer reads the same six columns on any family and engine:
#   arm       the linear predictor the effect enters (`psi`, `p`, `lambda`, ...;
#             `psi+p` for an effect shared across arms; NA when the fit does not
#             distinguish arms)
#   group     the grouping variable (`"species"` for a community fit's
#             per-species deviations)
#   level     the group level
#   term      the coefficient the effect shifts, `(Intercept)` for an intercept
#   estimate  the BLUP / posterior mean
#   std.error its posterior SD, NA where the fitter reports none
# =============================================================================

.TOBS_RANEF_COLS <- c("arm", "group", "level", "term", "estimate", "std.error")

.tobs_ranef_table <- function(estimate, level, term = "(Intercept)",
                              arm = NA_character_, group = NA_character_,
                              std.error = NA_real_) {
  n <- length(estimate)
  data.frame(arm       = rep_len(as.character(arm), n),
             group     = rep_len(as.character(group), n),
             level     = rep_len(as.character(level), n),
             term      = rep_len(as.character(term), n),
             estimate  = as.numeric(estimate),
             std.error = rep_len(as.numeric(std.error), n),
             stringsAsFactors = FALSE)
}

.tobs_ranef_empty <- function() .tobs_ranef_table(numeric(0), character(0))

# Stack per-term tables into one; an empty list is the zero-row table.
.tobs_ranef_stack <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(.tobs_ranef_empty())
  out <- do.call(rbind, unname(rows))
  rownames(out) <- NULL
  out
}

# The grouping variable's name: the term's grouping expression as written in
# the formula (`g`, `g:h`), else the fallback label.
.tobs_re_group_name <- function(group_expr, fallback) {
  if (is.null(group_expr)) return(as.character(fallback))
  paste(deparse(group_expr, width.cutoff = 500L), collapse = "")
}

# Arm name of a random-effect term from its process index into the fit's arms
# (`.tobs_fit_arms()`): empty = a visit-level term on the detection arm; several
# = one effect shared across those arms.
.tobs_ranef_arm <- function(object, process) {
  if (is.null(process)) return(NA_character_)
  arms <- tryCatch(names(.tobs_fit_arms(object)), error = function(e) character(0))
  idx <- if (length(process)) as.integer(process) else 2L
  if (!length(arms) || any(idx < 1L | idx > length(arms))) return(NA_character_)
  paste(arms[idx], collapse = "+")
}

# A per-term `re_effects` table (.tobs_re_effects_table()) in the ranef layout.
.tobs_ranef_from_effects <- function(object, tab) {
  .tobs_ranef_table(estimate  = tab$estimate, level = tab$level, term = tab$term,
                    arm       = .tobs_ranef_arm(object, attr(tab, "process")),
                    group     = .tobs_re_group_name(attr(tab, "group_expr"),
                                                    tab$group),
                    std.error = tab$std.error)
}

# A BLUP block a joint / sampled fit stores on `fit$re`: `blup` is a per-level
# vector (random intercept) or an [n_groups x n_coefs] matrix (random slopes),
# `blup_sd` its SD in the same shape, `levels` the group labels, `var` the
# grouping variable (`group_label` on a block that records only the label).
.tobs_ranef_from_re_block <- function(re) {
  bl  <- re[["blup"]]; bsd <- re[["blup_sd"]]
  if (is.null(bsd)) bsd <- NA_real_
  group <- re[["var"]] %||% re[["group_label"]] %||% NA_character_
  if (is.matrix(bl)) {
    cn  <- colnames(bl) %||% re[["coef_names"]] %||%
           paste0("coef", seq_len(ncol(bl)))
    lev <- re[["levels"]] %||% as.character(seq_len(nrow(bl)))
    .tobs_ranef_table(estimate = as.numeric(bl),
                      level = rep(lev, times = ncol(bl)),
                      term  = rep(cn, each = nrow(bl)),
                      arm = re[["arm"]] %||% NA_character_, group = group,
                      std.error = as.numeric(bsd))
  } else {
    lev <- re[["levels"]] %||% as.character(seq_along(bl))
    .tobs_ranef_table(estimate = as.numeric(bl), level = lev,
                      arm = re[["arm"]] %||% NA_character_, group = group,
                      std.error = as.numeric(bsd))
  }
}

# The engine's own ranef table (tulpa's `ranef.tulpa_fit()`: one `term` row per
# group level and coefficient, named `<group>[<level>]` / `<group>.<coef>[<level>]`
# from `re_layout`, with `estimate` and `sd`) in the ranef layout. The layout
# is walked in the engine's own order (levels outer, coefficients inner) so
# the names are never parsed; a table of another width keeps its names as
# levels.
.tobs_ranef_from_layout <- function(object, tab) {
  if (!is.data.frame(tab) || !nrow(tab)) return(.tobs_ranef_empty())
  rows <- lapply(object$re_layout %||% list(), function(rt) {
    cls <- rt$coef_labels %||% "(Intercept)"
    data.frame(group = rt$group_var %||% NA_character_,
               level = rep(as.character(rt$levels), each = length(cls)),
               term  = rep(cls, times = length(rt$levels)),
               stringsAsFactors = FALSE)
  })
  key <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(key) || nrow(key) != nrow(tab)) {
    key <- data.frame(group = NA_character_, level = as.character(tab$term),
                      term = "(Intercept)", stringsAsFactors = FALSE)
  }
  .tobs_ranef_table(estimate = tab$estimate, level = key$level, term = key$term,
                    group = key$group,
                    std.error = tab[["sd"]] %||% tab[["std.error"]] %||% NA_real_)
}
