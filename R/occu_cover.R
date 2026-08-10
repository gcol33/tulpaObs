# =============================================================================
# occu_cover.R - Joint occupancy-detection + cover hurdle family.
#
# A per-cell mixture combining MacKenzie single-season occupancy with a
# positive-cover observation on every detected visit. Plays the role of
# N-mixture, but for vegetation cover instead of counts: cell-level latent
# presence z_i ~ Bernoulli(psi_i), per-visit detection y_ij | z_i = 1
# ~ Bernoulli(p_ij), and per-visit cover y_pos_ij | y_ij = 1 ~ f_pos (beta or
# lognormal) on a third linear predictor. z marginalises out in closed form
# (two states), so the marginal log-likelihood is exact and the fit is a
# direct Laplace on the packed parameter vector
# (beta_occ, beta_p, beta_pos, log_dispersion).
#
# Per-cell likelihood:
#
#   any_det_i : L_i = psi_i * prod_j h_ij
#   no_det_i  : L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)
#
#   h_ij = (1 - p_ij) * 1{y_ij = 0}
#        + p_ij       * f_pos(y_pos_ij; eta_pos_ij, dispersion) * 1{y_ij = 1}
#
# Reduces to occu() when f_pos is degenerate and to the plot-level cover
# hurdle when J = 1 and p = 1. v1 covers the non-spatial Laplace path; a
# shared spatial field across the occ and cover arms (the analogue of
# cover()'s nested-Laplace joint engine) is v2.
#
# Files this touches:
#   R/obs_families.R    - occu_cover(response = ) constructor
#   R/tobs.R            - dispatch switch + .tobs_family_methods entry
#   tests/testthat/     - test-occu-cover.R recovery test
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Validate the cover response against the positive-arm density and reduce it to
# the gated form the likelihood reads. Cover is meaningful only where the visit
# is detected AND the value is observed; a detected visit with a missing cover
# carries the NA sentinel (kept in the occupancy / detection arms, dropped from
# the cover density by the same is.finite gate the C++ spec applies), an
# undetected visit carries 0. `y_pos_num` and `pos_mask` are conformable -- the
# dense [n_sites x max_visits] grid or the compact per-valid-visit vector -- so
# both binders share this.
#
# Beta requires the open unit interval and lognormal a positive value; the
# identity-Gaussian arm (gcol33/tulpaObs#112) lives on an unbounded real scale,
# where every finite value is admissible and the is.finite gate is the whole
# check.
.occu_cover_validate_pos_values <- function(y_pos_num, pos_mask, positive) {
  cover_obs <- pos_mask & is.finite(y_pos_num)
  if (identical(positive, "beta")) {
    if (any(cover_obs & (y_pos_num <= 0 | y_pos_num >= 1)))
      stop("Beta positive arm requires 0 < y_pos < 1 at every detected visit ",
           "with an observed cover; clip with pmin(pmax(y_pos, eps), 1 - eps).",
           call. = FALSE)
  } else if (!identical(positive, "gaussian")) {
    if (any(cover_obs & (y_pos_num <= 0)))
      stop("Lognormal positive arm requires y_pos > 0 at every detected visit ",
           "with an observed cover.", call. = FALSE)
  }
  y_pos_num[pos_mask & !cover_obs] <- NA_real_
  y_pos_num[!pos_mask]             <- 0
  y_pos_num
}

# Bind the joint occupancy-cover model. The psi predictor is cell-level
# (X_occ, n_sites rows). The detection and positive-cover predictors are
# visit-level (X_p_visit, X_pos_visit, n_sites * max_visits rows in
# site-major order). Visits with NA y are masked out of the likelihood; the
# binder zero-fills NA design rows so matrix algebra stays defined.
.tobs_build_occu_cover <- function(occ_formula, det_formula, pos_formula,
                                   data, y, y_pos, positive,
                                   det_visit_formula = NULL,
                                   det_visit_data   = NULL,
                                   pos_visit_formula = NULL,
                                   pos_visit_data    = NULL) {
  if (!is.matrix(y) || !is.matrix(y_pos)) {
    stop("y and y_pos must be matrices (n_sites x max_visits).", call. = FALSE)
  }
  if (!all(dim(y) == dim(y_pos))) {
    stop("y and y_pos must have identical dimensions.", call. = FALSE)
  }
  .tobs_check_site_count(nrow(y), nrow(data), "rows")

  n_sites    <- nrow(y)
  max_visits <- ncol(y)

  # Detection arm carries 0/1 values; NA = visit not observed.
  y_int <- matrix(as.integer(y), n_sites, max_visits)
  valid <- !is.na(y_int)
  if (any(y_int[valid] != 0L & y_int[valid] != 1L)) {
    stop("y must contain only 0, 1, or NA (binary detection per visit).",
         call. = FALSE)
  }
  y_int[!valid] <- 0L

  # Cover arm: meaningful only where detected AND the cover value is observed.
  # A detected visit with a missing cover (y_pos = NA) stays in the detection /
  # occupancy arms but drops out of the Beta/lognormal cover factor
  # (missing-at-random cover): the likelihood gates the cover density on
  # `y == 1 & is.finite(cover)`, the same semantic as the C++ std::isfinite gate.
  y_pos_num <- .occu_cover_validate_pos_values(
    matrix(as.numeric(y_pos), n_sites, max_visits),
    valid & (y_int == 1L), positive)

  # Reject structured terms in v1 (spatial sharing across arms is v2).
  .occu_cover_reject_structured(occ_formula, "occupancy")
  .occu_cover_reject_structured(det_formula, "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  # Site-level occupancy design.
  X_occ <- stats::model.matrix(occ_formula, data)

  # Site-level detection / positive design (intercept + any site-level covs).
  # The visit-level path adds visit-varying covariates on top.
  X_det_site <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)

  X_det_visit <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                           n_sites, max_visits, arm = "detection")
  X_pos_visit <- .tobs_build_visit_X(pos_visit_formula, pos_visit_data,
                                           n_sites, max_visits, arm = "positive cover")

  det_coef_names <- colnames(X_det_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_det_visit)) det_coef_names <- c(det_coef_names, colnames(X_det_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  structure(list(
    model_type  = "occu_cover",
    positive    = positive,
    y           = y_int,
    y_pos       = y_pos_num,
    valid       = valid,
    n_sites     = n_sites,
    max_visits  = max_visits,
    X_occ       = X_occ,
    X_det_site  = X_det_site,
    X_pos_site  = X_pos_site,
    X_det_visit = X_det_visit,
    X_pos_visit = X_pos_visit,
    formulas    = list(occ = occ_formula, det = det_formula, pos = pos_formula,
                       det_visit = det_visit_formula, pos_visit = pos_visit_formula),
    data        = data,
    process_info = list(
      list(name = "psi", p = ncol(X_occ),
           coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = length(det_coef_names),
           coef_names = det_coef_names, link = "logit"),
      list(name = "pos", p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}

# ---------------------------------------------------------------------------
# Compact (ragged) data binder
# ---------------------------------------------------------------------------

# Bind the joint occupancy-cover model from compact (ragged) inputs: one row per
# VALID visit instead of a padded [n_sites x max_visits] grid. Produces the same
# `tobs_model` the dense binder does, except the visit-level fields are already
# compacted -- `y_det_visit` / `y_pos_visit` / `site_of_visit` are length-V and
# `X_det_visit` / `X_pos_visit` are V-row designs -- and it carries `ragged =
# TRUE`. The joint-coupled arm builder reads these directly (it would otherwise
# flatten the dense grid to exactly this), so the engine sees an identical
# `responses` list. No dense `y` / `valid` matrices are built, so there is no
# per-site visit cap: memory is O(V) = total observations, not O(n_sites x
# max_visits). Scoped to the joint nested-Laplace path (the only consumer of the
# compacted arms); other engines / the cell-aggregated cover path read the dense
# fields and are gated off in .dispatch_occu_cover.
.tobs_build_occu_cover_ragged <- function(occ_formula, det_formula, pos_formula,
                                          data, y_ragged, y_pos_values, positive,
                                          det_visit_formula = NULL,
                                          det_visit_data   = NULL,
                                          pos_visit_formula = NULL,
                                          pos_visit_data    = NULL) {
  if (!inherits(y_ragged, "tobs_ragged")) {
    stop("ragged occu_cover binder expects a tobs_ragged `y`.", call. = FALSE)
  }
  n_sites        <- y_ragged$n_sites
  max_visits     <- y_ragged$max_visits
  site_of_visit  <- y_ragged$site
  n_visits_valid <- length(site_of_visit)
  .tobs_check_site_count(n_sites, nrow(data), "rows")

  # Detection arm carries 0/1 per valid visit.
  y_det_visit <- as.integer(y_ragged$values)
  if (anyNA(y_det_visit) || any(!(y_det_visit %in% c(0L, 1L)))) {
    stop("compact occu_cover: detection `y` must be 0/1 at every valid visit ",
         "(no NA -- a compact row is, by construction, a sampled visit).",
         call. = FALSE)
  }
  if (length(y_pos_values) != n_visits_valid) {
    stop(sprintf(paste0(
      "compact occu_cover: y_pos has %d rows but the detection response has %d ",
      "valid visits. The occurrence and cover tobs_data(compact = TRUE) calls ",
      "must use the same df / site / visit so they align."),
      length(y_pos_values), n_visits_valid), call. = FALSE)
  }

  # Cover meaningful only where detected AND observed. A detected visit with a
  # missing cover (NA) stays in the detection / occupancy arms but drops out of
  # the cover density (missing-at-random cover); the spec gates the cover term on
  # y_det == 1 & is.finite(cover), the same semantic as the C++ std::isfinite gate.
  y_pos_num <- .occu_cover_validate_pos_values(
    as.numeric(y_pos_values), y_det_visit == 1L, positive)

  .occu_cover_reject_structured(occ_formula, "occupancy")
  .occu_cover_reject_structured(det_formula, "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  # Site-level designs (n_sites rows) -- identical to the dense path.
  X_occ      <- stats::model.matrix(occ_formula, data)
  X_det_site <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)

  # Visit-level designs built directly on the V valid rows (max_per_unit = NULL
  # -> the compact signal that skips the padded-grid row check).
  X_det_visit <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                     n_sites, NULL, arm = "detection")
  X_pos_visit <- .tobs_build_visit_X(pos_visit_formula, pos_visit_data,
                                     n_sites, NULL, arm = "positive cover")

  det_coef_names <- colnames(X_det_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_det_visit)) det_coef_names <- c(det_coef_names, colnames(X_det_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  structure(list(
    model_type     = "occu_cover",
    positive       = positive,
    ragged         = TRUE,
    site_of_visit  = site_of_visit,
    y_det_visit    = y_det_visit,
    y_pos_visit    = y_pos_num,
    n_visits_valid = n_visits_valid,
    n_sites        = n_sites,
    max_visits     = max_visits,
    X_occ          = X_occ,
    X_det_site     = X_det_site,
    X_pos_site     = X_pos_site,
    X_det_visit    = X_det_visit,
    X_pos_visit    = X_pos_visit,
    formulas       = list(occ = occ_formula, det = det_formula, pos = pos_formula,
                          det_visit = det_visit_formula, pos_visit = pos_visit_formula),
    data           = data,
    process_info = list(
      list(name = "psi", p = ncol(X_occ),
           coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = length(det_coef_names),
           coef_names = det_coef_names, link = "logit"),
      list(name = "pos", p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}

# Reject formulas containing structured terms (bym2, icar, car, gp, etc.) on the
# non-spatial backends. A field shared across the psi and cover arms is fitted by
# the nested-Laplace joint engine, which resolves its own terms through
# .occu_cover_spatial_fields(); this guard is what keeps a structured term from
# being silently dropped on the routes that cannot carry one.
.occu_cover_reject_structured <- function(formula, arm) {
  if (is.null(formula)) return(invisible(NULL))
  labs <- attr(stats::terms(formula), "term.labels")
  structured <- c("bym2", "icar", "car", "car_proper", "gp", "spde",
                  "multiscale_gp", "re", "temporal", "svc", "latent", "copy")
  hits <- character(0)
  for (lab in labs) {
    fn <- tryCatch(as.character(as.call(parse(text = lab)[[1]])[[1]]),
                   error = function(e) NA_character_)
    if (!is.na(fn) && fn %in% structured) hits <- c(hits, fn)
  }
  if (length(hits) > 0L) {
    stop(sprintf(paste0(
      "occu_cover() v1 does not support structured terms (%s) on the %s arm. ",
      "Shared spatial / temporal / RE fields across the three arms is v2; ",
      "for now use a plain fixed-effects formula on each."),
      paste(unique(hits), collapse = ", "), arm), call. = FALSE)
  }
  invisible(NULL)
}

# Parse an observation-arm (detection / positive-cover) formula for random
# effects (gcol33/tulpaObs#102, #103). lme4 bars are desugared to re() first, so
# `(1 | g)`, `(x | g)`, `(x || g)`, `(0 + x | g)`, crossed `(1 | g) + (1 | h)`,
# and nested `(1 | g/h)` (-> re(g) + re(g:h)) all arrive as a list of re() terms.
# Returns NULL when the arm carries no random effect (the formula is then used
# unchanged downstream); otherwise a list with the fixed-effects formula (every
# re() term stripped, for the design build), the parsed term specs (one per re()
# term, in formula order), and `has_slope` (TRUE if any term is a random slope --
# the joint-engine slope blocks are gated on gcol33/tulpa#114). Each spec keeps
# its grouping expression, slope covariate, and intercept / correlated flags;
# group codes and the per-row design are resolved later against `data` / `visits`
# once the model's `valid` mask is built. Other structured terms (icar(), gp(),
# temporal(), ...) error; copy() is stripped separately.
.occu_cover_obs_re_term_spec <- function(label) {
  e    <- str2lang(label)
  args <- as.list(e)[-1L]
  group_expr <- if (!is.null(args$group)) args$group else args[[1L]]
  flag <- function(nm, default) {
    if (is.null(args[[nm]])) return(default)
    tryCatch(isTRUE(eval(args[[nm]])), error = function(...) default)
  }
  type <- tryCatch(if (!is.null(args$type)) eval(args$type) else "intercept",
                   error = function(...) "intercept")
  list(group_expr = group_expr, vars = all.vars(group_expr),
       type = type, covariate = args$covariate,
       intercept = flag("intercept", TRUE), correlated = flag("correlated", TRUE),
       term_label = label)
}

.occu_cover_obs_re_parse <- function(formula, arm) {
  if (is.null(formula)) return(NULL)
  f    <- .tobs_desugar_bars(formula)
  tt   <- stats::terms(f, keep.order = TRUE)
  labs <- attr(tt, "term.labels")
  # `re` is extracted here; `copy` is the cross-arm coupling, stripped separately
  # by .occu_cover_extract_pos_copies(); every other structured term is rejected.
  reg  <- setdiff(.tobs_term_names(), c("re", "copy"))
  re_labels <- character(0)
  for (lab in labs) {
    e    <- tryCatch(str2lang(lab), error = function(...) NULL)
    head <- if (is.call(e) && is.symbol(e[[1L]])) as.character(e[[1L]])
            else NA_character_
    if (identical(head, "re")) {
      re_labels <- c(re_labels, lab)
    } else if (!is.na(head) && head %in% reg) {
      stop(sprintf(paste0(
        "occu_cover(): structured term `%s()` is not supported on the %s arm; ",
        "random effects (`(1 | g)`, `(x | g)`, crossed / nested) are ",
        "(gcol33/tulpaObs#102, #103)."), head, arm), call. = FALSE)
    }
  }
  if (length(re_labels) == 0L) return(NULL)
  specs     <- lapply(re_labels, .occu_cover_obs_re_term_spec)
  fe_labels <- setdiff(labs, re_labels)
  fe <- stats::reformulate(
    termlabels = if (length(fe_labels)) fe_labels else "1",
    intercept  = as.logical(attr(tt, "intercept")))
  environment(fe) <- environment(formula)
  list(fe = fe, terms = specs,
       has_slope = any(vapply(specs, function(s) identical(s$type, "slope"),
                              logical(1))))
}

# Evaluate an obs-arm RE expression (a grouping factor or a slope covariate) to
# a flat value in site-major order (site 1 visits 1..max_visits, site 2, ...),
# matching `as.logical(t(model$valid))`. When all the expression's variables live
# in `data` (one row per site) the result is broadcast across each site's visits;
# when they live in `visits` (one row per site-visit) it is taken as-is. A
# matrix-valued expression (a cbind() of slope covariates) keeps its columns.
.occu_cover_obs_flat_eval <- function(expr, data, visit_df, n_sites, max_visits,
                                      arm) {
  vars   <- all.vars(expr)
  n_flat <- n_sites * max_visits
  if (length(vars) && all(vars %in% names(data))) {
    val <- eval(expr, envir = data)
    if (is.matrix(val)) {
      if (nrow(val) != n_sites)
        stop(sprintf("occu_cover(): a %s RE expression resolved to %d rows but there are %d sites.",
                     arm, nrow(val), n_sites), call. = FALSE)
      return(val[rep(seq_len(n_sites), each = max_visits), , drop = FALSE])
    }
    if (length(val) != n_sites)
      stop(sprintf("occu_cover(): a %s RE expression resolved to %d values but there are %d sites.",
                   arm, length(val), n_sites), call. = FALSE)
    return(rep(val, each = max_visits))
  }
  if (!is.null(visit_df) && length(vars) && all(vars %in% names(visit_df))) {
    val <- eval(expr, envir = visit_df)
    if (is.matrix(val)) {
      if (nrow(val) != n_flat)
        stop(sprintf("occu_cover(): a %s RE expression resolved to %d rows but expected %d (n_sites x max_visits).",
                     arm, nrow(val), n_flat), call. = FALSE)
      return(val)
    }
    if (length(val) != n_flat)
      stop(sprintf("occu_cover(): a %s RE expression resolved to %d values but expected %d (one per site-visit).",
                   arm, length(val), n_flat), call. = FALSE)
    return(val)
  }
  stop(sprintf(paste0(
    "occu_cover(): the %s random-effect variable(s) %s were not found in `data` ",
    "(one value per site) or `visits` (one per site-visit)."),
    arm, paste(vars, collapse = ", ")), call. = FALSE)
}

# Resolve the parsed obs-arm RE term specs to a per-term design list. Each term:
#   codes_flat : per-(site, visit) integer group codes, site-major, 0 for a
#                padded / dropped / unseen-level row (the engine's no-RE sentinel);
#                those rows fall out under the joint builder's `keep` subset.
#   levels     : sorted observed factor levels (BLUP naming + predict matching).
#   n_groups   : number of levels.
#   var        : grouping-variable label (predict matches newdata[[var]]).
#   type / correlated / has_intercept / n_coefs / coef_names : block shape.
#   Z          : per-row design weights [n_flat x n_coefs] for a slope term
#                (intercept column all-ones, slope columns the covariate values);
#                NULL for a plain random intercept (weight 1, the iid block).
# Factor levels are taken from the OBSERVED visits only, so an NA-padded or
# never-observed level adds no group.
.occu_cover_obs_re_design <- function(re_parse, data, visit_df, valid,
                                      n_sites, max_visits, arm) {
  .occu_cover_obs_re_terms(
    re_parse,
    function(expr) .occu_cover_obs_flat_eval(expr, data, visit_df,
                                             n_sites, max_visits, arm),
    as.logical(t(valid)), arm)
}

# Per-term design shared by the dense and compact routes. `eval_flat(expr)`
# returns one value (or matrix row) per design row in the arm's row order --
# site-major over the padded grid, or one row per valid visit -- and `observed`
# masks the rows carrying data (TRUE recycles to "every row" under the compact
# layout, where by construction every row is a sampled visit). Everything the
# two routes share -- level extraction from observed rows only, the 0 sentinel
# for padded / unseen-level rows, slope-matrix assembly and standardization --
# lives here; only the evaluator and the mask differ.
.occu_cover_obs_re_terms <- function(re_parse, eval_flat, observed, arm) {
  lapply(re_parse$terms, function(spec) {
    g_flat <- eval_flat(spec$group_expr)
    var <- paste(spec$vars, collapse = ":")
    lev <- sort(unique(as.character(g_flat[observed])))
    if (length(lev) < 2L) {
      stop(sprintf(paste0(
        "occu_cover(): the %s random-effect grouping `%s` has %d level(s) among ",
        "the observed visits; a random effect needs at least 2 groups."),
        arm, var, length(lev)), call. = FALSE)
    }
    codes <- match(as.character(g_flat), lev)
    codes[is.na(codes)] <- 0L

    if (identical(spec$type, "slope")) {
      cov_expr <- spec$covariate
      # A bare-string covariate (`re(g, covariate = "x")`) names columns; wrap it
      # as cbind(<symbols>) so the same flat evaluator builds the slope matrix.
      cov_chr <- tryCatch(eval(cov_expr), error = function(...) NULL)
      if (is.character(cov_chr)) {
        cov_expr <- as.call(c(list(as.name("cbind")), lapply(cov_chr, as.name)))
      }
      Xs <- as.matrix(eval_flat(cov_expr)); storage.mode(Xs) <- "double"
      if (is.null(colnames(Xs)))
        colnames(Xs) <- paste0("slope", seq_len(ncol(Xs)))
      # Standardize each slope covariate to unit SD (over the observed visits), so
      # the per-coefficient RE variance is on an O(1) scale the fixed Sigma grid
      # brackets regardless of the covariate's raw scale; a constant column (sd 0)
      # keeps scale 1. The fit runs on the scaled column and the reported slope
      # BLUP / SD are divided back by the scale, so the random slope is on the
      # covariate's natural units (correlation is scale-free). Mirrors the
      # fixed-effect design autoscaling.
      slope_scale <- apply(Xs, 2L, function(col) {
        s <- stats::sd(col[observed]); if (!is.finite(s) || s <= 0) 1 else s
      })
      Xs <- sweep(Xs, 2L, slope_scale, "/")
      has_int <- isTRUE(spec$intercept)
      Z <- if (has_int) cbind(`(Intercept)` = 1, Xs) else Xs
      coef_names  <- colnames(Z)
      n_coefs     <- ncol(Z)
      correlated  <- isTRUE(spec$correlated) && n_coefs > 1L
      coef_scales <- if (has_int) c(1, slope_scale) else slope_scale
    } else {
      Z <- NULL; coef_names <- "(Intercept)"; n_coefs <- 1L
      correlated <- FALSE; has_int <- TRUE; coef_scales <- 1
    }
    list(codes_flat = as.integer(codes), levels = lev, n_groups = length(lev),
         var = var, type = spec$type, n_coefs = n_coefs,
         coef_names = coef_names, correlated = correlated,
         has_intercept = has_int, Z = Z, coef_scales = as.numeric(coef_scales),
         term_label = spec$term_label)
  })
}

# Compact (ragged) counterpart of .occu_cover_obs_flat_eval: evaluate an RE
# expression to one value per VALID visit (length V, in site_of_visit order). A
# site-level variable (found in `data`) is broadcast to its visits via
# site_of_visit; a visit-level variable is read from the V-row `visit_df`
# directly. No padded grid, so no n_sites * max_visits layout.
.occu_cover_obs_flat_eval_ragged <- function(expr, data, visit_df, site_of_visit,
                                             arm) {
  vars <- all.vars(expr)
  n_sites <- if (is.null(data)) 0L else nrow(data)
  if (length(vars) && !is.null(data) && all(vars %in% names(data))) {
    val <- eval(expr, envir = data)
    if (is.matrix(val)) {
      if (nrow(val) != n_sites)
        stop(sprintf("occu_cover(): a %s RE expression resolved to %d rows but there are %d sites.",
                     arm, nrow(val), n_sites), call. = FALSE)
      return(val[site_of_visit, , drop = FALSE])
    }
    if (length(val) != n_sites)
      stop(sprintf("occu_cover(): a %s RE expression resolved to %d values but there are %d sites.",
                   arm, length(val), n_sites), call. = FALSE)
    return(val[site_of_visit])
  }
  if (!is.null(visit_df) && length(vars) && all(vars %in% names(visit_df))) {
    return(eval(expr, envir = visit_df))            # already one row per valid visit
  }
  stop(sprintf(paste0(
    "occu_cover(): the %s random-effect variable(s) %s were not found in `data` ",
    "(one value per site) or `visits` (one per site-visit)."),
    arm, paste(vars, collapse = ", ")), call. = FALSE)
}

# Compact (ragged) counterpart of .occu_cover_obs_re_design: the per-term design
# built over the V valid visits directly (every row observed, so no valid mask
# and no 0-padding). `codes_flat` is length V in site_of_visit order, which is
# exactly the order the joint-coupled arm builder reads under the ragged branch
# (keep = 1..V), so the identical list shape feeds the same engine path.
.occu_cover_obs_re_design_ragged <- function(re_parse, data, visit_df,
                                             site_of_visit, arm) {
  .occu_cover_obs_re_terms(
    re_parse,
    function(expr) .occu_cover_obs_flat_eval_ragged(expr, data, visit_df,
                                                    site_of_visit, arm),
    TRUE, arm)
}

# Resolve the parsed observation-arm RE terms against the built model and attach
# the per-term design lists (`model$re_det` / `model$re_pos`). The grouping codes
# need the model's `valid` mask (dense) or `site_of_visit` (compact), so this runs
# after the model build; every engine that hosts an observation-arm RE -- the
# grid-integrated joint fit and the sampler -- reads the SAME designs from here.
# The positive-cover RE is per visit, so it needs per-visit cover: a
# cell-aggregated cover arm has one row per unit and no per-visit grouping.
.occu_cover_attach_obs_re <- function(model, det_re_parse, pos_re_parse,
                                      data, det_visits, pos_visits) {
  design <- function(re_parse, visits, arm) {
    if (isTRUE(model$ragged))
      .occu_cover_obs_re_design_ragged(re_parse, data, visits,
                                       model$site_of_visit, arm)
    else
      .occu_cover_obs_re_design(re_parse, data, visits, model$valid,
                                model$n_sites, model$max_visits, arm)
  }
  if (!is.null(det_re_parse))
    model$re_det <- design(det_re_parse, det_visits, "detection")
  if (!is.null(pos_re_parse)) {
    if (!identical(model$cover_aggregate %||% "none", "none")) {
      stop("occu_cover(): a random effect on the positive-cover arm needs ",
           "per-visit cover (cover_aggregate = \"none\"); it cannot map onto ",
           "cell-aggregated cover rows (one per unit).", call. = FALSE)
    }
    model$re_pos <- design(pos_re_parse, pos_visits, "positive cover")
  }
  model
}

# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------

# Negative log-posterior at packed parameter vector
#   par = c(beta_occ, beta_p, beta_pos, log_dispersion)
# evaluated against the bound model. Gaussian prior with diagonal (mean, prec)
# aligned with par; flat prior when pprec == 0.
.tobs_occu_cover_nlp <- function(par, model, pmean, pprec) {
  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p

  bo  <- par[seq_len(p_occ)]
  bp  <- par[p_occ + seq_len(p_p)]
  bpos<- par[p_occ + p_p + seq_len(p_pos)]
  log_disp <- par[length(par)]

  eta <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  ll <- sum(.occu_cover_site_ll(model, eta$psi, eta$p_mat, eta$ep_mat, log_disp))
  penalty <- 0.5 * sum(pprec * (par - pmean)^2)
  -ll + penalty
}

# Build the three arm linear predictors from an arm-split coefficient triple.
# `bo` is the occupancy (psi) coefficient vector; `bp` / `bpos` are the
# detection / cover coefficient vectors, each packed as the site-level block
# followed by the optional visit-level block exactly as the fitter stacks them.
# Returns `psi` (length n_sites), `p_mat` (the per-visit detection probability,
# [n_sites x max_visits]), and `ep_mat` (the per-visit cover linear predictor on
# its link scale, [n_sites x max_visits]). The eta does not depend on the
# response, so the single-species fit and the community per-species marginal
# share one builder (single source of truth for the occu_cover predictors).
.occu_cover_eta_from_par <- function(model, bo, bp, bpos) {
  cl <- .tobs_clamp_eta
  n_sites    <- model$n_sites
  max_visits <- model$max_visits

  psi <- stats::plogis(cl(as.numeric(model$X_occ %*% bo)))

  bp_site <- bp[seq_len(ncol(model$X_det_site))]
  bp_visit <- if (!is.null(model$X_det_visit)) {
    bp[ncol(model$X_det_site) + seq_len(ncol(model$X_det_visit))]
  } else numeric(0)
  bpos_site <- bpos[seq_len(ncol(model$X_pos_site))]
  bpos_visit <- if (!is.null(model$X_pos_visit)) {
    bpos[ncol(model$X_pos_site) + seq_len(ncol(model$X_pos_visit))]
  } else numeric(0)

  # X_*_site is n_sites x ?, X_*_visit is (n_sites * max_visits) x ? in
  # site-major order. Broadcast the site-level eta across visits.
  eta_p_site <- as.numeric(model$X_det_site %*% bp_site)
  p_mat <- matrix(eta_p_site, n_sites, max_visits)
  if (length(bp_visit)) {
    eta_p_visit <- as.numeric(model$X_det_visit %*% bp_visit)
    p_mat <- p_mat + matrix(eta_p_visit, n_sites, max_visits, byrow = TRUE)
  }
  p_mat <- stats::plogis(cl(p_mat))

  eta_pos_site <- as.numeric(model$X_pos_site %*% bpos_site)
  ep_mat <- matrix(eta_pos_site, n_sites, max_visits)
  if (length(bpos_visit)) {
    eta_pos_visit <- as.numeric(model$X_pos_visit %*% bpos_visit)
    ep_mat <- ep_mat + matrix(eta_pos_visit, n_sites, max_visits, byrow = TRUE)
  }

  list(psi = psi, p_mat = p_mat, ep_mat = ep_mat)
}

# Detected occupancy units and their per-unit detected-visit cover values. The
# single source of truth for "which units carry a cover observation and what
# covers they hold", shared by the joint-coupled arm builder (which collapses
# these to one aggregated / latent pos-arm row per unit) and the pointwise
# log-likelihood (which must score the cover term at the same granularity the
# fitter optimised, gcol33/tulpaObs#34). `pos_site` indexes the occupancy units
# with at least one detection; `vals[[k]]` is that unit's detected covers.
.occu_cover_unit_cover <- function(model) {
  # Cover-observed visits (detected AND finite cover); a detected visit with a
  # missing cover (NA) contributes no cover term. A unit enters `pos_site` only
  # if it has at least one observed cover.
  obs_mat  <- model$valid & (model$y == 1L) & is.finite(model$y_pos)
  pos_site <- which(rowSums(obs_mat) > 0L)
  vals <- lapply(pos_site, function(i) as.numeric(model$y_pos[i, obs_mat[i, ]]))
  list(pos_site = pos_site, vals = vals)
}

# Positive-arm log-density of cover value(s) `y` at cover predictor `eta`
# (link scale) and dispersion `disp` (lognormal residual SD or beta precision).
# Vectorised over y / eta (and matrices), so the per-visit and the per-unit
# aggregated cover terms read one formula. Beta clamps the predictor before the
# logistic; lognormal uses the raw predictor (matching the historical kernels).
# Integer code for the positive-arm policy, shared with the C++ dispatch in
# src/occu_coupling_shared.h (pos_log_density / pos_grad_eta / pos_grad_logdisp)
# and the simulate draws: 0 = lognormal, 3 = beta, 4 = gaussian.
.occu_cover_pos_code <- function(positive) {
  switch(positive, beta = 3L, gaussian = 4L, lognormal = 0L, 0L)
}

# Positive-arm log-density. Single source of truth with the fit kernel
# src/occu_coupling_shared.h::pos_log_density: no eta clamp (the coupling kernel
# does not clamp -- at |eta| <= the converged mode the +-30 clamp never bit, and
# dropping it makes the WAIC / LOO pointwise density the number the model was fit
# with) and .tobs_log_safe on every log so the density is finite at the cover
# boundary (cover exactly 0 or 1) instead of -Inf (gcol33/tulpaObs#133).
.occu_cover_pos_logdens <- function(y, eta, disp, positive) {
  if (identical(positive, "beta")) {
    mu <- stats::plogis(eta)
    a  <- mu * disp
    b  <- (1 - mu) * disp
    lgamma(disp) - lgamma(a) - lgamma(b) +
      (a - 1) * .tobs_log_safe(y) + (b - 1) * .tobs_log_safe(1 - y)
  } else if (identical(positive, "gaussian")) {
    # Identity-Gaussian arm (gcol33/tulpaObs#112): residual on the raw response,
    # no change-of-variable Jacobian. mu = eta.
    -.tobs_log_safe(disp) - 0.5 * log(2 * pi) -
      0.5 * ((y - eta) / disp)^2
  } else {
    # lognormal: residual on log-cover, with the -log(y) Jacobian.
    ly <- .tobs_log_safe(y)
    -ly - .tobs_log_safe(disp) - 0.5 * log(2 * pi) -
      0.5 * ((ly - eta) / disp)^2
  }
}

# Closed-form per-unit lognormal latent-cover marginal log M_i (the compound-
# symmetry integral over the per-unit cover RE u_i ~ N(0, sigma_u^2) with fixed
# within-unit residual SD `disp2`). Mirrors src/occu_cover_latent.h::LognormalLatent
# exactly: Sigma = a I + b 11', a = disp2^2, b = sigma_u^2, plus the lognormal
# change-of-variables Jacobian -sum log y. `eta` is the unit-level predictor.
.occu_cover_latent_lognormal_logm <- function(vals, eta, disp2, sigma_u) {
  a <- disp2^2
  b <- sigma_u^2
  vapply(seq_along(vals), function(i) {
    v  <- vals[[i]]; m <- length(v); ly <- log(v)
    t1 <- sum(ly); t2 <- sum(ly * ly)
    denom <- a + m * b
    s1  <- t1 - m * eta[i]
    s2c <- t2 - 2 * eta[i] * t1 + m * eta[i]^2
    quad   <- if (a > 0) (s2c - (b / denom) * s1 * s1) / a else 0
    logdet <- (m - 1) * log(a) + log(denom)
    -0.5 * m * log(2 * pi) - 0.5 * logdet - 0.5 * quad - t1
  }, numeric(1))
}

# Probabilist Gauss-Hermite nodes / weights (weight exp(-x^2/2)/sqrt(2 pi),
# weights sum to 1) by Golub-Welsch on the symmetric Jacobi matrix. Dependency-
# free; integrates E_{N(0,1)}[g] ~ sum_k w_k g(z_k).
.gauss_hermite_prob <- function(n) {
  if (n <= 1L) return(list(nodes = 0, weights = 1))
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i)
  J[cbind(i + 1L, i)] <- sqrt(i)
  e   <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = (e$vectors[1L, ])[ord]^2)
}

# Per-unit beta latent-cover marginal log M_i = log integral of
# prod_j Beta(y_ij | sigmoid(eta + u), phi) * N(u; 0, sigma_u^2) du, by
# Gauss-Hermite against the cover-RE prior. Same marginal as
# src/occu_cover_latent.h::BetaLatent (non-adaptive quadrature of the same
# integral). `eta` is the unit-level predictor.
.occu_cover_latent_beta_logm <- function(vals, eta, phi, sigma_u, n_quad) {
  gh <- .gauss_hermite_prob(max(as.integer(n_quad), 15L))
  z  <- gh$nodes
  lw <- log(gh$weights)
  vapply(seq_along(vals), function(i) {
    v  <- vals[[i]]
    lt <- vapply(seq_along(z), function(k) {
      ell <- sum(.occu_cover_pos_logdens(v, eta[i] + sigma_u * z[k], phi, "beta"))
      lw[k] + ell
    }, numeric(1))
    mx <- max(lt)
    mx + log(sum(exp(lt - mx)))
  }, numeric(1))
}

# Per-unit cover contribution to the marginal log-likelihood (length n_sites,
# zero for units with no detection). For `cover_aggregate = "none"` this is the
# per-visit sum of the positive-arm density at detected visits; for "mean" /
# "median" it is one density at the per-unit aggregated cover; for "latent" it is
# the per-unit cover-RE marginal. `ep_mat` is the [n_sites x max_visits] cover
# predictor; under aggregation the cover design is unit-level so the predictor is
# constant across a unit's visits (column 1 is the unit value).
.occu_cover_cover_term <- function(model, ep_mat, log_disp, units = NULL) {
  n_sites  <- model$n_sites
  positive <- model$positive %||% "lognormal"
  mode     <- model$cover_aggregate %||% "none"

  if (identical(mode, "none")) {
    # Cover density at detected visits with an observed cover (missing-at-random
    # cover drops out); the NA-cover cells score 0.
    pos_mask <- model$valid & (model$y == 1L) & is.finite(model$y_pos)
    log_f_pos <- matrix(0, n_sites, model$max_visits)
    if (any(pos_mask)) {
      log_f_pos[pos_mask] <- .occu_cover_pos_logdens(
        model$y_pos[pos_mask], ep_mat[pos_mask], exp(log_disp), positive)
    }
    return(rowSums(log_f_pos))
  }

  if (is.null(units)) units <- .occu_cover_unit_cover(model)
  out <- numeric(n_sites)
  ps  <- units$pos_site
  if (length(ps) == 0L) return(out)
  eta <- ep_mat[ps, 1L]
  if (identical(mode, "latent")) {
    sigma_u <- exp(log_disp)
    disp2   <- model$cover_latent_disp2
    out[ps] <- if (identical(positive, "beta")) {
      .occu_cover_latent_beta_logm(units$vals, eta, disp2, sigma_u,
                                   model$cover_latent_nquad %||% 15L)
    } else {
      # gaussian has no latent cover-aggregate variant; the dispatcher rejects
      # cover_aggregate = "latent" with a gaussian arm before reaching here.
      .occu_cover_latent_lognormal_logm(units$vals, eta, disp2, sigma_u)
    }
  } else {
    aggfun <- if (identical(mode, "median")) stats::median else mean
    yv  <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
    out[ps] <- .occu_cover_pos_logdens(yv, eta, exp(log_disp), positive)
  }
  out
}

# Per-site marginal log-likelihood (latent occupancy state z integrated out in
# closed form over its two states), returned as a length-`n_sites` vector. The
# inputs are the per-cell occupancy probability `psi`, the per-visit detection
# probability matrix `p_mat` [n_sites x max_visits], the per-visit cover
# linear-predictor matrix `ep_mat` [n_sites x max_visits], and the scalar
# `log_dispersion`. This is the single source of truth shared by the fit's
# negative-log-posterior and the WAIC / PSIS-LOO pointwise log-likelihood.
.occu_cover_site_ll <- function(model, psi, p_mat, ep_mat, log_disp,
                                 units = NULL) {
  valid <- model$valid
  y     <- model$y

  log_p   <- ifelse(valid, log(p_mat),     0)
  log_1mp <- ifelse(valid, log(1 - p_mat), 0)

  # Detection mixture under z = 1, then the cover term at the granularity the
  # fitter optimised (per-visit / aggregated / latent, gcol33/tulpaObs#34). The
  # cover term is non-zero only for units with a detection, matching the
  # any-detection branch below.
  log_h_det  <- ifelse(valid, ifelse(y == 1L, log_p, log_1mp), 0)
  cover_term <- .occu_cover_cover_term(model, ep_mat, log_disp, units)

  any_det <- rowSums(y * valid, na.rm = FALSE) > 0
  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))

  det_ll <- log_psi + rowSums(log_h_det) + cover_term
  # No detection: psi * prod(1-p) + (1-psi). Logsumexp form for stability.
  ln_a <- log_psi   + rowSums(log_1mp)
  ln_b <- log_1mpsi
  nodet_ll <- .tobs_logsumexp2(ln_a, ln_b)

  ifelse(any_det, det_ll, nodet_ll)
}


# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_occu_cover <- function(model,
                                  method   = c("laplace"),
                                  priors   = NULL,
                                  max.iter = 200L,
                                  tol      = 1e-6,
                                  verbose  = TRUE,
                                  sigma.beta = 5,
                                  ...) {
  method <- match.arg(method)

  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_par   <- p_occ + p_p + p_pos + 1L  # +1 for log_dispersion

  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos"
  )

  start <- numeric(n_par)
  names(start) <- par_names

  # Warm starts from a separate-fit baseline. Occurrence intercept from
  # empirical detection-any rate; detection intercept at logit(0.5); cover
  # intercept at the marginal positive-mean on its arm's link scale;
  # dispersion at a modestly broad value.
  any_det <- rowSums(model$y * model$valid) > 0
  det_rate <- max(mean(any_det), 1e-3)
  start[1L] <- stats::qlogis(min(max(det_rate, 1e-3), 1 - 1e-3))

  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  pos_vals <- pos_vals[is.finite(pos_vals)]
  if (length(pos_vals) > 0L) {
    if (identical(model$positive, "beta")) {
      start[p_occ + p_p + 1L] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
      start[n_par]            <- log(10)   # phi ~ 10 = moderate beta concentration
    } else if (identical(model$positive, "gaussian")) {
      # Identity-Gaussian arm (gcol33/tulpaObs#112): the response is raw (may be
      # negative), so seed the intercept / sigma on the natural scale, not log.
      start[p_occ + p_p + 1L] <- mean(pos_vals)
      start[n_par]            <- log(stats::sd(pos_vals) + 0.1)
    } else {
      start[p_occ + p_p + 1L] <- mean(log(pos_vals))
      start[n_par]            <- log(stats::sd(log(pos_vals)) + 0.1)
    }
  } else {
    start[n_par] <- if (identical(model$positive, "beta")) log(10) else log(0.4)
  }

  # Gaussian prior aligned with par. Default sigma.beta on the betas
  # (a weakly-informative N(0, sigma.beta^2)); dispersion stays flat.
  pmean <- numeric(n_par)
  pprec <- numeric(n_par)
  if (isTRUE(is.null(priors)) || !isFALSE(priors)) {
    beta_idx <- seq_len(n_par - 1L)
    pprec[beta_idx] <- 1 / (sigma.beta^2)
  }

  opt <- stats::optim(start, .tobs_occu_cover_nlp,
                       model = model, pmean = pmean, pprec = pprec,
                       method = "BFGS", hessian = TRUE,
                       control = list(maxit = max.iter, reltol = tol,
                                      trace = if (isTRUE(verbose)) 1L else 0L))

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V)) {
    warning("occu_cover: observed-Fisher Hessian not invertible; SEs unreliable.",
            call. = FALSE)
    V <- matrix(NA_real_, n_par, n_par)
  }
  se <- sqrt(pmax(diag(V), 0))

  means <- opt$par
  names(means) <- par_names
  names(se)    <- par_names
  dimnames(V)  <- list(par_names, par_names)

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = se,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    positive     = model$positive,
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter    = opt$counts[1L])
  )), class = c("tobs_fit", "tulpa_fit"))
}


