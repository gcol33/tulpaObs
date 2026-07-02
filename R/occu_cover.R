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

  # Cover arm: meaningful only where y == 1. Zero-fill the rest so matrix
  # algebra stays defined; the likelihood gates evaluation by `y == 1`.
  y_pos_num <- matrix(as.numeric(y_pos), n_sites, max_visits)
  pos_mask  <- valid & (y_int == 1L)
  if (any(!is.finite(y_pos_num[pos_mask]))) {
    stop("y_pos must be finite at every detected visit (y == 1).",
         call. = FALSE)
  }
  if (identical(positive, "beta")) {
    bad <- pos_mask & (y_pos_num <= 0 | y_pos_num >= 1)
    if (any(bad)) {
      stop("Beta positive arm requires 0 < y_pos < 1 at every detected ",
           "visit; clip with pmin(pmax(y_pos, eps), 1 - eps).", call. = FALSE)
    }
  } else {
    bad <- pos_mask & (y_pos_num <= 0)
    if (any(bad)) {
      stop("Lognormal positive arm requires y_pos > 0 at every detected ",
           "visit.", call. = FALSE)
    }
  }
  y_pos_num[!pos_mask] <- 0

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

  # Cover meaningful only where detected; zero-fill the rest (the spec gates the
  # positive density on y_det == 1, exactly as the dense binder does).
  pos_mask  <- y_det_visit == 1L
  y_pos_num <- as.numeric(y_pos_values)
  if (any(!is.finite(y_pos_num[pos_mask]))) {
    stop("y_pos must be finite at every detected visit (y == 1).", call. = FALSE)
  }
  if (identical(positive, "beta")) {
    bad <- pos_mask & (y_pos_num <= 0 | y_pos_num >= 1)
    if (any(bad)) {
      stop("Beta positive arm requires 0 < y_pos < 1 at every detected visit; ",
           "clip with pmin(pmax(y_pos, eps), 1 - eps).", call. = FALSE)
    }
  } else {
    bad <- pos_mask & (y_pos_num <= 0)
    if (any(bad)) {
      stop("Lognormal positive arm requires y_pos > 0 at every detected visit.",
           call. = FALSE)
    }
  }
  y_pos_num[!pos_mask] <- 0

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

# Reject formulas containing structured terms (bym2, icar, car, gp, etc.).
# v1 of occu_cover is non-spatial; spatial sharing across the three arms is v2.
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
  valid_flat <- as.logical(t(valid))
  lapply(re_parse$terms, function(spec) {
    g_flat <- .occu_cover_obs_flat_eval(spec$group_expr, data, visit_df,
                                        n_sites, max_visits, arm)
    var <- paste(spec$vars, collapse = ":")
    lev <- sort(unique(as.character(g_flat[valid_flat])))
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
      Xs <- .occu_cover_obs_flat_eval(cov_expr, data, visit_df,
                                      n_sites, max_visits, arm)
      Xs <- as.matrix(Xs); storage.mode(Xs) <- "double"
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
        s <- stats::sd(col[valid_flat]); if (!is.finite(s) || s <= 0) 1 else s
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
  lapply(re_parse$terms, function(spec) {
    g <- .occu_cover_obs_flat_eval_ragged(spec$group_expr, data, visit_df,
                                          site_of_visit, arm)
    var <- paste(spec$vars, collapse = ":")
    lev <- sort(unique(as.character(g)))
    if (length(lev) < 2L) {
      stop(sprintf(paste0(
        "occu_cover(): the %s random-effect grouping `%s` has %d level(s) among ",
        "the observed visits; a random effect needs at least 2 groups."),
        arm, var, length(lev)), call. = FALSE)
    }
    codes <- match(as.character(g), lev)
    codes[is.na(codes)] <- 0L

    if (identical(spec$type, "slope")) {
      cov_expr <- spec$covariate
      cov_chr <- tryCatch(eval(cov_expr), error = function(...) NULL)
      if (is.character(cov_chr))
        cov_expr <- as.call(c(list(as.name("cbind")), lapply(cov_chr, as.name)))
      Xs <- .occu_cover_obs_flat_eval_ragged(cov_expr, data, visit_df,
                                             site_of_visit, arm)
      Xs <- as.matrix(Xs); storage.mode(Xs) <- "double"
      if (is.null(colnames(Xs)))
        colnames(Xs) <- paste0("slope", seq_len(ncol(Xs)))
      slope_scale <- apply(Xs, 2L, function(col) {
        s <- stats::sd(col); if (!is.finite(s) || s <= 0) 1 else s
      })
      Xs <- sweep(Xs, 2L, slope_scale, "/")
      has_int <- isTRUE(spec$intercept)
      Z <- if (has_int) cbind(`(Intercept)` = 1, Xs) else Xs
      coef_names  <- colnames(Z); n_coefs <- ncol(Z)
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
  det_mat  <- model$valid & (model$y == 1L)
  pos_site <- which(rowSums(det_mat) > 0L)
  vals <- lapply(pos_site, function(i) as.numeric(model$y_pos[i, det_mat[i, ]]))
  list(pos_site = pos_site, vals = vals)
}

# Positive-arm log-density of cover value(s) `y` at cover predictor `eta`
# (link scale) and dispersion `disp` (lognormal residual SD or beta precision).
# Vectorised over y / eta (and matrices), so the per-visit and the per-unit
# aggregated cover terms read one formula. Beta clamps the predictor before the
# logistic; lognormal uses the raw predictor (matching the historical kernels).
.occu_cover_pos_logdens <- function(y, eta, disp, is_beta) {
  if (is_beta) {
    mu <- stats::plogis(.tobs_clamp_eta(eta))
    a  <- mu * disp
    b  <- (1 - mu) * disp
    lgamma(disp) - lgamma(a) - lgamma(b) +
      (a - 1) * log(y) + (b - 1) * log(1 - y)
  } else {
    -log(y) - log(disp) - 0.5 * log(2 * pi) -
      0.5 * ((log(y) - eta) / disp)^2
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
      ell <- sum(.occu_cover_pos_logdens(v, eta[i] + sigma_u * z[k], phi, TRUE))
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
  n_sites <- model$n_sites
  is_beta <- identical(model$positive, "beta")
  mode    <- model$cover_aggregate %||% "none"

  if (identical(mode, "none")) {
    pos_mask <- model$valid & (model$y == 1L)
    dens <- .occu_cover_pos_logdens(model$y_pos, ep_mat, exp(log_disp), is_beta)
    log_f_pos <- matrix(0, n_sites, model$max_visits)
    log_f_pos[pos_mask] <- dens[pos_mask]
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
    out[ps] <- if (is_beta) {
      .occu_cover_latent_beta_logm(units$vals, eta, disp2, sigma_u,
                                   model$cover_latent_nquad %||% 15L)
    } else {
      .occu_cover_latent_lognormal_logm(units$vals, eta, disp2, sigma_u)
    }
  } else {
    aggfun <- if (identical(mode, "median")) stats::median else mean
    yv  <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
    out[ps] <- .occu_cover_pos_logdens(yv, eta, exp(log_disp), is_beta)
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
  if (length(pos_vals) > 0L) {
    if (identical(model$positive, "beta")) {
      start[p_occ + p_p + 1L] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
      start[n_par]            <- log(10)   # phi ~ 10 = moderate beta concentration
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
  draws <- .occu_cover_rmvn(n_draws, means, V)
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

# Draw from MVN via Cholesky; fall back to independent normals if not PD.
.occu_cover_rmvn <- function(n, mu, sigma) {
  p <- length(mu)
  if (any(!is.finite(sigma))) {
    return(matrix(rep(mu, each = n), n, p, byrow = FALSE))
  }
  L <- tryCatch(chol(sigma), error = function(e) NULL)
  z <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(L)) {
    sds <- sqrt(pmax(diag(sigma), 1e-8))
    return(sweep(z * rep(sds, each = n), 2L, mu, "+"))
  }
  sweep(z %*% L, 2L, mu, "+")
}


# ---------------------------------------------------------------------------
# Formula-native cross-arm coupling (INLA-style copy())
#
# The occurrence arm carries a spatial field; the positive (cover) arm declares
# it carries a scaled copy of that field with a copy() selector, the DAG edge
# u_occ -> cover placed in the formula:
#
#   copy(spatial(), alpha = grid(g))   the unique occurrence spatial effect,
#                                      one amplitude g over every block
#   copy(spatial(cell_idx), ...)       disambiguate by grouping variable when
#                                      the occurrence arm has several spatials
#   copy(spatial(), terms = list(intercept = grid(g0), time.sc = grid(g1)))
#                                      a per-block amplitude
#   copy("occ_space", ...)             explicit-name reference (lower-level),
#                                      requires spatial(..., name = "occ_space")
#
# No name is needed in the common case: spatial() selects the occurrence arm's
# spatial effect structurally. The engine still reads the coupling amplitude axes
# off `control$alpha.grid` / `control$alpha.grid.trend`, so the formula copy() is
# translated into those axes here and the downstream fit is unchanged:
#
#   whole-field amplitude g            -> alpha.grid = g, alpha.grid.trend = g
#   terms = list(intercept=, time.sc=) -> alpha.grid = ., alpha.grid.trend = .
#
# `alpha = grid(g)` integrates over g; a scalar fixes it (a length-1 grid). The
# block layout follows .occu_cover_spatial_fields(): the unweighted intercept
# field is block 1, weighted trend field(s) block 2+. Decoupling an arm is
# structural -- write spatial() (an own field) or omit copy() -- not a magic
# alpha of 0; 0 is only ever one value you could place in a grid.
# ---------------------------------------------------------------------------

# Parse copy() terms off the positive formula, returning the stripped
# fixed-effects positive formula plus the list of tobs_copy specs. The copy
# special is the only structured term allowed on the positive arm; any other is
# left in place for .occu_cover_reject_structured() to reject.
.occu_cover_extract_pos_copies <- function(pos_formula) {
  if (is.null(pos_formula)) return(list(formula = pos_formula, copies = list()))
  parsed <- .tobs_parse_formula(pos_formula, data = NULL)
  copies <- Filter(function(t) inherits(t, "tobs_copy"), parsed$terms)
  list(formula = parsed$fe_formula, copies = copies)
}

# Spatial-field constructors that declare a NEW latent field (unlike copy(), which
# reuses a named one). A term with one of these heads is a field; placement in an
# arm's formula puts the field on that arm.
.occu_cover_field_ctors <- c("spatial", "icar", "bym2", "car", "car_proper")

# Placement is the canonical way to put a field on an arm: a spatial-field term
# written in the detection or positive formula declares a field ON that arm. Lift
# such terms onto the occurrence formula, tagged with their arm via `to`, so the
# arm-generic resolver (.occu_cover_spatial_fields) sees every arm's fields in one
# place. copy() and RE terms are handled separately and left untouched. Returns the
# augmented occurrence formula and the stripped detection / positive formulas.
.occu_cover_lift_arm_fields <- function(occ_formula, det_formula, pos_formula) {
  arm_extra <- list()

  strip_arm <- function(arm_formula, arm) {
    if (is.null(arm_formula)) return(arm_formula)
    tt   <- stats::terms(arm_formula, keep.order = TRUE)
    labs <- attr(tt, "term.labels")
    keep <- character(0)
    for (lab in labs) {
      e    <- tryCatch(str2lang(lab), error = function(...) NULL)
      head <- if (is.call(e) && is.symbol(e[[1L]])) as.character(e[[1L]]) else NA_character_
      if (!is.na(head) && head %in% .occu_cover_field_ctors) {
        # The arm is fixed by placement. An explicit `to` that disagrees is an
        # error; a matching one is redundant but allowed.
        if (!is.null(e$to)) {
          to_val <- tryCatch(as.character(eval(e$to)), error = function(...) NULL)
          if (length(to_val) && !all(to_val == arm))
            stop(sprintf(paste0(
              "occu_cover(): a spatial field written in the %s formula is on the ",
              "%s arm by placement; drop the conflicting `to =`."), arm, arm),
              call. = FALSE)
        }
        e$to <- arm
        arm_extra[[length(arm_extra) + 1L]] <<- e
      } else if (!is.na(head) && head %in% c("|", "||")) {
        # An lme4 RE bar. terms() strips the parentheses off `(1 | g)` down to the
        # label `1 | g`; reformulate() would rebuild it as `... + 1 | g`, which R
        # re-parses as `(... + 1) | g` -- no longer an RE bar. Restore the parens
        # so the downstream RE parse (.occu_cover_obs_re_parse) still sees a bar.
        keep <- c(keep, sprintf("(%s)", lab))
      } else {
        keep <- c(keep, lab)
      }
    }
    fe <- stats::reformulate(
      termlabels = if (length(keep)) keep else "1",
      intercept  = as.logical(attr(tt, "intercept")))
    environment(fe) <- environment(arm_formula)
    fe
  }

  det2 <- strip_arm(det_formula, "detection")
  pos2 <- strip_arm(pos_formula, "positive")

  if (length(arm_extra)) {
    occ_tt   <- stats::terms(occ_formula)
    occ_labs <- attr(occ_tt, "term.labels")
    extra    <- vapply(arm_extra, function(e) paste(deparse(e), collapse = ""),
                       character(1))
    occ2 <- stats::reformulate(termlabels = c(occ_labs, extra),
                               intercept  = as.logical(attr(occ_tt, "intercept")))
    environment(occ2) <- environment(occ_formula)
  } else {
    occ2 <- occ_formula
  }
  list(occ = occ2, det = det2, pos = pos2)
}

# Map the positive arm's copy() specs onto the coupling-amplitude grids the
# joint_coupled fitter reads (control$alpha.grid for the intercept block,
# control$alpha.grid.trend for the trend block). `spatial_info` carries the
# resolved fields (block 1 = intercept, block 2+ = weighted trend), each with a
# `field_name` and a `component` label. Returns the updated control list.
#
# On the formula-native path (a named occupancy field, or any copy() present)
# the amplitude axes come ENTIRELY from copy(): a field block with no copy() is
# pinned at alpha = 0 (decoupled). On the back-compat path (no name, no copy())
# control$alpha.grid / .trend are left untouched, so old fits are byte-identical.
.occu_cover_apply_copy_coupling <- function(copies, spatial_info, control) {
  has_control_alpha <- any(c("alpha.grid", "alpha.grid.trend") %in% names(control))
  if (has_control_alpha && length(copies) > 0L) {
    stop("occu_cover(): set the cross-arm coupling with copy() in the positive ",
         "formula OR control$alpha.grid[.trend], not both.", call. = FALSE)
  }
  # control$alpha.grid is the low-level amplitude knob: when set (and no copy())
  # the engine reads the grids as given.
  if (has_control_alpha) return(control)

  if (is.null(spatial_info)) {
    if (length(copies) > 0L) {
      stop("occu_cover(): copy() needs a spatial field on the occurrence ",
           "formula, e.g. spatial(~ 1 || cell, graph = adj).", call. = FALSE)
    }
    return(control)
  }

  # Coupling is formula-native and explicit: a copy() carries the occurrence
  # spatial field onto the cover arm with the amplitude it names; a block with no
  # copy() is decoupled (alpha pinned 0), the field rides occupancy only. There
  # is no implicit default coupling.
  #
  # Component labels of the resolved field blocks. Block 1 is the intercept
  # field; blocks 2+ are weighted trend fields. "trend" is an alias for the
  # single trend block, the column name (e.g. "time.sc") names it explicitly.
  has_trend  <- length(spatial_info$fields) > 1L
  components <- vapply(spatial_info$fields,
                       function(f) f$component %||% NA_character_, character(1))
  node       <- spatial_info$group_var

  # Default to decoupled: every block pinned at alpha = 0. A copy() then sets the
  # amplitude axis on the block(s) it names.
  alpha_int   <- 0
  alpha_trend <- if (has_trend) 0 else NULL

  # Apply one (component, amplitude) assignment, returning the canonical block
  # role ("intercept" / "trend") it resolved to. `comp = NULL` is the whole
  # field. `cp_label` names the copy() in any error.
  apply_component <- function(comp, g, cp_label) {
    if (is.null(comp)) {
      alpha_int <<- g %||% .occu_cover_default_alpha_grid()
      if (has_trend) alpha_trend <<- g %||% .occu_cover_default_alpha_grid()
      return(c("intercept", if (has_trend) "trend"))
    }
    if (identical(comp, "intercept")) {
      alpha_int <<- g %||% .occu_cover_default_alpha_grid()
      return("intercept")
    }
    if (identical(comp, "trend") ||
        (has_trend && comp %in% stats::na.omit(components[-1L]))) {
      if (!has_trend) {
        stop(sprintf(paste0(
          "%s: the spatial field has no trend component (it is a single ",
          "intercept field)."), cp_label), call. = FALSE)
      }
      alpha_trend <<- g %||% .occu_cover_default_alpha_grid()
      return("trend")
    }
    avail <- paste0("\"", stats::na.omit(components), "\"", collapse = ", ")
    stop(sprintf(paste0(
      "%s: unknown field component \"%s\". Available component(s): %s, or the ",
      "whole field."), cp_label, comp, avail), call. = FALSE)
  }

  for (cp in copies) {
    cp_label <- sprintf("copy(%s)", cp$id)

    # The positive arm couples the occurrence spatial field through a selector
    # (spatial() / spatial(<grouping_var>)); a string reference is not a coupling
    # selector here.
    if (is.null(cp$selector_type)) {
      stop(sprintf(paste0(
        "%s: select the occurrence spatial field with copy(spatial()) or ",
        "copy(spatial(%s)), not a string."), cp_label, node %||% "<grouping_var>"),
        call. = FALSE)
    }
    if (!is.null(cp$selector_group)) {
      if (is.null(node)) {
        stop(sprintf(paste0(
          "%s: the occurrence spatial field has no named grouping variable; ",
          "use copy(spatial())."), cp_label), call. = FALSE)
      }
      if (!identical(cp$selector_group, node)) {
        stop(sprintf(paste0(
          "%s: no spatial effect grouped on \"%s\"; the occurrence spatial ",
          "field is on \"%s\". Use copy(spatial(%s)) or copy(spatial())."),
          cp_label, cp$selector_group, node, node), call. = FALSE)
      }
    }

    # Per-component grids (terms = list(...)) must address every field block, so
    # no block is silently left at alpha = 0. The whole-field / single-component
    # forms set the blocks they name; any unnamed block stays decoupled.
    if (!is.null(cp$copy_terms)) {
      covered <- character(0)
      for (k in names(cp$copy_terms)) {
        res <- cp$copy_terms[[k]]
        g   <- if (isTRUE(is.na(res$integrate))) NULL else res$grid
        covered <- union(covered, apply_component(k, g, cp_label))
      }
      required <- c("intercept", if (has_trend) "trend")
      missing_blocks <- setdiff(required, covered)
      if (length(missing_blocks)) {
        blocks <- paste0("\"", stats::na.omit(components), "\"", collapse = ", ")
        stop(sprintf(paste0(
          "%s: terms = must give an amplitude for every field block; %s left ",
          "unaddressed. Field blocks: %s."), cp_label,
          paste0("\"", missing_blocks, "\"", collapse = ", "), blocks),
          call. = FALSE)
      }
    } else {
      g <- if (is.na(cp$alpha_integrate)) NULL else cp$alpha_grid
      apply_component(cp$component, g, cp_label)
    }
  }

  control[["alpha.grid"]] <- alpha_int
  if (has_trend) control[["alpha.grid.trend"]] <- alpha_trend
  control
}

# The engine's default coupling grid (a copy() with no alpha = falls back to it,
# matching the joint_coupled fitter's own default). Single source of truth with
# .tobs_fit_occu_cover_joint_coupled().
.occu_cover_default_alpha_grid <- function() {
  c(0, exp(seq(log(0.1), log(3), length.out = 5)))
}


# ---------------------------------------------------------------------------
# Dispatcher (wired into tobs.R's switch)
# ---------------------------------------------------------------------------

.dispatch_occu_cover <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)

  if (is.null(detection)) {
    stop("occu_cover() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_cover() requires `y` (N x J detection-history matrix).",
         call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("occu_cover() requires `y_pos` (N x J positive-cover matrix; ",
         "values used only where y == 1).", call. = FALSE)
  }

  # Compact (ragged) input: `y` is a tobs_ragged carrier (tobs_data(compact =
  # TRUE)). It feeds the joint nested-Laplace engine one valid visit at a time,
  # with no padded grid and so no per-site visit cap. `y_pos` may be the paired
  # ragged carrier or a bare length-V vector aligned to `y`. Scoped to the
  # spatial joint path (the only consumer of the compacted arms); every other
  # configuration that would read the dense `y` / `valid` grid is gated off below
  # with a clear error rather than a silent dense rebuild.
  ragged <- inherits(y, "tobs_ragged")
  if (ragged) {
    y_pos_arg <- dots$y_pos
    if (inherits(y_pos_arg, "tobs_ragged")) {
      if (!identical(y_pos_arg$site, y$site) || !identical(y_pos_arg$visit, y$visit))
        stop("occu_cover(): compact `y` and `y_pos` are not aligned (different ",
             "site / visit order). Build both with tobs_data(compact = TRUE) on ",
             "the same df / site / visit.", call. = FALSE)
      y_pos_values <- y_pos_arg$values
    } else {
      y_pos_values <- y_pos_arg
    }
  }

  pos_formula <- dots$positive

  # Placement -> arm: a spatial-field term written in the positive (or detection)
  # formula is lifted onto the occurrence formula tagged with its arm, so the
  # arm-generic spatial resolver sees every arm's fields together. Done before the
  # RE parse and the positive-defaults-to-detection fallback, so each downstream
  # step sees a field-free arm formula and a shared formula is not lifted twice.
  # copy() stays on its arm's formula (it is a reference, not a new field).
  #
  # `to =` is retired: an arm is chosen by placement (write the field in that
  # arm's formula) and a field is shared across arms with copy(). Reject an
  # explicit user `to =` on the raw arm formulas, before placement sets the
  # internal `to` on the lifted calls.
  .tobs_reject_user_to(formula,   "occurrence")
  .tobs_reject_user_to(detection, "detection")
  .tobs_reject_user_to(pos_formula, "positive")

  lifted      <- .occu_cover_lift_arm_fields(formula, detection, pos_formula)
  formula     <- lifted$occ
  detection   <- lifted$det
  pos_formula <- lifted$pos

  if (is.null(pos_formula)) pos_formula <- detection

  # Observation-arm random intercept (gcol33/tulpaObs#102): a `(1 | g)` / `re(g)`
  # on the detection or positive-cover formula adds an iid RE latent block on
  # that arm of the joint nested-Laplace fit. Parse it off FIRST -- before the
  # copy() extraction and every design build -- so each downstream consumer sees
  # a clean fixed-effects formula; the grouping is resolved to per-visit codes
  # once the model (and its `valid` mask) is built. The block rides the joint
  # engine, so it needs method = "nested_laplace" (the non-spatial laplace / nuts
  # paths are coefficient-marginal fits with no latent block) -- error early
  # rather than silently drop the RE under another engine.
  det_re_parse <- .occu_cover_obs_re_parse(detection,   "detection")
  pos_re_parse <- .occu_cover_obs_re_parse(pos_formula, "positive cover")
  has_obs_re   <- !is.null(det_re_parse) || !is.null(pos_re_parse)
  if (has_obs_re && !identical(engine, "nested_laplace")) {
    stop("occu_cover(): a random effect on the detection / positive-cover arm ",
         "needs method = \"nested_laplace\" (the joint nested-Laplace engine ",
         "carries the RE as a latent block); got method = \"", engine, "\".",
         call. = FALSE)
  }
  if (ragged && identical(engine, "nuts")) {
    stop("occu_cover(): compact (ragged) input feeds the joint nested-Laplace ",
         "engine; method = \"nuts\" reads the dense detection grid. Use ",
         "method = \"nested_laplace\", or build the data densely for NUTS.",
         call. = FALSE)
  }
  # Random intercepts (crossed, nested) AND random slopes are supported on the
  # detection / positive-cover arms (gcol33/tulpaObs#103): an intercept rides one
  # `iid` block per term, an uncorrelated slope one weighted `iid` block per
  # coefficient, a correlated slope one multivariate free-Sigma `miid` block.
  # The slope blocks need tulpa's joint engine >= 0.0.39 (gcol33/tulpa#114), which
  # the DESCRIPTION Imports floor enforces.
  if (!is.null(det_re_parse)) detection   <- det_re_parse$fe
  if (!is.null(pos_re_parse)) pos_formula <- pos_re_parse$fe

  # Formula-native cross-arm coupling: pull any copy() term off the positive
  # formula and keep the stripped fixed-effects design. The copy() specs are
  # translated into the engine's coupling-amplitude grids once the occupancy
  # field blocks are resolved (below). Stripping here lets every downstream
  # consumer (NUTS branch, design build, structured-term rejection) see a clean
  # fixed-effects positive formula.
  pos_copy   <- .occu_cover_extract_pos_copies(pos_formula)
  pos_copies <- pos_copy$copies
  pos_formula <- pos_copy$formula

  # Spatial NUTS path (gcol33/tulpaObs#74): a car_proper() term on the psi formula
  # under method = "nuts" samples a FIXED-HYPER non-centered coupled proper-CAR
  # field jointly with the coefficient marginal (rather than grid-integrating it,
  # as nested_laplace does). Detected separately from the grid-integrated
  # icar/bym2 fields below, because the proper-CAR field is a NUTS-only structure.
  if (identical(engine, "nuts")) {
    nuts_sp <- .occu_cover_nuts_spatial_term(formula, data)
    if (!is.null(nuts_sp)) {
      .occu_cover_reject_structured(detection,   "detection")
      .occu_cover_reject_structured(pos_formula, "positive cover")
      vd_det  <- .normalize_visits(visits, detection,
                                   n_sites = nrow(y), max_visits = ncol(y))
      vd_pos  <- .normalize_visits(visits, pos_formula,
                                   n_sites = nrow(y), max_visits = ncol(y))
      model_sp <- .tobs_build_occu_cover(
        occ_formula = nuts_sp$fe, det_formula = vd_det$det_formula,
        pos_formula = vd_pos$det_formula, data = data, y = y,
        y_pos = dots$y_pos, positive = family$params$positive,
        det_visit_formula = vd_det$det_visit_formula,
        det_visit_data    = vd_det$visits,
        pos_visit_formula = vd_pos$det_visit_formula,
        pos_visit_data    = vd_pos$visits)
      model_sp$cover_aggregate <- "none"
      # Resolve the site -> field-node map (group_var lets sites > cells).
      sp_graph <- nuts_sp$spatial$graph
      n_cells_f <- nrow(sp_graph)
      gv <- nuts_sp$group_var
      if (!is.null(gv)) {
        if (!gv %in% names(data))
          stop(sprintf("occu_cover() group_var '%s' is not a column of data.",
                       gv), call. = FALSE)
        site_cell <- as.integer(data[[gv]])
        if (length(site_cell) != model_sp$n_sites || anyNA(site_cell) ||
            min(site_cell) < 1L || max(site_cell) > n_cells_f)
          stop(sprintf(paste0(
            "occu_cover() group_var '%s' must be an integer cell index in ",
            "1..%d, one per site (%d sites)."), gv, n_cells_f, model_sp$n_sites),
            call. = FALSE)
      } else {
        if (model_sp$n_sites != n_cells_f)
          stop(sprintf(paste0(
            "occu_cover() NUTS spatial: %d sites but the graph has %d nodes. ",
            "Map sites to cells with group_var on the car_proper() term, or ",
            "match the site count to the graph."),
            model_sp$n_sites, n_cells_f), call. = FALSE)
        site_cell <- seq_len(model_sp$n_sites)
      }
      model_sp$site_cell <- site_cell
      model_sp$n_cells   <- n_cells_f
      return(do.call(.tobs_fit_occu_cover_nuts_spatial,
                     c(list(model = model_sp, spatial = nuts_sp$spatial,
                            priors = priors), control)))
    }
  }

  # Detect the coupled spatial field(s) on the psi formula. The spatial path is
  # the joint nested-Laplace engine (shared field(s) across the psi and cover
  # arms); the non-spatial path is plain Laplace on the exact two-state
  # marginal. A weighted areal term adds a second coupled (SVC) field.
  spatial_info <- .occu_cover_spatial_fields(formula, data)
  has_spatial  <- !is.null(spatial_info)

  # Translate the positive arm's copy() spec(s) into the engine coupling grids
  # now that the occupancy field blocks are resolved. On the formula-native path
  # (named field and/or copy()) this sets control$alpha.grid[.trend]; on the
  # back-compat path it is a no-op, so control-driven fits are unchanged. A
  # copy() with no named field, or a named field outside the spatial path, is an
  # error surfaced here.
  if (length(pos_copies) > 0L && !has_spatial) {
    stop("occu_cover(): copy() on the positive arm needs a spatial field on the ",
         "occurrence formula (e.g. spatial(~ 1 || cell, graph = adj, name = ",
         "\"occ_space\")) under method = \"nested_laplace\".", call. = FALSE)
  }
  control <- .occu_cover_apply_copy_coupling(pos_copies, spatial_info, control)

  # Resolve cover aggregation (tulpaObs#33). NULL (unset) -> "mean" on the
  # shared-field spatial path (so the cover arm contributes at the cell scale and
  # does not outweigh occupancy on the shared field), "none" (per-visit) on the
  # non-spatial path (no shared field to over-weight). `agg_explicit` records
  # whether the user set it: an explicit mean / median on an unsupported
  # configuration errors, whereas the bare default quietly falls back to
  # per-visit cover so a plain visit-level fit keeps working.
  #
  # Aggregated cover is a per-cell observation, so its positive design must be
  # cell-level (resolved from the cell `data`). A `positive` formula that
  # references a visit-level covariate (a name carried in `visits`) is a
  # per-visit design and cannot be aggregated: an explicit request errors, the
  # bare default falls back to per-visit cover.
  visit_cov_names <- if (is.null(visits)) character(0)
                     else if (is.data.frame(visits) || is.list(visits)) names(visits)
                     else character(0)
  pos_is_visit_level <- length(intersect(all.vars(pos_formula),
                                          visit_cov_names)) > 0L

  agg_explicit    <- !is.null(family$params$cover_aggregate)
  cover_aggregate <- family$params$cover_aggregate %||%
                     (if (has_spatial) "mean" else "none")
  if (!has_spatial && cover_aggregate != "none") {
    stop(sprintf(paste0(
      "occu_cover(cover_aggregate = \"%s\") aggregates the cover arm on the ",
      "shared-field spatial path (method = \"nested_laplace\"); the non-spatial ",
      "laplace fit uses per-visit cover (cover_aggregate = \"none\")."),
      cover_aggregate), call. = FALSE)
  }
  if (cover_aggregate != "none" && pos_is_visit_level) {
    if (agg_explicit) {
      stop(sprintf(paste0(
        "occu_cover() cell-aggregated cover (cover_aggregate = \"%s\") needs a ",
        "cell-level positive design, but the `positive` formula references the ",
        "visit-level covariate(s) %s (carried in `visits`). Use a cell-level ",
        "positive covariate (a column of `data`), or cover_aggregate = \"none\" ",
        "for per-visit cover."), cover_aggregate,
        paste(intersect(all.vars(pos_formula), visit_cov_names),
              collapse = ", ")), call. = FALSE)
    }
    cover_aggregate <- "none"
  }
  # An arm-specific cover field (gcol33/tulpaObs#110) is scored per detected visit
  # (its node/weight index the pos-arm visit rows), so it needs per-visit cover;
  # an explicit aggregation errors, the bare default falls back to per-visit.
  if (has_spatial && !is.null(spatial_info$pos_armspec) &&
      cover_aggregate != "none") {
    if (agg_explicit) {
      stop(sprintf(paste0(
        "occu_cover(): an arm-specific cover field (to = \"positive\") uses ",
        "per-visit cover (cover_aggregate = \"none\"); it cannot map onto ",
        "cell-aggregated cover rows. Got cover_aggregate = \"%s\"."),
        cover_aggregate), call. = FALSE)
    }
    cover_aggregate <- "none"
  }

  if (has_spatial && engine == "laplace") {
    stop("occu_cover() found a spatial term (icar/bym2) in the psi formula ",
         "but method = \"laplace\" is non-spatial. Use method = ",
         "\"nested_laplace\" for the spatial v2 path.", call. = FALSE)
  }
  if (has_spatial && engine == "nuts") {
    # A car_proper() term would already have routed to the spatial NUTS fitter
    # above; reaching here means an intrinsic icar()/bym2() field, whose flat
    # field-mean direction needs the grid-integrated nested-Laplace path (or the
    # sampled-field community route) -- it is not the fixed-hyper NUTS structure.
    stop("occu_cover() with method = \"nuts\" samples a FIXED-HYPER proper-CAR ",
         "shared field; the intrinsic icar()/bym2() field on the psi formula ",
         "has a flat field-mean direction needing a sum-to-zero reparameterisation ",
         "for NUTS -- use method = \"nested_laplace\" (the shared coupled field is ",
         "grid-integrated), car_proper() for the NUTS shared field, or ",
         "ms_occu_cover() + icar() for a sampled shared field. (gcol33/tulpaObs#74)",
         call. = FALSE)
  }
  if (!has_spatial && engine == "nested_laplace") {
    stop("occu_cover() with method = \"nested_laplace\" requires a spatial ",
         "term (icar() or bym2()) on the psi formula.", call. = FALSE)
  }
  if (ragged && !has_spatial) {
    stop("occu_cover(): compact (ragged) input is implemented for the joint ",
         "nested-Laplace path (a spatial term on the occurrence formula). For a ",
         "non-spatial fit, build the data densely (tobs_data() without ",
         "compact = TRUE).", call. = FALSE)
  }
  if (ragged && cover_aggregate != "none") {
    stop("occu_cover(): compact (ragged) input uses per-visit cover ",
         "(cover_aggregate = \"none\"); the cell-aggregated cover path reads the ",
         "dense detection grid. Pass cover_aggregate = \"none\", or build densely.",
         call. = FALSE)
  }

  fe_formula <- if (has_spatial) spatial_info$fe else formula

  # Detection / cover arms never carry a spatial term (the shared field is
  # on the latent state z, not on the observation process). Other structured
  # terms (re, temporal, ...) are not supported on any arm in v1/v2.
  .occu_cover_reject_structured(detection,   "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  # Visit-design normalization + model build. The compact (ragged) path builds
  # the visit designs on the V valid rows directly (no padded grid); the dense
  # path flattens the [n_sites x max_visits] grid. Both converge to the same
  # `tobs_model` shape (the compact one carries `ragged = TRUE` and pre-compacted
  # visit-level fields), so the joint-coupled arm builder downstream is the same.
  if (ragged) {
    n_v    <- y$n_visits
    vd_det <- .normalize_visits_ragged(visits, detection, n_visits_valid = n_v)
    vd_pos <- .normalize_visits_ragged(visits, pos_formula, n_visits_valid = n_v)
    model  <- .tobs_build_occu_cover_ragged(
      occ_formula       = fe_formula,
      det_formula       = vd_det$det_formula,
      pos_formula       = vd_pos$det_formula,
      data              = data,
      y_ragged          = y,
      y_pos_values      = y_pos_values,
      positive          = family$params$positive,
      det_visit_formula = vd_det$det_visit_formula,
      det_visit_data    = vd_det$visits,
      pos_visit_formula = vd_pos$det_visit_formula,
      pos_visit_data    = vd_pos$visits)
  } else {
    vd_det <- .normalize_visits(visits, detection,
                                n_sites = nrow(y), max_visits = ncol(y))
    # Positive design. Per-visit cover reads the visit-level positive formula from
    # `visits`; cell-aggregated cover reads a cell-level positive design directly
    # from `data` (one value per occupancy unit) and carries no visit-level term.
    if (cover_aggregate == "none") {
      vd_pos            <- .normalize_visits(visits, pos_formula,
                                             n_sites = nrow(y), max_visits = ncol(y))
      pos_site_formula  <- vd_pos$det_formula
      pos_visit_formula <- vd_pos$det_visit_formula
      pos_visit_data    <- vd_pos$visits
    } else {
      pos_site_formula  <- pos_formula
      pos_visit_formula <- NULL
      pos_visit_data    <- NULL
    }

    model <- .tobs_build_occu_cover(
      occ_formula      = fe_formula,
      det_formula      = vd_det$det_formula,
      pos_formula      = pos_site_formula,
      data             = data,
      y                = y,
      y_pos            = dots$y_pos,
      positive         = family$params$positive,
      det_visit_formula = vd_det$det_visit_formula,
      det_visit_data    = vd_det$visits,
      pos_visit_formula = pos_visit_formula,
      pos_visit_data    = pos_visit_data
    )
  }

  model$cover_aggregate <- cover_aggregate

  if (has_spatial) {
    fields      <- spatial_info$fields
    base_graph  <- fields[[1L]]$graph

    # Resolve the site -> field-node map. With group_var the occupancy units
    # (sites, one per row of `data` / `y`) map onto fewer field nodes (cells),
    # so the same cell field is shared across that cell's sites (e.g. cell-year
    # sites sharing one cell). Without group_var the two coincide 1:1.
    n_cells_field <- nrow(base_graph)
    gv <- spatial_info$group_var
    if (!is.null(gv)) {
      if (!gv %in% names(data)) {
        stop(sprintf("occu_cover() group_var '%s' is not a column of data.", gv),
             call. = FALSE)
      }
      site_cell <- as.integer(data[[gv]])
      if (length(site_cell) != model$n_sites || anyNA(site_cell) ||
          min(site_cell) < 1L || max(site_cell) > n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() group_var '%s' must be an integer cell index in 1..%d, ",
          "one per site (%d sites)."), gv, n_cells_field, model$n_sites),
          call. = FALSE)
      }
    } else {
      if (model$n_sites != n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() spatial: %d sites but the graph has %d nodes. Map sites ",
          "to cells with group_var = \"<col>\" on the icar()/bym2() term (e.g. ",
          "site = cell-year), or match the site count to the graph."),
          model$n_sites, n_cells_field), call. = FALSE)
      }
      site_cell <- seq_len(model$n_sites)
    }
    model$site_cell <- site_cell
    model$n_cells   <- n_cells_field

    # Optional per-group random intercept on the occupancy arm, layered on the
    # shared field (gcol33/tulpaObs#56). The grouping is per occupancy unit (one
    # code per site / data row); validate its length and carry it to the fitter.
    re_spec <- spatial_info$re
    if (!is.null(re_spec)) {
      if (length(re_spec$group_idx) != model$n_sites) {
        stop(sprintf(paste0(
          "occu_cover() spatial + RE: the random-effect grouping has %d codes ",
          "but there are %d occupancy units (sites)."),
          length(re_spec$group_idx), model$n_sites), call. = FALSE)
      }
    }

    # Observation-arm random effects (gcol33/tulpaObs#102, #103): resolve each
    # detection / positive-cover RE term to its per-(site, visit) design (group
    # codes + slope weights) now that the model carries its `valid` mask, and
    # attach the per-term LIST to the model. The fitter reads model$re_det /
    # model$re_pos (one entry per crossed / nested term) and the joint-coupled
    # builder subsets each term's codes by the same `keep` used for the arm rows.
    # The positive-cover RE is per visit, so it needs per-visit cover
    # (cover_aggregate = "none"); a cell-aggregated cover arm has one row per unit
    # and no per-visit grouping.
    if (!is.null(det_re_parse)) {
      model$re_det <- if (isTRUE(model$ragged))
        .occu_cover_obs_re_design_ragged(det_re_parse, data, vd_det$visits,
                                         model$site_of_visit, "detection")
      else
        .occu_cover_obs_re_design(det_re_parse, data, vd_det$visits, model$valid,
                                  model$n_sites, model$max_visits, "detection")
    }
    if (!is.null(pos_re_parse)) {
      if (cover_aggregate != "none") {
        stop("occu_cover(): a random effect on the positive-cover arm needs ",
             "per-visit cover (cover_aggregate = \"none\"); it cannot map onto ",
             "cell-aggregated cover rows (one per unit).", call. = FALSE)
      }
      model$re_pos <- if (isTRUE(model$ragged))
        .occu_cover_obs_re_design_ragged(pos_re_parse, data, vd_pos$visits,
                                         model$site_of_visit, "positive cover")
      else
        .occu_cover_obs_re_design(pos_re_parse, data, vd_pos$visits, model$valid,
                                  model$n_sites, model$max_visits, "positive cover")
    }

    # joint_coupled (3-arm nested-Laplace via tulpa's cell_coupling spec) is the
    # default: outer-grid integration over (sigma, alpha [, sigma_trend,
    # alpha_trend]) with inner Newton driven by the occu_cover_{lognormal,beta}
    # cell-coupling spec. 150-300x faster than v3 at N=100 and reliably completes
    # at N=200+ where v3 trips on a missing-value compare in its outer BFGS. v3
    # pure-R nested-Laplace and v2's joint Laplace stay reachable via
    # control$engine = "v3_nested" / "v2_joint" as debug escape hatches; both
    # take only the single intercept field.
    correlated <- isTRUE(spatial_info$correlated)
    engine_pick <- control[["engine"]] %||% "joint"
    control[["engine"]] <- NULL
    # The coupling now lives in the formula (copy()), so the default engine is
    # "joint"; "joint_coupled" is the deprecated alias for the same fitter.
    if (identical(engine_pick, "joint_coupled")) {
      message("occu_cover(): control$engine = \"joint_coupled\" is deprecated; ",
              "use \"joint\" (the cross-arm coupling now lives in the positive ",
              "formula via copy()).")
      engine_pick <- "joint"
    }
    if (correlated && engine_pick %in% c("v2_joint", "v3_nested")) {
      stop(sprintf(paste0(
        "occu_cover(): a correlated spatial bar (`|`, free-Sigma MCAR) needs ",
        "the default joint_coupled engine; the \"%s\" escape hatch couples a ",
        "single shared field only."), engine_pick), call. = FALSE)
    }
    if (engine_pick %in% c("v2_joint", "v3_nested")) {
      if (length(spatial_info$armspec)) {
        stop(sprintf(paste0(
          "occu_cover() an arm-specific field (to = \"positive\" / ",
          "\"detection\") needs the default joint_coupled engine; the \"%s\" ",
          "escape hatch couples a single shared field only."), engine_pick),
          call. = FALSE)
      }
      if (!is.null(re_spec) || !is.null(model$re_det) || !is.null(model$re_pos)) {
        stop(sprintf(paste0(
          "occu_cover() per-group RE needs the default joint_coupled engine; ",
          "the \"%s\" escape hatch has no RE block."),
          engine_pick), call. = FALSE)
      }
      # The v2/v3 escape hatches model per-visit cover only; cell-aggregated
      # cover is a joint_coupled feature. An explicit request errors; the bare
      # default falls back to per-visit on these engines.
      if (cover_aggregate != "none") {
        if (agg_explicit) {
          stop(sprintf(paste0(
            "occu_cover() cell-aggregated cover (cover_aggregate = \"%s\") is ",
            "wired on the default joint_coupled engine; the \"%s\" escape hatch ",
            "models per-visit cover only."), cover_aggregate, engine_pick),
            call. = FALSE)
        }
        model$cover_aggregate <- "none"
      }
      if (length(fields) > 1L) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" couples a single shared field; ",
          "weighted SVC field(s) need the default joint_coupled engine."),
          engine_pick), call. = FALSE)
      }
      if (!is.null(gv)) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" binds the field 1:1 to sites and does ",
          "not support group_var; use the default joint_coupled engine."),
          engine_pick), call. = FALSE)
      }
      fit_args <- c(list(model = model, adj = base_graph, priors = priors),
                    control)
      fitter <- if (engine_pick == "v2_joint") .tobs_fit_occu_cover_spatial
                else .tobs_fit_occu_cover_nested
      return(do.call(fitter, fit_args))
    }
    fit_args <- c(list(model = model, fields = fields, priors = priors,
                       re_spec = re_spec, correlated = correlated,
                       pos_armspec = spatial_info$armspec[["pos"]],
                       det_armspec = spatial_info$armspec[["p"]]),
                  control)
    return(do.call(.tobs_fit_occu_cover_joint_coupled, fit_args))
  }

  # Non-spatial NUTS: sample the exact two-state coefficient marginal (the
  # in-tree FullGradFn), warm-started at the Laplace mode. Other non-spatial
  # routes (only "laplace" here) fit the direct Laplace optim.
  if (identical(engine, "nuts")) {
    return(do.call(.tobs_fit_occu_cover_nuts,
                   c(list(model = model, priors = priors), control)))
  }

  fit_args <- list(model = model, method = engine, priors = priors)
  fit_args <- c(fit_args, control)
  do.call(.tobs_fit_occu_cover, fit_args)
}


# ---------------------------------------------------------------------------
# Pointwise log-likelihood (WAIC / PSIS-LOO) -- gcol33/tulpaObs#26
# ---------------------------------------------------------------------------

# Pointwise log-likelihood [n_draws x n_sites] for an occu_cover() fit: the
# per-site marginal log-likelihood (latent occupancy state integrated out)
# evaluated at each posterior draw. The spatial nested-Laplace fit samples the
# grid-integrated joint (betas + shared field) via the joint substrate; the
# non-spatial Laplace fit reuses the stored coefficient draws (no field). Both
# feed the same per-draw site-likelihood accumulation.
# Per-draw arm coefficients + dispersion + shared-field contributions for an
# occu_cover() fit. Joint path: sample the grid-integrated posterior and
# accumulate each block's field on the occupancy / cover arms. Non-spatial path:
# read the stored coefficient draws (no field). Returns a list consumed by both
# the pointwise log-likelihood and the posterior-mean plug-in.
.tobs_occu_cover_components <- function(object, n.draws = 1000L) {
  model   <- object$model
  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_det   <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_sites <- model$n_sites

  if (!is.null(.tobs_joint_fit(object))) {
    bundle <- .tobs_joint_draws(object, n = n.draws)
    S      <- bundle$n
    b_occ  <- bundle$b$occ
    b_det  <- bundle$b$det
    b_pos  <- bundle$b$pos
    disp   <- bundle$disp
    # Each site reads its field node (cell) via site_cell; the field draws
    # blk$z carry one column per node, so blk$z[, units] broadcasts the shared
    # cell field across that cell's sites. Identity when no group_var.
    units  <- model$site_cell %||% seq_len(n_sites)
    wfun   <- function(nm) {
      if (!nm %in% names(model$data)) {
        stop("occu_cover WAIC: trend-field weight column '", nm,
             "' is not in the fitted data.", call. = FALSE)
      }
      as.numeric(model$data[[nm]])
    }
    field_occ <- matrix(0, n_sites, S)
    field_pos <- matrix(0, n_sites, S)
    for (blk in bundle$blocks) {
      z_unit <- blk$z[, units, drop = FALSE]   # [S x n_sites]
      w <- if (is.null(blk$weight)) rep(1, n_sites) else wfun(blk$weight)
      field_occ <- field_occ + t(z_unit * blk$amp_occ) * w
      field_pos <- field_pos + t(z_unit * blk$amp_pos) * w
    }
  } else {
    draws <- object$draws
    if (is.null(draws) || !is.matrix(draws)) {
      stop("occu_cover WAIC: the fit carries no posterior draw matrix.",
           call. = FALSE)
    }
    if (!is.null(n.draws) && n.draws < nrow(draws)) {
      draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
    }
    b_occ <- draws[, seq_len(p_occ), drop = FALSE]
    b_det <- draws[, p_occ + seq_len(p_det), drop = FALSE]
    b_pos <- draws[, p_occ + p_det + seq_len(p_pos), drop = FALSE]
    disp  <- exp(draws[, ncol(draws)])         # trailing log-dispersion column
    S     <- nrow(draws)
    fld <- .tobs_occu_cover_v3_field(object, n_sites, S)
    field_occ <- fld$field_occ
    field_pos <- fld$field_pos
  }
  list(b_occ = b_occ, b_det = b_det, b_pos = b_pos, disp = disp,
       field_occ = field_occ, field_pos = field_pos)
}

# Per-cell field draws for the v3 nested-Laplace occu_cover spatial path, which
# stores no joint object to sample but DOES carry the per-cell field posterior in
# `field_table` (z_mean / z_sd) plus the copy coefficient `alpha` (and, for a
# spatially-varying trend, `trend_field_table` / `alpha_trend` weighted by the
# trend covariate). The per-site pointwise log-likelihood for site i depends only
# on the field at that site's own cell, so per-cell marginal Gaussian draws
# N(z_mean, z_sd) give the exact per-observation marginal -- the cross-cell
# covariance never enters a single site's term. Returns the [n_sites x S]
# occupancy / cover field contributions (zero matrices when no field table is
# present, recovering the non-spatial behaviour).
.tobs_occu_cover_v3_field <- function(object, n_sites, S) {
  zero <- matrix(0, n_sites, S)
  ft <- object$field_table
  if (is.null(ft) || is.null(ft$z_mean)) {
    return(list(field_occ = zero, field_pos = zero))
  }
  model    <- object$model
  site_cell <- model$site_cell %||% seq_len(n_sites)
  # Per-cell field draws (one column per draw), mapped to sites by site_cell.
  draw_cell_field <- function(tab) {
    n_cell <- nrow(tab)
    z <- matrix(stats::rnorm(n_cell * S, tab$z_mean, tab$z_sd), n_cell, S)
    z[site_cell, , drop = FALSE]                # [n_sites x S]
  }
  alpha <- object$means[["alpha"]] %||% object$spatial$alpha_mean %||% 0
  f_occ <- draw_cell_field(ft)
  field_occ <- f_occ
  field_pos <- alpha * f_occ

  # Spatially-varying trend field, weighted per cell by its trend covariate.
  tt <- object$trend_field_table
  if (!is.null(tt) && !is.null(tt$z_mean)) {
    alpha_tr <- object$means[["alpha_trend"]] %||% 0
    wname <- object$trend_weight %||% object$trend_weights[[1L]]
    w <- if (!is.null(wname) && wname %in% names(model$data)) {
      as.numeric(model$data[[wname]])
    } else rep(1, n_sites)
    f_tr <- draw_cell_field(tt) * w             # [n_sites x S], row-scaled
    field_occ <- field_occ + f_tr
    field_pos <- field_pos + alpha_tr * f_tr
  }
  list(field_occ = field_occ, field_pos = field_pos)
}

.tobs_ploglik_occu_cover <- function(object, n.draws = 1000L, n.threads = 1L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  .occu_cover_ploglik_core(object$model, c0$b_occ, c0$b_det, c0$b_pos,
                           c0$disp, c0$field_occ, c0$field_pos,
                           n_threads = n.threads)
}

# Pointwise log-likelihood at the posterior mean (length n_sites): the per-site
# marginal evaluated at the mean arm coefficients, mean dispersion, and mean
# shared field -- the same core driven by a one-draw mean.
.tobs_occu_cover_loglik_at_mean <- function(object, n.draws = 1000L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  as.numeric(.occu_cover_ploglik_core(
    object$model,
    matrix(colMeans(c0$b_occ), nrow = 1L),
    matrix(colMeans(c0$b_det), nrow = 1L),
    matrix(colMeans(c0$b_pos), nrow = 1L),
    mean(c0$disp),
    matrix(rowMeans(c0$field_occ), ncol = 1L),
    matrix(rowMeans(c0$field_pos), ncol = 1L)
  ))
}

# Accumulate the [S x n_sites] pointwise log-likelihood from per-draw arm
# coefficients (`b_occ` / `b_det` / `b_pos`, each [S x p_arm]), per-draw
# dispersion `disp` [S], and per-arm shared-field contributions `field_occ` /
# `field_pos` [n_sites x S] (zero matrices for a non-spatial fit). The detection
# and cover arms split into site-level and visit-level coefficient blocks exactly
# as the fitter packs them.
# Per-draw linear-predictor blocks for an occu_cover() fit: the occupancy
# `eta_psi_all` [n_sites x S] (with the shared field), and the site-level +
# optional visit-level detection / cover blocks the fitter packs. Split out so
# the pointwise log-likelihood, the posterior predictive check, and the PIT all
# build the per-draw psi / p / cover predictors from one source.
.occu_cover_eta_components <- function(model, b_occ, b_det, b_pos,
                                       field_occ, field_pos) {
  p_det_site <- ncol(model$X_det_site)
  p_pos_site <- ncol(model$X_pos_site)
  has_det_visit <- !is.null(model$X_det_visit)
  has_pos_visit <- !is.null(model$X_pos_visit)
  list(
    eta_psi_all = tcrossprod(model$X_occ, b_occ) + field_occ,
    eta_p_site_all = tcrossprod(model$X_det_site,
                                b_det[, seq_len(p_det_site), drop = FALSE]),
    eta_pos_site_all = tcrossprod(model$X_pos_site,
                                  b_pos[, seq_len(p_pos_site), drop = FALSE]) +
                       field_pos,
    eta_p_visit_all = if (has_det_visit) {
      tcrossprod(model$X_det_visit,
                 b_det[, p_det_site + seq_len(ncol(model$X_det_visit)),
                       drop = FALSE])
    } else NULL,
    eta_pos_visit_all = if (has_pos_visit) {
      tcrossprod(model$X_pos_visit,
                 b_pos[, p_pos_site + seq_len(ncol(model$X_pos_visit)),
                       drop = FALSE])
    } else NULL,
    has_det_visit = has_det_visit, has_pos_visit = has_pos_visit
  )
}

# Draw `d`'s occupancy predictor and the [n_sites x max_visits] detection /
# cover linear predictors, folding the visit-level block in site-major order.
.occu_cover_draw_eta <- function(comp, d, n_sites, max_visits) {
  p_eta <- matrix(comp$eta_p_site_all[, d], n_sites, max_visits)
  if (comp$has_det_visit) {
    p_eta <- p_eta + matrix(comp$eta_p_visit_all[, d], n_sites, max_visits,
                            byrow = TRUE)
  }
  ep_mat <- matrix(comp$eta_pos_site_all[, d], n_sites, max_visits)
  if (comp$has_pos_visit) {
    ep_mat <- ep_mat + matrix(comp$eta_pos_visit_all[, d], n_sites, max_visits,
                              byrow = TRUE)
  }
  list(psi_eta = comp$eta_psi_all[, d], p_eta = p_eta, ep_mat = ep_mat)
}

# Compact (ragged) counterpart of .occu_cover_draw_eta: the detection / cover
# predictors are length-V (one per valid visit), built by broadcasting the
# site-level block onto each visit's site and adding the per-visit block --
# never a padded [n_sites x max_visits] matrix.
.occu_cover_draw_eta_ragged <- function(comp, d, site_of_visit) {
  p_eta <- comp$eta_p_site_all[site_of_visit, d]
  if (comp$has_det_visit) p_eta <- p_eta + comp$eta_p_visit_all[, d]
  ep <- comp$eta_pos_site_all[site_of_visit, d]
  if (comp$has_pos_visit) ep <- ep + comp$eta_pos_visit_all[, d]
  list(psi_eta = comp$eta_psi_all[, d], p_eta = p_eta, ep = ep)
}

# Compact (ragged) counterpart of .occu_cover_site_ll: the per-site occu-cover
# marginal from length-V per-visit predictors grouped by `model$site_of_visit`,
# returned as a length-n_sites vector. Algebraically identical to the dense
# version (same MacKenzie mixture, same per-visit Beta/lognormal cover term);
# rowsum() accumulates each site's visits instead of a matrix rowSums. Scoped to
# cover_aggregate = "none" (the only mode the ragged path builds).
.occu_cover_site_ll_ragged <- function(model, psi, p_vec, ep_vec, log_disp) {
  site    <- model$site_of_visit
  y       <- model$y_det_visit
  n_sites <- model$n_sites
  g       <- factor(site, levels = seq_len(n_sites))

  log_p   <- log(p_vec)
  log_1mp <- log(1 - p_vec)
  log_h_det <- ifelse(y == 1L, log_p, log_1mp)

  is_beta <- identical(model$positive %||% "lognormal", "beta")
  pos_ld  <- numeric(length(y))
  det     <- y == 1L
  if (any(det)) {
    pos_ld[det] <- .occu_cover_pos_logdens(model$y_pos_visit[det], ep_vec[det],
                                           exp(log_disp), is_beta)
  }

  sum_log_h  <- as.numeric(rowsum(log_h_det, g))
  sum_log1mp <- as.numeric(rowsum(log_1mp,   g))
  cover_term <- as.numeric(rowsum(pos_ld,    g))
  n_det      <- as.numeric(rowsum(as.numeric(y), g))

  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))
  det_ll    <- log_psi + sum_log_h + cover_term
  nodet_ll  <- .tobs_logsumexp2(log_psi + sum_log1mp, log_1mpsi)
  ifelse(n_det > 0, det_ll, nodet_ll)
}

# Available system RAM in bytes, or NA if it cannot be read without a hard
# dependency. Linux (the LiSC server) exposes /proc/meminfo MemAvailable; the
# optional `ps` package covers the other platforms. NA -> the caller falls back
# to a fixed default budget.
.occu_cover_free_ram <- function() {
  v <- tryCatch({
    if (file.exists("/proc/meminfo")) {
      ln <- grep("^MemAvailable:", readLines("/proc/meminfo", n = 60L), value = TRUE)
      if (length(ln)) as.numeric(sub("\\D+(\\d+).*", "\\1", ln[1L])) * 1024 else NA_real_
    } else if (requireNamespace("ps", quietly = TRUE)) {
      as.numeric(ps::ps_system_memory()[["avail"]])
    } else NA_real_
  }, error = function(e) NA_real_)
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}

# Draw-chunk size for the WAIC pointwise log-likelihood. The heavy transient is
# the two [n_plots x chunk] per-visit eta matrices (detection + cover), so bound
# `2 * n_plots * chunk * 8` bytes to a fraction of free RAM (a 4 GB default when
# RAM cannot be probed). WAIC is a sum over draws, so chunking is exact. Returns
# a chunk in 1..n_draws.
.occu_cover_waic_chunk <- function(n_plots, n_draws, frac = 0.4) {
  free   <- .occu_cover_free_ram()
  budget <- if (is.finite(free)) frac * free else 4e9
  per_draw <- 2 * max(as.numeric(n_plots), 1) * 8
  max(1L, min(as.integer(n_draws), as.integer(budget / per_draw)))
}

# Flatten a dense (padded [n_sites x max_visits]) no-aggregation occu_cover model
# to the ragged one-row-per-valid-visit form the C++ pointwise kernel consumes.
# The dense visit designs are site-major with n_sites * max_visits rows (cell
# (i, v) at row (i - 1) * max_visits + v), which is exactly the column-major
# position of that cell in t(valid); so the valid-cell selector indexes both the
# response and the visit designs. Cells are enumerated site-major, visit
# ascending -- the same order the dense rowSums accumulates -- so the kernel's
# per-site sums are byte-identical to .occu_cover_site_ll on the dense grid.
.occu_cover_dense_ragged <- function(model) {
  mv  <- model$max_visits
  sel <- which(t(model$valid))               # site-major, visit-ascending
  v_idx <- ((sel - 1L) %% mv) + 1L
  s_idx <- ((sel - 1L) %/% mv) + 1L
  cell  <- cbind(s_idx, v_idx)
  list(
    site_of_visit = as.integer(s_idx),
    y_det_visit   = as.integer(model$y[cell]),
    y_pos_visit   = as.numeric(model$y_pos[cell]),
    X_det_visit   = if (!is.null(model$X_det_visit))
                      model$X_det_visit[sel, , drop = FALSE] else NULL,
    X_pos_visit   = if (!is.null(model$X_pos_visit))
                      model$X_pos_visit[sel, , drop = FALSE] else NULL,
    V = length(sel)
  )
}

# `chunk` (draws per block) defaults to a memory-adaptive size; the [n_plots x
# chunk] eta matrices are the WAIC's memory peak, and processing draws in blocks
# keeps that bounded while the returned [S x n_sites] pointwise log-likelihood is
# byte-identical to the unchunked result (each draw's row depends only on that
# draw).
.occu_cover_ploglik_core <- function(model, b_occ, b_det, b_pos, disp,
                                     field_occ, field_pos, chunk = NULL,
                                     n_threads = 1L) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  S  <- nrow(b_occ)
  cl <- .tobs_clamp_eta
  is_ragged <- isTRUE(model$ragged)
  mode <- model$cover_aggregate %||% "none"

  # No-aggregation path (every ragged fit, and the dense grid without cover
  # aggregation): the C++ kernel mirrors .occu_cover_site_ll_ragged draw for
  # draw and parallelises over draws, with no draw-chunking (each draw's
  # per-visit predictors live in thread-private scratch, not the [V x n_draws]
  # transient .occu_cover_waic_chunk bounds). A dense grid is flattened to the
  # same one-row-per-valid-visit form (.occu_cover_dense_ragged), summed in the
  # same visit order as the dense rowSums, so both feed one kernel.
  if (identical(mode, "none")) {
    rg <- if (is_ragged) {
      list(site_of_visit = as.integer(model$site_of_visit),
           y_det_visit   = as.integer(model$y_det_visit),
           y_pos_visit   = as.numeric(model$y_pos_visit),
           X_det_visit   = model$X_det_visit, X_pos_visit = model$X_pos_visit,
           V = length(model$site_of_visit))
    } else .occu_cover_dense_ragged(model)
    empty_v <- function(m) if (is.null(m)) matrix(0, rg$V, 0L) else m
    return(cpp_occu_cover_ploglik_ragged(
      X_occ = model$X_occ, X_det_site = model$X_det_site,
      X_pos_site = model$X_pos_site,
      X_det_visit = empty_v(rg$X_det_visit),
      X_pos_visit = empty_v(rg$X_pos_visit),
      site_of_visit = rg$site_of_visit,
      y_det_visit   = rg$y_det_visit,
      y_pos_visit   = rg$y_pos_visit,
      b_occ = b_occ, b_det = b_det, b_pos = b_pos, disp = disp,
      field_occ = field_occ, field_pos = field_pos,
      is_beta   = identical(model$positive %||% "lognormal", "beta"),
      eta_bound = .TOBS_ETA_BOUND,
      n_threads = max(1L, as.integer(n_threads))))
  }

  # Aggregated (mean / median / latent) paths stay in R, draw-chunked
  # to bound the [n_plots x chunk] per-visit eta transient. Detected-unit cover
  # values are draw-invariant, so resolve them once (gcol33/tulpaObs#34).
  units <- if (identical(mode, "none")) NULL else .occu_cover_unit_cover(model)

  n_plots <- max(
    if (!is.null(model$X_det_visit)) nrow(model$X_det_visit) else 0L,
    if (!is.null(model$X_pos_visit)) nrow(model$X_pos_visit) else 0L,
    n_sites)
  if (is.null(chunk)) chunk <- .occu_cover_waic_chunk(n_plots, S)

  ll <- matrix(0, S, n_sites)
  for (st in seq.int(1L, S, by = chunk)) {
    idx  <- st:min(st + chunk - 1L, S)
    comp <- .occu_cover_eta_components(model,
              b_occ[idx, , drop = FALSE], b_det[idx, , drop = FALSE],
              b_pos[idx, , drop = FALSE],
              field_occ[, idx, drop = FALSE], field_pos[, idx, drop = FALSE])
    for (j in seq_along(idx)) {
      d <- idx[j]
      de    <- .occu_cover_draw_eta(comp, j, n_sites, max_visits)
      psi   <- stats::plogis(cl(de$psi_eta))
      p_mat <- stats::plogis(cl(de$p_eta))
      ll[d, ] <- .occu_cover_site_ll(model, psi, p_mat, de$ep_mat,
                                     log(disp[d]), units = units)
    }
  }
  ll
}


# ---------------------------------------------------------------------------
# Posterior predictive check + PIT (gcol33/tulpaObs#27)
# ---------------------------------------------------------------------------

# Posterior-predictive cover discrepancy for one draw, at the granularity the
# fitter optimised (gcol33/tulpaObs#34). Returns the observed and replicated
# positive-part discrepancy (`obs` / `rep`) the PPC adds to the detection term.
# `ep_mat` is the draw's [n_sites x max_visits] cover predictor; `disp` its
# dispersion (beta precision / lognormal residual SD, or the cover-RE SD under
# "latent"); `units` the per-unit detected covers from .occu_cover_unit_cover().
# Mode branches mirror .occu_cover_cover_term: "none" scores the per-visit cells;
# "mean" / "median" one aggregated cover per detected unit at the unit predictor
# and dispersion the fit held; "latent" the per-unit covers replicated through
# the shared cover RE (u_i ~ N(0, sigma_u^2)) at the within-unit residual disp2.
.occu_cover_ppc_cover <- function(model, ep_mat, disp, units, is_beta, stat_fn,
                                  cl) {
  draw_pos <- function(eta, d) {
    if (is_beta) {
      mu <- stats::plogis(cl(eta))
      stats::rbeta(length(eta), mu * d, (1 - mu) * d)
    } else {
      exp(stats::rnorm(length(eta), eta, d))
    }
  }
  mean_pos <- function(eta, d) {
    if (is_beta) stats::plogis(cl(eta)) else exp(cl(eta) + d^2 / 2)
  }
  mode <- model$cover_aggregate %||% "none"

  if (identical(mode, "none")) {
    pos_mask <- model$valid & (model$y == 1L)
    if (!any(pos_mask)) return(c(obs = 0, rep = 0))
    Epos <- mean_pos(ep_mat, disp)
    yrep <- matrix(draw_pos(as.vector(ep_mat), disp), nrow(ep_mat), ncol(ep_mat))
    return(c(obs = stat_fn(model$y_pos[pos_mask], Epos[pos_mask]),
             rep = stat_fn(yrep[pos_mask],        Epos[pos_mask])))
  }

  ps <- units$pos_site
  if (length(ps) == 0L) return(c(obs = 0, rep = 0))
  eta <- ep_mat[ps, 1L]                     # unit-level cover predictor

  if (identical(mode, "latent")) {
    disp2   <- model$cover_latent_disp2
    m       <- lengths(units$vals)
    unit_of <- rep(seq_along(ps), m)
    v_all   <- unlist(units$vals, use.names = FALSE)
    eta_all <- eta[unit_of]
    u_all   <- stats::rnorm(length(ps), 0, disp)[unit_of]
    e_all   <- mean_pos(eta_all, disp2)
    return(c(obs = stat_fn(v_all, e_all),
             rep = stat_fn(draw_pos(eta_all + u_all, disp2), e_all)))
  }

  aggfun <- if (identical(mode, "median")) stats::median else mean
  yv   <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
  Epos <- mean_pos(eta, disp)
  c(obs = stat_fn(yv, Epos), rep = stat_fn(draw_pos(eta, disp), Epos))
}

# Posterior predictive check for an occu_cover() fit. Per draw, the latent
# occupancy z is sampled from its full conditional given the detection history
# (the spOccupancy construction), detection replicates from Bernoulli(z p), and
# the cover replicate is built at the granularity the fit used (per-visit for
# cover_aggregate = "none", one aggregated cover per detected unit for "mean" /
# "median", and the shared cover-RE marginal for "latent"; gcol33/tulpaObs#34).
# The discrepancy is a Freeman-Tukey (or chi-squared) sum over the detection
# cells plus the positive-part term, returning a Bayesian p-value.
.tobs_ppc_occu_cover <- function(object,
                                 fit.stat = c("freeman-tukey", "chi-squared"),
                                 n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  model    <- object$model
  positive <- model$positive %||% "lognormal"
  is_beta  <- identical(positive, "beta")
  c0   <- .tobs_occu_cover_components(object, n.samples)
  comp <- .occu_cover_eta_components(model, c0$b_occ, c0$b_det, c0$b_pos,
                                     c0$field_occ, c0$field_pos)
  disp <- c0$disp
  S <- nrow(c0$b_occ)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  y <- model$y; valid <- model$valid
  cl <- .tobs_clamp_eta
  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(o, e) sum((sqrt(o) - sqrt(e))^2, na.rm = TRUE)
  } else {
    function(o, e) sum((o - e)^2 / (e + 1e-10), na.rm = TRUE)
  }
  any_det <- rowSums(y * valid, na.rm = TRUE) > 0
  n_valid <- rowSums(valid)
  mode <- model$cover_aggregate %||% "none"

  # No-aggregation path: the per-draw simulation (latent z, detection replicate,
  # cover replicate) runs in cpp_occu_cover_ppc, which draws from R's RNG stream
  # in the SAME order as the former R loop, so under a fixed seed the discrepancy
  # is byte-identical. Build the per-draw predictors (deterministic) here.
  if (identical(mode, "none")) {
    psi_all <- matrix(0, n_sites, S)
    p_all   <- matrix(0, n_sites, S * max_visits)
    ep_all  <- matrix(0, n_sites, S * max_visits)
    for (s in seq_len(S)) {
      de   <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
      cols <- (s - 1L) * max_visits + seq_len(max_visits)
      psi_all[, s]   <- stats::plogis(cl(de$psi_eta))
      p_all[, cols]  <- stats::plogis(cl(de$p_eta))
      ep_all[, cols] <- de$ep_mat
    }
    vint <- valid; storage.mode(vint) <- "integer"
    yint <- y;     storage.mode(yint) <- "integer"
    r <- cpp_occu_cover_ppc(psi_all, p_all, ep_all, yint, model$y_pos, vint,
                            as.integer(any_det), as.integer(n_valid), disp,
                            is_beta, identical(fit.stat, "freeman-tukey"))
    return(list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
                bayesian.p = mean(r$fit.y.rep > r$fit.y)))
  }

  # Aggregated (mean / median / latent) cover discrepancy: the detection
  # replicate plus the aggregated / latent cover replicate run in
  # cpp_occu_cover_ppc_agg with matched RNG order (byte-identical). The observed
  # aggregates / detected covers are draw-invariant and gathered once.
  units <- .occu_cover_unit_cover(model)
  psi_all <- matrix(0, n_sites, S)
  p_all   <- matrix(0, n_sites, S * max_visits)
  ep_all  <- matrix(0, n_sites, S * max_visits)
  for (s in seq_len(S)) {
    de   <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
    cols <- (s - 1L) * max_visits + seq_len(max_visits)
    psi_all[, s]   <- stats::plogis(cl(de$psi_eta))
    p_all[, cols]  <- stats::plogis(cl(de$p_eta))
    ep_all[, cols] <- de$ep_mat
  }
  vint <- valid; storage.mode(vint) <- "integer"
  yint <- y;     storage.mode(yint) <- "integer"
  ps0  <- as.integer(units$pos_site - 1L)
  if (identical(mode, "latent")) {
    mode_code <- 2L
    vals_flat <- as.numeric(unlist(units$vals, use.names = FALSE))
    unit_off  <- as.integer(c(0L, cumsum(lengths(units$vals))))
    yv <- numeric(0); disp2 <- model$cover_latent_disp2 %||% 0
  } else {
    mode_code <- 1L
    aggfun <- if (identical(mode, "median")) stats::median else mean
    yv <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
    vals_flat <- numeric(0); unit_off <- 0L; disp2 <- 0
  }
  r <- cpp_occu_cover_ppc_agg(psi_all, p_all, ep_all, yint, vint,
    as.integer(any_det), as.integer(n_valid), disp, mode_code, ps0,
    as.numeric(yv), vals_flat, unit_off, disp2, is_beta,
    identical(fit.stat, "freeman-tukey"))
  list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
       bayesian.p = mean(r$fit.y.rep > r$fit.y))
}

# Per-draw CDF limits for the occu_cover() per-site detection summary
# (any-detection vs all-zero), marginalized over the latent occupancy state with
# the shared field folded in per site. Returns the [S x n_sites] lower / upper
# CDF limits of the ordered detected / non-detected outcome -- the randomized-PIT
# building block shared by the posterior PIT and the LOO-PIT.
.occu_cover_pit_cdf_limits <- function(object, n.samples) {
  model <- object$model
  c0    <- .tobs_occu_cover_components(object, n.samples)
  valid <- model$valid; y <- model$y
  any_det <- as.integer(rowSums(y * valid, na.rm = TRUE) > 0)
  vint <- valid; storage.mode(vint) <- "integer"
  empty_v <- function(m) if (is.null(m))
    matrix(0, model$n_sites * model$max_visits, 0L) else m
  # The per-draw detection-summary CDF limits are deterministic; the former R
  # loop now runs in cpp_occu_cover_cdf_limits, parallel over draws.
  cpp_occu_cover_cdf_limits(model$X_occ, model$X_det_site,
                            empty_v(model$X_det_visit), c0$b_occ, c0$b_det,
                            c0$field_occ, vint, any_det, 1L)
}

# Randomized PIT for an occu_cover() fit, on the per-site detection summary
# (any-detection vs all-zero) marginalized over the latent occupancy state, with
# the shared field projected per site. The detected / non-detected outcome is the
# ordered event; the left and right CDF limits feed the engine's randomized PIT.
.tobs_pit_occu_cover <- function(object, n.samples = 250) {
  lim <- .occu_cover_pit_cdf_limits(object, n.samples)
  tulpa::tulpa_pit(lim$cdf_upper, cdf_lower = lim$cdf_lower)
}

# Leave-one-out randomized PIT for an occu_cover() fit -- the INLA `cpo$pit`
# analogue. Per site, the draws are reweighted by their PSIS leave-one-out
# importance weights (so the predictive distribution excludes that site's own
# contribution), then the LOO-weighted CDF limits feed the randomized PIT. The
# loglik matrix supplying the weights is the field-folded pointwise log-
# likelihood, so the LOO-PIT inherits the full-model predictor. A site whose PSIS
# weights are degenerate (k unavailable) falls back to the equal-weight posterior
# CDF for that site.
.tobs_loo_pit_occu_cover <- function(object, n.samples = 1000L, ll = NULL) {
  lim <- .occu_cover_pit_cdf_limits(object, n.samples)
  if (is.null(ll)) ll <- .tobs_ploglik_occu_cover(object, n.samples)
  .tobs_loo_pit_from_limits(ll, lim$cdf_lower, lim$cdf_upper)
}

# LOO-weighted randomized PIT from a pointwise log-likelihood matrix `ll`
# [S x N] and the per-draw CDF limits `Fl` / `Fu` [S x N] of the ordered outcome.
# For each observation the PSIS leave-one-out weights w_is (from tulpa_psis on
# -ll[, i]) reweight the per-draw CDF limits to the LOO predictive limits
# Fl_loo_i = sum_s w_is Fl[s, i], Fu_loo_i likewise, then a single uniform draw
# placed between them gives the randomized LOO-PIT. The single source of truth
# behind every family's LOO-PIT so the construction matches INLA's cpo$pit.
.tobs_loo_pit_from_limits <- function(ll, Fl, Fu) {
  # Per-observation PSIS leave-one-out weighting of the CDF limits + a uniform
  # jitter, batched in tulpa's cpp_psis_loo_pit (PSIS columns parallel, the runif
  # in index order), so it is byte-identical to the former per-column R loop.
  tail_len <- getFromNamespace(".psis_tail_len", "tulpa")(nrow(ll), NULL)
  getFromNamespace("cpp_psis_loo_pit", "tulpa")(
    ll, Fl, Fu, as.integer(tail_len), 1L)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate joint occupancy-detection + cover data
#'
#' Per-cell mixture: latent z_i ~ Bernoulli(psi_i), per-visit detection
#' y_ij | z_i = 1 ~ Bernoulli(p_ij), per-visit cover y_pos_ij | y_ij = 1
#' drawn from `positive` (`"beta"` or `"lognormal"`) on the cover-arm linear
#' predictor. Used by the recovery test and as the generator for the
#' joint occupancy + cover hurdle family (see [occu_cover()]).
#'
#' When `adj` is supplied (a square adjacency matrix), an ICAR field
#' `f[1..N]` is drawn from `MVN(0, Q^-)` (with sum-to-zero constraint),
#' and the linear predictors become
#'
#'     eta_psi_i = X_psi[i, ] %*% beta_occ + sigma * f[i]
#'     eta_pos_ij = X_pos[i, ] %*% beta_pos + alpha * sigma * f[i]
#'
#' matching the v2 nested-Laplace fit's parameterisation.
#'
#' @param N Number of sites (cells).
#' @param J Number of visits per site.
#' @param n_occ_covs,n_det_covs,n_pos_covs Number of covariates on each arm
#'   (drawn IID standard normal).
#' @param beta_occ,beta_p,beta_pos Coefficient vectors c(intercept, slopes).
#'   Defaults pick weakly-informative values: psi intercept at logit(0.4),
#'   p intercept at logit(0.5), cover intercept on the appropriate link.
#' @param response `"beta"` or `"lognormal"`.
#' @param phi Beta precision when `positive = "beta"` (default 30).
#' @param sigma_pos Lognormal residual SD when `positive = "lognormal"`
#'   (default 0.4).
#' @param adj Optional N x N adjacency matrix. When supplied, generates the
#'   shared ICAR field; when NULL, the simulator is non-spatial (matches v1).
#' @param sigma Spatial field amplitude (used only when `adj` is supplied).
#' @param alpha Cover-arm scaling on the shared field (used only when `adj`
#'   is supplied). 1.0 = arms see the field identically; positive = same sign,
#'   negative = opposite.
#' @param trend Logical; when `TRUE` (and `adj` is supplied) a SECOND shared
#'   ICAR field `f2` (a spatially-varying temporal trend) is generated on the
#'   same graph, weighted by a per-cell covariate `time` drawn IID standard
#'   normal. The trend enters the occupancy and cover predictors as
#'   `sigma_trend * time_i * f2[i]` (occupancy) and
#'   `alpha_trend * sigma_trend * time_i * f2[i]` (cover); the detection
#'   predictor is unaffected. The `time` covariate is per-cell and broadcast
#'   to every visit of that cell.
#' @param sigma_trend Trend-field amplitude (used only when `trend = TRUE`).
#' @param alpha_trend Cover-arm scaling on the trend field (used only when
#'   `trend = TRUE`).
#' @param pos_field Logical; when `TRUE` (and `adj` is supplied) draw an
#'   INDEPENDENT areal field on the cover (positive) arm only -- an intercept
#'   field plus a time-weighted trend field, each unrelated to the occupancy
#'   field and with no alpha copy (gcol33/tulpaObs#110). Adds a `time` column to
#'   the returned `data` and reports `g0` / `g1` (the two fields) and their SDs in
#'   `truth`. Fit by placing `spatial(~ 1 + time || cell, graph = adj)` in the
#'   `positive` formula.
#' @param sigma_pos_int Cover-arm intercept-field SD (used only when
#'   `pos_field = TRUE`).
#' @param sigma_pos_trend Cover-arm trend-field SD (used only when
#'   `pos_field = TRUE`).
#' @param re_det_groups Optional integer `>= 2`: the number of levels of a
#'   per-visit detection random intercept (a `habitat` factor on `visit_data`,
#'   levels `hab1..K`), drawn `b_g ~ N(0, sigma_re_p^2)` and centred. `NULL`
#'   (default) adds no detection RE. Recover it with
#'   `detection = ~ ... + (1 | habitat)`.
#' @param sigma_re_p SD of the `re_det_groups` random intercept (default 0.7).
#' @param re_det Optional named list of FURTHER per-visit detection random
#'   effects, for crossed / nested / slope designs. Each element
#'   `list(K =, sigma =, prefix =, nested_in =, slope_cov =, sigma_slope =, rho =)`
#'   adds a factor column (levels `<prefix>1..K`). Without `slope_cov` it is a
#'   centred `N(0, sigma^2)` random intercept; `nested_in = "<name>"` nests its
#'   codes within a previously listed grouping (matching `(1 | parent/child)`),
#'   otherwise crossed. With `slope_cov = "<column>"` (a per-visit covariate,
#'   generated `N(0, slope_sd)` if absent, `slope_sd` default 1) it is a random
#'   slope: a slope-only uncorrelated block when `rho` is unset, or a correlated
#'   intercept + slope block (covariance from `sigma` / `sigma_slope` / `rho`)
#'   when `rho` is given.
#'   Truth is returned in `truth$re_det[[name]]` (named by the level label a fit
#'   reconstructs): `b` / `b_slope` BLUP vectors, or the `B` BLUP matrix plus
#'   `s0` / `s1` / `rho` for a correlated slope.
#' @param seed Optional integer seed.
#' @return A list with `y` (N x J detection matrix), `y_pos` (N x J cover
#'   matrix, NA where not detected), `data` (per-site covariate frame, gaining
#'   a `time` column when `trend = TRUE`), `visit_data` (per-visit covariate
#'   frame, N*J rows in site-major order), and `truth` (the coefficients,
#'   dispersion, and field(s) if generated; `f2`, `sigma_trend`, `alpha_trend`,
#'   and `time` when `trend = TRUE`).
#' @export
simulate_occu_cover <- function(N             = 200L,
                                 J             = 4L,
                                 n_occ_covs    = 1L,
                                 n_det_covs    = 1L,
                                 n_pos_covs    = 1L,
                                 beta_occ      = NULL,
                                 beta_p        = NULL,
                                 beta_pos      = NULL,
                                 positive      = c("lognormal", "beta"),
                                 phi           = 30,
                                 sigma_pos     = 0.4,
                                 adj           = NULL,
                                 sigma         = 0.6,
                                 alpha         = 1.0,
                                 trend         = FALSE,
                                 sigma_trend   = 0.6,
                                 alpha_trend   = 1.0,
                                 pos_field       = FALSE,
                                 sigma_pos_int   = 0.5,
                                 sigma_pos_trend = 0.6,
                                 det_field       = FALSE,
                                 sigma_p_int     = 0.5,
                                 sigma_p_trend   = 0.6,
                                 re_det_groups = NULL,
                                 sigma_re_p    = 0.7,
                                 re_det        = NULL,
                                 seed          = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N); J <- as.integer(J)

  if (is.null(beta_occ)) beta_occ <- c(stats::qlogis(0.4), stats::runif(n_occ_covs, -0.5, 0.5))
  if (is.null(beta_p))   beta_p   <- c(0, stats::runif(n_det_covs, -0.5, 0.5))
  if (is.null(beta_pos)) {
    pos_int <- if (positive == "beta") stats::qlogis(0.3) else log(0.1)
    beta_pos <- c(pos_int, stats::runif(n_pos_covs, -0.5, 0.5))
  }

  # Optional shared ICAR field(s). Draw each f as N(0, Q^-) via the
  # eigendecomposition of Q; the constant (null) component is dropped, giving a
  # zero-mean draw on the constrained space, then divide by sqrt(scale_q) so the
  # field has geo-mean marginal variance 1 (the Sorbye-Rue convention; `sigma *
  # f` then has geo-mean marginal SD sigma, matching INLA's `scale.model = TRUE`
  # and the fitter's parameterisation).
  f  <- numeric(N)
  f2 <- numeric(N)
  g0 <- numeric(N)   # arm-specific cover intercept field (gcol33/tulpaObs#110)
  g1 <- numeric(N)   # arm-specific cover trend field
  h0 <- numeric(N)   # arm-specific detection intercept field
  h1 <- numeric(N)   # arm-specific detection trend field
  time_cov <- numeric(N)
  if (!is.null(adj)) {
    if (!is.matrix(adj) || nrow(adj) != N || ncol(adj) != N) {
      stop("adj must be an N x N adjacency matrix.", call. = FALSE)
    }
    Q       <- .occu_cover_icar_Q(adj)
    scale_q <- .occu_cover_icar_scale(adj)
    eig <- eigen(Q, symmetric = TRUE)
    keep <- eig$values > 1e-8
    draw_field <- function() {
      z_white <- stats::rnorm(sum(keep))
      fk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                         (z_white / sqrt(eig$values[keep])))
      fk <- fk - mean(fk)
      fk / sqrt(scale_q)
    }
    f <- draw_field()
    # A time covariate is needed by either the shared trend field or the
    # arm-specific cover trend field.
    if (isTRUE(trend) || isTRUE(pos_field) || isTRUE(det_field)) {
      time_cov <- as.numeric(scale(stats::rnorm(N)))
    }
    if (isTRUE(trend)) f2 <- draw_field()
    # Arm-specific cover field(s) (gcol33/tulpaObs#110): an INDEPENDENT cover-arm
    # intercept field g0 and time-weighted trend field g1, each unrelated to the
    # occupancy field f. They enter the cover linear predictor only (no psi
    # contribution, no alpha copy), so delta_cover_cond carries a spatial
    # structure the occupancy field's alpha copy cannot express.
    if (isTRUE(pos_field)) {
      g0 <- draw_field()
      g1 <- draw_field()
    }
    # Arm-specific detection field(s): an INDEPENDENT detection-arm intercept field
    # h0 and time-weighted trend field h1, unrelated to the occupancy and cover
    # fields. They enter the detection linear predictor only (no psi / cover
    # contribution, no copy), so p varies spatially on its own.
    if (isTRUE(det_field)) {
      h0 <- draw_field()
      h1 <- draw_field()
    }
  }

  # Site-level covariates (psi predictor).
  occ_covs <- data.frame(matrix(stats::rnorm(N * n_occ_covs), N, n_occ_covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n_occ_covs))
  X_occ <- stats::model.matrix(~ ., occ_covs)
  eta_psi <- as.vector(X_occ %*% beta_occ) + sigma * f
  if (!is.null(adj) && isTRUE(trend)) {
    eta_psi <- eta_psi + sigma_trend * time_cov * f2
  }
  psi <- stats::plogis(eta_psi)
  z_state <- stats::rbinom(N, 1L, psi)

  # Visit-level covariates (p and cover predictors). Same draw used for both
  # arms, mirroring how `tobs_data()`'s `det.covs` matrices feed both formulas.
  det_covs <- data.frame(matrix(stats::rnorm(N * J * n_det_covs), N * J, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  pos_covs <- data.frame(matrix(stats::rnorm(N * J * n_pos_covs), N * J, n_pos_covs))
  names(pos_covs) <- paste0("pos_cov", seq_len(n_pos_covs))
  visit_data <- cbind(det_covs, pos_covs)

  X_p   <- stats::model.matrix(~ ., det_covs)
  X_pos <- stats::model.matrix(~ ., pos_covs)
  eta_p   <- as.vector(X_p   %*% beta_p)
  eta_pos_base <- as.vector(X_pos %*% beta_pos)

  # Optional per-visit detection random effects (gcol33/tulpaObs#102, #103). Each
  # grouping is a categorical visit-level factor (e.g. an EUNIS habitat class)
  # with a random intercept b_g ~ N(0, sigma^2) on the detection linear
  # predictor; the factor rides `visit_data`, so a fit reads it via `visits` and
  # `detection = ~ ... + (1 | <factor>)`. `re_det_groups` / `sigma_re_p` set the
  # first ("habitat") grouping (back-compat); `re_det` is a named list of further
  # groupings, each `list(K =, sigma =, prefix =, nested_in =)`, for CROSSED
  # (`nested_in = NULL`) or NESTED (`nested_in = "<parent>"`, sub-codes nested
  # within the parent factor's codes -- matching `(1 | parent/child)`) designs.
  # BLUPs are centred so the detection intercept stays identified. Truth BLUPs are
  # stored NAMED by the level label a fit reconstructs (the interaction label for
  # a nested grouping), so recovery checks align by name, not factor sort order.
  grp_specs <- list()
  if (!is.null(re_det_groups)) {
    grp_specs[["habitat"]] <- list(K = as.integer(re_det_groups),
                                   sigma = sigma_re_p, prefix = "hab",
                                   nested_in = NULL)
  }
  if (!is.null(re_det)) {
    for (nm in names(re_det)) {
      s <- re_det[[nm]]
      grp_specs[[nm]] <- list(K = as.integer(s$K), sigma = s$sigma %||% 0.7,
                              prefix = s$prefix %||% nm, nested_in = s$nested_in,
                              # Random-slope fields: `slope_cov` names a per-visit
                              # covariate column (generated N(0, 1) if absent). With
                              # `rho` set it is a correlated intercept + slope block
                              # (Sigma from sigma / sigma_slope / rho); without, a
                              # slope-only uncorrelated block.
                              slope_cov = s$slope_cov,
                              sigma_slope = s$sigma_slope,
                              rho = s$rho, slope_sd = s$slope_sd)
    }
  }
  re_codes  <- list()    # within-grouping per-visit code (for nesting)
  re_truth  <- list()
  b_p_re <- NULL; re_det_levels <- NULL
  for (nm in names(grp_specs)) {
    s <- grp_specs[[nm]]
    if (s$K < 2L)
      stop(sprintf("re_det grouping '%s' needs K >= 2.", nm), call. = FALSE)
    sub <- sample.int(s$K, N * J, replace = TRUE)
    re_codes[[nm]] <- sub
    visit_data[[nm]] <- factor(paste0(s$prefix, sub),
                               levels = paste0(s$prefix, seq_len(s$K)))
    if (is.null(s$nested_in)) {
      code   <- sub
      labels <- paste0(s$prefix, seq_len(s$K))
    } else {
      parent <- grp_specs[[s$nested_in]]
      pcode  <- re_codes[[s$nested_in]]
      code   <- (pcode - 1L) * s$K + sub                 # unique (parent, sub)
      pc     <- rep(seq_len(parent$K), each = s$K)
      sc     <- rep(seq_len(s$K),      times = parent$K)
      labels <- paste0(parent$prefix, pc, ".", s$prefix, sc)  # interaction label
    }
    if (is.null(s$slope_cov)) {
      # Random intercept.
      b <- stats::rnorm(length(labels), 0, s$sigma); b <- b - mean(b)
      names(b) <- labels
      eta_p <- eta_p + b[code]
      re_truth[[nm]] <- list(b = b, sigma = s$sigma, levels = labels,
                             kind = "intercept")
      if (identical(nm, "habitat")) {                    # back-compat truth slots
        b_p_re <- unname(b); re_det_levels <- labels
      }
    } else {
      # Random slope on a per-visit covariate (generated N(0, slope_sd) if
      # absent; `slope_sd` != 1 exercises the slope-covariate standardization).
      if (is.null(visit_data[[s$slope_cov]]))
        visit_data[[s$slope_cov]] <- stats::rnorm(N * J, 0, s$slope_sd %||% 1)
      xv <- as.numeric(visit_data[[s$slope_cov]])
      if (is.null(s$rho)) {
        # Slope-only uncorrelated block (0 + x | g).
        b1 <- stats::rnorm(length(labels), 0, s$sigma); b1 <- b1 - mean(b1)
        names(b1) <- labels
        eta_p <- eta_p + xv * b1[code]
        re_truth[[nm]] <- list(b_slope = b1, sigma = s$sigma, levels = labels,
                               cov = s$slope_cov, kind = "slope")
      } else {
        # Correlated intercept + slope block (1 + x | g): (b0, b1) ~ N(0, Sigma).
        s0 <- s$sigma; s1 <- s$sigma_slope %||% s$sigma; rho <- s$rho
        Sig <- matrix(c(s0^2, rho * s0 * s1, rho * s0 * s1, s1^2), 2L, 2L)
        L   <- t(chol(Sig))
        B2  <- matrix(stats::rnorm(2L * length(labels)), length(labels), 2L) %*% t(L)
        B2  <- sweep(B2, 2L, colMeans(B2))               # centre each coef
        rownames(B2) <- labels; colnames(B2) <- c("(Intercept)", s$slope_cov)
        eta_p <- eta_p + B2[code, 1L] + xv * B2[code, 2L]
        re_truth[[nm]] <- list(B = B2, s0 = s0, s1 = s1, rho = rho,
                               levels = labels, cov = s$slope_cov, kind = "corr")
      }
    }
  }

  y     <- matrix(0L, N, J)
  y_pos <- matrix(NA_real_, N, J)

  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      idx <- (i - 1L) * J + j
      if (z_state[i] == 1L) {
        eta_p_ij <- eta_p[idx]
        if (!is.null(adj) && isTRUE(det_field)) {
          eta_p_ij <- eta_p_ij + sigma_p_int * h0[i] +
                      sigma_p_trend * time_cov[i] * h1[i]
        }
        p_ij <- stats::plogis(eta_p_ij)
        d <- stats::rbinom(1L, 1L, p_ij)
        y[i, j] <- d
        if (d == 1L) {
          eta_pos_ij <- eta_pos_base[idx] + alpha * sigma * f[i]
          if (!is.null(adj) && isTRUE(trend)) {
            eta_pos_ij <- eta_pos_ij + alpha_trend * sigma_trend * time_cov[i] * f2[i]
          }
          if (!is.null(adj) && isTRUE(pos_field)) {
            eta_pos_ij <- eta_pos_ij +
              sigma_pos_int * g0[i] + sigma_pos_trend * time_cov[i] * g1[i]
          }
          if (positive == "beta") {
            mu <- stats::plogis(eta_pos_ij)
            y_pos[i, j] <- stats::rbeta(1L, mu * phi, (1 - mu) * phi)
          } else {
            y_pos[i, j] <- exp(stats::rnorm(1L, eta_pos_ij, sigma_pos))
          }
        }
      }
    }
  }

  occ_out <- occ_covs
  has_trend    <- !is.null(adj) && isTRUE(trend)
  has_posfield <- !is.null(adj) && isTRUE(pos_field)
  has_detfield <- !is.null(adj) && isTRUE(det_field)
  if (has_trend || has_posfield || has_detfield) occ_out$time <- time_cov
  # The arm-specific field bars index their graph node by a `cell` column.
  if (has_posfield || has_detfield) occ_out$cell <- seq_len(N)

  list(
    y          = y,
    y_pos      = y_pos,
    data       = occ_out,
    visit_data = visit_data,
    adj        = adj,
    truth      = list(
      beta_occ    = beta_occ,
      beta_p      = beta_p,
      beta_pos    = beta_pos,
      psi         = psi,
      z           = z_state,
      positive    = positive,
      phi         = if (positive == "beta")      phi       else NA_real_,
      sigma_pos   = if (positive == "lognormal") sigma_pos else NA_real_,
      f           = f,
      sigma       = if (!is.null(adj)) sigma else NA_real_,
      alpha       = if (!is.null(adj)) alpha else NA_real_,
      f2          = if (has_trend) f2          else NULL,
      time        = if (has_trend || has_posfield) time_cov else NULL,
      sigma_trend = if (has_trend) sigma_trend else NA_real_,
      alpha_trend = if (has_trend) alpha_trend else NA_real_,
      g0              = if (has_posfield) g0              else NULL,
      g1              = if (has_posfield) g1              else NULL,
      sigma_pos_int   = if (has_posfield) sigma_pos_int   else NA_real_,
      sigma_pos_trend = if (has_posfield) sigma_pos_trend else NA_real_,
      h0              = if (has_detfield) h0              else NULL,
      h1              = if (has_detfield) h1              else NULL,
      sigma_p_int     = if (has_detfield) sigma_p_int     else NA_real_,
      sigma_p_trend   = if (has_detfield) sigma_p_trend   else NA_real_,
      sigma_re_p  = if (!is.null(re_det_groups)) sigma_re_p else NA_real_,
      b_p_re      = b_p_re,
      re_det_levels = re_det_levels,
      re_det      = if (length(re_truth)) re_truth else NULL
    )
  )
}
