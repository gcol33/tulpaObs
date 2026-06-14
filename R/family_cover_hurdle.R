# =============================================================================
# family_cover_hurdle.R — Vegetation cover hurdle on the tulpa Laplace backend
#
# Phase 1a: lognormal-positive variant via two independent tulpa_laplace()
# calls (binomial on occur, gaussian on log(cover[occur==1])).
# Phase 1c: joint shared-field model via tulpa_nested_laplace_joint(),
# lognormal-positive on BYM2/ICAR/CAR_proper.
# Phase 1d: beta-positive on the joint engine. phi is profiled (pre-fit on
# the positive subset alone via tulpa_laplace_beta(), then plugged into the
# joint as fixed dispersion). Mirrors the sigma_pos handling for lognormal.
# Full posterior integration over phi is scheduled for Phase 3.
# =============================================================================


# ---------------------------------------------------------------------------
# Dispatcher (called from tobs())
# ---------------------------------------------------------------------------

.dispatch_cover <- function(formula, data, family, detection, y,
                            visits, engine, priors, control,
                            approx = "gaussian_laplace",
                            correction = "none", ...) {
  positive <- family$params$positive
  if (!positive %in% c("lognormal", "beta")) {
    stop("cover(positive = '", positive, "') is not supported. ",
         "Use 'lognormal' or 'beta'.", call. = FALSE)
  }
  # gibbs/mi are rejected centrally by the per-family method registry
  # (.tobs_family_methods), so `correction` is always "none" here.
  if (!is.null(detection)) {
    stop("`cover()` does not use a detection formula ",
         "(replicates = 'single'). Drop the `detection` argument.",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("`cover()` requires `y` (a length-N numeric vector of cover ",
         "in [0, 1]).", call. = FALSE)
  }
  # A spatially varying trend is model structure and so lives in the formula,
  # as a second weighted areal term (gcol33/tulpaObs#59). `control$trend` is
  # removed: control carries fitting behaviour only. `[[` (exact), never `$`.
  if (!is.null(control[["trend"]])) {
    stop("control$trend is no longer supported for cover hurdle models.\n",
         "Declare spatially varying trends directly in the formula, e.g.\n\n",
         "  ~ time.sc +\n",
         "    icar(graph = adj, group_var = \"cell_idx\") +\n",
         "    icar(graph = adj, weight = time.sc, group_var = \"cell_idx\")",
         call. = FALSE)
  }
  enc      <- encode_cover_hurdle(formula, data, y, positive = positive)
  temporal <- enc$temporal
  re       <- enc$re

  # NUTS: the non-spatial sampler over the exact two-arm coefficient marginal.
  # Any structured term (areal field, weighted trend, correlated / arm-specific
  # bar, temporal, re) is integrated on the nested-Laplace outer grid, not
  # sampled here, so reject it with a pointer rather than dropping it silently.
  if (identical(engine, "nuts")) {
    has_struct <- !is.null(enc$spatial_spec) || !is.null(enc$trend) ||
                  !is.null(enc$mcar) || !is.null(enc$armspec) ||
                  !is.null(temporal) || (!is.null(re) && length(re) > 0L)
    if (has_struct) {
      stop("cover() NUTS is the non-spatial sampler: a spatial / temporal / re ",
           "term in the formula is not yet wired for method = 'nuts'. Use ",
           "method = 'nested_laplace' for structured terms.", call. = FALSE)
    }
    return(.tobs_fit_cover_nuts_dispatch(formula, data, y, positive, family,
                                         priors, control))
  }

  has_multi <- !is.null(temporal) || (!is.null(re) && length(re) > 0L)
  if (has_multi && !identical(engine, "nested_laplace")) {
    stop("temporal()/re() terms in a cover() formula require ",
         "method = 'nested_laplace' or 'nested_laplace_sla' (got engine '",
         engine, "'). The single-Laplace path is fixed-effects + spatial only.",
         call. = FALSE)
  }
  if (!is.null(enc$trend) && !identical(engine, "nested_laplace")) {
    stop("a weighted areal trend term (icar(..., weight = )) in a cover() ",
         "formula requires method = 'nested_laplace' or 'nested_laplace_sla' ",
         "(got engine '", engine, "'). The single-Laplace path is ",
         "fixed-effects + a single intercept field only.", call. = FALSE)
  }
  if (!is.null(enc$mcar) && !identical(engine, "nested_laplace")) {
    stop("a correlated spatial bar (single `|`) in a cover() formula requires ",
         "method = 'nested_laplace' or 'nested_laplace_sla' (got engine '",
         engine, "'). The correlated (MCAR) coefficient fields integrate their ",
         "cross-covariance on the outer nested-Laplace grid.", call. = FALSE)
  }
  if (!is.null(enc$armspec) && !identical(engine, "nested_laplace")) {
    stop("an arm-specific spatial bar (single-arm `to`) in a cover() formula ",
         "requires method = 'nested_laplace' or 'nested_laplace_sla' (got ",
         "engine '", engine, "'). The separate per-arm latent fields integrate ",
         "each field's precision on the outer nested-Laplace grid ",
         "(gcol33/tulpaObs#65).", call. = FALSE)
  }

  if (identical(engine, "nested_laplace")) {
    return(decode_cover_hurdle_joint(
      fit_cover_hurdle_joint_nested(enc, data, positive, control,
                                    temporal = temporal, re = re,
                                    priors = priors),
      enc, family, approx = approx
    ))
  }

  fits <- fit_cover_hurdle(enc, positive, engine, priors, control)
  decode_cover_hurdle(fits, enc, family, approx = approx)
}


# ---------------------------------------------------------------------------
# Encode
# ---------------------------------------------------------------------------

#' Encode cover-hurdle data for the two-Laplace fit
#'
#' Splits `y` into a binomial occurrence indicator and a positive-cover
#' subset, builds design matrices for each arm using the fixed-effects part
#' of `formula`, and extracts any structured terms it carried (an areal
#' `icar()`/`bym2()` spatial field, plus `temporal()` / `re()` blocks). The
#' formula is parsed against the NA-dropped observations so the structured
#' index codes align with both arms.
#'
#' For `positive = "lognormal"` the positive arm's response is
#' `log(y[occur == 1])`. For `positive = "beta"` it is `y[occur == 1]` on
#' the natural (0, 1) scale; an additional eps-clip is applied so the
#' Laplace engine does not see exact 1's introduced upstream.
#'
#' @param formula State-process formula (no LHS); used for both occurrence
#'   and positive-cover arms.
#' @param data Data frame with `nrow(data) == length(y)`.
#' @param y Length-N numeric vector of cover in `[0, 1]`. NAs are dropped
#'   from both arms (treated as missing, not as zero cover).
#' @param positive One of `"lognormal"` or `"beta"`.
#' @return A list with: `occ_data`, `pos_data`, `spatial_spec` (a
#'   `tulpa_spatial` built from the unweighted areal formula term, or NULL),
#'   `trend` (per-observation weight + label from a weighted areal term, or
#'   NULL), `temporal` and `re` (structured terms from the formula, or NULL),
#'   `N`, `idx_pos` (row indices of the positive subset within `data`),
#'   `formula` (the fixed-effects formula), `positive`.
#' @keywords internal
encode_cover_hurdle <- function(formula, data, y,
                                positive = c("lognormal", "beta"),
                                autoscale = TRUE) {
  positive <- match.arg(positive)
  if (!is.numeric(y)) stop("`y` must be numeric.", call. = FALSE)
  if (length(y) != nrow(data)) {
    stop(sprintf("length(y) (%d) must equal nrow(data) (%d).",
                 length(y), nrow(data)), call. = FALSE)
  }
  rng <- range(y, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    stop("`y` must be in [0, 1] (got range [", rng[1], ", ", rng[2], "]).",
         call. = FALSE)
  }

  obs_keep <- !is.na(y)
  y_obs    <- y[obs_keep]
  data_obs <- data[obs_keep, , drop = FALSE]
  occur    <- as.integer(y_obs > 0)

  # Parse structured terms against the NA-dropped observations so re()/
  # temporal() index codes align with both hurdle arms. An areal spatial
  # term is converted to the tulpa_spatial spec the cover engine consumes.
  cover_struct <- .encode_cover_terms(formula, data_obs)
  fe_formula    <- cover_struct$fe

  X_occ_natural <- stats::model.matrix(fe_formula, data_obs)

  is_pos <- occur == 1L
  data_pos <- data_obs[is_pos, , drop = FALSE]
  y_pos    <- y_obs[is_pos]
  if (positive == "lognormal") {
    y_pos_resp <- log(y_pos)
  } else {
    # Beta arm needs y strictly in (0, 1). Cap at 1 - 1e-6; lower bound is
    # already guaranteed by occur == 1 + the range check above.
    y_pos_resp <- pmin(y_pos, 1 - 1e-6)
  }
  X_pos_natural <- stats::model.matrix(fe_formula, data_pos)

  # Autoscale numeric columns of each arm's design matrix so the optimizer
  # sees well-conditioned predictors (gcol33/tulpaObs#9). Each arm gets its
  # own scaling parameters; betas / SEs are transformed back to natural
  # scale by `decode_cover_hurdle*()`. Pass `autoscale = FALSE` to disable
  # — used by internal tests that probe `.loglik_cover_*` against a known
  # natural-scale truth, where any centering would shift the maximum.
  if (isTRUE(autoscale)) {
    occ_scaled <- .autoscale_design(X_occ_natural)
    pos_scaled <- .autoscale_design(X_pos_natural)
    X_occ      <- occ_scaled$X
    X_pos      <- pos_scaled$X
    scale_occ  <- occ_scaled$scale
    scale_pos  <- pos_scaled$scale
  } else {
    X_occ      <- X_occ_natural
    X_pos      <- X_pos_natural
    scale_occ  <- .scale_meta(X_occ_natural)
    scale_pos  <- .scale_meta(X_pos_natural)
    scale_occ$cols <- integer(0); scale_occ$means <- numeric(0); scale_occ$sds <- numeric(0)
    scale_pos$cols <- integer(0); scale_pos$means <- numeric(0); scale_pos$sds <- numeric(0)
  }

  list(
    occ_data = list(y = occur, n_trials = rep(1L, length(occur)), X = X_occ),
    pos_data = list(y = y_pos_resp, X = X_pos),
    spatial_spec = cover_struct$spatial,
    trend        = cover_struct$trend,
    mcar         = cover_struct$mcar,
    armspec      = cover_struct$armspec,
    temporal     = cover_struct$temporal,
    re           = cover_struct$re,
    N            = length(occur),
    idx_pos      = which(is_pos),
    formula      = fe_formula,
    positive     = positive,
    obs_keep     = obs_keep,
    scale_occ    = scale_occ,
    scale_pos    = scale_pos
  )
}

# Enforce the shared-latent `to` gate: a single shared latent field copied across
# arms requires BOTH cover arms on `to` (presence + positive). This still gates
# the correlated (`|` / MCAR) bar, whose semantics are copy-only -- an
# arm-specific correlated field is not defined (the cross-field Sigma lives within
# one arm; without a copy there is no cross-arm transfer to estimate). The
# INDEPENDENT (`||`) single-arm bar IS wired, as a separate per-arm latent
# (gcol33/tulpaObs#65), and is routed past this gate in `.encode_cover_terms`.
.cover_bar_check_to <- function(spec) {
  to <- spec$to %||% .tobs_cover_arms
  if (!setequal(to, .tobs_cover_arms)) {
    stop(sprintf(paste0(
      "spatial(<bar> with `|`, to = %s): a correlated (MCAR) coefficient field ",
      "is copy-only and requires both cover arms (to = c(\"presence\", ",
      "\"positive\")). For an arm-specific separate latent use the INDEPENDENT ",
      "spelling `||` with a single-arm `to` (gcol33/tulpaObs#65)."),
      paste0("c(", paste0("\"", to, "\"", collapse = ", "), ")")),
      call. = FALSE)
  }
  invisible(spec)
}

# Tag the cover-arm label a bar's `to` selects onto the joint-engine arm name:
# presence -> "occ", positive -> "pos". The cover hurdle's two response arms are
# named occ/pos internally (the binomial occurrence arm and the positive-cover
# arm); the formula `to =` labels are presence/positive (gcol33/tulpaObs#61).
.cover_arm_to_slot <- function(arm) {
  switch(arm,
         presence = "occ",
         positive = "pos",
         stop(sprintf("internal: unknown cover arm label '%s'.", arm),
              call. = FALSE))
}

# Desugar a captured INDEPENDENT (`||`) varying-coefficient spatial bar
# (gcol33/tulpaObs#61) into the intercept + per-covariate trend `tobs_spatial`
# terms the cover machinery already consumes. The expanded terms are plain
# icar/bym2/car/car_proper specs identical to the two-term form, so the bar
# desugars to exactly the existing #59 coupled path. A correlated (`|`) bar is
# routed earlier in `.encode_cover_terms` to `.cover_build_mcar_spec()`.
.cover_desugar_spatial_bar <- function(spec, data_obs) {
  .cover_bar_check_to(spec)
  .tobs_expand_spatial_bar(spec, data_obs)
}

# Build the correlated separable-MCAR field spec (gcol33/tulpaObs#64) for the
# cover hurdle from a captured correlated bar (single `|`). The bar's design
# columns (intercept + covariates) become the p coupled areal fields sharing a
# free cross-covariance Sigma (x) Q^-1; both arms (presence + positive) see the
# same fields and the whole correlated field is copied onto the positive arm
# with one estimated amplitude alpha (the cross-arm transfer). The within-arm
# covariance Sigma (the relationship AMONG the fields, e.g. does a high-baseline
# cell trend up?) is integrated over the outer CCD grid in log-Cholesky coords.
#
# `data_obs` is the NA-dropped data. The intercept field's per-observation
# weight is all-ones; a covariate field's is the design column. Returns a
# `tulpa_spatial`-shaped list with the engine-facing pieces over the OCCURRENCE
# arm (cell index + per-field weights of length nrow(data_obs)); the positive
# arm is sliced from these by `enc$idx_pos` in fit_cover_hurdle_joint_nested,
# mirroring how the `||` trend weight is subset.
.cover_build_mcar_spec <- function(spec, data_obs) {
  .cover_bar_check_to(spec)
  if (!spec$type %in% c("icar", "car")) {
    stop(sprintf(paste0(
      "spatial(<bar> with `|`): a correlated (MCAR) coefficient field uses the ",
      "intrinsic CAR (icar); model = \"%s\" is not supported. Use ",
      "model = \"icar\" (the default) or the independent spelling `||`."),
      spec$type), call. = FALSE)
  }
  specs <- tulpa::tulpa_bar_field_specs(spec$bar_formula, data_obs)
  node  <- attr(specs, "node")
  .tobs_validate_bar_node(node, spec$graph, data_obs)
  # Replicate over the `by` levels (gcol33/tulpaObs#82): the separable-MCAR field
  # is built over I_L (x) Q -- one disjoint correlated (intercept, slope) field
  # per level -- with the cross-field Sigma (x) Q^-1 shared across levels (no
  # `by` is the identity). The copy onto the positive arm carries the whole
  # replicated field at the one estimated amplitude alpha, unchanged.
  rg      <- .tobs_bar_resolve_graph(spec, data_obs, node)
  graph   <- rg$graph
  n_nodes <- nrow(graph)
  idx_occ <- rg$idx
  csr <- adjacency_to_csr(graph)

  # tulpa_bar_field_specs() returns weight = NULL for the intercept column (it is
  # the unweighted all-ones field); the MCAR factory needs an explicit per-row
  # design column for every field, so materialize ones for the intercept.
  n_obs <- nrow(data_obs)
  field_weight_occ <- lapply(specs, function(col) {
    if (isTRUE(col$is_intercept) || is.null(col$weight)) rep(1.0, n_obs)
    else as.numeric(col$weight)
  })
  field_names <- vapply(specs, function(col)
    paste(node, col$column_name, sep = "."), character(1))

  list(
    type            = "mcar",
    graph           = graph,
    n_spatial_units = as.integer(n_nodes),
    n_fields        = length(specs),
    adj_row_ptr     = as.integer(csr$row_ptr),
    adj_col_idx     = as.integer(csr$col_idx),
    n_neighbors     = as.integer(csr$n_neighbors),
    idx_occ         = as.integer(idx_occ),
    field_weight_occ = field_weight_occ,
    field_names     = field_names,
    by              = rg$by
  )
}

# Parse a cover() formula against the NA-dropped observations: return the
# fixed-effects formula plus the structured terms it carried, split by kind.
# The areal terms (icar/bym2/car/car_proper) split by their `weight`:
#   * an unweighted areal term is the shared intercept field, converted to the
#     tulpa_spatial spec the engine consumes (`spatial`);
#   * a weighted areal term (`icar(graph = adj, weight = col, group_var = ...)`)
#     is the spatially-varying TREND field -- the formula-DSL spelling of the
#     coupled second besag block (gcol33/tulpaObs#59). Its per-observation weight
#     `col` and label come back in `trend`.
# A varying-coefficient bar (`spatial(~ 1 + w || node, graph = adj, to = ...)`,
# gcol33/tulpaObs#61) is the compact single-term spelling: it desugars in place
# to the intercept field (its `1` column) plus a weight-scaled trend field per
# covariate column, all on the bar's node index, so the two forms feed the same
# machinery. The shared `to = c("presence", "positive")` path is the only one
# wired here; `|` and arm-specific `to` are gated below.
# svc()/latent() are not meaningful for the cover hurdle and are rejected.
.encode_cover_terms <- function(formula, data_obs) {
  bind <- .tobs_bind_formulas(list(state = formula), data_obs)
  spatial_specs <- list(); temporal <- NULL; re <- list(); mcar <- NULL
  armspec <- list()
  for (t in bind$terms) {
    spec <- t$spec
    if (inherits(spec, "tobs_spatial") && isTRUE(spec$is_bar) &&
        isTRUE(spec$correlated)) {
      # Correlated bar (single `|`): one separable-MCAR block over the bar's
      # design columns sharing a free Sigma, copied onto the positive arm
      # (gcol33/tulpaObs#64). Distinct from the `||` (independent) desugaring.
      if (!is.null(mcar)) {
        stop("cover(): only one correlated spatial bar (single `|`) is ",
             "supported.", call. = FALSE)
      }
      mcar <- .cover_build_mcar_spec(spec, data_obs)
    } else if (inherits(spec, "tobs_spatial") && isTRUE(spec$is_bar) &&
               length(spec$to %||% .tobs_cover_arms) == 1L) {
      # Arm-specific INDEPENDENT bar (single-arm `to`, gcol33/tulpaObs#65): a
      # separate per-arm latent field on ONLY that arm, with its own precision
      # and NO cross-arm copy. Distinct from the shared `||` (both-arm) desugar,
      # which copies a presence-anchored field onto the positive arm. Each
      # single-arm bar is collected as a self-describing field block; the fitter
      # places each on its arm via a 0-sentinel spatial_idx on the other arm.
      armspec[[length(armspec) + 1L]] <- .tobs_armspecific_bar_fields(spec, data_obs)
    } else if (inherits(spec, "tobs_spatial") && isTRUE(spec$is_bar)) {
      expanded <- .cover_desugar_spatial_bar(spec, data_obs)
      spatial_specs <- c(spatial_specs, expanded)
    } else if (inherits(spec, "tobs_spatial")) {
      spatial_specs[[length(spatial_specs) + 1L]] <- spec
    } else if (inherits(spec, "tobs_temporal")) {
      if (!is.null(temporal)) {
        stop("cover(): only one temporal term is supported.", call. = FALSE)
      }
      temporal <- spec
    } else if (inherits(spec, "tobs_re")) {
      re[[length(re) + 1L]] <- spec
    } else {
      stop(sprintf("cover() does not support `%s` terms in the formula.",
                   spec$label %||% class(spec)[1]), call. = FALSE)
    }
  }
  # Soft guard (gcol33/tulpaObs#62): a bare `| / ||` RE bar whose grouping factor
  # is also an areal term's graph-node group_var is the engine-bar-idiom papercut
  # -- the bar is fitted as an IID random effect, not a spatial field. RE bars are
  # legitimate, so this informs (message) rather than rejecting; it is silent when
  # the bar's factor is unrelated to any spatial term.
  .tobs_cover_bar_re_guard(formula, spatial_specs)
  if (!is.null(mcar) && length(spatial_specs) > 0L) {
    stop("cover(): a correlated spatial bar (single `|`) is the whole spatial ",
         "structure; it cannot be combined with other areal terms in the same ",
         "formula. Put the intercept and slope fields in the one bar ",
         "(~ 1 + w | node).", call. = FALSE)
  }
  if (length(armspec) > 0L) {
    # Arm-specific separate latents (gcol33/tulpaObs#65) are an independent
    # spatial structure: each is its own per-arm block with its own precision and
    # no cross-arm copy. They do not compose with the shared/copied intercept +
    # trend machinery, the correlated MCAR copy, or temporal()/re() blocks in the
    # same fit (those are coupled structures; mixing would silently re-introduce a
    # cross-arm transfer the user opted out of). Reject the combination with a
    # pointer rather than ignoring one half.
    if (length(spatial_specs) > 0L || !is.null(mcar)) {
      stop("cover(): arm-specific spatial fields (single-arm `to`) are a ",
           "separate per-arm structure with no cross-arm copy; they cannot be ",
           "combined with a shared field (both-arm bar or icar()/bym2()), a ",
           "correlated `|` bar, or a weighted trend term in the same formula. ",
           "Use only single-arm spatial() bars, or only the shared form.",
           call. = FALSE)
    }
    # At most one field per arm in the first ship: two presence-only (or two
    # positive-only) bars would be two independent fields on the same arm, which
    # the joint driver carries but the cover-arm field reconstruction does not yet
    # disambiguate. Each arm takes one separate latent.
    arms_used <- vapply(armspec, function(a) a$arm, character(1))
    if (anyDuplicated(arms_used)) {
      stop("cover(): each arm-specific spatial field must target a distinct arm ",
           "(got two on the same arm). Combine multiple coefficient fields into ",
           "one bar, e.g. spatial(~ 1 + w || cell, graph = adj, to = \"",
           arms_used[anyDuplicated(arms_used)], "\").", call. = FALSE)
    }
  }
  fields <- .cover_resolve_spatial_fields(spatial_specs, data_obs)
  list(fe = bind$fe$state, spatial = fields$spatial, trend = fields$trend,
       temporal = temporal, re = if (length(re)) re else NULL, mcar = mcar,
       armspec = if (length(armspec)) armspec else NULL)
}

# Partition the cover() formula's areal terms into the shared intercept field
# and the optional spatially-varying trend field (gcol33/tulpaObs#59). An
# unweighted areal term is the intercept; a weighted areal term
# (`icar(..., weight = col)`) is the trend -- the second coupled besag block
# that `control$trend` used to introduce. Both spellings of the weighted term --
# bare `icar(..., weight = )` and the umbrella `spatial(model = "icar",
# weight = )` -- resolve to the same `tobs_spatial` term and so to the same
# trend block.
#
# Returns list(spatial, trend):
#   * spatial: the tulpa_spatial spec for the intercept field, or NULL.
#   * trend:   list(w_occ [N], label) of the per-observation trend weight over
#              `data_obs`, or NULL when no weighted term is present.
.cover_resolve_spatial_fields <- function(specs, data_obs) {
  if (length(specs) == 0L) return(list(spatial = NULL, trend = NULL))

  weighted   <- vapply(specs, function(s) !is.null(s$weight), logical(1))
  base_specs <- specs[!weighted]
  wt_specs   <- specs[weighted]

  if (length(base_specs) == 0L) {
    stop("cover() spatial requires one unweighted intercept field ",
         "(e.g. icar(graph = adj)); only weighted trend field(s) were given. ",
         "Add the bare areal term first.", call. = FALSE)
  }
  if (length(base_specs) > 1L) {
    stop("cover() supports exactly one unweighted intercept field; got ",
         length(base_specs), ". A spatially-varying trend must be a weighted ",
         "areal term, e.g. icar(graph = adj, weight = time.sc).", call. = FALSE)
  }
  if (length(wt_specs) > 1L) {
    stop("cover() supports a single weighted trend field; got ",
         length(wt_specs), ".", call. = FALSE)
  }

  base_spec <- base_specs[[1L]]
  spatial   <- .tobs_term_to_tulpa_spatial(base_spec)

  trend <- NULL
  if (length(wt_specs) == 1L) {
    ws <- wt_specs[[1L]]
    # The trend field shares the intercept field's areal graph and group_var;
    # the cell-coupling engine carries one node per cell, so the two blocks
    # differ only in the per-observation weight.
    if (!identical(dim(ws$graph), dim(base_spec$graph)) ||
        !all(ws$graph == base_spec$graph)) {
      stop("cover() trend field must share the same areal graph as the ",
           "intercept field (same nodes / adjacency).", call. = FALSE)
    }
    if (!identical(ws$group_var, base_spec$group_var)) {
      stop("cover() trend field must share the intercept field's group_var ",
           "(or both name none).", call. = FALSE)
    }
    w_occ <- as.numeric(ws$weight)
    if (length(w_occ) != nrow(data_obs)) {
      stop(sprintf(paste0(
        "cover() trend weight has length %d but the data has %d ",
        "observations; supply it as a per-observation covariate."),
        length(w_occ), nrow(data_obs)), call. = FALSE)
    }
    if (anyNA(w_occ) || !all(is.finite(w_occ))) {
      stop("cover() trend weight must be a finite numeric covariate.",
           call. = FALSE)
    }
    trend <- list(w_occ = w_occ, label = ws$weight_label %||% "trend")
  }

  list(spatial = spatial, trend = trend)
}


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

#' Fit the two arms of a cover hurdle
#'
#' Two independent `tulpa::tulpa_laplace()` calls. The joint shared-field
#' fit is Phase 1c. For `positive = "lognormal"` the positive arm is a
#' Gaussian fit on `log(cover)` with sigma estimated post-hoc as the
#' residual standard error. For `positive = "beta"` the positive arm uses
#' `tulpa::tulpa_laplace_beta()` which estimates the precision `phi` via
#' an outer 1-D optimisation and weights the Hessian accordingly.
#'
#' @param enc Output of [encode_cover_hurdle()].
#' @param positive `"lognormal"` or `"beta"` (taken from `enc$positive`).
#' @param engine `"laplace"` (default) or `"nested_laplace"`. The latter is
#'   routed through [fit_cover_hurdle_joint_nested()].
#' @param priors Optional [cover_priors()] object (or a coercible list /
#'   `FALSE`). Adds a weakly-informative fixed-effect penalty on both arms
#'   (occurrence + positive, beta or lognormal); `NULL` / `FALSE` fit
#'   unpenalised. Rejected with a spatial formula (the spatial solver carries
#'   its own prior).
#' @param control List with optional `max.iter`, `tol`, `n.threads`.
#' @return List with `m_occ`, `m_pos`, `positive`, `pos_fit_n`, `pos_fit_p`,
#'   plus one of `sigma_pos` (lognormal) or `phi_pos` (beta).
#' @keywords internal
fit_cover_hurdle <- function(enc, positive = enc$positive,
                             engine = "laplace",
                             priors = NULL, control = list()) {
  if (!engine %in% c("laplace", "auto")) {
    stop("cover() currently supports only method = 'laplace'/'laplace_sla' ",
         "or 'nested_laplace'/'nested_laplace_sla' (got engine '", engine,
         "'). nuts lands in later phases.", call. = FALSE)
  }
  max_iter  <- control$max.iter  %||% 100L
  tol       <- control$tol       %||% 1e-6
  n_threads <- control$n.threads %||% 1L

  # Opt-in fixed-effect priors (cover_priors()): the same quadratic beta_prior
  # tulpa_laplace() applies on the occupancy path, specified on natural-scale
  # coefficients and applied on the autoscaled design (occupancy convention).
  # NULL = unpenalised. Both arms are penalisable -- the occurrence and
  # lognormal arms through tulpa_laplace(), the beta arm through
  # tulpa_laplace_beta()'s beta_prior. Spatial cover formulas still reject the
  # prior (the spatial solver carries its own).
  cprior <- .resolve_cover_priors(priors)
  if (!is.null(cprior) && !is.null(enc$spatial_spec)) {
    stop("cover priors are not supported with a spatial term in the formula ",
         "(the spatial Laplace solver carries its own fixed-effect prior). ",
         "Drop the prior or the spatial term.", call. = FALSE)
  }
  occ_bp <- .cover_arm_prior(cprior, "occ", colnames(enc$occ_data$X))
  pos_bp <- .cover_arm_prior(cprior, "pos", colnames(enc$pos_data$X))

  m_occ <- tulpa::tulpa_laplace(
    y        = enc$occ_data$y,
    n_trials = enc$occ_data$n_trials,
    X        = enc$occ_data$X,
    family   = "binomial",
    spatial  = enc$spatial_spec,
    max_iter = max_iter, tol = tol, n_threads = n_threads,
    beta_prior = occ_bp
  )

  if (length(enc$pos_data$y) < ncol(enc$pos_data$X) + 1L) {
    stop("Too few positive-cover sites (", length(enc$pos_data$y),
         ") for the requested formula (", ncol(enc$pos_data$X),
         " coefficients). Need at least ncol(X) + 1.", call. = FALSE)
  }

  n_pos <- length(enc$pos_data$y)
  p_pos <- ncol(enc$pos_data$X)

  if (positive == "lognormal") {
    m_pos <- tulpa::tulpa_laplace(
      y        = enc$pos_data$y,
      n_trials = rep(1L, n_pos),
      X        = enc$pos_data$X,
      family   = "gaussian",
      spatial  = enc$spatial_spec,
      max_iter = max_iter, tol = tol, n_threads = n_threads,
      beta_prior = pos_bp
    )
    # Gaussian Laplace runs with phi = 1; estimate residual SD post-hoc.
    beta_pos <- m_pos$mode[seq_len(p_pos)]
    eta_pos  <- as.numeric(enc$pos_data$X %*% beta_pos)
    resid    <- enc$pos_data$y - eta_pos
    sigma_pos <- sqrt(sum(resid^2) / max(n_pos - p_pos, 1L))
    return(list(
      m_occ     = m_occ,
      m_pos     = m_pos,
      positive  = "lognormal",
      sigma_pos = sigma_pos,
      pos_fit_n = n_pos,
      pos_fit_p = p_pos
    ))
  }

  # positive == "beta": the beta arm is penalised via the beta solver's own
  # beta_prior (gcol33/tulpa, tulpa_laplace_beta gained beta_prior); the
  # occurrence arm is penalised by occ_bp above.
  m_pos <- tulpa::tulpa_laplace_beta(
    y         = enc$pos_data$y,
    X         = enc$pos_data$X,
    spatial   = enc$spatial_spec,
    max_iter  = max_iter, tol = tol, n_threads = n_threads,
    beta_prior = pos_bp
  )
  list(
    m_occ     = m_occ,
    m_pos     = m_pos,
    positive  = "beta",
    phi_pos   = m_pos$phi,
    pos_fit_n = n_pos,
    pos_fit_p = p_pos
  )
}


# ---------------------------------------------------------------------------
# Decode
# ---------------------------------------------------------------------------

#' Decode the two-arm fit into a cover_fit object
#'
#' Extracts beta vectors and SEs for each arm. SEs are scaled to match
#' each arm's dispersion convention:
#'
#' * lognormal arm: `tulpa_laplace(family = "gaussian")` computes the
#'   Hessian assuming phi = 1, so SEs are rescaled by `sigma_pos^2`.
#' * beta arm: `tulpa_laplace_beta()` already weights the Hessian by phi
#'   (Fisher information), so SEs are returned at scale 1.
#'
#' Under an SLA method (`method = "laplace_sla"` / `"nested_laplace_sla"`), the
#' cover-hurdle SLA gamma is
#' computed via [`.sla_compute_cover_hurdle()`]: a per-arm 5-point FD of
#' the *original* Bernoulli / Beta / Lognormal log-likelihood against the
#' arm's `solve(H_beta)` Sigma (raw Hessian — no Louis correction needed
#' here because both arms run real likelihoods at the mode, not the
#' pseudo-binomial M-step encoding). Per-arm pseudo-draws are then
#' resampled from skew-normals fit by moment-matching `(beta_arm,
#' se_arm, gamma_arm)`.
#'
#' @keywords internal
decode_cover_hurdle <- function(fits, enc, family,
                                approx = "gaussian_laplace") {
  p_occ_n <- ncol(enc$occ_data$X)
  p_pos_n <- ncol(enc$pos_data$X)

  # Modes / SEs come back from the engine in the *scaled* design's
  # parameterization. Transform back to the user-facing natural scale via
  # the (mean, sd) cache stashed by `encode_cover_hurdle()`. The covariance
  # transform `V_nat = T %*% V_sc %*% t(T)` is the right object because the
  # Hessian-based vcov has informative off-diagonals between intercept and
  # slope; treating it as diagonal here would inflate the intercept SE.
  beta_occ_sc <- fits$m_occ$mode[seq_len(p_occ_n)]
  beta_pos_sc <- fits$m_pos$mode[seq_len(p_pos_n)]

  V_occ_sc <- if (!is.null(fits$m_occ$H_beta)) {
    tryCatch(solve(fits$m_occ$H_beta), error = function(e) NULL)
  } else NULL
  pos_vcov_scale <- if (fits$positive == "lognormal") fits$sigma_pos^2 else 1
  V_pos_sc <- if (!is.null(fits$m_pos$H_beta)) {
    tryCatch(pos_vcov_scale * solve(fits$m_pos$H_beta), error = function(e) NULL)
  } else NULL

  beta_occ <- .unscale_beta_vec(beta_occ_sc, enc$scale_occ)
  beta_pos <- .unscale_beta_vec(beta_pos_sc, enc$scale_pos)
  names(beta_occ) <- colnames(enc$occ_data$X)
  names(beta_pos) <- colnames(enc$pos_data$X)

  V_occ <- .unscale_vcov_block(V_occ_sc, enc$scale_occ)
  V_pos <- .unscale_vcov_block(V_pos_sc, enc$scale_pos)
  se_occ <- if (is.null(V_occ)) rep(NA_real_, p_occ_n) else
    sqrt(pmax(diag(as.matrix(V_occ)), 0))
  se_pos <- if (is.null(V_pos)) rep(NA_real_, p_pos_n) else
    sqrt(pmax(diag(as.matrix(V_pos)), 0))
  if (length(se_occ)) names(se_occ) <- names(beta_occ)
  if (length(se_pos)) names(se_pos) <- names(beta_pos)

  hyperpar <- list(
    occ = .extract_spatial_hyperpar(fits$m_occ, enc$spatial_spec),
    pos = .extract_spatial_hyperpar(fits$m_pos, enc$spatial_spec)
  )
  if (fits$positive == "lognormal") {
    hyperpar$sigma_pos <- fits$sigma_pos
  } else {
    hyperpar$phi_pos <- fits$phi_pos
  }

  # Simplified-Laplace gamma + skew-normal pseudo-draws per arm.
  skew_occ <- NULL
  skew_pos <- NULL
  draws_occ <- NULL
  draws_pos <- NULL
  sla_status <- "off"
  if (identical(approx, "simplified_laplace")) {
    sla_res <- .sla_compute_cover_hurdle(fits, enc, fits$positive)
    sla_draws <- .sla_build_cover_hurdle_draws(
      beta_occ, se_occ, beta_pos, se_pos, sla_res,
      V_occ = V_occ, V_pos = V_pos
    )
    draws_occ <- sla_draws$draws_occ
    draws_pos <- sla_draws$draws_pos
    sla_status <- sla_draws$sla_status
    if (isTRUE(sla_res$valid)) {
      skew_occ <- sla_res$gamma_occ
      skew_pos <- sla_res$gamma_pos
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
      sigma_pos    = if (fits$positive == "lognormal") fits$sigma_pos else NA_real_,
      sigma_pos_sd = NA_real_,
      phi_pos      = if (fits$positive == "beta")      fits$phi_pos    else NA_real_,
      phi_pos_sd   = NA_real_,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
      converged    = isTRUE(fits$m_occ$converged) && isTRUE(fits$m_pos$converged),
      # Unified convergence record, the same list every other family stores
      # (gcol33/tulpaObs#88), so a mixed-family QC pass reads one accessor
      # (`convergence(fit)` / `fit$convergence$converged`) across occu /
      # occu_cover / cover. The top-level `converged` is kept for glance() and
      # back-compat; `sla_status` carries the simplified-Laplace marginal code.
      convergence  = list(
        converged  = isTRUE(fits$m_occ$converged) && isTRUE(fits$m_pos$converged),
        n_iter     = fits$m_occ$n_iter %||% NA_integer_,
        sla_status = sla_status),
      log_marginal = c(occ = fits$m_occ$log_marginal,
                       pos = fits$m_pos$log_marginal),
      skew_occ     = skew_occ,
      skew_pos     = skew_pos,
      draws_occ    = draws_occ,
      draws_pos    = draws_pos,
      sla_status   = sla_status
    ),
    class = c("cover_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


# ---------------------------------------------------------------------------
# Pointwise log-likelihood (for WAIC / PSIS-LOO / tobs_stack)
# ---------------------------------------------------------------------------

# Pointwise log-likelihood [n_draws x N] for a cover hurdle fit. The hurdle is
# fully observed (nothing to marginalize): per site, an occurrence Bernoulli
# plus, when y > 0, the positive-part density. Per-arm betas are sampled from
# the Gaussian-Laplace posterior in the optimizer's (scaled) parameterization
# and paired with the stored scaled designs, so eta matches the fitted model
# without needing the raw data frame. Dispersion (sigma_pos / phi_pos) is held
# at its fitted value; SLA skew is not applied to the LOO marginals.
#
# Only the separate-Laplace path (method = "laplace" / "laplace_sla") carries
# the per-arm mode + Hessian this needs; the nested-joint path errors.
# Shared hurdle log-likelihood kernel: given per-draw linear predictors
# `eta_occ` [S x N] (occurrence) and `eta_pos` [S x N_pos] (positive part),
# dispersion `disp` (scalar or length-S), and the response encoding, accumulate
# the [S x N] pointwise hurdle log-likelihood -- log(1 - p) at absent sites,
# log(p) + positive-part log-density at occupied sites. The single source of
# truth for both the separate-Laplace and nested-joint cover paths and for the
# posterior-mean plug-in (S = 1).
.tobs_cover_hurdle_ll <- function(eta_occ, eta_pos, disp, occur, y_pos, idx_pos,
                                  positive) {
  S <- nrow(eta_occ); N <- ncol(eta_occ)
  sd_disp <- if (length(disp) == 1L) rep(disp, S) else disp
  log_p   <- .tobs_log_p(eta_occ)
  log_1mp <- .tobs_log_1mp(eta_occ)

  ll <- matrix(0, S, N)
  absent <- occur == 0L
  if (any(absent)) ll[, absent] <- log_1mp[, absent, drop = FALSE]

  pos_col <- match(seq_len(N), idx_pos)   # eta_pos column per site (NA if absent)
  for (i in which(occur == 1L)) {
    j <- pos_col[i]
    if (positive == "lognormal") {
      # y_pos = log(y); density of natural-scale y is the Gaussian on log(y)
      # times the Jacobian 1/y, i.e. dnorm(log y, eta, sigma, log) - log(y).
      dens <- stats::dnorm(y_pos[j], mean = eta_pos[, j], sd = sd_disp,
                           log = TRUE) - y_pos[j]
    } else {
      mu   <- stats::plogis(eta_pos[, j])
      dens <- stats::dbeta(y_pos[j], mu * sd_disp, (1 - mu) * sd_disp,
                           log = TRUE)
    }
    ll[, i] <- log_p[, i] + dens
  }
  ll
}

# Per-draw cover linear predictors [S x N] / [S x N_pos]. Separate-Laplace path:
# sample each arm's Gaussian-Laplace posterior at its scaled design. Nested-joint
# path: sample the grid-integrated joint and project the shared field at each
# observation's spatial unit (`spi_full` / `spi_pos`); the dispersion is then the
# per-draw grid value. Returns list(eta_occ, eta_pos, disp).
.tobs_cover_eta_draws <- function(object, n.draws = 1000L) {
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  if (!is.null(.tobs_joint_fit(object))) {
    bundle  <- .tobs_joint_draws(object, n = n.draws)
    if (isTRUE(object$armspecific)) {
      # Arm-specific separate latents store no node map / weight on the bundle
      # block (gcol33/tulpaObs#95); the pointwise-loglik consumer runs over the
      # fit's observations, so it rebuilds the per-arm per-observation node map
      # and covariate-weight lookup from `armspec_blocks` and hands them to
      # .tobs_joint_arm_eta exactly as the shared-field path passes spi_* / wfun.
      u_occ  <- .tobs_armspec_obs_units(object, 1L, nrow(enc$occ_data$X))
      u_pos  <- .tobs_armspec_obs_units(object, 2L, nrow(enc$pos_data$X))
      wf_occ <- .tobs_armspec_obs_wfun(object, 1L)
      wf_pos <- .tobs_armspec_obs_wfun(object, 2L)
      eta_occ <- t(.tobs_joint_arm_eta(bundle, enc$occ_data$X, "occ", u_occ, wf_occ))
      eta_pos <- t(.tobs_joint_arm_eta(bundle, enc$pos_data$X, "pos", u_pos, wf_pos))
      return(list(eta_occ = eta_occ, eta_pos = eta_pos, disp = bundle$disp))
    }
    spi_full <- object$spi_full
    spi_pos  <- object$spi_pos
    if (is.null(spi_full) || is.null(spi_pos)) {
      stop("Pointwise log-likelihood for the nested-joint cover() fit needs the ",
           "per-observation spatial-unit index (`spi_full` / `spi_pos`); refit ",
           "with the current tulpaObs so they are stored on the fit.",
           call. = FALSE)
    }
    # A coupled trend field carries a per-observation weight; resolve it per arm
    # via wfun (the engine's svc_weight replayed at predict time). NULL when the
    # fit has no trend field, in which case .tobs_joint_arm_eta never calls it.
    w_occ_fun <- if (!is.null(object$trend_w_occ))
                   function(col) as.numeric(object$trend_w_occ) else NULL
    w_pos_fun <- if (!is.null(object$trend_w_pos))
                   function(col) as.numeric(object$trend_w_pos) else NULL
    eta_occ <- t(.tobs_joint_arm_eta(bundle, enc$occ_data$X, "occ", spi_full,
                                     wfun = w_occ_fun))
    eta_pos <- t(.tobs_joint_arm_eta(bundle, enc$pos_data$X, "pos", spi_pos,
                                     wfun = w_pos_fun))
    return(list(eta_occ = eta_occ, eta_pos = eta_pos, disp = bundle$disp))
  }
  if (!is.null(object$nuts) && is.matrix(object$draws)) {
    # NUTS path: project the exact posterior coefficient draws through the
    # natural-scale presence / positive designs. The dispersion is the per-draw
    # exp(log_disp) trailing column, so the score marginalizes the calibrated
    # NUTS posterior (the point of the sampler) rather than the Laplace mode.
    draws  <- object$draws
    if (!is.null(n.draws) && n.draws < nrow(draws)) {
      draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
    }
    p_pres <- ncol(enc$occ_data$X); p_pos <- ncol(enc$pos_data$X)
    B_occ <- draws[, seq_len(p_pres), drop = FALSE]
    B_pos <- draws[, p_pres + seq_len(p_pos), drop = FALSE]
    disp  <- exp(draws[, ncol(draws)])
    return(list(eta_occ = B_occ %*% t(enc$occ_data$X),
                eta_pos = B_pos %*% t(enc$pos_data$X), disp = disp))
  }
  if (is.null(enc) || is.null(object$occ$mode) || is.null(object$occ$H_beta) ||
      is.null(object$pos$mode) || is.null(object$pos$H_beta)) {
    stop("Pointwise log-likelihood for cover() is implemented for the ",
         "separate-Laplace path (method = 'laplace' / 'laplace_sla') and the ",
         "nested-joint shared-field path (method = 'nested_laplace'); this fit ",
         "carries neither a per-arm mode + Hessian nor a joint object.",
         call. = FALSE)
  }
  X_occ <- enc$occ_data$X
  X_pos <- enc$pos_data$X
  p_occ <- ncol(X_occ); p_pos <- ncol(X_pos)
  mode_occ   <- object$occ$mode[seq_len(p_occ)]
  mode_pos   <- object$pos$mode[seq_len(p_pos)]
  pos_vscale <- if (positive == "lognormal") (object$sigma_pos %||% 1)^2 else 1
  V_occ <- tryCatch(solve(object$occ$H_beta), error = function(e) NULL)
  V_pos <- tryCatch(pos_vscale * solve(object$pos$H_beta), error = function(e) NULL)
  S <- as.integer(n.draws)
  B_occ <- .tobs_mvn_draws(mode_occ, V_occ, S)   # [S x p_occ]
  B_pos <- .tobs_mvn_draws(mode_pos, V_pos, S)   # [S x p_pos]
  disp  <- if (positive == "lognormal") object$sigma_pos else object$phi_pos
  list(eta_occ = B_occ %*% t(X_occ), eta_pos = B_pos %*% t(X_pos), disp = disp)
}

# Pointwise log-likelihood [n_draws x N] for a cover hurdle fit (separate-Laplace
# or nested-joint shared-field), assembled from the per-draw linear predictors
# through the shared hurdle kernel.
.tobs_ploglik_cover <- function(object, n.draws = 1000L) {
  enc <- object$encoding
  e   <- .tobs_cover_eta_draws(object, n.draws)
  .tobs_cover_hurdle_ll(e$eta_occ, e$eta_pos, e$disp, enc$occ_data$y,
                        enc$pos_data$y, enc$idx_pos,
                        object$positive %||% "lognormal")
}

# Pointwise log-likelihood at the posterior mean of the parameters (length N):
# the hurdle kernel evaluated at the mean linear predictors and mean dispersion.
.tobs_cover_loglik_at_mean <- function(object, n.draws = 1000L) {
  enc <- object$encoding
  e   <- .tobs_cover_eta_draws(object, n.draws)
  mean_eta_occ <- matrix(colMeans(e$eta_occ), nrow = 1L)
  mean_eta_pos <- matrix(colMeans(e$eta_pos), nrow = 1L)
  as.numeric(.tobs_cover_hurdle_ll(
    mean_eta_occ, mean_eta_pos, mean(e$disp), enc$occ_data$y,
    enc$pos_data$y, enc$idx_pos, object$positive %||% "lognormal"
  ))
}


# ---------------------------------------------------------------------------
# Posterior predictive check + PIT for the cover hurdle (gcol33/tulpaObs#27)
# ---------------------------------------------------------------------------

# Randomized PIT for a cover() hurdle fit (length N). The predictive CDF mixes a
# point mass 1 - p at the structural zero with p * F_pos on the positive part
# (lognormal / beta CDF at the fitted per-draw predictor), projecting the shared
# field per observation for the nested-joint fit. Absent sites use the left /
# right limits [0, 1 - p] around the zero mass; occupied sites are continuous
# (F = 1 - p + p F_pos), so the engine's randomized PIT is degenerate there.
.tobs_pit_cover <- function(object, n.samples = 250) {
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  e <- .tobs_cover_eta_draws(object, n.draws = n.samples)
  eta_occ <- e$eta_occ; eta_pos <- e$eta_pos
  S <- nrow(eta_occ); N <- ncol(eta_occ)
  occur <- enc$occ_data$y; y_pos <- enc$pos_data$y; idx_pos <- enc$idx_pos
  sd_disp <- if (length(e$disp) == 1L) rep(e$disp, S) else e$disp

  p      <- stats::plogis(eta_occ)            # [S x N]
  one_mp <- 1 - p
  Fl <- matrix(0, S, N)
  Fu <- one_mp                                # absent: F in [0, 1 - p]

  pos_col <- match(seq_len(N), idx_pos)
  for (i in which(occur == 1L)) {
    j <- pos_col[i]
    if (positive == "lognormal") {
      Fpos <- stats::pnorm((y_pos[j] - eta_pos[, j]) / sd_disp)
    } else {
      mu   <- stats::plogis(eta_pos[, j])
      Fpos <- stats::pbeta(y_pos[j], mu * sd_disp, (1 - mu) * sd_disp)
    }
    val <- one_mp[, i] + p[, i] * Fpos        # continuous -> no randomization
    Fu[, i] <- val
    Fl[, i] <- val
  }
  tulpa::tulpa_pit(Fu, cdf_lower = Fl)
}

# Posterior predictive check for a cover() hurdle fit. Per draw, occurrence
# replicates from Bernoulli(p) and cover replicates from the fitted positive
# part; the discrepancy is a Freeman-Tukey (or chi-squared) sum on the
# occurrence arm (all N sites, against p) plus the positive arm (occupied
# subset, against the positive-part mean), returning a Bayesian p-value.
.tobs_ppc_cover <- function(object,
                            fit.stat = c("freeman-tukey", "chi-squared"),
                            n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  e <- .tobs_cover_eta_draws(object, n.draws = n.samples)
  eta_occ <- e$eta_occ; eta_pos <- e$eta_pos
  S <- nrow(eta_occ); N <- ncol(eta_occ)
  occur <- enc$occ_data$y; y_pos <- enc$pos_data$y
  sd_disp <- if (length(e$disp) == 1L) rep(e$disp, S) else e$disp
  # Observed cover on the natural scale at occupied sites (positive subset).
  y_pos_nat <- if (positive == "lognormal") exp(y_pos) else y_pos
  n_pos <- length(y_pos_nat)

  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(o, ex) sum((sqrt(o) - sqrt(ex))^2, na.rm = TRUE)
  } else {
    function(o, ex) sum((o - ex)^2 / (ex + 1e-10), na.rm = TRUE)
  }
  p_all <- stats::plogis(eta_occ)

  fit_y <- fit_rep <- numeric(S)
  for (s in seq_len(S)) {
    p_s     <- p_all[s, ]
    occ_rep <- stats::rbinom(N, 1, p_s)
    occ_obs <- stat_fn(occur, p_s)
    occ_rp  <- stat_fn(occ_rep, p_s)

    pos_obs <- pos_rp <- 0
    if (n_pos > 0L) {
      if (positive == "lognormal") {
        mu_log <- eta_pos[s, ]; sg <- sd_disp[s]
        Epos     <- exp(mu_log + sg^2 / 2)
        ypos_rep <- exp(stats::rnorm(n_pos, mu_log, sg))
      } else {
        mu   <- stats::plogis(eta_pos[s, ]); phi <- sd_disp[s]
        Epos <- mu
        ypos_rep <- stats::rbeta(n_pos, mu * phi, (1 - mu) * phi)
      }
      pos_obs <- stat_fn(y_pos_nat, Epos)
      pos_rp  <- stat_fn(ypos_rep, Epos)
    }
    fit_y[s]   <- occ_obs + pos_obs
    fit_rep[s] <- occ_rp + pos_rp
  }
  list(fit.y = fit_y, fit.y.rep = fit_rep,
       bayesian.p = mean(fit_rep > fit_y))
}
# mode (point mass) if V is unavailable, and jitters a near-singular V.
.tobs_mvn_draws <- function(mu, V, S) {
  p <- length(mu)
  if (is.null(V)) return(matrix(mu, S, p, byrow = TRUE))
  R <- tryCatch(chol(V), error = function(e) {
    jit <- 1e-8 * (mean(diag(V)) + 1e-8)
    chol(V + diag(jit, p))
  })
  Z <- matrix(stats::rnorm(S * p), S, p)
  matrix(mu, S, p, byrow = TRUE) + Z %*% R
}


# ---------------------------------------------------------------------------
# Predict
# ---------------------------------------------------------------------------

#' Predict cover from a cover_fit
#'
#' Occurrence probability is always `p = plogis(X * beta_occ)`. The
#' conditional positive cover `mu` depends on the positive-part family:
#'
#' * `positive = "lognormal"`: `mu = exp(eta_pos + sigma_pos^2 / 2)`
#'   (lognormal back-transform on log-cover).
#' * `positive = "beta"`: `mu = plogis(eta_pos)` (mean of the beta on
#'   the natural cover scale with logit link).
#'
#' Expected cover is `E[y] = p * mu` under both positive parts.
#'
#' The separate-Laplace fit (`method = "laplace"`) returns a fixed-effects-only
#' numeric vector. The nested-Laplace shared-field fit
#' (`method = "nested_laplace"`) instead projects the shared occupancy-cover
#' field and returns a `tobs_prediction` of posterior draws -- the same tidy /
#' `change` contract as [predict.tobs_fit()] for `occu_cover()`: pass
#' `type = "change"` with `times = c(t1, t2)` and `time_col` for a per-cell
#' delta map. Each prediction unit is a row of `newdata` (or a `cell` column
#' indexing the field cells), and every quantity is marginalized per draw over
#' the grid-integrated joint posterior.
#'
#' @param object A `cover_fit`.
#' @param newdata A data frame of covariates matching the original formula. For
#'   the nested-Laplace fit, one row per spatial unit (or a `cell` column).
#' @param type Separate-Laplace fit: one of `"expected"`, `"occupancy"`,
#'   `"conditional"`. Nested-Laplace fit: `"occurrence"`, `"cover_cond"`,
#'   `"cover_exp"`, or `"change"` (the legacy aliases are accepted and mapped).
#' @param include_RE Ignored for the separate-Laplace fit (no spatial
#'   projection); the nested-Laplace fit always projects the shared field.
#' @param times,time_col,level,nsim,draws Nested-Laplace fit only: `times =
#'   c(t1, t2)` and `time_col` drive the `"change"` map; `level` is the credible
#'   level, `nsim` the draw count, `draws` whether to attach the draw matrices.
#' @param ... Unused.
#' @return Separate-Laplace fit: a numeric vector. Nested-Laplace fit: a
#'   `tobs_prediction`.
#' @export
predict.cover_fit <- function(object, newdata = NULL,
                                     type = NULL, include_RE = FALSE,
                                     times = NULL, time_col = NULL,
                                     level = 0.95, nsim = 1000L, draws = TRUE,
                                     ...) {
  # Nested-Laplace shared-field fit: route through the unified joint predict
  # substrate (gcol33/tulpaObs#23). Map the legacy fixed-effects type names onto
  # the joint vocabulary so old calls keep working.
  if (!is.null(.tobs_joint_fit(object))) {
    if (is.null(type)) type <- "occurrence"
    type <- switch(type,
                   expected    = "cover_exp",
                   occupancy   = "occurrence",
                   conditional = "cover_cond",
                   type)
    return(.tobs_predict_joint(object, newdata = newdata, type = type,
                               times = times, level = level, nsim = nsim,
                               draws = draws, time_col = time_col))
  }

  if (is.null(type)) type <- "expected"
  type <- match.arg(type, c("expected", "occupancy", "conditional"))
  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` is required.", call. = FALSE)
  }
  if (isTRUE(include_RE) && !is.null(object$encoding$spatial_spec)) {
    message("predict.cover_fit(): spatial RE projection at new ",
            "locations is not implemented for the separate-Laplace fit; ",
            "returning fixed-effects-only predictions.")
  }

  X <- stats::model.matrix(object$encoding$formula, newdata)
  if (ncol(X) != length(object$beta_occ)) {
    stop("Design-matrix column count (", ncol(X), ") does not match the ",
         "fitted model (", length(object$beta_occ), "). Check `newdata`.",
         call. = FALSE)
  }

  eta_occ <- as.numeric(X %*% object$beta_occ)
  eta_pos <- as.numeric(X %*% object$beta_pos)
  p  <- stats::plogis(eta_occ)
  positive <- object$positive %||% "lognormal"
  mu <- if (positive == "beta") {
    stats::plogis(eta_pos)
  } else {
    exp(eta_pos + object$sigma_pos^2 / 2)
  }

  switch(
    type,
    expected    = p * mu,
    occupancy   = p,
    conditional = mu
  )
}


# ---------------------------------------------------------------------------
# Print / summary
# ---------------------------------------------------------------------------

#' @export
print.cover_fit <- function(x, ...) {
  positive <- x$positive %||% "lognormal"
  cat(sprintf("<cover_fit (%s positive part)>\n", positive))
  cat(sprintf("  N total      : %d\n", x$n_total))
  cat(sprintf("  N positive   : %d (%.1f%%)\n",
              x$n_positive, 100 * x$n_positive / x$n_total))
  if (positive == "lognormal") {
    cat(sprintf("  sigma_pos    : %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos      : %.4f\n", x$phi_pos))
  }
  cat(sprintf("  converged    : %s\n",
              if (isTRUE(x$converged)) "yes" else "no"))
  if (!is.null(x$sla_status) && !identical(x$sla_status, "off")) {
    cat(sprintf("  marginals    : %s\n", x$sla_status))
  }
  cat("\nPresence (binomial logit):\n")
  print(.coef_table(x$beta_occ, x$se_occ))
  pos_header <- if (positive == "beta") {
    "Positive (beta, logit link, on y > 0):"
  } else {
    "Positive (Gaussian on log y > 0):"
  }
  cat("\n", pos_header, "\n", sep = "")
  print(.coef_table(x$beta_pos, x$se_pos))
  invisible(x)
}

#' @export
summary.cover_fit <- function(object, ...) {
  # NUTS fit: return the per-parameter posterior table (mean / sd / quantiles
  # plus the cross-chain Rhat / ESS the convergence list carries), matching the
  # generic NUTS summary surface so the sampler diagnostics are visible.
  if (!is.null(object$nuts) && !is.null(object$draws)) {
    return(.tobs_cover_nuts_summary(object))
  }
  # Arm labels (gcol33/tulpaObs#61): the two hurdle arms are `presence`
  # (the y > 0 Bernoulli arm) and `positive` (the y | y > 0 arm). The `to =`
  # argument of a spatial() bar validates against these labels, so summary()
  # prints the same names (formula label == output label).
  out <- list(
    family       = object$family,
    positive     = object$positive %||% "lognormal",
    n_total      = object$n_total,
    n_positive   = object$n_positive,
    sigma_pos    = object$sigma_pos,
    phi_pos      = object$phi_pos,
    converged    = object$converged,
    presence     = .coef_table(object$beta_occ, object$se_occ),
    positive_arm = .coef_table(object$beta_pos, object$se_pos),
    log_marginal = object$log_marginal,
    hyperpar     = object$hyperpar
  )
  class(out) <- "summary.cover_fit"
  out
}

#' @export
print.summary.cover_fit <- function(x, ...) {
  cat("Cover hurdle fit summary\n")
  cat(sprintf("  positive part: %s\n", x$positive))
  cat(sprintf("  N total = %d, N positive = %d\n", x$n_total, x$n_positive))
  if (x$positive == "lognormal") {
    cat(sprintf("  sigma_pos = %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos   = %.4f\n", x$phi_pos))
  }
  cat(sprintf("  log marginal: occ = %.3f, pos = %.3f\n",
              x$log_marginal["occ"], x$log_marginal["pos"]))
  cat("\nPresence:\n"); print(x$presence)
  pos_header <- if (x$positive == "beta") "Positive (beta, logit):" else "Positive (Gaussian):"
  cat("\n", pos_header, "\n", sep = ""); print(x$positive_arm)
  invisible(x)
}


# ---------------------------------------------------------------------------
# Joint nested-Laplace fit (Phase 1c lognormal, Phase 1d beta)
# ---------------------------------------------------------------------------

# Collapse the occurrence (binomial) arm to its exact sufficient statistic.
# Observations agreeing on the occurrence design row AND every per-observation
# component of the linear predictor -- the spatial cell, any temporal / RE block
# index, any SVC weight -- are exchangeable Bernoulli trials: replacing them with
# one Binomial row (n = count, y = successes) leaves the log-likelihood, gradient
# and Hessian pointwise unchanged (the engine's binomial kernel carries no
# combinatorial constant) while cutting the row count. `keys` is the named list
# of those per-observation vectors; the X-row plus every key forms the grouping
# key. Returns the aggregated `y`, `n` (n_trials), `X`, the per-group
# representative of each key, and `rep_i` (the representative row per group).
.cover_aggregate_occ <- function(y, X, keys) {
  parts <- c(as.data.frame(X, stringsAsFactors = FALSE), keys)
  gid   <- as.integer(factor(do.call(paste, c(parts, list(sep = "\r")))))
  ord   <- order(gid)
  rep_i <- ord[!duplicated(gid[ord])]               # first row per group, group order
  list(
    y     = as.numeric(rowsum(as.numeric(y), gid)),   # rowsum orders by sorted gid
    n     = as.integer(tabulate(gid, nbins = max(gid))),
    X     = X[rep_i, , drop = FALSE],
    keys  = lapply(keys, function(v) v[rep_i]),
    rep_i = rep_i
  )
}

# Aggregate the occurrence arm in place against `keys` and reset its trial /
# RE bookkeeping. Returns the updated arm plus the per-group key representatives
# for the caller to scatter back onto the latent blocks.
.cover_apply_occ_agg <- function(arm_occ, keys) {
  og <- .cover_aggregate_occ(arm_occ$y, arm_occ$X, keys)
  arm_occ$y        <- og$y
  arm_occ$n_trials <- og$n
  arm_occ$X        <- og$X
  arm_occ$re_idx   <- rep(0, length(og$y))
  list(arm_occ = arm_occ, keys = og$keys)
}

# Collapse the positive (beta) arm to its exact grouped sufficient statistics
# (tulpaObs#49). Beta is not a count family, so there is no single-row collapse:
# a group of plots sharing the positive design row AND every per-observation
# component of its linear predictor (cell, trend weight, RE/time index) are
# exchangeable Beta(mu*phi, (1-mu)*phi) draws, and the beta log-density is linear
# in log(y) and log(1-y). One row carrying (n = count, slog_y = sum log(y),
# slog_1my = sum log(1-y)) therefore leaves the log-likelihood, gradient and
# (Fisher) Hessian pointwise unchanged -- the per-arm `n_trials`, `slog_y`,
# `slog_1my` are read by tulpa's built-in beta spec (gcol33/tulpa). `keys` is the
# named list of per-observation latent components; the X-row plus every key forms
# the grouping key.
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

# Build one areal latent block restricted to a SINGLE arm with NO cross-arm copy
# (gcol33/tulpaObs#65). The block carries the field's own precision axis (tau for
# icar/car/car_proper, sigma + rho for bym2) integrated on the outer grid; the
# OTHER arm's per-arm spatial_idx is the all-zero sentinel, which the joint
# multi-block scatter reads as "this arm's rows do not see this block" (the
# `l_b > 0` guard in nested_laplace_joint_multi.h), so the field contributes to
# exactly one arm. This is the no-copy assembler; the shared/copied intercept +
# trend path uses the copy spec instead (one assembler, copy off here, copy on
# there). `slot` is the active arm (1 = occ, 2 = pos); `n_occ` / `n_pos` size the
# sentinel; `idx_active` is the per-obs node code on the active arm.
# `svc_weight` (NULL for the intercept field) is the per-obs design column,
# placed on the active arm with a zero sentinel weight on the inactive arm.
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
  # control$sigma.grid (translated to tau for the intrinsic backends).
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
    block$tau_grid <- sort(1.0 / sigma_grid^2)
    if (block$type == "car_proper") {
      block$rho_car_grid <- as.numeric(control$rho.car.grid %||%
                                       c(0.5, 0.8, 0.95, 0.99))
    }
  } else {
    # bym2 fits as a non-copied block, but its unit-variance field is the
    # rho-mixed phi + theta pair; the arm-specific draw projection reads phi only
    # (the per-arm field reconstruction for the bym2 mix is not wired), so a bym2
    # arm-specific field would mis-project at predict / WAIC time. Restrict to the
    # intrinsic ICAR / proper CAR backends, whose unit field IS phi.
    stop(sprintf(paste0(
      "cover() arm-specific field (%s): areal type '%s' is not supported. ",
      "Arm-specific separate fields (single-arm `to`) use icar / car / ",
      "car_proper (the bym2 phi+theta mix is deferred)."),
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

# Fit the cover hurdle with arm-specific separate spatial latent field(s)
# (gcol33/tulpaObs#65): one or more single-arm `||` bars, each a per-arm areal
# field with its own precision and NO cross-arm copy. Each bar's design columns
# (intercept + covariates) become independent non-copied areal blocks placed on
# that bar's single arm via the 0-sentinel spatial_idx on the other arm; each
# field's precision is integrated on the outer nested-Laplace grid. Separate
# single-arm bars (one to = "presence", one to = "positive") are independent
# per-arm fields with no coupling between them. Same output contract as
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

  # Positive-arm dispersion grid (same regime as the single-field path).
  if (positive == "lognormal") {
    pos_family   <- "gaussian"
    sigma_hat    <- .prefit_lognormal_sigma(enc, control)
    phi_hat      <- sigma_hat
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(sigma_hat / 3), log(sigma_hat * 3), length.out = 7))
  } else {
    pos_family   <- "beta"
    phi_hat      <- 1.0
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(2), log(300), length.out = 7))
  }

  arm_occ <- list(
    y = as.numeric(enc$occ_data$y), n_trials = enc$occ_data$n_trials,
    X = enc$occ_data$X, re_idx = rep(0, N), n_re_groups = 0L,
    sigma_re = 1.0, family = "binomial", phi = 1.0
  )
  arm_pos <- list(
    y = as.numeric(enc$pos_data$y), n_trials = rep(1L, N_pos),
    X = enc$pos_data$X, re_idx = rep(0, N_pos), n_re_groups = 0L,
    sigma_re = 1.0, family = pos_family, phi = phi_hat
  )

  # Opt-in fixed-effect priors (cover_priors()), as on the single-field path.
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
        label = label)
    }
  }

  joint_control <- list(
    max_iter  = control$max.iter  %||% 50L,
    tol       = control$tol       %||% 1e-6,
    n_threads = control$n.threads %||% 1L,
    n_threads_outer = control$n.threads.outer %||% 1L,
    store_Q   = TRUE,
    hessian   = control$hessian   %||% (if (positive == "beta") "fisher" else "lm"),
    adaptive_grid             = control$adaptive.grid             %||% TRUE,
    adaptive_grid_edge_thresh = control$adaptive.grid.edge.thresh %||% 0.02,
    adaptive_grid_max_passes  = control$adaptive.grid.max.passes  %||% 1L,
    progress          = control[["progress"]]      %||% TRUE,
    progress.every    = control$progress.every    %||% 0L,
    progress.throttle = control$progress.throttle %||% 2,
    progress.file     = control$progress.file     %||% "",
    checkpoint        = control$checkpoint,
    # Each field carries 1 (icar/car) or 2 (bym2/car_proper) latent axes; with
    # the pos-arm phi axis the dense outer tensor grows fast, so the mode-centred
    # CCD is the default for >= 3 axes (forwarded; the engine declines back to
    # the tensor grid on a ridge). Override via control$integration.
    integration       = control$integration %||% "ccd"
  )

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = blocks,         # list-of-blocks -> multi-block path
    copy      = NULL,           # NO copy: each field on one arm only
    phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
    prior_sigma = control$prior.sigma,
    control = joint_control
  )

  # Shared per-arm beta post-processing (identical to the single-field path).
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
  se_occ <- sqrt(pmax(0, var_of_means_occ + mean_of_var_occ))
  se_pos <- sqrt(pmax(0, var_of_means_pos + mean_of_var_pos))

  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive == "lognormal") {
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
# onto the positive arm (gcol33/tulpaObs#64). The p design columns of the bar
# become p coupled areal fields with a free cross-covariance Sigma (x) Q^-1
# (within-arm covariance among the fields, integrated over the outer CCD in
# log-Cholesky coordinates); the whole correlated field is copied onto the
# positive arm with one estimated amplitude alpha (the cross-arm transfer). Same
# output contract as fit_cover_hurdle_joint_nested so decode_cover_hurdle_joint
# consumes it unchanged.
.fit_cover_hurdle_joint_mcar <- function(enc, data, positive = enc$positive,
                                         control = list(), priors = NULL) {
  mc    <- enc$mcar
  p     <- mc$n_fields
  N     <- enc$N
  N_pos <- length(enc$pos_data$y)
  idx_pos <- enc$idx_pos

  # Positive-arm dispersion grid (same regime as the single-field path).
  if (positive == "lognormal") {
    pos_family   <- "gaussian"
    sigma_hat    <- .prefit_lognormal_sigma(enc, control)
    phi_hat      <- sigma_hat
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(sigma_hat / 3), log(sigma_hat * 3), length.out = 7))
  } else {
    pos_family   <- "beta"
    phi_hat      <- 1.0
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(2), log(300), length.out = 7))
  }

  arm_occ <- list(
    y = as.numeric(enc$occ_data$y), n_trials = enc$occ_data$n_trials,
    X = enc$occ_data$X, re_idx = rep(0, N), n_re_groups = 0L,
    sigma_re = 1.0, family = "binomial", phi = 1.0
  )
  arm_pos <- list(
    y = as.numeric(enc$pos_data$y), n_trials = rep(1L, N_pos),
    X = enc$pos_data$X, re_idx = rep(0, N_pos), n_re_groups = 0L,
    sigma_re = 1.0, family = pos_family, phi = phi_hat
  )

  # Opt-in fixed-effect priors (cover_priors()), as on the single-field path.
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

  # MCAR block: per-arm cell index (occ, pos) and per-field per-arm design
  # weight (occ, pos). The positive arm slices the occ weights / cell index by
  # `idx_pos`, exactly as the `||` trend weight is subset.
  idx_occ <- mc$idx_occ
  idx_pos_cell <- idx_occ[idx_pos]
  field_weight <- lapply(mc$field_weight_occ, function(w_occ)
    list(as.numeric(w_occ), as.numeric(w_occ[idx_pos])))

  mcar_block <- list(
    type            = "mcar",
    n_spatial_units = mc$n_spatial_units,
    n_fields        = as.integer(p),
    adj_row_ptr     = mc$adj_row_ptr,
    adj_col_idx     = mc$adj_col_idx,
    n_neighbors     = mc$n_neighbors,
    spatial_idx     = list(as.integer(idx_occ), as.integer(idx_pos_cell)),
    field_weight    = field_weight
  )

  alpha_grid <- control$alpha.grid %||%
    c(0, exp(seq(log(0.1), log(3), length.out = 5)))
  copy_spec <- list(arm = "pos", block = 1L, alpha_grid = as.numeric(alpha_grid))

  joint_control <- list(
    max_iter  = control$max.iter  %||% 50L,
    tol       = control$tol       %||% 1e-6,
    n_threads = control$n.threads %||% 1L,
    n_threads_outer = control$n.threads.outer %||% 1L,
    store_Q   = TRUE,
    hessian   = control$hessian   %||% (if (positive == "beta") "fisher" else "lm"),
    adaptive_grid             = control$adaptive.grid             %||% TRUE,
    adaptive_grid_edge_thresh = control$adaptive.grid.edge.thresh %||% 0.02,
    adaptive_grid_max_passes  = control$adaptive.grid.max.passes  %||% 1L,
    progress          = control[["progress"]]      %||% TRUE,
    progress.every    = control$progress.every    %||% 0L,
    progress.throttle = control$progress.throttle %||% 2,
    progress.file     = control$progress.file     %||% "",
    checkpoint        = control$checkpoint,
    # The MCAR block carries p(p+1)/2 + 1 latent axes (log-Cholesky + alpha),
    # so the outer grid uses the mode-centred CCD by default (the same recipe
    # the single-arm MCAR uses); the dense tensor would blow up.
    integration       = control$integration %||% "ccd"
  )

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = list(mcar_block),
    copy      = copy_spec,
    phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
    prior_sigma = control$prior.sigma,
    prior_alpha = control$prior.alpha,
    control = joint_control
  )

  # Shared per-arm beta post-processing (identical to the single-field path).
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
  se_occ <- sqrt(pmax(0, var_of_means_occ + mean_of_var_occ))
  se_pos <- sqrt(pmax(0, var_of_means_pos + mean_of_var_pos))

  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive == "lognormal") {
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
  alpha_mu <- as.numeric(fit$theta_mean[["b1.alpha"]] %||% NA_real_)
  alpha_sd <- as.numeric(fit$theta_sd[["b1.alpha"]]   %||% NA_real_)

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
#'   `sigma_grid`, `rho_grid`, `rho_car_grid`, `sigma_pos_grid`,
#'   `phi_init`, `phi_bounds` (the last two are forwarded to the beta
#'   pre-fit when `positive = "beta"`). For ICAR / CAR_proper backends
#'   `tau_grid` is also accepted and translated to `sigma_grid` as
#'   `sigma = 1 / sqrt(tau)`. The cover-arm field amplitude rides its
#'   own `sigma_pos_grid` axis. Regularizing hyperpriors on the joint
#'   (sigma, alpha) axes can be set via `prior_sigma` (donor amplitude)
#'   and `prior_alpha` (copy coefficient) — each a length-2 list
#'   `list(family, params)` matching tulpa's `prior_sigma` / `prior_alpha`
#'   args. The prior on alpha directly regularizes the copy scalar at small
#'   `n_pos`, replacing the per-arm `prior_sigma_pos` of the pre-reparam
#'   API.
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
  if (!positive %in% c("lognormal", "beta")) {
    stop("method = 'nested_laplace'/'nested_laplace_sla' for cover() supports positive = ",
         "'lognormal' or 'beta'. Got '", positive, "'.", call. = FALSE)
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
  # Replicated CAR (gcol33/tulpaObs#82): a `by` factor on the shared bar replicates
  # the field across the factor's levels. Build I_L (x) Q and offset each
  # observation's node into its level's copy (tulpa::tulpa_bar_field_replicate),
  # then resolve the field over the replicated graph so its precision is the
  # block-diagonal Kronecker and the field hyperparameters are shared across
  # levels (one sigma[, rho_car]). The coupled trend block copies the same
  # replicated structure and the copy onto the positive arm carries the whole
  # replicated field at the one estimated alpha, so `by` composes with the shared
  # field, the trend, and the cross-arm copy unchanged. No `by` is the identity.
  # Injecting the offset index as the field's node column lets the SAME
  # prior_from_spec build the precision and resolve the index over the replicated
  # graph in one pass.
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
  # `col` and copied onto the positive arm with its own alpha axis
  # (gcol33/tulpaObs#59). This is the analogue of the INLA joint model's
  # `f(cell.slope, time, model = "besag") + f(cell.slope.ab, time, copy =)`
  # spatially-varying time trend (two coupled besag fields, two copy
  # coefficients). The per-observation weight enters each arm's field
  # contribution via the engine's per-block `svc_weight`. `enc$trend` carries
  # the per-observation weight over `data_obs`; subset it to the positive arm.
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

  arm_occ <- list(
    y           = as.numeric(enc$occ_data$y),
    n_trials    = enc$occ_data$n_trials,
    X           = enc$occ_data$X,
    spatial_idx = as.integer(spi_full),
    re_idx      = rep(0, N),
    n_re_groups = 0L,
    sigma_re    = 1.0,
    family      = "binomial",
    phi         = 1.0
  )

  # Positive-arm dispersion. Both regimes integrate phi on the outer joint
  # hyperparameter grid; `arm_pos$phi` is a placeholder overridden per grid
  # point by the joint engine.
  #   * lognormal: gaussian phi is the noise SD. The non-spatial residual SD
  #     from `.prefit_lognormal_sigma()` is an upper bound on the truth (it
  #     absorbs the alpha-mediated field variance), so we use it as the
  #     *centre* of a 7-point log-spaced grid spanning [sigma_hat/3,
  #     sigma_hat*3]. Neighbour-ratio ~1.44 keeps the inner Laplace
  #     warm-starts close enough that adaptive densification rarely fires.
  #     Override via `control$phi.grid`.
  #   * beta:      phi is integrated on the outer joint hyperparameter grid
  #     (tulpaObs#7). 7 log-spaced points span 2..300 (neighbour-ratio
  #     ~2.4); the joint engine's mode-tracked interior densification
  #     (gcol33/tulpa#19 follow-up) adds 1-2 midpoint cells around the peak
  #     when adjacent grid levels carry density above the edge threshold,
  #     so the *effective* phi resolution near the peak matches the
  #     previous 13-point default while the baseline cell count drops ~46%.
  #     History: fixed 13 was set when adaptive_grid was off-by-default;
  #     5 points gave ~18% mean bias, 9 points ~12% under the static grid.
  #     Refinement now closes that gap dynamically.
  if (positive == "lognormal") {
    pos_family   <- "gaussian"
    sigma_hat    <- .prefit_lognormal_sigma(enc, control)
    phi_hat      <- sigma_hat
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(sigma_hat / 3), log(sigma_hat * 3), length.out = 7))
  } else {
    pos_family   <- "beta"
    phi_hat      <- 1.0
    phi_grid_pos <- control$phi.grid %||%
      exp(seq(log(2), log(300), length.out = 7))
  }

  arm_pos <- list(
    y           = as.numeric(enc$pos_data$y),
    n_trials    = rep(1L, N_pos),
    X           = enc$pos_data$X,
    spatial_idx = as.integer(spi_pos),
    re_idx      = rep(0, N_pos),
    n_re_groups = 0L,
    sigma_re    = 1.0,
    family      = pos_family,
    phi         = phi_hat
  )

  # Opt-in fixed-effect priors (cover_priors()). The joint engine reads a
  # per-arm `beta_prior_mean` / `beta_prior_prec` on each response and replaces
  # its uniform weak default with the quadratic penalty. Mirrors the separate-
  # Laplace path (.cover_arm_prior -> tulpa_laplace beta_prior): the natural-
  # scale numbers are applied at face value to the (autoscaled) design, where
  # every predictor is O(1), so a weakly-informative sd is a sensible ridge.
  # Attaching here lets the priors ride through the aggregation/scatter steps
  # (they mutate rows, not the design columns the prior keys on). Precisions are
  # floored at the engine's own weak default so an Inf-sd bucket reproduces the
  # pre-existing weak ridge rather than dropping the diagonal. NULL / FALSE /
  # "none" leave both arms unpenalised.
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

  # Strip the per-obs spatial_idx (tulpa_nested_laplace_joint takes it per
  # arm) and the legacy rho_bounds field (joint car_proper uses rho_car_grid).
  # Forward control-grid overrides per backend.
  #
  # gcol33/tulpa#18: the engine now parameterizes the joint outer grid as
  # (sigma_occ, sigma_pos) instead of (sigma, alpha). For ICAR / CAR_proper
  # the spatial field is unit-precision and the donor-arm amplitude lives
  # on `prior$sigma_grid` (legacy `tau_grid` from prior_from_spec is
  # translated below). The cover-arm field amplitude lives on
  # `copy$sigma_pos_grid`. alpha is recovered post-hoc as
  # sigma_pos / sigma_occ from the joint posterior.
  prior_for_joint <- prior
  prior_for_joint$spatial_idx <- NULL
  prior_for_joint$rho_bounds  <- NULL
  # Translate any legacy tau_grid that prior_from_spec attached to an
  # ICAR / CAR_proper prior — the joint engine takes sigma = 1/sqrt(tau)
  # as its donor-amplitude axis under the new parameterization.
  if (!is.null(prior_for_joint$tau_grid) &&
      is.null(prior_for_joint$sigma_grid)) {
    prior_for_joint$sigma_grid <- 1.0 / sqrt(as.numeric(prior_for_joint$tau_grid))
    prior_for_joint$tau_grid   <- NULL
  } else if (!is.null(prior_for_joint$tau_grid)) {
    prior_for_joint$tau_grid <- NULL
  }
  if (!is.null(control$sigma.grid))   prior_for_joint$sigma_grid   <- control$sigma.grid
  if (!is.null(control$rho.grid))     prior_for_joint$rho_grid     <- control$rho.grid
  if (!is.null(control$tau.grid)) {
    prior_for_joint$sigma_grid <- 1.0 / sqrt(as.numeric(control$tau.grid))
  }
  if (!is.null(control$rho.car.grid)) prior_for_joint$rho_car_grid <- control$rho.car.grid

  if (!is.null(control$sigma.pos.grid)) {
    sigma_pos_grid <- as.numeric(control$sigma.pos.grid)
  } else {
    sigma_donor <- prior_for_joint$sigma_grid %||%
      exp(seq(log(0.1), log(3), length.out = 5))
    sigma_pos_grid <- as.numeric(sigma_donor)
  }

  # Direct (sigma, alpha) copy axis: the cover arm sees the shared field at
  # amplitude alpha * sigma_donor. The single-block path declares alpha on the
  # pos arm via field_coef (the engine takes no single-block `copy`); the
  # multi-block branches carry their own per-block alpha grids.
  alpha_grid <- control$alpha.grid %||%
    c(0, exp(seq(log(0.1), log(3), length.out = 5L)))

  # Outer joint-grid integration controls, shared by the multi-block and
  # single-block dispatch. The dense outer tensor (sigma_occ x [rho] x
  # sigma_pos x phi_pos) concentrates almost all posterior mass on a handful
  # of cells, but the inner latent mode moves substantially across the grid,
  # so the cheap-pass prune is OFF by default: the full-grid solve is the
  # correct default (gcol33/tulpaObs#20). The rank-safe speed path is the
  # adaptive grid (`adaptive_grid = TRUE`): it brackets the mode with FULL
  # inner solves and densifies near it, so it never approximates the
  # marginal and cannot drop the true mode. The cheap-pass prune is available
  # opt-in via control$prune = TRUE; it is now rank-faithful (a neighbour-
  # warm-started lattice sweep) and gated (a safety check falls back to the
  # full grid if the screen's ranking looks unreliable), but the correct
  # full grid remains the default. Override via control$prune /
  # control$prune.tol.
  joint_control <- list(
    max_iter  = control$max.iter  %||% 50L,
    tol       = control$tol       %||% 1e-6,
    n_threads = control$n.threads %||% 1L,
    # Outer-grid parallelism (gcol33/tulpa#46 lever 2). The sparse joint driver
    # dispatches grid cells across n_threads_outer threads, each with its own
    # replicated cell-solve state (the engine clamps the count if the replicas
    # would be too large). Default 1 (serial outer, prior behaviour). Preferred
    # over inner per-observation threads on many-core hardware for the cover
    # hurdle's large outer grid, where the expensive mode-region cells dominate.
    n_threads_outer = control$n.threads.outer %||% 1L,
    store_Q   = TRUE,
    # Inner-Newton curvature (gcol33/tulpa#46). The beta positive arm's observed
    # mixture Hessian is indefinite away from the mode, so observed-curvature
    # Newton steps stall and the inner Newton hits max.iter in every grid cell.
    # Expected/Fisher curvature is PSD by construction and converges in ~12
    # steps. The final mode-pass always re-factorizes with the observed Hessian,
    # so the reported SEs, log_det and grid weights are unchanged. The lognormal
    # arm is exactly quadratic (one inner step), so observed curvature is already
    # optimal -> keep "lm".
    hessian   = control$hessian   %||% (if (positive == "beta") "fisher" else "lm"),
    prune     = control$prune     %||% FALSE,
    prune_tol = control$prune.tol %||% 1e-4,
    adaptive_grid             = control$adaptive.grid             %||% TRUE,
    adaptive_grid_edge_thresh = control$adaptive.grid.edge.thresh %||% 0.02,
    adaptive_grid_max_passes  = control$adaptive.grid.max.passes  %||% 1L,
    # Outer-grid progress + ETA (gcol33/tulpa#45, tulpaObs#43). Two independent
    # channels, both ON by default:
    #   * `progress` gates the Rcout console line -- the progress bar. ON by
    #     default (NOT tied to `verbose`); set control$progress = FALSE to
    #     silence it. This fit runs for hours, so the bar earns its place.
    #   * `progress.file` writes the ETA to disk and is emitted whenever it is
    #     non-empty, INDEPENDENT of `progress`/`verbose` -- it is the only
    #     channel that survives a detached Start-Process stdout buffer, and a
    #     detached fit is exactly when it is the sole liveness signal
    #     (gcol33/tulpa#53). The engine builds GridProgress when either channel
    #     is wanted, so a quiet console with a heartbeat file still reports.
    # `[[` (exact) not `$`: `control$progress` prefix-matches `progress.file`.
    progress          = control[["progress"]]      %||% TRUE,
    progress.every    = control$progress.every    %||% 0L,
    progress.throttle = control$progress.throttle %||% 2,
    progress.file     = control$progress.file     %||% "",
    # Grid-cell checkpoint/resume (gcol33/tulpa#50). A full-field cover-hurdle
    # fit runs for hours; `control$checkpoint = list(path =, resume =)` makes
    # the outer grid append each completed cell to `path` and a resume run load
    # the finished cells and solve only the rest, so a killed/rebooted fit
    # resumes instead of restarting. Forwarded verbatim to the engine.
    checkpoint        = control$checkpoint,
    # Outer-grid node layout (gcol33/tulpa#61, tulpaObs#31). On the coupled-trend
    # multi-block path (>= 3 latent axes: intercept + trend sigma/alpha) "ccd"
    # places a central composite design over the latent axes and crosses the
    # pos-arm phi tensor on top; "grid" forces the dense tensor. Forwarded so a
    # two-field trend fit can request CCD; NULL falls through to the engine
    # default.
    integration       = control$integration
  )

  # Exact sufficient-statistic reduction of the occurrence (binomial) arm,
  # default ON (tulpaObs#48). The collapse is pointwise exact -- observations
  # sharing the occurrence design row AND every per-observation latent component
  # (cell, trend weight, RE/time index) are exchangeable Bernoulli trials, so one
  # Binomial row (n = count, y = successes) leaves the log-likelihood, gradient
  # and Hessian unchanged. Multi-seed parameter recovery on the aggregated path
  # holds against simulated truth (test-cover-hurdle-aggregate-recovery.R), so
  # the reduction is the default; set control$aggregate.occ = FALSE for the
  # full per-plot occurrence arm. `[[` (exact), never `$` (prefix-matching).
  do_agg_occ <- !isFALSE(control[["aggregate.occ"]])

  # Exact grouped sufficient-statistic reduction of the positive (beta) arm,
  # default ON for the beta arm (tulpaObs#49). Beta has no single-row collapse,
  # so plots sharing the positive design row AND every per-observation latent
  # component are collapsed to one row carrying (n, sum log y, sum log(1-y));
  # tulpa's built-in beta spec reads those sufficient statistics. Byte-identical
  # to the full per-plot beta arm (test-cover-hurdle-aggregate-pos.R), with the
  # both-arms-aggregated default behind a multi-seed parameter-recovery suite
  # (test-cover-hurdle-aggregate-recovery.R), so the reduction is the default;
  # set control$aggregate.pos = FALSE for the full per-plot positive arm. The
  # collapse is beta-only -- a lognormal positive arm would need its own
  # (n, sum, sum-of-squares) statistics, so an EXPLICIT aggregate.pos = TRUE
  # errors there rather than silently no-op, while the default leaves a non-beta
  # arm untouched. `[[` (exact), never `$` (prefix-matching).
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
  # joint engine. Copy semantics remain on the spatial block (sigma_occ /
  # sigma_pos), other blocks are shared identically across the two arms
  # (no per-arm scale).
  if (has_trend) {
    # Coupled trend path: two shared areal blocks on the same graph -- block 1
    # the unweighted intercept field, block 2 the per-observation-weighted SVC
    # field -- each copied onto the positive arm with its own alpha axis. The
    # engine's multi-block driver (list-valued prior + list-valued copy) carries
    # the per-block svc_weight and per-block alpha; the returned theta_grid axes
    # are b<k>.sigma / b<k>.alpha (gcol33/tulpaObs#15 on the cover hurdle).
    base_block <- prior_for_joint
    if (is.null(base_block$sigma_grid)) {
      base_block$sigma_grid <- exp(seq(log(0.1), log(3), length.out = 5))
    }
    if (tolower(base_block$type) == "bym2" && is.null(base_block$rho_grid)) {
      base_block$rho_grid <- c(0.25, 0.5, 0.75)
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
      c(0, exp(seq(log(0.1), log(3), length.out = 5)))
    alpha_grid_trend <- control$alpha.grid.trend %||% alpha_grid_base

    prior_coupled <- list(base_block, trend_block)
    copy_coupled  <- list(
      list(arm = "pos", block = 1L, alpha_grid = as.numeric(alpha_grid_base)),
      list(arm = "pos", block = 2L, alpha_grid = as.numeric(alpha_grid_trend))
    )
    arm_occ$spatial_idx <- NULL
    arm_pos$spatial_idx <- NULL
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = prior_coupled,
      copy      = copy_coupled,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma = control$prior.sigma,
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
      sigma_pos_grid = sigma_pos_grid
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
      arm_occ$spatial_idx <- as.integer(ag$keys$idx1)
    }
    if (do_agg_pos) {
      agp     <- .cover_apply_pos_agg(arm_pos, list(idx1 = arm_pos$spatial_idx))
      arm_pos <- agp$arm_pos
      arm_pos$spatial_idx <- as.integer(agp$keys$idx1)
    }
    # Adaptive grid forwarding. Defaults match the joint engine's defaults
    # (`adaptive_grid = TRUE`, threshold 0.02, one pass) and triggered the
    # under-coverage fix in INLAabun D3 — see gcol33/tulpaObs#8. Pass
    # `control$adaptive.grid = FALSE` to recover the legacy fixed-grid
    # behaviour for reproducibility checks.
    arm_pos$field_coef <- list(name = "alpha", grid = as.numeric(alpha_grid))
    fit <- tulpa::tulpa_nested_laplace_joint(
      responses = list(occ = arm_occ, pos = arm_pos),
      prior     = prior_for_joint,
      phi_grid  = if (!is.null(phi_grid_pos)) list(pos = phi_grid_pos) else NULL,
      prior_sigma = control$prior.sigma,
      prior_alpha = control$prior.alpha,
      control = joint_control
    )
  }

  # Posterior-weighted mean / SE for the per-arm beta blocks.
  layout <- fit$arm_layout
  p_occ  <- layout$p[1]
  p_pos  <- layout$p[2]
  bocc_idx <- layout$beta_start[1] + seq_len(p_occ)
  bpos_idx <- layout$beta_start[2] + seq_len(p_pos)

  # The engine returns per-cell modes and per-cell precision blocks in the
  # *scaled* design's parameterization (see `encode_cover_hurdle()` and
  # gcol33/tulpaObs#9). Transform per cell, then aggregate to natural-scale
  # posterior moments. Doing it cell-by-cell on the full constrained vcov
  # block preserves the intercept's cross-covariance contribution; a
  # diag-only approach would underestimate the intercept SE.
  scale_occ <- enc$scale_occ %||% .scale_meta(enc$occ_data$X)
  scale_pos <- enc$scale_pos %||% .scale_meta(enc$pos_data$X)
  T_occ <- .scale_transform(scale_occ)
  T_pos <- .scale_transform(scale_pos)

  modes_occ_sc <- fit$modes[, bocc_idx, drop = FALSE]
  modes_pos_sc <- fit$modes[, bpos_idx, drop = FALSE]
  # Per-cell transform: modes_nat[k, ] = T %*% modes_sc[k, ]
  modes_occ <- modes_occ_sc %*% t(T_occ)
  modes_pos <- modes_pos_sc %*% t(T_pos)

  beta_occ <- as.numeric(crossprod(fit$weights, modes_occ))
  beta_pos <- as.numeric(crossprod(fit$weights, modes_pos))

  # Var-of-means + Mean-of-Var, both in natural scale.
  var_of_means_occ <- as.numeric(crossprod(fit$weights, modes_occ^2)) - beta_occ^2
  var_of_means_pos <- as.numeric(crossprod(fit$weights, modes_pos^2)) - beta_pos^2

  inner_blocks <- .joint_inner_vcov_block(fit, c(bocc_idx, bpos_idx))
  if (is.null(inner_blocks)) {
    mean_of_var_occ <- rep(0, p_occ)
    mean_of_var_pos <- rep(0, p_pos)
  } else {
    occ_rows <- seq_along(bocc_idx)
    pos_rows <- length(bocc_idx) + seq_along(bpos_idx)
    n_grid_eff <- length(inner_blocks)
    diag_occ <- matrix(0, n_grid_eff, p_occ)
    diag_pos <- matrix(0, n_grid_eff, p_pos)
    for (k in seq_len(n_grid_eff)) {
      V_block <- inner_blocks[[k]]
      if (is.null(V_block)) next
      V_occ_sc <- V_block[occ_rows, occ_rows, drop = FALSE]
      V_pos_sc <- V_block[pos_rows, pos_rows, drop = FALSE]
      V_occ_nat <- T_occ %*% V_occ_sc %*% t(T_occ)
      V_pos_nat <- T_pos %*% V_pos_sc %*% t(T_pos)
      diag_occ[k, ] <- pmax(diag(V_occ_nat), 0)
      diag_pos[k, ] <- pmax(diag(V_pos_nat), 0)
    }
    w_eff <- fit$weights[seq_len(n_grid_eff)]
    mean_of_var_occ <- as.numeric(crossprod(w_eff, diag_occ))
    mean_of_var_pos <- as.numeric(crossprod(w_eff, diag_pos))
  }

  se_occ <- sqrt(pmax(0, var_of_means_occ + mean_of_var_occ))
  se_pos <- sqrt(pmax(0, var_of_means_pos + mean_of_var_pos))

  # Dispersion summary on the positive arm. Both regimes integrate the
  # dispersion scalar on the outer joint hyperparameter grid; read the
  # posterior mean and SD from the engine's `theta_mean` / `theta_sd`. Those
  # are computed against the phi-axis marginal (foreign-axis slice cells
  # filtered out by `.joint_recalibrate_axis_moments`) with Laplace-at-mode
  # SD at the modal cell (gcol33/tulpa#20), so they are grid-spacing-
  # independent. Hand-rolling `sum(weights * theta_grid^2) - mean^2` against
  # `theta_grid[, "phi_pos"]` underestimates SD on sharply peaked axes and
  # additionally collapses on slice cells that pin phi at the modal value
  # while varying other axes -- that's the legacy pattern tulpa#20/#21 were
  # added to replace.
  #
  # The phi axis carries the gaussian residual SD for lognormal and the
  # beta precision for beta; surface under the respective slot names.
  phi_mu <- as.numeric(fit$theta_mean[["phi_pos"]])
  phi_sd <- as.numeric(fit$theta_sd[["phi_pos"]])
  if (positive == "lognormal") {
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
  if (fits$positive == "lognormal") {
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
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = length(enc$idx_pos),
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
    class = c("cover_fit", "tobs_fit", "tulpa_fit")
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
