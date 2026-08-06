# =============================================================================
# family_cover_hurdle.R — Vegetation cover hurdle on the tulpa Laplace backend
#
# A cover response is a hurdle: a Bernoulli occurrence indicator plus a
# positive-magnitude arm fit on the occupied subset. Three routes fit it.
#
# laplace: two independent tulpa_laplace() calls, one per arm. The arms share
# no latent structure, so each is fit at its own mode and the positive-arm
# dispersion comes from that fit (residual SD for the gaussian arms, the beta
# solver's own phi for the beta arms).
#
# nested_laplace: the joint shared-field model via
# tulpa_nested_laplace_joint(), on a bym2 / icar / car / car_proper graph. The
# positive-arm dispersion becomes a hyperparameter axis integrated on the outer
# grid, with a non-spatial pre-fit supplying the grid centre.
#
# nuts: the non-spatial sampler over the exact two-arm coefficient marginal; a
# structured term is rejected with a pointer.
# =============================================================================


# ---------------------------------------------------------------------------
# Dispatcher (called from tobs())
# ---------------------------------------------------------------------------

.dispatch_cover <- function(formula, data, family, detection, y,
                            visits, engine, priors, control,
                            approx = "gaussian_laplace",
                            correction = "none", ...) {
  positive <- family$params$positive
  if (!positive %in% c("lognormal", "lognormal_trunc", "beta", "beta_oi",
                       "ordinal", "gaussian")) {
    stop("cover(response = '", positive, "') is not supported. ",
         "Use 'lognormal', 'lognormal_trunc', 'beta', 'beta_oi', 'ordinal', ",
         "or 'gaussian'.",
         call. = FALSE)
  }
  # The ordinal (interval-censored Gaussian) and truncated-lognormal
  # (upper-truncated Gaussian) positive arms are wired only on the joint
  # nested-Laplace engine: their extra per-observation bound (interval
  # (lower, upper] for ordinal, the truncation ceiling for lognormal_trunc) is
  # consumed by the joint solver's built-in family, which the single-Laplace /
  # NUTS paths do not carry.
  if (positive %in% c("ordinal", "lognormal_trunc") &&
      !identical(engine, "nested_laplace")) {
    stop("cover(response = \"", positive, "\") requires method = ",
         "'nested_laplace' / 'nested_laplace_sla' (got engine '", engine, "'). ",
         "The ", positive, " arm carries a per-observation bound integrated on ",
         "the joint outer grid; the single-Laplace and NUTS paths are not wired ",
         "for it.", call. = FALSE)
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
  # Per-arm formulas (arm = formula): tobs()'s `presence` and `positive` formula
  # args (folded into `...`) give each hurdle arm its own fixed effects. `positive`
  # here is the formula, distinct from the family's positive-arm likelihood
  # (family$params$positive). Absent both, the single `formula` is shared.
  .cover_dots <- list(...)
  enc      <- encode_cover_hurdle(formula, data, y, positive = positive,
                                  breaks = family$params$breaks,
                                  presence_formula = .cover_dots$presence,
                                  positive_formula = .cover_dots$positive)
  temporal <- enc$temporal
  re       <- enc$re

  # A copy() in a per-arm positive formula sets the cross-arm coupling-amplitude
  # grid(s) the presence field is transferred onto the positive arm with. The
  # whole-field form (copy(spatial()) / copy(spatial(), alpha = )) sets one grid
  # on both the intercept (alpha.grid) and trend (alpha.grid.trend) blocks; a
  # bare copy(spatial()) with no alpha leaves control untouched so the fitter's
  # default grid applies (byte-identical to the shared to = both spelling). The
  # per-component form copy(spatial(), terms = list(intercept = , trend = ))
  # sets each block's grid independently, so the intercept and trend decouple --
  # mirroring occu_cover()'s per-component copy grammar.
  if (!is.null(enc$copy_alpha) || !is.null(enc$copy_terms)) {
    if (any(c("alpha.grid", "alpha.grid.trend") %in% names(control))) {
      stop("cover(): set the cross-arm coupling with copy() in the positive ",
           "formula OR control$alpha.grid[.trend], not both.", call. = FALSE)
    }
  }
  if (!is.null(enc$copy_terms)) {
    control <- .cover_apply_copy_terms(enc$copy_terms, enc$trend, control)
  } else if (!is.null(enc$copy_alpha)) {
    control$alpha.grid       <- as.numeric(enc$copy_alpha)
    control$alpha.grid.trend <- as.numeric(enc$copy_alpha)
  }

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
#' @param positive Positive-arm likelihood: one of `"lognormal"`,
#'   `"lognormal_trunc"`, `"beta"`, `"beta_oi"`, `"ordinal"`, `"gaussian"`.
#' @param breaks Ordinal cut points, required when `positive = "ordinal"` and
#'   ignored otherwise.
#' @param autoscale Logical (default `TRUE`); rescale the positive arm away
#'   from the boundary so the Laplace engine never sees an exact 0 or 1.
#' @param presence_formula Optional formula overriding `formula` on the
#'   occurrence arm. `NULL` uses `formula` for both arms.
#' @param positive_formula Optional formula overriding `formula` on the
#'   positive-cover arm. `NULL` uses `formula` for both arms.
#' @return A list with: `occ_data`, `pos_data`, `spatial_spec` (a
#'   `tulpa_spatial` built from the unweighted areal formula term, or NULL),
#'   `trend` (per-observation weight + label from a weighted areal term, or
#'   NULL), `temporal` and `re` (structured terms from the formula, or NULL),
#'   `N`, `idx_pos` (row indices of the positive subset within `data`),
#'   `formula` (the fixed-effects formula), `positive`.
#' @keywords internal
encode_cover_hurdle <- function(formula, data, y,
                                positive = c("lognormal", "lognormal_trunc",
                                             "beta", "beta_oi", "ordinal",
                                             "gaussian"),
                                breaks = NULL,
                                autoscale = TRUE,
                                presence_formula = NULL,
                                positive_formula = NULL) {
  positive <- match.arg(positive)
  if (!is.numeric(y)) stop("`y` must be numeric.", call. = FALSE)
  .tobs_check_site_count(length(y), nrow(data), "values")
  # The identity-Gaussian arm (gcol33/tulpaObs#112) is the delta-normal hurdle
  # for a response on a real, unbounded scale: absence is the exact 0 sentinel,
  # presence is any nonzero magnitude (which may be negative). The bounded-cover
  # arms (beta / lognormal / ordinal) keep the [0, 1] cover-fraction contract and
  # the y > 0 presence rule.
  if (!identical(positive, "gaussian")) {
    rng <- range(y, na.rm = TRUE)
    if (rng[1] < 0 || rng[2] > 1) {
      stop("`y` must be in [0, 1] (got range [", rng[1], ", ", rng[2], "]).",
           call. = FALSE)
    }
  }

  obs_keep <- !is.na(y)
  y_obs    <- y[obs_keep]
  data_obs <- data[obs_keep, , drop = FALSE]
  occur    <- if (identical(positive, "gaussian")) as.integer(y_obs != 0)
              else                                 as.integer(y_obs > 0)

  # Per-arm formulas (arm = formula): `presence` and `positive` carry their own
  # fixed effects, so the two arms get independent designs. The single `formula`
  # (shared across arms) remains the back-compat spelling. An arm is chosen by
  # placement -- write a field in that arm's formula -- and shared across arms
  # with copy(). When per-arm formulas are absent this branch is skipped, so the
  # shared path is untouched.
  copy_alpha <- NULL
  copy_terms_alpha <- NULL
  per_arm <- !is.null(presence_formula) || !is.null(positive_formula)
  if (per_arm) {
    if (is.null(presence_formula) || is.null(positive_formula)) {
      stop("cover(): give BOTH `presence` and `positive` per-arm formulas, or a ",
           "single shared formula (not one arm only).", call. = FALSE)
    }
    # Placement -> arm: split each arm formula into its fixed effects and its
    # spatial-field calls. Each field call is evaluated into its spec and tagged
    # with the arm (spec$to) directly; the tagged specs feed .encode_cover_specs
    # (a field in `positive` becomes an arm-specific positive field, indexed onto
    # the positive rows by the fitter). Fixed effects give each arm its design.
    #
    # copy() is the canonical shared-field spelling: it lives in the `positive`
    # formula and couples the presence field onto the positive arm. Each arm's
    # copy() calls come back from .cover_lift_arm_fields as unevaluated language
    # objects (so a data-dependent term elsewhere in the arm formula is never
    # forced); a copy() on the `presence` formula is rejected -- the anchor arm
    # needs no copy.
    lift_occ <- .cover_lift_arm_fields(presence_formula, "presence")
    lift_pos <- .cover_lift_arm_fields(positive_formula, "positive")
    if (length(lift_occ$copies)) {
      stop("cover(): copy() belongs in the `positive` formula; it copies the ",
           "presence field onto the positive arm.", call. = FALSE)
    }
    occ_arm <- "presence"
    if (length(lift_pos$copies)) {
      reg_env    <- list2env(.tobs_terms, parent = environment(positive_formula))
      pos_copies <- lapply(lift_pos$copies, function(cl) eval(cl, envir = reg_env))
      promo <- .cover_promote_copied_fields(pos_copies, lift_occ$fields)
      lift_occ$fields <- promo$fields
      occ_arm          <- promo$arm
      copy_alpha       <- promo$alpha
      copy_terms_alpha <- promo$alpha_terms
    }
    fe_occ_formula <- lift_occ$fe
    fe_pos_formula <- lift_pos$fe
    arm_specs <- c(
      lapply(lift_occ$fields, function(e) .tobs_eval_arm_field(
        e, occ_arm, data_obs, environment(presence_formula))),
      lapply(lift_pos$fields, function(e) .tobs_eval_arm_field(
        e, "positive", data_obs, environment(positive_formula)))
    )
    # Only the bar form carries an arm-specific / copied field; a plain areal
    # constructor placed in an arm formula has no such spelling.
    for (s in arm_specs) {
      if (!isTRUE(s$is_bar)) {
        stop(sprintf(paste0(
          "cover(): a field placed in a per-arm formula must use the bar form ",
          "spatial(~ 1 + w || cell, graph = adj); got %s()."),
          s$type %||% class(s)[1L]), call. = FALSE)
      }
    }
    if (length(arm_specs)) {
      cover_struct <- .encode_cover_specs(arm_specs, data_obs)
    } else {
      cover_struct <- list(fe = NULL, spatial = NULL, trend = NULL,
                           temporal = NULL, re = NULL, mcar = NULL, armspec = NULL)
    }
  } else {
    # Parse structured terms against the NA-dropped observations so re()/
    # temporal() index codes align with both hurdle arms. An areal spatial
    # term is converted to the tulpa_spatial spec the cover engine consumes.
    cover_struct   <- .encode_cover_terms(formula, data_obs)
    fe_occ_formula <- cover_struct$fe
    fe_pos_formula <- cover_struct$fe
  }
  fe_formula    <- fe_occ_formula

  X_occ_natural <- stats::model.matrix(fe_occ_formula, data_obs)

  # One-inflated Beta (gcol33/tulpaObs#108): the presence arm still models y > 0,
  # but plots recorded at full cover (y = 1) are a genuine point mass, not a near-1
  # continuous value. With a constant inflation probability the likelihood
  # factorizes -- pi is the share of positive plots at the ceiling, and the
  # interior Beta is fit on the (0, 1) plots only -- so the Beta arm drops the
  # ceiling rows and `oi` carries pi for the decode / predict.
  oi <- NULL
  is_pos <- occur == 1L
  if (positive == "beta_oi") {
    tol         <- 1e-6
    at_ceiling  <- is_pos & (y_obs >= 1 - tol)
    n_pos_total <- sum(is_pos)
    n_ceiling   <- sum(at_ceiling)
    is_pos      <- is_pos & (y_obs < 1 - tol)   # Beta arm = interior positives
    if (sum(is_pos) < 2L) {
      stop("cover(response = \"beta_oi\"): fewer than 2 interior (0 < cover < 1) ",
           "plots after setting the ceiling (cover = 1) aside. The interior Beta ",
           "is not identifiable; use positive = \"beta\" or check the response ",
           "scale.", call. = FALSE)
    }
    pi_hat <- n_ceiling / n_pos_total
    oi <- list(n_ceiling = n_ceiling, n_positive = n_pos_total, pi_one = pi_hat,
               pi_one_sd = sqrt(pi_hat * (1 - pi_hat) / n_pos_total))
  }
  data_pos <- data_obs[is_pos, , drop = FALSE]
  y_pos    <- y_obs[is_pos]
  pos_lower <- NULL
  pos_upper <- NULL
  pos_class <- NULL
  if (positive %in% c("lognormal", "lognormal_trunc")) {
    # Both fit a Gaussian on log-cover; lognormal_trunc additionally truncates the
    # latent log-cover at log(1) = 0 (cover <= 1), the truncation ceiling carried
    # in pos_data$trunc_upper below.
    y_pos_resp <- log(y_pos)
  } else if (positive == "gaussian") {
    # Identity-Gaussian arm (gcol33/tulpaObs#112): a plain Gaussian on the raw
    # positive response (already on a real, unbounded scale), with no log
    # transform and no change-of-variable Jacobian.
    y_pos_resp <- y_pos
  } else if (positive == "ordinal") {
    # Interval-censored Gaussian on log-cover with fixed Braun-Blanquet
    # thresholds. Each positive plot's cover (a class midpoint, or an
    # aggregate of several records, 1 - prod(1 - cover)) is censored to the
    # ordinal class band it falls in. `breaks` are the interior boundaries on
    # the (0, 1) cover-fraction scale; class k = (brk[k-1], brk[k]] with open
    # outer classes (lower = -Inf for the lowest, upper = +Inf for the highest).
    # The latent linear predictor lives on log-cover, so the bounds are the log
    # of the class boundaries -- the same scale the lognormal arm's eta uses.
    if (is.null(breaks)) {
      stop("encode_cover_hurdle(positive = \"ordinal\") requires `breaks`.",
           call. = FALSE)
    }
    lthr      <- log(as.numeric(breaks))            # length K-1 (interior, log)
    K1        <- length(breaks) + 1L                # number of classes
    # Lower-closed bands [brks[k-1], brks[k]) so a class representative that sits
    # exactly on its lower boundary maps to its OWN class. The MOTIVATE myscale
    # rep for class "a" is 0.03 = the (1.5-3] | (3-5] cut at 3% -- with upper-
    # closed bands it would fall one class too low. (Open/closed at the boundary
    # is measure-zero for the continuous latent likelihood; only the class
    # ASSIGNMENT of the discrete representatives is affected.)
    pos_class <- findInterval(y_pos, breaks) + 1L   # 1..K1
    pos_lower <- ifelse(pos_class == 1L,  -Inf, lthr[pmax(pos_class - 1L, 1L)])
    pos_upper <- ifelse(pos_class == K1,   Inf, lthr[pmin(pos_class, length(lthr))])
    # The arm's nominal response is the log-cover; the engine's interval family
    # reads (lower, upper) and ignores this, but downstream prefit / summaries
    # use it as the per-plot point value.
    y_pos_resp <- log(y_pos)
  } else {
    # Beta arm needs y strictly in (0, 1). Cap at 1 - 1e-6; lower bound is
    # already guaranteed by occur == 1 + the range check above.
    y_pos_resp <- pmin(y_pos, 1 - 1e-6)
  }
  X_pos_natural <- stats::model.matrix(fe_pos_formula, data_pos)

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

  pos_data <- list(y = y_pos_resp, X = X_pos)
  if (positive == "ordinal") {
    pos_data$lower <- as.numeric(pos_lower)
    pos_data$upper <- as.numeric(pos_upper)
    pos_data$class <- as.integer(pos_class)
  } else if (positive == "lognormal_trunc") {
    # Upper-truncate the latent log-cover at log(1) = 0 (cover in (0, 1]); one
    # ceiling per positive plot, on the same log-cover (predictor) scale as eta.
    pos_data$trunc_upper <- rep(0, length(y_pos_resp))
  }

  list(
    occ_data = list(y = occur, n_trials = rep(1L, length(occur)), X = X_occ),
    pos_data = pos_data,
    spatial_spec = cover_struct$spatial,
    trend        = cover_struct$trend,
    mcar         = cover_struct$mcar,
    armspec      = cover_struct$armspec,
    temporal     = cover_struct$temporal,
    re           = cover_struct$re,
    N            = length(occur),
    idx_pos      = which(is_pos),
    formula      = fe_formula,
    fe_occ       = fe_occ_formula,
    fe_pos       = fe_pos_formula,
    per_arm      = per_arm,
    copy_alpha   = copy_alpha,
    copy_terms   = copy_terms_alpha,
    positive     = positive,
    breaks       = if (positive == "ordinal") as.numeric(breaks) else NULL,
    oi           = oi,
    obs_keep     = obs_keep,
    scale_occ    = scale_occ,
    scale_pos    = scale_pos
  )
}

# Validate the placement tag on a bar field. `spec$to` is set by the lift, not by
# the user, and must resolve to a subset of the cover arms (presence, positive). A
# single shared latent copied across arms carries both; a correlated (`|` / MCAR)
# bar is copy-only across both arms -- an arm-specific correlated field is not
# defined (the cross-field Sigma lives within one arm, so without a copy there is
# no cross-arm transfer to estimate). The independent (`||`) single-arm bar is a
# separate per-arm latent (gcol33/tulpaObs#65), routed past this check in
# `.encode_cover_terms`.
.cover_bar_check_to <- function(spec) {
  to  <- spec$to %||% .tobs_cover_arms
  bad <- setdiff(to, .tobs_cover_arms)
  if (length(bad) > 0L || length(to) < 1L) {
    stop(sprintf(paste0(
      "internal: a cover bar field has placement `to` = %s; expected ",
      "\"presence\", \"positive\", or both."),
      paste0("c(", paste0("\"", to, "\"", collapse = ", "), ")")),
      call. = FALSE)
  }
  # Both cover arms (the copy-only correlated field, gcol33/tulpaObs#64) OR a
  # single arm (a free-Sigma correlated field on that arm alone, no cross-arm
  # copy, gcol33/tulpaObs#109) are both supported.
  invisible(spec)
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
    to              = spec$to %||% .tobs_cover_arms,
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
# Split a per-arm cover formula into its fixed-effects formula and its spatial-
# field terms (tagged with the arm via `to =`, so a field placed in `presence` /
# `positive` lands on that arm through the same machinery a shared `to =` bar
# uses). lme4 bars are desugared first so `(1 | g)` reads as re(); temporal() /
# re() / other structured terms are not routed per-arm yet (declare them on the
# shared `formula`).
.cover_lift_arm_fields <- function(formula, arm) {
  if (is.null(formula)) return(list(fe = NULL, fields = list()))
  tt   <- stats::terms(.tobs_desugar_bars(formula), keep.order = TRUE)
  labs <- attr(tt, "term.labels")
  field_ctors <- c("spatial", "icar", "bym2", "car", "car_proper")
  keep <- character(0); fields <- list(); copies <- list()
  for (lab in labs) {
    e    <- tryCatch(str2lang(lab), error = function(...) NULL)
    head <- if (is.call(e) && is.symbol(e[[1L]])) as.character(e[[1L]]) else NA_character_
    if (!is.na(head) && head %in% field_ctors) {
      # The arm is fixed by placement (this formula's arm); the field call is
      # kept unevaluated and tagged with the arm later, on its evaluated spec.
      fields[[length(fields) + 1L]] <- e
    } else if (!is.na(head) && identical(head, "copy")) {
      # copy() is not a fixed effect and not a placed field: it references the
      # presence field and couples it onto this arm. Collect the call (evaluated
      # into a tobs_copy spec by the caller) and keep it out of the
      # fixed-effects formula.
      copies[[length(copies) + 1L]] <- e
    } else if (!is.na(head) && head %in% .tobs_term_names()) {
      stop(sprintf(paste0(
        "cover(): the per-arm `%s` formula supports fixed effects and spatial() ",
        "fields; `%s()` is not routed per-arm. Declare it on the shared `formula`."),
        arm, head), call. = FALSE)
    } else {
      keep <- c(keep, lab)
    }
  }
  fe <- stats::reformulate(if (length(keep)) keep else "1",
                           intercept = as.logical(attr(tt, "intercept")))
  environment(fe) <- environment(formula)
  list(fe = fe, fields = fields, copies = copies)
}

# Promote the presence-arm spatial field(s) that a copy() selects to the shared
# (both-arm) arm tag, so encode routes them through the presence-anchored,
# positive-coupled machinery. copy() is the canonical shared-field spelling for
# per-arm cover formulas: the field is placed on the `presence` formula and
# copy(spatial()) in the `positive` formula couples it across, mirroring
# occu_cover(). Returns the presence field calls, the both-arm tag to set on
# their specs, and the coupling amplitude grid (NULL = the fitter's default,
# estimated on the standard alpha grid).
.cover_promote_copied_fields <- function(copies, occ_fields) {
  if (length(copies) > 1L) {
    stop("cover(): one copy() per fit is supported (it couples the whole ",
         "presence field onto the positive arm). Per-component coupling is not ",
         "yet wired for cover().", call. = FALSE)
  }
  cp <- copies[[1L]]
  if (is.null(cp$selector_type)) {
    stop("cover(): copy() must select the presence spatial field, e.g. ",
         "copy(spatial()); a field-name string is not a coupling selector here.",
         call. = FALSE)
  }
  if (length(occ_fields) == 0L) {
    stop("cover(): copy() needs a spatial field on the `presence` formula to ",
         "copy onto the positive arm, e.g. ",
         "presence = ~ x + spatial(~ 1 || cell, graph = adj).", call. = FALSE)
  }
  # Per-component coupling: copy(terms = list(<component> = ...)) couples each
  # field block (intercept, weighted trend) with its own amplitude grid, so the
  # two can decouple. Resolved here into a named list of grids (NULL = the
  # fitter's default grid); the trend-block-existence check runs downstream in
  # .dispatch_cover, where the encoded field is known. The whole-field form
  # (no terms) keeps its single amplitude, applied to both blocks.
  alpha_terms <- if (!is.null(cp$copy_terms)) {
    lapply(cp$copy_terms, function(res)
      if (isTRUE(is.na(res$integrate))) NULL else res$grid)
  } else NULL
  alpha_grid <- if (!is.null(cp$copy_terms)) NULL
                else if (is.na(cp$alpha_integrate)) NULL
                else cp$alpha_grid
  list(fields = occ_fields, arm = c("presence", "positive"),
       alpha = alpha_grid, alpha_terms = alpha_terms)
}

# Map a copy(terms = list(...)) per-component amplitude spec onto the fitter's
# per-block coupling grids: the presence field's intercept block reads
# control$alpha.grid, the weighted-trend block reads control$alpha.grid.trend
# (see the joint-coupled cover fitter: base_block = block 1, trend_block =
# block 2). `trend` is the resolved weighted-trend field (NULL when the presence
# field is a single intercept). Every existing block must be named, so no block
# is silently left at the default -- mirroring .occu_cover_apply_copy_coupling().
.cover_apply_copy_terms <- function(copy_terms, trend, control) {
  has_trend   <- !is.null(trend)
  trend_label <- if (has_trend) trend$label else NULL

  grids <- list()   # canonical role ("intercept" / "trend") -> amplitude grid
  for (k in names(copy_terms)) {
    role <- if (identical(k, "intercept")) {
      "intercept"
    } else if (identical(k, "trend") ||
               (!is.null(trend_label) && identical(k, trend_label))) {
      "trend"
    } else {
      NA_character_
    }
    if (is.na(role)) {
      avail <- if (has_trend) {
        sprintf("\"intercept\", \"trend\" (or \"%s\")", trend_label)
      } else "\"intercept\""
      stop(sprintf(paste0(
        "cover(): copy(terms = ): unknown field component \"%s\". ",
        "Name %s."), k, avail), call. = FALSE)
    }
    if (identical(role, "trend") && !has_trend) {
      stop("cover(): copy(terms = ) names a \"trend\" component, but the ",
           "presence field has no weighted trend block (it is a single ",
           "intercept field). Drop the trend component or add a weighted areal ",
           "term, e.g. spatial(~ 1 + time.sc || cell, graph = adj).",
           call. = FALSE)
    }
    grids[[role]] <- copy_terms[[k]]
  }

  required <- c("intercept", if (has_trend) "trend")
  missing_roles <- setdiff(required, names(grids))
  if (length(missing_roles)) {
    stop(sprintf(paste0(
      "cover(): copy(terms = ) must give an amplitude for every field block; ",
      "%s left unaddressed. Field blocks: %s."),
      paste0("\"", missing_roles, "\"", collapse = ", "),
      paste0("\"", required, "\"", collapse = ", ")), call. = FALSE)
  }

  # NULL grid = the fitter's default (assigning NULL drops the control element).
  control$alpha.grid <- if (is.null(grids$intercept)) NULL
                        else as.numeric(grids$intercept)
  if (has_trend) {
    control$alpha.grid.trend <- if (is.null(grids$trend)) NULL
                                else as.numeric(grids$trend)
  }
  control
}

# Parse a cover() formula against the NA-dropped observations.
.encode_cover_terms <- function(formula, data_obs) {
  bind   <- .tobs_bind_formulas(list(state = formula), data_obs)
  specs  <- lapply(bind$terms, function(t) t$spec)
  enc    <- .encode_cover_specs(specs, data_obs, re_guard_formula = formula)
  enc$fe <- bind$fe$state
  enc
}

# Route a list of parsed cover structured-term specs into the encoded field
# blocks. Split out from .encode_cover_terms() so the per-arm placement path can
# feed specs it has already evaluated and tagged with `spec$to` (the arm),
# instead of round-tripping them through a deparsed formula. `re_guard_formula`
# is the source formula for the bare-bar RE soft guard; the placement path passes
# NULL (its per-arm formulas reject re() upstream, so there is nothing to guard).
.encode_cover_specs <- function(specs, data_obs, re_guard_formula = NULL) {
  spatial_specs <- list(); temporal <- NULL; re <- list(); mcar <- NULL
  armspec <- list()
  for (spec in specs) {
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
  if (!is.null(re_guard_formula)) {
    .tobs_cover_bar_re_guard(re_guard_formula, spatial_specs)
  }
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
           "(got two on the same arm). Combine the coefficient fields into one ",
           "bar in that arm's formula, e.g. ",
           "positive = ~ x + spatial(~ 1 + w || cell, graph = adj).",
           call. = FALSE)
    }
  }
  fields <- .cover_resolve_spatial_fields(spatial_specs, data_obs)
  list(fe = NULL, spatial = fields$spatial, trend = fields$trend,
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
#' Two independent `tulpa::tulpa_laplace()` calls. For
#' `positive = "lognormal"` the positive arm is a Gaussian fit on
#' `log(cover)` with sigma estimated post-hoc as the residual standard
#' error. For `positive = "beta"` the positive arm uses
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
    stop("cover() supports method = 'laplace'/'laplace_sla', ",
         "'nested_laplace'/'nested_laplace_sla' and 'nuts' (got engine '",
         engine, "'). This fitter takes the two-arm laplace path only.",
         call. = FALSE)
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

  if (positive %in% c("lognormal", "gaussian")) {
    # Both are location-Gaussian arms with dispersion sigma_pos; they differ only
    # in whether the response was log-transformed in encode_cover_hurdle()
    # (lognormal on log-cover, gaussian on the raw response, #112). From here the
    # fit machinery is identical.
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
      positive  = positive,
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
    positive  = positive,           # "beta" or "beta_oi" (interior Beta)
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
  pos_vcov_scale <- if (fits$positive %in% c("lognormal", "gaussian"))
                      fits$sigma_pos^2 else 1
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
  if (fits$positive %in% c("lognormal", "gaussian")) {
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
      sigma_pos    = if (fits$positive %in% c("lognormal", "gaussian"))
                       fits$sigma_pos else NA_real_,
      sigma_pos_sd = NA_real_,
      phi_pos      = if (fits$positive %in% c("beta", "beta_oi")) fits$phi_pos else NA_real_,
      phi_pos_sd   = NA_real_,
      pi_one       = enc$oi$pi_one    %||% NA_real_,
      pi_one_sd    = enc$oi$pi_one_sd %||% NA_real_,
      n_ceiling    = enc$oi$n_ceiling %||% NA_integer_,
      hyperpar     = hyperpar,
      encoding     = enc,
      family       = family,
      n_total      = enc$N,
      n_positive   = enc$oi$n_positive %||% length(enc$idx_pos),
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
    class = c("cover_fit", "tobs_multiarm_fit", "tobs_fit", "tulpa_fit")
  )
  out
}


