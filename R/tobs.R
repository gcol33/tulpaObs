# =============================================================================
# tobs.R — Unified entry point for tulpaObs latent-state observation models
#
# The public API. Takes a state-process formula, data, a family object, and
# (depending on family) detection formula + response. Dispatches to the
# internal engine path for the family's `name` slug.
# =============================================================================


#' Fit a hierarchical latent-state observation model
#'
#' Unified entry point for occupancy, abundance, distance, removal, and related
#' models that share the latent-state-plus-imperfect-detection generative
#' template. The specific model is chosen via the `family` argument; the engine
#' (Laplace / nested Laplace / NUTS) via `engine`.
#'
#' @param formula state-process formula, e.g. `~ elev + forest`. For occupancy
#'   this is the occupancy probability formula; for N-mixture the abundance
#'   formula; for the cover hurdle the latent-presence formula.
#' @param data data frame of site-level covariates with `nrow(data) ==
#'   nrow(y)`.
#' @param family a `tobs_family` object (see [obs_family()] and the concrete
#'   constructors [occu()], [abun()], [cover()], ...).
#' @param detection detection-process formula, e.g. `~ observer + effort`.
#'   Family-dependent: required for [occu()] and [abun()], ignored for
#'   [jsdm()] and (currently) [cover()].
#' @param y response. Shape depends on family:
#'   * [occu()] — N x J detection-history matrix.
#'   * [abun()] — N x J integer count matrix.
#'   * `ms_*` — S x N x J array.
#'   * [cover()] — length-N vector of cover proportions in \[0, 1\].
#' @param visit_data optional visit-level detection covariates. Accepts
#'   either:
#'   * a named list of `[n_sites, max_visits]` matrices (the shape returned by
#'     `tobs_data()` in `det.covs`) — flattened internally to a long data
#'     frame in site-major order;
#'   * a data frame with `nrow(y) * ncol(y)` rows in site-major order.
#'
#'   When `visit_data` is provided without a `"formula"` attribute, the
#'   `detection` argument is interpreted as the visit-level detection
#'   formula and the site-level detection design matrix is an intercept
#'   only. To split visit-level and site-level detection covariates
#'   (e.g. visit-level effort plus site-level observer category), pass a
#'   long data frame with `attr(visit_data, "formula") = ~ effort` and use
#'   `detection = ~ observer` for the site-level terms.
#' @param spatial optional spatial spec from `tulpa` or `tulpaMesh` (e.g.
#'   [tulpa::spatial_bym2()], [tulpa::spatial_nngp()], or a [tulpaMesh::mesh]
#'   object), or one of the convenience helpers.
#' @param temporal optional temporal spec (e.g. [tulpa::temporal_ar1()]).
#' @param re optional random-effects spec from [tobs_re()] or
#'   [tobs_community_re()], or a list of such specs. Can be combined with
#'   `spatial` and `temporal` in the same fit. Under `engine = "nested_laplace"`
#'   each `re` term with `model = "iid"` becomes an IID latent block in the
#'   multi-block prior; for other engines `re` is forwarded to the underlying
#'   fitter as the per-observation RE structure.
#' @param engine inference engine: `"auto"`, `"laplace"`, `"nested_laplace"`,
#'   or `"nuts"`. `"auto"` chooses the family's `default_engine`.
#' @param priors optional prior specification. For occupancy families fit
#'   with `engine = "laplace"`, pass a list or [occu_priors()] object to
#'   set weakly-informative quadratic priors on the fixed-effect
#'   coefficients (defaults pull the detection intercept toward
#'   `p = 0.5` and break the psi-p identifiability ridge at small `J`).
#'   Pass `priors = FALSE` to disable the default prior and recover the
#'   historical unpenalised MAP behavior. For NUTS, this is forwarded
#'   to the underlying tulpa engine.
#' @param control list of low-level controls (`n_threads`, `max_iter`, `tol`,
#'   etc.).
#' @param ... family-specific named arguments forwarded to the underlying
#'   engine builder.
#'
#' @return An object of class `c("tobs_fit", "<family>_fit", "tulpa_fit")`.
#'
#' @examples
#' \dontrun{
#' # Single-season occupancy
#' fit <- tobs(
#'   formula   = ~ elev,
#'   data      = sites,
#'   family    = occu(),
#'   detection = ~ effort,
#'   y         = y_matrix
#' )
#' }
#'
#' @export
tobs <- function(formula,
                 data,
                 family,
                 detection  = NULL,
                 y          = NULL,
                 visit_data = NULL,
                 spatial    = NULL,
                 temporal   = NULL,
                 re         = NULL,
                 engine     = c("auto", "laplace", "nested_laplace", "nuts"),
                 approx     = c("gaussian_laplace", "simplified_laplace"),
                 priors     = NULL,
                 control    = list(),
                 ...) {

  approx <- match.arg(approx)

  if (missing(family)) {
    stop(
      "`family` is required. Use one of: occu(), dyn_occu(), ms_occu(), ",
      "int_occu(), jsdm(), abun(), cover(), ... See `?obs_family` for the ",
      "full list.",
      call. = FALSE
    )
  }
  if (!inherits(family, "tobs_family")) {
    stop(
      "`family` must be a tobs_family object (e.g. occu() or abun()), ",
      "not a ", paste(class(family), collapse = "/"), ".",
      call. = FALSE
    )
  }

  re <- .normalize_re_arg(re)

  engine <- match.arg(engine)
  if (engine == "auto") engine <- family$default_engine

  if (family$status == "planned") {
    .stop_planned_family(family)
  }

  dispatch <- switch(
    family$name,
    occu     = .dispatch_occu,
    dyn_occu = .dispatch_dyn_occu,
    ms_occu  = .dispatch_ms_occu,
    int_occu = .dispatch_int_occu,
    jsdm     = .dispatch_jsdm,
    cover    = .dispatch_cover,
    stop(sprintf(
      "Internal error: family %q has status 'working' but no dispatcher.",
      family$name
    ), call. = FALSE)
  )

  fit <- dispatch(
    formula    = formula,
    data       = data,
    family     = family,
    detection  = detection,
    y          = y,
    visit_data = visit_data,
    spatial    = spatial,
    temporal   = temporal,
    re         = re,
    engine     = engine,
    approx     = approx,
    priors     = priors,
    control    = control,
    ...
  )

  if (!inherits(fit, "tobs_fit")) {
    class(fit) <- c("tobs_fit", class(fit))
  }
  attr(fit, "tobs_family") <- family
  fit
}


# ---------------------------------------------------------------------------
# Internal dispatchers — one per working family. Each wraps the internal
# model builder + fitter for that family's structure.
# ---------------------------------------------------------------------------

.dispatch_occu <- function(formula, data, family, detection, y, visit_data,
                           spatial, temporal, re, engine, priors, control,
                           approx = "gaussian_laplace", ...) {
  if (is.null(detection)) {
    stop("occu() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu() requires `y` (an N x J detection-history matrix).", call. = FALSE)
  }
  vd <- .normalize_visit_data(visit_data, detection, n_sites = nrow(y),
                              max_visits = ncol(y))
  model <- .tobs_build_model(
    occ_formula        = formula,
    det_formula        = vd$det_formula,
    data               = data,
    y                  = y,
    det_visit_formula  = vd$det_visit_formula,
    det_visit_data     = vd$visit_data
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal, re = re,
         method = .map_engine(engine, family = "occu"), priors = priors,
         approx = approx),
    control
  ))
}

.dispatch_dyn_occu <- function(formula, data, family, detection, y, visit_data,
                               spatial, temporal, re, engine, priors, control,
                               approx = "gaussian_laplace", ...) {
  dots <- list(...)
  if (is.null(dots$col_formula)) {
    stop("dyn_occu() requires a `col_formula = ~ ...` argument.", call. = FALSE)
  }
  if (is.null(dots$ext_formula)) {
    stop("dyn_occu() requires an `ext_formula = ~ ...` argument.", call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula  = formula,
    det_formula  = detection,
    data         = data,
    y            = y,
    col_formula  = dots$col_formula,
    ext_formula  = dots$ext_formula
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal, re = re,
         method = .map_engine(engine, family = "dyn_occu"), priors = priors,
         approx = approx),
    control
  ))
}

.dispatch_ms_occu <- function(formula, data, family, detection, y, visit_data,
                              spatial, temporal, re, engine, priors, control,
                              approx = "gaussian_laplace", ...) {
  dots <- list(...)
  if (is.null(dots$species)) {
    stop("ms_occu() requires a `species` argument.", call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    species     = dots$species
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal, re = re,
         method = .map_engine(engine, family = "ms_occu"), priors = priors,
         approx = approx),
    control
  ))
}

.dispatch_int_occu <- function(formula, data, family, detection, y, visit_data,
                               spatial, temporal, re, engine, priors, control,
                               approx = "gaussian_laplace", ...) {
  model <- .tobs_build_model(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    integrated  = TRUE
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal, re = re,
         method = .map_engine(engine, family = "int_occu"), priors = priors,
         approx = approx),
    control
  ))
}

.dispatch_jsdm <- function(formula, data, family, detection, y, visit_data,
                           spatial, temporal, re, engine, priors, control,
                           approx = "gaussian_laplace", ...) {
  dots <- list(...)
  model <- .tobs_build_model(
    occ_formula = formula,
    data        = data,
    y           = y,
    jsdm        = TRUE,
    species     = dots$species
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal, re = re,
         method = .map_engine(engine, family = "jsdm"), priors = priors,
         approx = approx),
    control
  ))
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.map_engine <- function(engine, family = NULL) {
  # Engine name translation between the tobs vocabulary and what the underlying
  # fitter currently understands. Single-season occupancy (`family = "occu"`)
  # has a real nested-Laplace path that routes through tulpa's multi-block
  # nested-Laplace engine via `.tobs_em_nested_laplace()`; other families
  # still fall back to single-Laplace with a NOTE because their dispatch
  # hasn't been wired through yet.
  if (engine == "nested_laplace") {
    if (identical(family, "occu")) return("nested_laplace")
    message(
      "tobs(): nested_laplace is currently wired only for single-season ",
      "occupancy (`family = occu()`); falling back to single-Laplace for ",
      "family '", family %||% "(unspecified)", "'."
    )
    return("laplace")
  }
  switch(
    engine,
    laplace = "laplace",
    nuts    = "nuts",
    engine
  )
}

# Coerce the user's `re` argument into a list-of-tobs_re or NULL. Accepts:
#   * NULL                    -> NULL
#   * a single tobs_re object -> list(re)
#   * a list of tobs_re       -> as-is (validated element-wise)
.normalize_re_arg <- function(re) {
  if (is.null(re)) return(NULL)
  if (inherits(re, "tobs_re")) return(list(re))
  if (is.list(re)) {
    if (length(re) == 0L) return(NULL)
    bad <- !vapply(re, inherits, logical(1), what = "tobs_re")
    if (any(bad)) {
      stop("`re` must be a tobs_re object or a list of tobs_re objects; ",
           "element(s) ", paste(which(bad), collapse = ", "),
           " are not tobs_re.", call. = FALSE)
    }
    return(re)
  }
  stop("`re` must be a tobs_re object, a list of tobs_re, or NULL; got ",
       paste(class(re), collapse = "/"), ".", call. = FALSE)
}

# Normalize the user's `visit_data` argument to the shape `.tobs_build_single`
# expects: a long data frame with `n_sites * max_visits` rows in site-major
# order (row `r` corresponds to site `(r-1) %/% max_visits + 1`, visit
# `(r-1) %% max_visits + 1`), plus the formula to apply to it.
#
# Accepts:
#   * NULL                                  -> NULL, no visit-level path
#   * named list of [n_sites, max_visits]   -> flatten to long DF, treat
#       matrices (the `tobs_data()` shape)     `detection` as visit-level
#                                              (intercept dropped, site-level
#                                              X_det is intercept-only)
#   * data.frame with N*J rows + a          -> existing dual-formula behavior:
#       "formula" attribute                    site-level `detection` against
#                                              `data`, attr formula against
#                                              `visit_data`
#   * data.frame with N*J rows, no formula  -> treat `detection` as visit-level
#       attribute                              (intercept dropped); site-level
#                                              X_det is intercept-only
#
# Returns a list with:
#   visit_data         — long data frame (or NULL)
#   det_visit_formula  — formula applied to visit_data (or NULL)
#   det_formula        — formula applied to site-level `data`
.normalize_visit_data <- function(visit_data, detection,
                                  n_sites, max_visits) {
  if (is.null(visit_data)) {
    return(list(visit_data = NULL,
                det_visit_formula = NULL,
                det_formula = detection))
  }

  expected_rows <- n_sites * max_visits

  # Case 1: list of [n_sites, max_visits] matrices (tobs_data() output)
  if (is.list(visit_data) && !is.data.frame(visit_data)) {
    nms <- names(visit_data)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("`visit_data` (list of matrices) must be a named list; ",
           "names become the column names of the flattened frame.",
           call. = FALSE)
    }
    bad <- vapply(visit_data, function(m) {
      !is.matrix(m) || nrow(m) != n_sites || ncol(m) != max_visits
    }, logical(1))
    if (any(bad)) {
      stop(sprintf(
        "`visit_data` elements must be [%d x %d] matrices matching y; ",
        n_sites, max_visits),
        sprintf("element(s) %s have wrong shape.",
                paste(nms[bad], collapse = ", ")),
        call. = FALSE)
    }
    flat <- as.data.frame(
      lapply(visit_data, function(m) as.vector(t(m))),
      stringsAsFactors = FALSE
    )
    return(list(
      visit_data = flat,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  # Case 2 / 3: data frame
  if (is.data.frame(visit_data)) {
    if (nrow(visit_data) != expected_rows) {
      stop(sprintf(
        "`visit_data` (data frame) must have %d rows (n_sites * max_visits); ",
        expected_rows),
        sprintf("got %d.", nrow(visit_data)),
        call. = FALSE)
    }
    attached <- attr(visit_data, "formula")
    if (!is.null(attached)) {
      # Dual-formula power-user mode: detection stays site-level
      return(list(
        visit_data = visit_data,
        det_visit_formula = attached,
        det_formula = detection
      ))
    }
    return(list(
      visit_data = visit_data,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  stop("`visit_data` must be NULL, a named list of [n_sites x max_visits] ",
       "matrices, or a long data frame with n_sites * max_visits rows; ",
       "got ", paste(class(visit_data), collapse = "/"), ".",
       call. = FALSE)
}

# Drop the intercept term from a formula (returns `~ . - 1`-style update).
# Preserves the LHS if any (none of our detection formulas have one).
.drop_intercept <- function(f) {
  stats::update(f, ~ . - 1)
}

.stop_planned_family <- function(family) {
  phase <- switch(
    family$name,
    cover    = "Phase 1 (cover hurdle, beta variant)",
    abun     = "Phase 2 (N-mixture)",
    ms_abun  = "Phase 2 (N-mixture)",
    dyn_abun = "Phase 3 (open N-mixture)",
    distance = "Phase 4 (distance / removal / FP)",
    removal  = "Phase 4 (distance / removal / FP)",
    fp_occu  = "Phase 4 (distance / removal / FP)",
    "a future phase"
  )
  stop(
    sprintf(
      "Family `%s` (%s) is planned but not yet implemented. Scheduled: %s. ",
      family$name, family$class_long, phase
    ),
    "See PLAN_tulpaObs.md for the rollout.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# Print method for fits
# ---------------------------------------------------------------------------

#' Print method for tobs_fit
#' @param x a `tobs_fit` object.
#' @param ... forwarded to underlying print methods.
#' @return `x`, invisibly.
#' @export
print.tobs_fit <- function(x, ...) {
  fam <- attr(x, "tobs_family")
  if (!is.null(fam)) {
    cat(sprintf("<tobs_fit: %s>\n", fam$class_long))
    cat(sprintf("  family         : %s (status: %s)\n", fam$name, fam$status))
    cat(sprintf("  default engine : %s\n", fam$default_engine))
    cat("\n")
  }
  model <- x$model
  if (!is.null(model)) {
    if (model$model_type == "single") {
      cat(sprintf("  Sites: %d, Max visits: %d\n", model$n_sites, model$max_visits))
    } else if (model$model_type == "dynamic") {
      cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
                  model$n_sites, model$n_seasons, model$max_visits))
    } else if (model$model_type == "community") {
      cat(sprintf("  Sites: %d, Species: %d\n", model$n_sites, model$n_species))
    } else if (model$model_type == "integrated") {
      cat(sprintf("  Sites: %d, Sources: %d\n", model$n_sites, model$n_sources))
    } else if (model$model_type == "jsdm") {
      cat(sprintf("  Sites: %d, Species: %d\n", model$n_sites, model$n_species))
    }
  }
  if (!is.null(x$n_samples)) {
    cat(sprintf("  Samples: %d", x$n_samples))
    if (!is.null(x$epsilon) && !is.na(x$epsilon)) {
      cat(sprintf(", step size: %.4f", x$epsilon))
    }
    cat("\n")
  }
  if (!is.null(x$sla_status) && !identical(x$sla_status, "off")) {
    cat(sprintf("  Marginals: %s\n", x$sla_status))
    clipped <- attr(x$draws, "sla_clipped")
    fallback <- attr(x$draws, "sla_fallback")
    if (length(clipped) > 0) {
      cat(sprintf("    skew clipped to +/-0.95 for: %s\n",
                  paste(clipped, collapse = ", ")))
    }
    if (length(fallback) > 0) {
      cat(sprintf("    fell back to Gaussian for: %s\n",
                  paste(fallback, collapse = ", ")))
    }
  }
  if (!is.null(x$divergent) && sum(x$divergent) > 0) {
    cat(sprintf("  WARNING: %d divergent transitions\n", sum(x$divergent)))
  }
  if (!is.null(x$intercepts)) {
    cat("\n")
    for (nm in names(x$intercepts)) {
      label <- switch(nm,
        psi  = "Mean occupancy (intercept)",
        psi1 = "Mean initial occupancy (intercept)",
        p    = "Mean detection (intercept)",
        gamma   = "Mean colonization (intercept)",
        epsilon = "Mean extinction (intercept)"
      )
      if (!is.null(label)) {
        cat(sprintf("%s: %.3f\n", label, x$intercepts[[nm]]))
      }
    }
  }
  invisible(x)
}
