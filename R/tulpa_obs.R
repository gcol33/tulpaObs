# =============================================================================
# tulpa_obs.R — Unified entry point for TulpaObs latent-state observation models
#
# This is the new public API. During the Phase 0 transition it is a thin
# dispatcher that calls the existing occu() / occu_fit() pipeline for the
# `occ()` family and raises informative errors for not-yet-implemented
# families. The existing entry points (occu, occu_fit) keep working
# unchanged.
# =============================================================================


#' Fit a hierarchical latent-state observation model
#'
#' Unified entry point for occupancy, abundance, distance, removal, and
#' related models that share the latent-state-plus-imperfect-detection
#' generative template. The specific model is chosen via the `family`
#' argument; the engine (Laplace / nested Laplace / NUTS) via `engine`.
#'
#' This function is the canonical TulpaObs API. The package-specific
#' entry points (`occu()`, `occu_fit()`) remain exported for backward
#' compatibility and are called internally during the Phase 0 transition.
#'
#' @section Phase 0 status:
#'
#' Only `family = occ()` is wired to a real engine in v0.1. Other families
#' raise informative errors pointing to the planned phase in
#' `PLAN_TulpaObs.md`.
#'
#' @param formula state-process formula, e.g. `~ elev + forest`. For
#'   occupancy this is the occupancy probability formula; for N-mixture
#'   the abundance formula; for cover hurdle the latent-presence formula.
#' @param data data frame of site-level covariates with `nrow(data) ==
#'   nrow(y)`.
#' @param family a `tulpa_obs_family` object (see [obs_family()] and the
#'   concrete constructors [occ()], [nmixture()], [cover_hurdle()], ...).
#' @param detection detection-process formula, e.g. `~ observer + effort`.
#'   Family-dependent: required for `occ()` and `nmixture()`, ignored for
#'   `jsdm()` and (currently) `cover_hurdle()`.
#' @param y response. Shape depends on family:
#'   * `occ()` — N x J detection-history matrix.
#'   * `nmixture()` — N x J integer count matrix.
#'   * `multispecies_*` — S x N x J array.
#'   * `cover_hurdle()` — length-N vector of cover proportions in \[0, 1\].
#' @param visit_data optional data frame of visit-level detection
#'   covariates with `nrow(visit_data) == nrow(y) * ncol(y)`.
#' @param spatial optional spatial spec from `tulpa` or `tulpaMesh`
#'   (e.g. [tulpa::spatial_bym2()], [tulpa::spatial_nngp()], or a
#'   [tulpaMesh::mesh] object).
#' @param temporal optional temporal spec (e.g. [tulpa::temporal_ar1()]).
#' @param engine inference engine: `"auto"`, `"laplace"`,
#'   `"nested_laplace"`, or `"nuts"`. `"auto"` chooses the family's
#'   `default_engine`.
#' @param priors optional [tulpa::tulpa_priors()] object.
#' @param control list of low-level controls (`n_threads`, `max_iter`,
#'   `tol`, etc.).
#' @param ... family-specific named arguments forwarded to the underlying
#'   engine builder.
#'
#' @return An object of class `c("tulpa_obs_fit", "<family>_fit", "tulpa_fit")`.
#'
#' @examples
#' \dontrun{
#' # Single-season occupancy (currently dispatches to occu())
#' fit <- tulpa_obs(
#'   formula   = ~ elev,
#'   data      = sites,
#'   family    = occ(),
#'   detection = ~ effort,
#'   y         = y_matrix
#' )
#' }
#'
#' @export
tulpa_obs <- function(formula,
                      data,
                      family,
                      detection  = NULL,
                      y          = NULL,
                      visit_data = NULL,
                      spatial    = NULL,
                      temporal   = NULL,
                      engine     = c("auto", "laplace", "nested_laplace", "nuts"),
                      priors     = NULL,
                      control    = list(),
                      ...) {

  if (missing(family)) {
    stop(
      "`family` is required. Use one of: occ(), dynamic_occ(), ",
      "multispecies_occ(), integrated_occ(), jsdm(), nmixture(), ",
      "cover_hurdle(), ... See `?obs_family` for the full list.",
      call. = FALSE
    )
  }
  if (!inherits(family, "tulpa_obs_family")) {
    stop(
      "`family` must be a tulpa_obs_family object (e.g. occ() or nmixture()), ",
      "not a ", paste(class(family), collapse = "/"), ".",
      call. = FALSE
    )
  }

  engine <- match.arg(engine)
  if (engine == "auto") engine <- family$default_engine

  if (family$status == "planned") {
    .stop_planned_family(family)
  }

  dispatch <- switch(
    family$name,
    occ              = .dispatch_occ,
    dynamic_occ      = .dispatch_dynamic_occ,
    multispecies_occ = .dispatch_multispecies_occ,
    integrated_occ   = .dispatch_integrated_occ,
    jsdm             = .dispatch_jsdm,
    cover_hurdle     = .dispatch_cover_hurdle,
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
    engine     = engine,
    priors     = priors,
    control    = control,
    ...
  )

  if (!inherits(fit, "tulpa_obs_fit")) {
    class(fit) <- c("tulpa_obs_fit", class(fit))
  }
  attr(fit, "tulpa_obs_family") <- family
  fit
}


# ---------------------------------------------------------------------------
# Internal dispatchers — one per working family.
#
# Each thin-wraps the existing tulpaOcc builder during the Phase 0
# transition. Phase 1+ will replace these with direct engine calls.
# ---------------------------------------------------------------------------

.dispatch_occ <- function(formula, data, family, detection, y, visit_data,
                          spatial, temporal, engine, priors, control, ...) {
  if (is.null(detection)) {
    stop("occ() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occ() requires `y` (an N x J detection-history matrix).", call. = FALSE)
  }
  if (!is.null(spatial) || !is.null(temporal)) {
    stop(
      "spatial / temporal arguments via tulpa_obs() are not yet wired for ",
      "family = occ(). Use occu() with the existing spatial helpers ",
      "(occu_icar, occu_bym2, occu_gp, occu_spde) and occu_temporal() for now; ",
      "the unified hook lands in Phase 1.",
      call. = FALSE
    )
  }
  model <- occu(
    occ_formula        = formula,
    det_formula        = detection,
    data               = data,
    y                  = y,
    det_visit_formula  = if (!is.null(visit_data)) attr(visit_data, "formula") else NULL,
    det_visit_data     = visit_data
  )
  do.call(occu_fit, c(
    list(model = model, engine = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_dynamic_occ <- function(formula, data, family, detection, y, visit_data,
                                  spatial, temporal, engine, priors, control, ...) {
  dots <- list(...)
  if (is.null(dots$col_formula)) {
    stop("dynamic_occ() requires a `col_formula = ~ ...` argument.", call. = FALSE)
  }
  if (is.null(dots$ext_formula)) {
    stop("dynamic_occ() requires an `ext_formula = ~ ...` argument.", call. = FALSE)
  }
  model <- occu(
    occ_formula  = formula,
    det_formula  = detection,
    data         = data,
    y            = y,
    col_formula  = dots$col_formula,
    ext_formula  = dots$ext_formula
  )
  do.call(occu_fit, c(
    list(model = model, engine = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_multispecies_occ <- function(formula, data, family, detection, y, visit_data,
                                       spatial, temporal, engine, priors, control, ...) {
  dots <- list(...)
  if (is.null(dots$species)) {
    stop("multispecies_occ() requires a `species` argument.", call. = FALSE)
  }
  model <- occu(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    species     = dots$species
  )
  do.call(occu_fit, c(
    list(model = model, engine = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_integrated_occ <- function(formula, data, family, detection, y, visit_data,
                                     spatial, temporal, engine, priors, control, ...) {
  model <- occu(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    integrated  = TRUE
  )
  do.call(occu_fit, c(
    list(model = model, engine = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_jsdm <- function(formula, data, family, detection, y, visit_data,
                           spatial, temporal, engine, priors, control, ...) {
  dots <- list(...)
  model <- occu(
    occ_formula = formula,
    data        = data,
    y           = y,
    jsdm        = TRUE,
    species     = dots$species
  )
  do.call(occu_fit, c(
    list(model = model, engine = .map_engine(engine), priors = priors),
    control
  ))
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.map_engine <- function(engine) {
  # Engine name translation between the TulpaObs vocabulary and what the
  # underlying occu_fit() currently understands. Phase 0 maps both
  # "laplace" and "nested_laplace" to occu_fit's "laplace" mode and
  # surfaces a NOTE when the user asked for nested_laplace specifically.
  switch(
    engine,
    laplace        = "laplace",
    nested_laplace = {
      message(
        "tulpa_obs(): nested_laplace requested but Phase 0 dispatches to ",
        "the single-Laplace engine for family = occ(). Nested Laplace ",
        "hookup lands in Phase 1."
      )
      "laplace"
    },
    nuts           = "nuts",
    engine
  )
}

.stop_planned_family <- function(family) {
  phase <- switch(
    family$name,
    cover_hurdle      = "Phase 1 (cover hurdle)",
    nmixture          = "Phase 2 (N-mixture)",
    multispecies_nmix = "Phase 2 (N-mixture)",
    dynamic_nmix      = "Phase 3 (open N-mixture)",
    distance          = "Phase 4 (distance / removal / FP)",
    removal           = "Phase 4 (distance / removal / FP)",
    false_positive    = "Phase 4 (distance / removal / FP)",
    "a future phase"
  )
  stop(
    sprintf(
      "Family `%s` (%s) is planned but not yet implemented. Scheduled: %s. ",
      family$name, family$class_long, phase
    ),
    "See PLAN_TulpaObs.md for the rollout. ",
    "In the meantime, use the existing entry points where they exist ",
    "(occu() for occupancy variants).",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# Print method for fits
# ---------------------------------------------------------------------------

#' Print method for tulpa_obs_fit
#' @param x a `tulpa_obs_fit` object.
#' @param ... forwarded to underlying print methods.
#' @return `x`, invisibly.
#' @export
print.tulpa_obs_fit <- function(x, ...) {
  fam <- attr(x, "tulpa_obs_family")
  if (!is.null(fam)) {
    cat(sprintf("<tulpa_obs_fit: %s>\n", fam$class_long))
    cat(sprintf("  family         : %s (status: %s)\n", fam$name, fam$status))
    cat(sprintf("  default engine : %s\n", fam$default_engine))
    cat("\n")
  }
  # Strip the tulpa_obs_fit class so the underlying print method runs.
  cls <- setdiff(class(x), "tulpa_obs_fit")
  class(x) <- cls
  NextMethod("print", x, ...)
  invisible(x)
}
