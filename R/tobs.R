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
#' @param visit_data optional data frame of visit-level detection covariates
#'   with `nrow(visit_data) == nrow(y) * ncol(y)`.
#' @param spatial optional spatial spec from `tulpa` or `tulpaMesh` (e.g.
#'   [tulpa::spatial_bym2()], [tulpa::spatial_nngp()], or a [tulpaMesh::mesh]
#'   object), or one of the convenience helpers.
#' @param temporal optional temporal spec (e.g. [tulpa::temporal_ar1()]).
#' @param engine inference engine: `"auto"`, `"laplace"`, `"nested_laplace"`,
#'   or `"nuts"`. `"auto"` chooses the family's `default_engine`.
#' @param priors optional [tulpa::tulpa_priors()] object.
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
                 engine     = c("auto", "laplace", "nested_laplace", "nuts"),
                 priors     = NULL,
                 control    = list(),
                 ...) {

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
    engine     = engine,
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
                           spatial, temporal, engine, priors, control, ...) {
  if (is.null(detection)) {
    stop("occu() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu() requires `y` (an N x J detection-history matrix).", call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula        = formula,
    det_formula        = detection,
    data               = data,
    y                  = y,
    det_visit_formula  = if (!is.null(visit_data)) attr(visit_data, "formula") else NULL,
    det_visit_data     = visit_data
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal,
         method = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_dyn_occu <- function(formula, data, family, detection, y, visit_data,
                               spatial, temporal, engine, priors, control, ...) {
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
    list(model = model, spatial = spatial, temporal = temporal,
         method = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_ms_occu <- function(formula, data, family, detection, y, visit_data,
                              spatial, temporal, engine, priors, control, ...) {
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
    list(model = model, spatial = spatial, temporal = temporal,
         method = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_int_occu <- function(formula, data, family, detection, y, visit_data,
                               spatial, temporal, engine, priors, control, ...) {
  model <- .tobs_build_model(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    integrated  = TRUE
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal,
         method = .map_engine(engine), priors = priors),
    control
  ))
}

.dispatch_jsdm <- function(formula, data, family, detection, y, visit_data,
                           spatial, temporal, engine, priors, control, ...) {
  dots <- list(...)
  model <- .tobs_build_model(
    occ_formula = formula,
    data        = data,
    y           = y,
    jsdm        = TRUE,
    species     = dots$species
  )
  do.call(.tobs_fit_model, c(
    list(model = model, spatial = spatial, temporal = temporal,
         method = .map_engine(engine), priors = priors),
    control
  ))
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.map_engine <- function(engine) {
  # Engine name translation between the tobs vocabulary and what the underlying
  # fitter currently understands. Phase 0 maps both "laplace" and
  # "nested_laplace" to the single-Laplace mode and surfaces a NOTE when the
  # user asked for nested_laplace specifically.
  switch(
    engine,
    laplace        = "laplace",
    nested_laplace = {
      message(
        "tobs(): nested_laplace requested but Phase 0 dispatches to the ",
        "single-Laplace engine for occupancy families. Nested Laplace ",
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
