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
#'
#'   Structured effects are written as terms inside the formula, the way
#'   `lme4`, `mgcv`, and `INLA` do: spatial fields `icar(graph = adj)`,
#'   `bym2(graph = adj)`, `gp(lon, lat)`, `spde(lon, lat)`; random effects
#'   `re(group)`; temporal fields `temporal(time)`; spatially varying
#'   coefficients `svc(lon, lat, indices = ...)`; community latent factors
#'   `latent(k)`. A term enters whichever linear predictor it is written in
#'   (occupancy `formula` or `detection`). To share one realization across
#'   both predictors, tag the term with `id = "u"` and write `copy("u")` in
#'   the other formula.
#'
#'   Random effects also accept `lme4` bar syntax as shorthand for `re()`:
#'   `(1 | g)` is `re(g)`, `(x | g)` is a correlated random intercept and
#'   slope `re(g, type = "slope", covariate = x)`, and `(x || g)` drops the
#'   correlation (`correlated = FALSE`). Use `re()` directly for AR1/RW
#'   structures or other options.
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
#' @param visits optional visit-level detection covariates. Accepts
#'   either:
#'   * a named list of `[n_sites, max_visits]` matrices (the shape returned by
#'     `tobs_data()` in `det.covs`) — flattened internally to a long data
#'     frame in site-major order;
#'   * a data frame with `nrow(y) * ncol(y)` rows in site-major order.
#'
#'   When `visits` is provided without a `"formula"` attribute, the
#'   `detection` argument is interpreted as the visit-level detection
#'   formula and the site-level detection design matrix is an intercept
#'   only. To split visit-level and site-level detection covariates
#'   (e.g. visit-level effort plus site-level observer category), pass a
#'   long data frame with `attr(visits, "formula") = ~ effort` and use
#'   `detection = ~ observer` for the site-level terms.
#' @param method inference route, naming a fully-specified path rather than a
#'   pair of orthogonal knobs:
#'   * `"auto"` — the family's default route (see `default_engine`).
#'   * `"laplace"` — EM + Laplace with Gaussian marginals (fast default).
#'   * `"laplace_sla"` — Laplace with skew-corrected (simplified-Laplace)
#'     marginals.
#'   * `"laplace_gibbs"` / `"laplace_mi"` — Laplace with a post-EM Gibbs /
#'     multiple-imputation correction. These run the unpenalised EM, so the
#'     weakly-informative fixed-effect prior is disabled unless you pass
#'     `priors` explicitly.
#'   * `"nested_laplace"` — multi-block nested Laplace (single-season
#'     occupancy and cover-hurdle joint).
#'   * `"nested_laplace_sla"` — nested Laplace with skew-corrected marginals.
#'   * `"nuts"` — HMC / NUTS sampler (every structure; reports Rhat / ESS).
#' @param priors optional prior specification. For occupancy families fit
#'   with a Laplace method (`method = "laplace"`, `"laplace_sla"`,
#'   `"nested_laplace"`), pass a list or [occu_priors()] object to set
#'   weakly-informative quadratic priors on the fixed-effect coefficients
#'   (defaults pull the detection intercept toward `p = 0.5` and break the
#'   psi-p identifiability ridge at small `J`). Pass `priors = FALSE` to
#'   disable the default prior and recover the unpenalised MAP. The
#'   `"laplace_gibbs"` / `"laplace_mi"` routes set `priors = FALSE`
#'   automatically (the correction is defined on the unpenalised EM; see
#'   gcol33/tulpa#27). For NUTS, this is forwarded to the underlying tulpa
#'   engine.
#' @param control list of low-level engine controls. Names follow the
#'   dotted-separator convention. Sampler controls (`method = "nuts"`):
#'   * `n.iter` — total iterations per chain, including warmup (default 2000).
#'   * `n.warmup` — warmup / adaptation iterations per chain (default 1000).
#'   * `n.thin` — keep every `n.thin`-th post-warmup draw (default 1).
#'   * `n.chains` — number of chains, run with offset seeds and pooled
#'     (default 1). Split-Rhat / bulk / tail ESS are reported on `$convergence`.
#'   * `n.threads` — chains to run in parallel (default 1, sequential). Values
#'     `> 1` use a PSOCK cluster and require tulpaObs to be installed.
#'   * `adapt.delta` — target acceptance probability (default 0.8).
#'   * `max.treedepth` — NUTS maximum tree depth (default 10).
#'   * `seed` — base RNG seed; chain `c` uses `seed + c - 1` (default 42).
#'     The resolved per-chain seeds are stored on `$seeds`.
#'
#'   Laplace controls (`method = "laplace"` / `"laplace_sla"` /
#'   `"nested_laplace"`): `max.iter`, `tol`, `damping`, `sigma.beta`,
#'   `sigma.re.scale`. Stochastic-correction controls (`"laplace_gibbs"` /
#'   `"laplace_mi"`): `n.gibbs` / `n.imputations` (Rubin-pooled draw count)
#'   and `seed` (stored on `$seed`).
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
                 visits     = NULL,
                 method     = c("auto", "laplace", "laplace_sla",
                                "laplace_gibbs", "laplace_mi",
                                "nested_laplace", "nested_laplace_sla", "nuts"),
                 priors     = NULL,
                 control    = list(),
                 ...) {

  method <- match.arg(method)

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

  route   <- .tobs_resolve_method(method, family)
  engine  <- route$engine
  approx  <- route$approx

  # Gibbs / MI corrections are defined on the unpenalised EM, so the
  # weakly-informative fixed-effect prior is switched off for those routes
  # unless the user asked for one explicitly. A penalized correction is
  # blocked upstream: tulpa's Laplace block fitter takes no beta prior, so
  # the correction phases cannot regularize (gcol33/tulpa#27).
  if (route$correction %in% c("gibbs", "mi") && is.null(priors)) {
    priors <- FALSE
  }

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
    visits     = visits,
    engine     = engine,
    approx     = approx,
    correction = route$correction,
    priors     = priors,
    control    = control,
    ...
  )

  if (!inherits(fit, "tobs_fit")) {
    class(fit) <- c("tobs_fit", class(fit))
  }
  attr(fit, "tobs_family") <- family
  # Record the resolved public route for provenance / reproducibility.
  fit$method <- method
  fit
}


# ---------------------------------------------------------------------------
# Internal dispatchers — one per working family. Each wraps the internal
# model builder + fitter for that family's structure.
# ---------------------------------------------------------------------------

.dispatch_occu <- function(formula, data, family, detection, y, visits,
                           engine, priors, control,
                           approx = "gaussian_laplace",
                           correction = "none", ...) {
  if (is.null(detection)) {
    stop("occu() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu() requires `y` (an N x J detection-history matrix).", call. = FALSE)
  }
  vd <- .normalize_visits(visits, detection, n_sites = nrow(y),
                          max_visits = ncol(y))
  model <- .tobs_build_model(
    occ_formula        = formula,
    det_formula        = vd$det_formula,
    data               = data,
    y                  = y,
    det_visit_formula  = vd$det_visit_formula,
    det_visit_data     = vd$visits
  )
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_dyn_occu <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
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
    list(model = model,
         method = .map_engine(engine, family = "dyn_occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_ms_occu <- function(formula, data, family, detection, y, visits,
                              engine, priors, control,
                              approx = "gaussian_laplace",
                              correction = "none", ...) {
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
    list(model = model,
         method = .map_engine(engine, family = "ms_occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_int_occu <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
  model <- .tobs_build_model(
    occ_formula = formula,
    det_formula = detection,
    data        = data,
    y           = y,
    integrated  = TRUE
  )
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "int_occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_jsdm <- function(formula, data, family, detection, y, visits,
                           engine, priors, control,
                           approx = "gaussian_laplace",
                           correction = "none", ...) {
  dots <- list(...)
  model <- .tobs_build_model(
    occ_formula = formula,
    data        = data,
    y           = y,
    jsdm        = TRUE,
    species     = dots$species
  )
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "jsdm"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Public `method` names are sugar over the orthogonal internal triple
# (engine, approx, correction). Each method names one fully-specified route;
# invalid cross-products (e.g. NUTS with an SLA marginal) simply have no name.
.tobs_method_table <- list(
  laplace            = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "none"),
  laplace_sla        = list(engine = "laplace",        approx = "simplified_laplace", correction = "none"),
  laplace_gibbs      = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "gibbs"),
  laplace_mi         = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "mi"),
  nested_laplace     = list(engine = "nested_laplace", approx = "gaussian_laplace",   correction = "none"),
  nested_laplace_sla = list(engine = "nested_laplace", approx = "simplified_laplace", correction = "none"),
  nuts               = list(engine = "nuts",           approx = "gaussian_laplace",   correction = "none")
)

# Resolve a public method name to the internal (engine, approx, correction)
# triple. `"auto"` maps the family's default engine to its base method.
.tobs_resolve_method <- function(method, family) {
  if (identical(method, "auto")) {
    method <- switch(
      family$default_engine,
      laplace        = "laplace",
      nested_laplace = "nested_laplace",
      nuts           = "nuts",
      stop(sprintf("Family '%s' has an unknown default_engine '%s'.",
                   family$name, family$default_engine), call. = FALSE)
    )
  }
  route <- .tobs_method_table[[method]]
  if (is.null(route)) {
    stop(sprintf("Unknown method '%s'. See `?tobs` for the route list.",
                 method), call. = FALSE)
  }
  route
}

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

# Normalize the user's `visits` argument to the shape `.tobs_build_single`
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
#                                              `visits`
#   * data.frame with N*J rows, no formula  -> treat `detection` as visit-level
#       attribute                              (intercept dropped); site-level
#                                              X_det is intercept-only
#
# Returns a list with:
#   visits             — long data frame (or NULL)
#   det_visit_formula  — formula applied to visits (or NULL)
#   det_formula        — formula applied to site-level `data`
.normalize_visits <- function(visits, detection,
                              n_sites, max_visits) {
  if (is.null(visits)) {
    return(list(visits = NULL,
                det_visit_formula = NULL,
                det_formula = detection))
  }

  expected_rows <- n_sites * max_visits

  # Case 1: list of [n_sites, max_visits] matrices (tobs_data() output)
  if (is.list(visits) && !is.data.frame(visits)) {
    nms <- names(visits)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("`visits` (list of matrices) must be a named list; ",
           "names become the column names of the flattened frame.",
           call. = FALSE)
    }
    bad <- vapply(visits, function(m) {
      !is.matrix(m) || nrow(m) != n_sites || ncol(m) != max_visits
    }, logical(1))
    if (any(bad)) {
      stop(sprintf(
        "`visits` elements must be [%d x %d] matrices matching y; ",
        n_sites, max_visits),
        sprintf("element(s) %s have wrong shape.",
                paste(nms[bad], collapse = ", ")),
        call. = FALSE)
    }
    flat <- as.data.frame(
      lapply(visits, function(m) as.vector(t(m))),
      stringsAsFactors = FALSE
    )
    return(list(
      visits = flat,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  # Case 2 / 3: data frame
  if (is.data.frame(visits)) {
    if (nrow(visits) != expected_rows) {
      stop(sprintf(
        "`visits` (data frame) must have %d rows (n_sites * max_visits); ",
        expected_rows),
        sprintf("got %d.", nrow(visits)),
        call. = FALSE)
    }
    attached <- attr(visits, "formula")
    if (!is.null(attached)) {
      # Dual-formula power-user mode: detection stays site-level
      return(list(
        visits = visits,
        det_visit_formula = attached,
        det_formula = detection
      ))
    }
    return(list(
      visits = visits,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  stop("`visits` must be NULL, a named list of [n_sites x max_visits] ",
       "matrices, or a long data frame with n_sites * max_visits rows; ",
       "got ", paste(class(visits), collapse = "/"), ".",
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
    if (!is.null(x$n_chains)) {
      cat(sprintf(" (%d chain%s%s)", x$n_chains,
                  if (x$n_chains > 1L) "s" else "",
                  if (!is.null(x$n_thin) && x$n_thin > 1L)
                    sprintf(", thin %d", x$n_thin) else ""))
    }
    if (!is.null(x$epsilon) && !is.na(x$epsilon)) {
      cat(sprintf(", step size: %.4f", x$epsilon))
    }
    cat("\n")
  }
  # Reproducibility: seeds for stochastic routes (NUTS chains / MI / Gibbs).
  if (!is.null(x$seeds)) {
    cat(sprintf("  Seeds: %s\n", paste(x$seeds, collapse = ", ")))
  } else if (!is.null(x$seed)) {
    cat(sprintf("  Seed: %d\n", x$seed))
  }
  if (!is.null(x$convergence)) {
    rh <- x$convergence$rhat
    eb <- x$convergence$ess_bulk
    if (any(is.finite(rh)) || any(is.finite(eb))) {
      cat(sprintf("  Convergence: max Rhat %.3f, min bulk ESS %.0f\n",
                  max(rh, na.rm = TRUE), min(eb, na.rm = TRUE)))
      if (any(rh > 1.01, na.rm = TRUE)) {
        cat("    WARNING: Rhat > 1.01 for some parameters; chains may not have ",
            "mixed. Increase n.iter / n.chains.\n", sep = "")
      }
    }
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
