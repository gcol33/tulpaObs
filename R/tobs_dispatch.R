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
  if (is.null(dots$colonization)) {
    stop("dyn_occu() requires a `colonization = ~ ...` argument.", call. = FALSE)
  }
  if (is.null(dots$extinction)) {
    stop("dyn_occu() requires an `extinction = ~ ...` argument.", call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula  = formula,
    det_formula  = detection,
    data         = data,
    y            = y,
    col_formula  = dots$colonization,
    ext_formula  = dots$extinction
  )
  # Season-varying colonization / extinction (a [n_sites x (T-1)] matrix column;
  # gcol33/tulpaObs#124) is wired for the Laplace-EM engines only. The C++ NUTS
  # forward reads one gamma / epsilon linear predictor per site, so an
  # interval-indexed rate is gated there with a pointer.
  if ((isTRUE(model$col_season_varying) || isTRUE(model$ext_season_varying)) &&
      identical(engine, "nuts")) {
    stop("dyn_occu(): season-varying colonization / extinction (a ",
         "[n_sites x (T-1)] matrix covariate) is not yet wired for ",
         "method = \"nuts\"; use method = \"laplace\" ",
         "(gcol33/tulpaObs#124).", call. = FALSE)
  }
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "dyn_occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_ms_count <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
  dots <- list(...)
  if (!is.null(detection)) {
    stop("ms_count() has no detection process; drop the `detection` formula.",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_count() requires `y` (an n_sites x n_species matrix, or a named ",
         "list of n_species count vectors).", call. = FALSE)
  }
  response <- family$params$response %||% "poisson"
  if (!is.null(dots$trials) && !identical(response, "binomial")) {
    stop("ms_count(): `trials` applies only to ms_count(response = ",
         "\"binomial\"); drop it for the ", response, " response.",
         call. = FALSE)
  }
  bind  <- .tobs_bind_formulas(list(mu = formula), data)
  model <- .tobs_build_ms_count(
    formula = bind$fe$mu, data = data, y = y, species = dots$species,
    response = response, trials = dots$trials, structured_terms = bind$terms)

  # Polya-Gamma Gibbs (method = "pg_gibbs", #126): the logistic responses
  # (binomial / bernoulli) admit PG augmentation, so the community coefficient
  # posterior is sampled exactly by a per-species conjugate Gibbs -- a calibrated
  # community-variance posterior. Non-spatial; the binder rejects a structured
  # term for pg_gibbs below via the shared structs check being bypassed here (v1
  # has no PG-spatial path).
  if (identical(engine, "pg_gibbs")) {
    if (!(response %in% c("bernoulli", "binomial"))) {
      stop("ms_count(): method = \"pg_gibbs\" applies to the logistic responses ",
           "(binomial; jsdm() is bernoulli). For ", response,
           " use method = \"laplace\" or \"nuts\".", call. = FALSE)
    }
    if (!is.null(bind$terms) && length(bind$terms) > 0L) {
      stop("ms_count(): method = \"pg_gibbs\" is the non-spatial community fit; ",
           "drop the structured term (PG-spatial is a follow-up).", call. = FALSE)
    }
    return(.tobs_fit_ms_count_pg_gibbs(
      model, priors = priors,
      sigma.beta = control[["sigma.beta"]] %||% 2.5,
      n.iter   = as.integer(control[["n.iter"]]   %||% 3000L),
      n.warmup = as.integer(control[["n.warmup"]] %||% 1500L),
      n.chains = max(as.integer(control[["n.chains"]] %||% 2L), 2L),
      n.thin   = as.integer(control[["n.thin"]]   %||% 1L),
      seed     = as.integer(control[["seed"]]     %||% 1L),
      verbose  = isTRUE(control[["verbose"]])))
  }

  # The binomial community response is wired for the non-spatial Laplace-EM only
  # (community svcPGBinom, gcol33/tulpaObs#125). NUTS needs a binomial family in
  # the in-tree C++ FullGradFn, and a shared field / latent factor needs the
  # binomial working callback in the block-coordinate driver -- both follow-ups.
  if (identical(response, "binomial")) {
    if (identical(engine, "nuts")) {
      stop("ms_count(response = \"binomial\"): NUTS is not yet wired for the ",
           "binomial community response; use method = \"laplace\" ",
           "(gcol33/tulpaObs#125).", call. = FALSE)
    }
    structs_b <- .tobs_structures_from_model(model)
    if (!is.null(structs_b$spatial) || !is.null(structs_b$latent)) {
      stop("ms_count(response = \"binomial\"): a shared areal field / latent ",
           "factor is not yet wired for the binomial community response; use ",
           "method = \"laplace\" without a structured term. The single-species ",
           "binomial areal field (svcPGBinom) is count(response = ",
           "\"binomial\") + icar() (gcol33/tulpaObs#125).", call. = FALSE)
    }
  }

  # A shared areal field (icar()) on the abundance formula routes to the
  # community-spatial fitter (the sfMsAbund analogue) under nested_laplace; a
  # non-areal structured term (temporal / re / svc / latent) or a non-icar field
  # is not yet wired -- error rather than silently drop (gcol33/tulpaObs#117).
  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc)) {
    stop("ms_count(): temporal / re / svc terms are not yet wired for the ",
         "community count family; a shared areal field icar() or a latent() ",
         "factor term are the structured terms supported (gcol33/tulpaObs#117).",
         call. = FALSE)
  }
  # A shared latent structure -- a shared areal field icar()/SVC bar (the
  # sfMsAbund / svcMsAbund analogue), latent() factors (lfMsAbund), or BOTH (the
  # spatial-factor case) -- routes to the unified community latent fitter. A
  # shared areal field needs method = "nested_laplace"; a factor-only model uses
  # method = "laplace". The fitter validates the field kind (icar only), the
  # one-node-per-site map, and a complete y.
  if (!is.null(structs$spatial) || !is.null(structs$latent)) {
    if (!is.null(structs$spatial)) {
      if (!identical(engine, "nested_laplace")) {
        stop("ms_count(): a shared areal field needs method = ",
             "\"nested_laplace\" (drop the icar() term for the non-spatial ",
             "community fit).", call. = FALSE)
      }
    } else if (identical(engine, "nested_laplace") || identical(engine, "nuts")) {
      stop("ms_count(): a latent() factor model uses method = \"laplace\" ",
           "(the block-coordinate community Laplace-EM).", call. = FALSE)
    }
    return(.tobs_fit_ms_count_latent(
      model, spatial = structs$spatial, latent = structs$latent,
      max.iter   = control[["max.iter"]] %||% 200L,
      tol        = control[["tol"]] %||% 1e-4,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      priors     = priors,
      max.outer  = control[["max.outer"]] %||% 25L,
      verbose    = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nested_laplace")) {
    stop("ms_count(): method = \"nested_laplace\" needs a shared areal field ",
         "icar() on the formula. For the non-spatial community fit use ",
         "method = \"laplace\".", call. = FALSE)
  }
  if (identical(engine, "nuts")) {
    # NUTS over the exact joint community count posterior (community means,
    # per-species deviations, community covariance) via the in-tree C++
    # FullGradFn, warm-started at the Laplace-EM mode (gcol33/tulpaObs#117).
    return(.tobs_fit_ms_count_nuts(
      model,
      sigma.beta    = control[["sigma.beta"]] %||% 10,
      sigma.logr    = control[["sigma.logr"]] %||% 1.5,
      n.iter        = as.integer(control[["n.iter"]]   %||% 1000L),
      n.warmup      = as.integer(control[["n.warmup"]] %||% 1000L),
      n.chains      = as.integer(control[["n.chains"]] %||% 1L),
      max.treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
      adapt.delta   = control[["adapt.delta"]] %||% 0.9,
      seed          = as.integer(control[["seed"]] %||% 1L),
      verbose       = isTRUE(control[["verbose"]])))
  }

  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_count, fit_args)
}

.dispatch_ms_occu <- function(formula, data, family, detection, y, visits,
                              engine, priors, control,
                              approx = "gaussian_laplace",
                              correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_occu() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_occu() requires `y` (a 3D array [n_sites x max_visits x ",
         "n_species] or a named list of detection-history matrices).",
         call. = FALSE)
  }

  # Resolve structured terms on the occupancy / detection formulas (the spatial
  # field on the occupancy arm routes to the areal community fitter; the FE part
  # is what model.matrix sees). Bind precomputes the spatial spec (graph / CSR /
  # n_units / type).
  bind <- .tobs_bind_formulas(list(psi = formula, p = detection), data)

  model <- .tobs_build_ms_occu(
    occ_formula = bind$fe$psi, det_formula = bind$fe$p,
    data = data, y = y, species = dots$species,
    structured_terms = bind$terms)
  structs <- .tobs_structures_from_model(model)

  # Any structured term other than a shared areal field or latent() factors is
  # not wired for the community occupancy family.
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc)) {
    stop("ms_occu(): temporal / re / svc terms are not wired for the community ",
         "occupancy family; a shared areal field (icar()/bym2()/car_proper()) ",
         "on the occupancy formula and latent() factors are the structured ",
         "terms supported.", call. = FALSE)
  }

  # A shared areal field (icar()/bym2()/car_proper()) on the occupancy formula
  # routes to the nested-Laplace community-occupancy fitter (the occupancy
  # analogue of sfMsNMix; tulpaObs#75). It must sit on the occupancy arm only.
  if (!is.null(structs$spatial)) {
    if (identical(engine, "nuts")) {
      stop("method = \"nuts\" for ms_occu() is the non-spatial community ",
           "occupancy sampler; a shared areal field uses ",
           "method = \"nested_laplace\".", call. = FALSE)
    }
    if (!identical(engine, "nested_laplace")) {
      stop("a shared areal field on the occupancy formula needs ",
           "method = \"nested_laplace\" for ms_occu().", call. = FALSE)
    }
    if (isTRUE(structs$spatial$shared[2L])) {
      stop("ms_occu() areal field sits on the occupancy arm only; a field on ",
           "the detection formula is not supported.", call. = FALSE)
    }
    # A varying-coefficient bar spatial(~ 1 + w || cell, graph) (the svcMsPGOcc
    # analogue, gcol33/tulpaObs#118), or latent() factors alongside the field
    # (the sfMsPGOcc analogue: a shared field plus per-species factor loadings,
    # gcol33/tulpaObs#119), route to the block-coordinate latent fitter. A plain
    # intercept field alone keeps the in-tree C++ community-spatial nested
    # Laplace-EM.
    if (isTRUE(structs$spatial$is_bar) || isTRUE(structs$spatial$is_multifield) ||
        !is.null(structs$latent)) {
      return(.tobs_fit_ms_occu_field(
        model, spatial = structs$spatial, latent = structs$latent,
        max.iter  = control[["max.iter"]] %||% 200L,
        tol       = control[["tol"]] %||% 1e-4,
        sigma.beta = control[["sigma.beta"]] %||% 5,
        priors    = priors,
        max.outer = control[["max.outer"]] %||% 20L,
        verbose   = isTRUE(control[["verbose"]])))
    }
    return(.tobs_fit_ms_occu_spatial(
      model, spatial = structs$spatial,
      max.iter = control[["max.iter"]] %||% 100L,
      verbose  = isTRUE(control[["verbose"]])))
  }

  # latent() factors with no shared field: the lfMsPGOcc analogue -- residual
  # species co-occurrence on the occupancy arm via Q per-site latent factors with
  # per-species loadings, by the block-coordinate community Laplace-EM
  # (gcol33/tulpaObs#119).
  if (!is.null(structs$latent)) {
    if (identical(engine, "nested_laplace") || identical(engine, "nuts")) {
      stop("ms_occu(): a latent() factor model uses method = \"laplace\" (the ",
           "block-coordinate community Laplace-EM).", call. = FALSE)
    }
    return(.tobs_fit_ms_occu_field(
      model, spatial = NULL, latent = structs$latent,
      max.iter  = control[["max.iter"]] %||% 200L,
      tol       = control[["tol"]] %||% 1e-4,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      priors    = priors,
      max.outer = control[["max.outer"]] %||% 20L,
      verbose   = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nested_laplace")) {
    stop("method = \"nested_laplace\" for ms_occu() needs a shared areal field ",
         "(icar()/bym2()/car_proper()) on the occupancy formula. For the ",
         "non-spatial community fit use method = \"laplace\".", call. = FALSE)
  }

  # Polya-Gamma Gibbs (method = "pg_gibbs", spOccupancy msPGOcc; tulpaObs#115,
  # #126): a hierarchical PG Gibbs over the exact community posterior -- per-
  # species PG-augmented conjugate coefficient updates with conjugate community
  # mean + Inverse-Gamma community variance draws. Gives a CALIBRATED community-
  # variance posterior (the Laplace-EM leaves those attenuated).
  if (identical(engine, "pg_gibbs")) {
    return(.tobs_fit_ms_occu_pg_gibbs(
      model, priors = priors,
      sigma.beta = control[["sigma.beta"]] %||% 2.5,
      n.iter   = as.integer(control[["n.iter"]]   %||% 3000L),
      n.warmup = as.integer(control[["n.warmup"]] %||% 1500L),
      n.chains = max(as.integer(control[["n.chains"]] %||% 2L), 2L),
      n.thin   = as.integer(control[["n.thin"]]   %||% 1L),
      seed     = as.integer(control[["seed"]]     %||% 1L),
      verbose  = isTRUE(control[["verbose"]])))
  }

  # NUTS (method = "nuts", tulpaObs#69): sample the exact joint posterior of the
  # non-spatial community single-season occupancy (community means, per-species
  # deviations, and the two independent per-arm community covariances) via the
  # in-tree C++ FullGradFn over the closed-form occupancy two-state per-(species,
  # site) marginal (R/ms_occu_nuts.R, src/ms_occu_nuts.cpp), warm-started at the
  # community Laplace-EM mode.
  if (identical(engine, "nuts")) {
    return(.tobs_fit_ms_occu_nuts(
      model,
      sigma.beta    = control[["sigma.beta"]] %||% 5,
      n.iter        = as.integer(control[["n.iter"]]   %||% 1000L),
      n.warmup      = as.integer(control[["n.warmup"]] %||% 1000L),
      n.chains      = as.integer(control[["n.chains"]] %||% 1L),
      max.treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
      adapt.delta   = control[["adapt.delta"]] %||% 0.9,
      seed          = as.integer(control[["seed"]] %||% 1L),
      verbose       = isTRUE(control[["verbose"]])))
  }

  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_occu, fit_args)
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
  if (!is.null(detection)) {
    stop("jsdm() has no detection process; drop the `detection` formula.",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("jsdm() requires `y` (an n_sites x n_species presence/absence matrix, ",
         "or a named list of n_species 0/1 vectors).", call. = FALSE)
  }
  # A JSDM is the community GLMM on an observed presence/absence response: the
  # same model class as ms_count() (per-species coefficients with a Gaussian
  # community covariance, no detection, no latent state) with a logit link. It
  # therefore shares the binder, the community Laplace-EM, the latent-structure
  # driver, the NUTS target, and every S3 method (gcol33/tulpaObs#121), which is
  # what makes latent() factors (lfJSDM) and a shared field + factors (sfJSDM)
  # fall out with no separate fitter.
  bind  <- .tobs_bind_formulas(list(mu = formula), data)
  model <- .tobs_build_ms_count(
    formula = bind$fe$mu, data = data, y = y, species = dots$species,
    response = "bernoulli", structured_terms = bind$terms)

  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc)) {
    stop("jsdm(): temporal / re / svc terms are not wired for the JSDM family; ",
         "a shared areal field (icar()/car_proper()/bym2()) and latent() ",
         "factors are the structured terms supported.", call. = FALSE)
  }

  # A shared areal field (sfJSDM), latent() factors (lfJSDM), or BOTH route to
  # the community latent driver. A field needs nested_laplace; a factor-only
  # model is the block-coordinate laplace.
  if (!is.null(structs$spatial) || !is.null(structs$latent)) {
    if (!is.null(structs$spatial)) {
      if (!identical(engine, "nested_laplace")) {
        stop("jsdm(): a shared areal field needs method = \"nested_laplace\" ",
             "(drop the areal term for the non-spatial JSDM).", call. = FALSE)
      }
    } else if (identical(engine, "nested_laplace") || identical(engine, "nuts")) {
      stop("jsdm(): a latent() factor model uses method = \"laplace\" (the ",
           "block-coordinate community Laplace-EM).", call. = FALSE)
    }
    return(.tobs_fit_ms_count_latent(
      model, spatial = structs$spatial, latent = structs$latent,
      max.iter   = control[["max.iter"]] %||% 200L,
      tol        = control[["tol"]] %||% 1e-4,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      priors     = priors,
      max.outer  = control[["max.outer"]] %||% 25L,
      verbose    = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nested_laplace")) {
    stop("jsdm(): method = \"nested_laplace\" needs a shared areal field ",
         "(icar()/car_proper()/bym2()) on the formula. For the non-spatial ",
         "JSDM use method = \"laplace\".", call. = FALSE)
  }
  # Polya-Gamma Gibbs (method = "pg_gibbs", spOccupancy msPGOcc-family; #126): the
  # community Bernoulli GLMM has no latent state, so the PG sampler is the pure
  # per-species conjugate update + community mean / Inverse-Gamma variance -- a
  # calibrated community-variance posterior (the Laplace-EM attenuates it).
  if (identical(engine, "pg_gibbs")) {
    return(.tobs_fit_ms_count_pg_gibbs(
      model, priors = priors,
      sigma.beta = control[["sigma.beta"]] %||% 2.5,
      n.iter   = as.integer(control[["n.iter"]]   %||% 3000L),
      n.warmup = as.integer(control[["n.warmup"]] %||% 1500L),
      n.chains = max(as.integer(control[["n.chains"]] %||% 2L), 2L),
      n.thin   = as.integer(control[["n.thin"]]   %||% 1L),
      seed     = as.integer(control[["seed"]]     %||% 1L),
      verbose  = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nuts")) {
    # Samples the exact joint community posterior (community means, per-species
    # deviations, community covariance) over the Bernoulli response via the
    # family-aware in-tree C++ FullGradFn, warm-started at the Laplace-EM mode.
    return(.tobs_fit_ms_count_nuts(
      model,
      sigma.beta    = control[["sigma.beta"]] %||% 10,
      sigma.logr    = control[["sigma.logr"]] %||% 1.5,
      n.iter        = as.integer(control[["n.iter"]]   %||% 1000L),
      n.warmup      = as.integer(control[["n.warmup"]] %||% 1000L),
      n.chains      = as.integer(control[["n.chains"]] %||% 1L),
      max.treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
      adapt.delta   = control[["adapt.delta"]] %||% 0.9,
      seed          = as.integer(control[["seed"]] %||% 1L),
      verbose       = isTRUE(control[["verbose"]])))
  }

  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_count, fit_args)
}

.dispatch_count <- function(formula, data, family, detection, y, visits,
                            engine, priors, control, trials = NULL,
                            approx = "gaussian_laplace",
                            correction = "none", ...) {
  if (!is.null(detection)) {
    stop("count() has no detection process; drop the `detection` formula.",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("count() requires `y` (a numeric vector, one value per site), or the ",
         "response on a two-sided `formula` left-hand side.", call. = FALSE)
  }
  response <- family$params$response %||% "poisson"
  if (!is.null(trials) && !identical(response, "binomial")) {
    stop("count(): `trials` applies only to count(response = \"binomial\"); ",
         "drop it for the ", response, " response.", call. = FALSE)
  }
  model <- .tobs_build_model(occ_formula = formula, data = data, y = y,
                             count = TRUE, count_response = response,
                             count_trials = trials)

  # A plain areal field -- icar()/car_proper() -- on the abundance formula routes
  # to nested-Laplace (the spAbund analogue): the field is a latent GMRF prior on
  # the count GLMM block, integrated over its hyperparameters. Every other
  # structured term (temporal / re / svc / latent) and the varying-coefficient /
  # bar, bym2, or group_var areal forms are not yet wired for the count family --
  # error with a pointer rather than silently drop them (gcol33/tulpaObs#117).
  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc) || !is.null(structs$latent)) {
    stop("count(): temporal / re / svc / latent terms are not yet wired for ",
         "the count family; a plain areal field (icar()/car_proper()) ",
         "is the only structured term supported (gcol33/tulpaObs#117).",
         call. = FALSE)
  }
  spatial_areal <- FALSE
  if (!is.null(structs$spatial)) {
    sp <- structs$spatial
    # A bar spatial(~ 1 + w || cell, graph) -- or the explicit two-term form --
    # carries an intercept field plus one varying-coefficient field per covariate
    # (the spAbundance svcAbund analogue, gcol33/tulpaObs#120). Resolve the
    # ordered field list and validate each field's OWN kind: the bar's `type` is
    # not the per-field kind. A plain intercept field resolves to a length-1 list,
    # so both forms share one validation path.
    sp_fields <- .tobs_resolve_occu_spatial_fields(sp, model)
    ftypes <- vapply(sp_fields, function(f) f$type %||% "unknown", character(1))
    # icar / car_proper only: their eta contribution is exactly w * f[cell]
    # (d_fac = 1), so the demeaned per-cell field reconstructs exactly from the
    # grid modes (the shared nested field summary). A bym2 field mixes a
    # structured (phi) and an unstructured (theta) component with
    # hyperparameter-dependent scales, so its per-cell field is not reconstructed
    # on this generic path; the improper non-intrinsic car() is likewise not
    # wired. Both are #117 follow-ups -- point to the supported fields rather
    # than return a fit with no field.
    if (!all(ftypes %in% c("icar", "car_proper"))) {
      stop(sprintf(paste0(
        "count(): an areal field on the count formula supports icar() or ",
        "car_proper(); got '%s'. bym2() (mixed structured/unstructured field) ",
        "and the improper car(), plus continuous-mesh spde()/gp(), are not yet ",
        "wired for the count family (gcol33/tulpaObs#117)."),
        paste(unique(setdiff(ftypes, c("icar", "car_proper"))),
              collapse = "' / '")), call. = FALSE)
    }
    # A bar spells its grouping as `|| cell`, so group_var is set even when the
    # field has one node per site (the identity map). Reject it only when it
    # actually AGGREGATES -- fewer field nodes than sites -- which the count
    # nested path does not yet reconstruct.
    n_nodes <- vapply(sp_fields, function(f)
      if (is.null(f$graph)) NA_integer_ else nrow(f$graph), integer(1))
    if (any(!is.na(n_nodes) & n_nodes < nrow(data))) {
      stop("count(): a spatial group_var mapping several sites to one field ",
           "node (sites > cells) is not yet wired for the count family; one ",
           "field node per site is required (gcol33/tulpaObs#117).",
           call. = FALSE)
    }
    # Areal count is Poisson- OR binomial-only. With one field node per site the
    # negbin size / gaussian residual variance and the latent field are
    # FUNDAMENTALLY not identified together, not merely awkward to fit: the field
    # (one free value per site) absorbs all extra-Poisson variation, so the
    # field-integrated marginal likelihood is monotone in the dispersion toward
    # the Poisson limit (verified: size -> Inf, residual variance -> 0). No
    # estimator -- outer loop OR a joint dispersion grid -- recovers the
    # dispersion in this design. Identification needs replication within a site
    # (an N-mixture, abun()) or more sites than field nodes (a group_var areal
    # term, not yet wired for count). Poisson has no dispersion parameter and is
    # cleanly identified. The BINOMIAL response is also cleanly identified: its
    # variance is pinned by the trial count n, so there is no free dispersion for
    # the field to absorb -- this is exactly spOccupancy's svcPGBinom, which fits
    # a per-node field even at trials = 1 (gcol33/tulpaObs#125).
    if (!response %in% c("poisson", "binomial")) {
      stop(sprintf(paste0(
        "count(response = \"%s\") with an areal field is not identifiable: with ",
        "one field node per site the %s and the latent field are confounded ",
        "(the field absorbs all overdispersion). Use a Poisson areal count -- ",
        "count() -- or drop the areal term for a non-spatial %s fit; for ",
        "overdispersion with a spatial signal use abun() (N-mixture, replicated ",
        "counts) instead (gcol33/tulpaObs#117)."),
        response,
        if (identical(response, "negbin")) "negbin size"
        else "gaussian residual variance",
        response), call. = FALSE)
    }
    spatial_areal <- TRUE
  }

  # Engine resolution mirrors jsdm(): an areal field needs nested_laplace; the
  # non-spatial GLMM needs laplace. Reject the mismatched pairing loudly.
  if (spatial_areal) {
    if (!identical(engine, "nested_laplace")) {
      stop("count(): a plain areal field on the formula needs ",
           "method = \"nested_laplace\". For a non-spatial count GLMM drop the ",
           "areal term (or use method = \"laplace\").", call. = FALSE)
    }
  } else if (identical(engine, "nested_laplace")) {
    stop("count(): method = \"nested_laplace\" needs a plain areal field ",
         "(icar()/car_proper()) on the formula. For a non-spatial count ",
         "GLMM use method = \"laplace\".", call. = FALSE)
  }

  fit_once <- function(phi) {
    model$count_phi <- phi
    do.call(.tobs_fit_model, c(
      list(model = model, method = engine, priors = priors,
           approx = approx, correction = correction),
      control))
  }

  # Poisson has no dispersion. For negbin (size) / gaussian (residual variance)
  # tulpa_laplace takes a FIXED phi per fit, so estimate it in an outer loop:
  # fit beta given phi, update phi from the fitted mean, refit, until log(phi)
  # converges. Gaussian with an identity link converges in one update (beta is
  # phi-free); negbin iterates (the NB weights depend on the size).
  p_mu   <- model$process_info[[1]]$p
  is_log <- identical(model$link, "log")
  yv     <- as.numeric(model$y_count)
  mu_of  <- function(fit) {
    # The dispersion loop (negbin size / gaussian variance) runs on the
    # non-spatial path only -- an areal count is Poisson-only, so there is never
    # a latent field to fold into the fitted mean here.
    beta <- fit$means[seq_len(p_mu)]
    eta  <- as.vector(model$X_occ %*% beta)
    if (is_log) exp(eta) else eta
  }
  update_phi <- function(fit) {
    mu <- mu_of(fit)
    if (identical(response, "gaussian")) {
      max(mean((yv - mu)^2), 1e-8)
    } else { # negbin size by profile MLE given the fitted mean
      nll <- function(r) -sum(stats::dnbinom(yv, size = r,
                                             mu = pmax(mu, 1e-8), log = TRUE))
      stats::optimize(nll, interval = c(1e-3, 1e4))$minimum
    }
  }

  if (response %in% c("poisson", "binomial")) {
    # No free dispersion (binomial variance is pinned by the trial count).
    fit <- fit_once(1.0)
  } else {
    phi <- 1.0
    fit <- fit_once(phi)
    for (it in seq_len(25L)) {
      phi_new <- update_phi(fit)
      if (abs(log(phi_new) - log(phi)) < 1e-4) { phi <- phi_new; break }
      phi <- phi_new
      fit <- fit_once(phi)
    }
    fit <- fit_once(phi)
    disp_name <- if (identical(response, "negbin")) "size" else "variance"
    fit$count_dispersion <- stats::setNames(list(response, phi),
                                            c("response", disp_name))
    fit$count_dispersion$phi <- phi
  }
  fit$model$response <- response
  fit$model$link     <- model$link
  fit
}

.dispatch_abun <- function(formula, data, family, detection, y, visits,
                           engine, priors, control,
                           approx = "gaussian_laplace",
                           correction = "none", ...) {
  if (is.null(detection)) {
    stop("abun() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("abun() requires `y` (an N x J integer count matrix).", call. = FALSE)
  }
  vd <- .normalize_visits(visits, detection, n_sites = nrow(y),
                          max_visits = ncol(y))
  # Zero-inflated N-mixture (zip / zinb, gcol33/tulpaObs#116) is the non-spatial
  # laplace path only in v1; gate here (before engine routing) so a nuts /
  # nested_laplace request errors with a pointer instead of silently dropping the
  # structural-zero component down a non-ZI sampler.
  if ((family$params$mixture %||% "poisson") %in% c("zip", "zinb") &&
      !identical(engine %||% "laplace", "laplace")) {
    stop(sprintf(paste0("abun(mixture = \"%s\") supports method = \"laplace\" ",
                        "only; got \"%s\"."), family$params$mixture, engine),
         call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula        = formula,
    det_formula        = vd$det_formula,
    data               = data,
    y                  = y,
    abundance          = TRUE,
    det_visit_formula  = vd$det_visit_formula,
    det_visit_data     = vd$visits
  )
  # K_max and the abundance mixture travel with the family object (abun(K_max =,
  # mixture =)); thread them into the fitter alongside the engine controls.
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "abun"), priors = priors,
         approx = approx, correction = correction,
         K.max = family$params$K_max, mixture = family$params$mixture),
    control
  ))
}

.dispatch_removal <- function(formula, data, family, detection, y, visits,
                              engine, priors, control,
                              approx = "gaussian_laplace",
                              correction = "none", ...) {
  if (is.null(detection)) {
    stop("removal() requires a `detection` formula (per-pass detection).",
         call. = FALSE)
  }
  if (is.null(y)) {
    stop("removal() requires `y` (an N x K integer matrix of per-pass ",
         "removals, passes in column order).", call. = FALSE)
  }
  vd <- .normalize_visits(visits, detection, n_sites = nrow(y),
                          max_visits = ncol(y))
  model <- .tobs_build_removal(
    abund_formula     = formula,
    det_formula       = vd$det_formula,
    data              = data,
    y                 = y,
    det_visit_formula = vd$det_visit_formula,
    det_visit_data    = vd$visits
  )
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "removal"), priors = priors,
         approx = approx, correction = correction,
         K.max = family$params$K_max, mixture = family$params$mixture),
    control
  ))
}

.dispatch_distance <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
  if (is.null(detection)) {
    stop("distance() requires a `detection` formula (the site-level log-sigma ",
         "detection-scale model, e.g. ~ habitat).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("distance() requires `y` (an n_sites x n_bins integer matrix of ",
         "per-distance-bin detected counts).", call. = FALSE)
  }
  cutpoints <- family$params$cutpoints
  if (is.null(cutpoints)) {
    stop("distance() requires `cutpoints` (the distance-bin edges, ",
         "0 = c_0 < ... < c_B); pass distance(cutpoints = ...).", call. = FALSE)
  }
  model <- .tobs_build_distance(
    abund_formula = formula, det_formula = detection, data = data, y = y,
    cutpoints = cutpoints, key = family$params$key,
    transect = family$params$transect, mixture = family$params$mixture)
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "distance"), priors = priors,
         approx = approx, correction = correction,
         K.max = family$params$K_max, mixture = family$params$mixture),
    control
  ))
}

.dispatch_dyn_abun <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("dyn_abun() requires a `detection` formula (the per-visit detection ",
         "model).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("dyn_abun() requires `y` (a 3D array [n_sites x max_visits x n_seasons] ",
         "of counts, or a list of per-season count matrices).", call. = FALSE)
  }
  model <- .tobs_build_dyn_abun(
    occ_formula = formula, det_formula = detection, data = data, y = y,
    omega_formula = dots$omega %||% ~1,
    gamma_formula = dots$gamma %||% ~1,
    mixture = family$params$mixture %||% "poisson", K_max = family$params$K_max)
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "dyn_abun"), priors = priors,
         approx = approx, correction = correction,
         K.max = family$params$K_max, mixture = family$params$mixture),
    control
  ))
}

.dispatch_fp_occu <- function(formula, data, family, detection, y, visits,
                              engine, priors, control,
                              approx = "gaussian_laplace",
                              correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("fp_occu() requires a `detection` formula (the true-detection p11 ",
         "model).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("fp_occu() requires `y` (an n_sites x J matrix of detection states in ",
         "{0, 1, 2}: 0 none, 1 ambiguous, 2 certain).", call. = FALSE)
  }
  model <- .tobs_build_fp_occu(
    occ_formula = formula, det_formula = detection, data = data, y = y,
    fp_formula = dots$p10 %||% ~1, b_formula = dots$certainty %||% ~1)
  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "fp_occu"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
}

.dispatch_ms_abun <- function(formula, data, family, detection, y, visits,
                              engine, priors, control,
                              approx = "gaussian_laplace",
                              correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_abun() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_abun() requires `y` (a 3D array [n_sites x max_visits x ",
         "n_species] or a named list of count matrices).", call. = FALSE)
  }
  if (is.null(dots$species)) {
    stop("ms_abun() requires a `species` argument (the species labels).",
         call. = FALSE)
  }
  model <- .tobs_build_ms_abun(
    abund_formula = formula, det_formula = detection,
    data = data, y = y, species = dots$species,
    det_visit_formula = dots$det_visit_formula,
    det_visit_data    = dots$det_visit_data)

  # NUTS (method = "nuts", tulpaObs#14): sample the exact joint posterior of the
  # non-spatial community N-mixture (community means, per-species deviations, and
  # community covariances) via the in-tree C++ FullGradFn over the closed-form
  # per-(species, site) marginal (R/ms_abun_nuts.R, src/ms_abun_nuts.cpp),
  # warm-started at the Laplace-EM mode. Spatial / temporal / RE terms are not yet
  # wired on the sampler.
  if (identical(engine, "nuts")) {
    structs <- .tobs_structures_from_model(model)
    if (!is.null(structs$temporal) || !is.null(structs$re) ||
        !is.null(structs$svc) || !is.null(structs$latent)) {
      stop("method = \"nuts\" for ms_abun() is the community N-mixture sampler ",
           "(optionally with a shared areal field on the abundance arm); a ",
           "temporal / random-effect / svc / latent term is not yet wired on the ",
           "sampler. Use method = \"laplace\".", call. = FALSE)
    }
    # A shared areal field on the abundance formula joins the community sampler as
    # a fixed-hyper non-centered field (proper-CAR only; tulpaObs#73). icar/bym2
    # stay on nested_laplace (their intrinsic field needs a sum-to-zero reparam).
    if (!is.null(structs$spatial)) {
      if (isTRUE(structs$spatial$shared[2L])) {
        stop("ms_abun() areal field sits on the abundance arm only; a field on ",
             "the detection formula is not supported.", call. = FALSE)
      }
      return(.tobs_fit_ms_abun_nuts_spatial(
        model, spatial = structs$spatial,
        mixture       = family$params$mixture %||% "poisson",
        K_max         = family$params$K_max,
        sigma.beta    = control[["sigma.beta"]] %||% 10,
        sigma.logr    = 1.5,
        n.iter        = as.integer(control[["n.iter"]]   %||% 1000L),
        n.warmup      = as.integer(control[["n.warmup"]] %||% 1000L),
        n.chains      = as.integer(control[["n.chains"]] %||% 1L),
        max.treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
        adapt.delta   = control[["adapt.delta"]] %||% 0.9,
        seed          = as.integer(control[["seed"]] %||% 1L),
        max.iter      = as.integer(control[["max.iter"]] %||% 100L),
        verbose       = isTRUE(control[["verbose"]])))
    }
    return(.tobs_fit_ms_abun_nuts(
      model, mixture = family$params$mixture %||% "poisson",
      K_max         = family$params$K_max,
      sigma.beta    = control[["sigma.beta"]] %||% 10,
      n.iter        = as.integer(control[["n.iter"]]   %||% 1000L),
      n.warmup      = as.integer(control[["n.warmup"]] %||% 1000L),
      n.chains      = as.integer(control[["n.chains"]] %||% 1L),
      max.treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
      adapt.delta   = control[["adapt.delta"]] %||% 0.9,
      seed          = as.integer(control[["seed"]] %||% 1L),
      verbose       = isTRUE(control[["verbose"]])))
  }

  # A spatial term on the abundance formula (icar() / bym2() / car_proper())
  # routes to the shared-field community N-mixture (gcol33/tulpaObs#12). The fit
  # is driven directly by the in-tree community-spatial C++ grid driver, so it
  # reads the dotted controls here (mirroring the non-spatial community path
  # below). control$inner_solver picks the inner method: "em" (default, the
  # closed-form Laplace-EM) or "newton" (the exact-Newton shared-field solve +
  # AGHQ community debias; areal Poisson only).
  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc)) {
    stop("ms_abun() supports a spatial term (icar()/bym2()/car_proper()/spde()) ",
         "and latent() factors on the abundance formula; temporal / re / svc ",
         "terms are not yet wired for the community N-mixture.", call. = FALSE)
  }
  # latent() factors -- residual species co-occurrence on the abundance arm via
  # Q per-site factors with per-species loadings (the lfMsNMix analogue), on
  # their own or alongside a shared field (the spatial-factor case) -- route to
  # the shared block-coordinate latent fitter (R/ms_abun_latent.R, over the
  # generic driver in R/community_latent.R). A plain shared field with NO factors
  # keeps its dedicated in-tree C++ nested Laplace-EM (tulpaObs#12) below: that
  # route is faster and already recovery-tested.
  if (!is.null(structs$latent)) {
    if (isTRUE(structs$latent$shared[2L])) {
      stop("ms_abun() latent() factors sit on the abundance arm only; a latent ",
           "term on the detection formula is not supported.", call. = FALSE)
    }
    if (!is.null(structs$spatial)) {
      if (!identical(engine, "nested_laplace")) {
        stop("ms_abun(): a shared field alongside latent() factors needs ",
             "method = \"nested_laplace\" (drop the field for the factor-only ",
             "community fit).", call. = FALSE)
      }
      if (isTRUE(structs$spatial$shared[2L])) {
        stop("ms_abun() areal field sits on the abundance arm only; a field on ",
             "the detection formula is not supported.", call. = FALSE)
      }
    } else if (identical(engine, "nested_laplace")) {
      stop("ms_abun(): a latent()-factor model with no shared field uses ",
           "method = \"laplace\" (the block-coordinate community Laplace-EM).",
           call. = FALSE)
    }
    return(.tobs_fit_ms_abun_latent(
      model, spatial = structs$spatial, latent = structs$latent,
      mixture    = family$params$mixture %||% "poisson",
      K_max      = family$params$K_max,
      max.iter   = control[["max.iter"]] %||% 100L,
      tol        = control[["tol"]] %||% 1e-4,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      priors     = priors,
      max.outer  = as.integer(control[["max.outer"]] %||% 25L),
      verbose    = isTRUE(control[["verbose"]])))
  }
  if (!is.null(structs$spatial)) {
    return(.tobs_fit_ms_nmix_spatial(
      model, spatial = structs$spatial,
      mixture      = family$params$mixture %||% "poisson",
      K_max        = family$params$K_max,
      max_iter     = control[["max.iter"]] %||% 100L,
      inner_solver = control[["inner_solver"]] %||% "em",
      n_quad       = as.integer(control[["n.quad"]] %||% 1L),
      lkj_eta      = control[["re.lkj"]] %||% 1.5,
      integration  = control[["integration"]] %||% "grid",
      verbose      = isTRUE(control[["verbose"]])))
  }

  # Community N-mixture is fit directly by tulpa over a shared native oracle (no
  # EM-around-INLA path), so it bypasses .tobs_fit_model and reads the dotted
  # controls here. `optimizer` selects the outer driver over that one oracle:
  #   - "em" (default): the fast Laplace-EM (block-Newton mode + closed-form
  #     covariance M-step). Production path -- profiling shows the FD-gradient
  #     joint optimizer dominates the residual runtime, so EM is the default.
  #   - "joint_fd": the finite-difference joint (theta, Sigma) optimizer. Opt-in
  #     for correctness / architecture validation, and the only driver that does
  #     the n.quad > 1 AGHQ variance-component debias. Slower than EM.
  #   - "joint_grad": reserved analytic-gradient extension; errors for now.
  # The AGHQ debias barely moves the community covariances for this family (each
  # species' count marginal is informative, so the per-group Laplace is already
  # accurate), so the EM default loses nothing in practice; n.quad is exposed via
  # joint_fd for sparse / rare-species regimes.
  optimizer <- control[["optimizer"]] %||% "em"
  n_quad    <- as.integer(control[["n.quad"]] %||% 1L)
  .tobs_fit_ms_nmix(
    model,
    mixture   = family$params$mixture %||% "poisson",
    K_max     = family$params$K_max,
    max_iter  = control[["max.iter"]] %||% 100L,
    optimizer = optimizer,
    n_quad    = n_quad,
    lkj_eta   = control[["re.lkj"]] %||% 1.5,
    verbose   = isTRUE(control[["verbose"]]))
}


.dispatch_ms_distance <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_distance() requires a `detection` formula (the log detection-scale ",
         "predictor).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_distance() requires `y` (a 3D array [n_sites x n_bins x ",
         "n_species] or a named list of per-bin count matrices).", call. = FALSE)
  }
  if (is.null(dots$species)) {
    stop("ms_distance() requires a `species` argument (the species labels).",
         call. = FALSE)
  }
  if (is.null(family$params$cutpoints)) {
    stop("ms_distance() requires `cutpoints` (the distance-bin edges) on the ",
         "family, e.g. family = ms_distance(cutpoints = c(0, 25, 50, 100)).",
         call. = FALSE)
  }
  model <- .tobs_build_ms_distance(
    abund_formula = formula, det_formula = detection,
    data = data, y = y, species = dots$species,
    cutpoints  = family$params$cutpoints,
    key        = family$params$key      %||% "halfnorm",
    transect   = family$params$transect %||% "line",
    mixture    = family$params$mixture  %||% "poisson",
    quad_order = as.integer(control[["quad.order"]] %||% 64L))

  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc)) {
    stop("ms_distance() supports a spatial term ",
         "(icar()/car_proper()/bym2()/spde()) and latent() factors on the ",
         "abundance formula; temporal / re / svc terms are not yet wired for ",
         "the community distance family.", call. = FALSE)
  }
  for (nm in c("spatial", "latent")) {
    if (!is.null(structs[[nm]]) && isTRUE(structs[[nm]]$shared[2L])) {
      stop(sprintf(paste0("ms_distance() %s() sits on the abundance arm only; ",
                          "a %s term on the detection formula is not ",
                          "supported."), nm, nm), call. = FALSE)
    }
  }
  # A shared field needs nested_laplace; factors alone are the plain
  # block-coordinate Laplace-EM, as for every other community family.
  if (!is.null(structs$spatial)) {
    if (!identical(engine, "nested_laplace")) {
      stop("ms_distance(): a shared field needs method = \"nested_laplace\" ",
           "(drop the field for the non-spatial community fit).", call. = FALSE)
    }
  } else if (identical(engine, "nested_laplace")) {
    stop("ms_distance(): method = \"nested_laplace\" needs a shared field ",
         "(icar()/car_proper()/bym2()/spde()) on the abundance formula. For the ",
         "non-spatial community fit use method = \"laplace\".", call. = FALSE)
  }
  .tobs_fit_ms_distance(
    model, spatial = structs$spatial, latent = structs$latent,
    mixture    = family$params$mixture %||% "poisson",
    K_max      = family$params$K_max,
    max.iter   = control[["max.iter"]] %||% 100L,
    tol        = control[["tol"]] %||% 1e-4,
    sigma.beta = control[["sigma.beta"]] %||% 5,
    priors     = priors,
    max.outer  = as.integer(control[["max.outer"]] %||% 25L),
    verbose    = isTRUE(control[["verbose"]]))
}


.dispatch_ms_occu_cover <- function(formula, data, family, detection, y, visits,
                                    engine, priors, control,
                                    approx = "gaussian_laplace",
                                    correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_occu_cover() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_occu_cover() requires `y` (a 3D array [n_sites x max_visits x ",
         "n_species] or a named list of detection-history matrices).",
         call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("ms_occu_cover() requires `y_pos` (a 3D array / list matching `y`; ",
         "values used only where y == 1).", call. = FALSE)
  }
  if (is.null(dots$species)) {
    stop("ms_occu_cover() requires a `species` argument (the species labels).",
         call. = FALSE)
  }

  pos_formula <- dots$positive %||% detection

  # A single icar() shared field on the occupancy arm routes to the K=1
  # reduced-rank spatial-factor community fitter (Laplace-EM): per-species
  # loadings on one ICAR field give correlated, range-shifted occupancy maps
  # the plain per-species RE cannot produce (gcol33/tulpa#67). Detected on the
  # user formulas before visit-normalization (which would move covariates off
  # the formula). The detector also gates the remaining structured terms
  # (errors), so the non-spatial path below is reached only when no arm carries
  # one.
  sp <- .tobs_ms_ocs_spatial_request(formula, detection, pos_formula, data)

  # Site / visit dimensions, robust to y being a 3D array or a list of matrices.
  if (is.list(y) && !is.array(y)) {
    n_sites <- nrow(y[[1L]]); max_visits <- ncol(y[[1L]])
  } else {
    n_sites <- dim(y)[1L]; max_visits <- dim(y)[2L]
  }

  # A cover-arm spatial factor carries an icar() term on the cover formula; strip
  # it to the fixed-effects part before visit-normalization (the field is wired
  # separately via the shared-graph builder, not as a cover covariate).
  pos_formula_eff <- if (!is.null(sp) && !is.null(sp$fe_pos)) sp$fe_pos else pos_formula

  vd_det <- .normalize_visits(visits, detection,
                              n_sites = n_sites, max_visits = max_visits)
  vd_pos <- .normalize_visits(visits, pos_formula_eff,
                              n_sites = n_sites, max_visits = max_visits)

  if (!is.null(sp)) {
    nf     <- control[["n.factors"]] %||% 1L
    auto_K <- is.character(nf) && identical(tolower(nf), "auto")
    use_nuts <- identical(engine, "nuts")
    model <- .tobs_build_ms_occu_cover_spatial(
      occ_formula      = sp$fe_occ,
      det_formula      = vd_det$det_formula,
      pos_formula      = vd_pos$det_formula,
      data             = data,
      y                = y,
      y_pos            = dots$y_pos,
      positive         = family$params$positive,
      species          = dots$species,
      adj              = sp$graph,
      K                = if (auto_K) 1L else as.integer(nf),
      cover_factor     = isTRUE(sp$cover_factor),
      field_type       = sp$field_type %||% "icar",
      det_visit_formula = vd_det$det_visit_formula,
      det_visit_data    = vd_det$visits,
      pos_visit_formula = vd_pos$det_visit_formula,
      pos_visit_data    = vd_pos$visits
    )
    if (auto_K && use_nuts) {
      stop("method = \"nuts\" needs an explicit n.factors (K); the auto-K ",
           "rank selection is a Laplace-evidence procedure. Fit the chosen K ",
           "with method = \"nuts\".", call. = FALSE)
    }
    if (auto_K) {
      # Choose the latent-factor rank K by the empirical-Bayes Laplace marginal
      # likelihood: latent-level pointwise criteria (held-out cells, WAIC) fail
      # because each ICAR field adds ~N effective latent parameters, so they
      # track field dimension rather than rank. .ms_ocs_select_K fits the
      # identified (constrained) K-ladder, integrates the field out so its prior
      # supplies the Occam penalty, and returns the argmax-evidence fit; the
      # per-K evidence table is attached as fit$spatial$K_selection.
      sel <- .ms_ocs_select_K(
        model,
        K.max      = as.integer(control[["n.factors.max"]] %||% 4L),
        sd_L       = control[["sd.load"]]    %||% 1.0,
        sigma.beta = control[["sigma.beta"]] %||% 5,
        max.em     = control[["max.iter"]]   %||% 30L,
        tol        = control[["tol"]]        %||% 1e-3,
        verbose    = isTRUE(control[["verbose"]]))
      model$K <- sel$K
      out <- build_ms_occu_cover_spatial_fit(model, sel$fit)
      out$spatial$K_selection <- sel$table
      return(out)
    }
    # K > 1 is identified only up to an orthogonal rotation, so the unconstrained
    # loading posterior is improper along that manifold; the identified
    # (lower-triangular) parameterisation gives well-posed loading uncertainty.
    # K = 1 carries no continuous rotation, so it stays unconstrained.
    constrain <- control[["constrain"]] %||% (model$K > 1L)
    if (use_nuts) {
      fit <- .tobs_fit_ms_occu_cover_spatial_nuts(
        model, sd_L = control[["sd.load"]] %||% 1.0,
        sigma.beta = control[["sigma.beta"]] %||% 5,
        constrain = constrain, control = control)
      out <- build_ms_occu_cover_spatial_fit(model, fit)
      # Surface the sampler diagnostics at the top level (the shared builder
      # stamps NA placeholders); the raw draws stay under fit$nuts.
      nd <- fit$nuts
      out$accept_prob <- nd$accept_prob
      out$divergent   <- nd$divergent
      out$treedepth   <- nd$treedepth
      out$epsilon     <- nd$epsilon
      out$n_samples   <- nrow(nd$draws)
      out$n_chains    <- nd$n_chains
      out$max_rhat    <- nd$max_rhat
      out$min_ess     <- nd$min_ess
      out$nuts        <- nd
      return(out)
    }
    fit <- .tobs_fit_ms_occu_cover_spatial(
      model,
      sd_L       = control[["sd.load"]]   %||% 1.0,
      max.em     = control[["max.iter"]]  %||% 30L,
      tol        = control[["tol"]]       %||% 1e-3,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      constrain  = constrain,
      verbose    = isTRUE(control[["verbose"]])
    )
    return(build_ms_occu_cover_spatial_fit(model, fit))
  }
  if (identical(engine, "nuts")) {
    stop("method = \"nuts\" for ms_occu_cover() requires a shared spatial-factor ",
         "term (icar()/car_proper()/bym2() on the occupancy arm); the ",
         "non-spatial community fit is Laplace only.", call. = FALSE)
  }

  model <- .tobs_build_ms_occu_cover(
    occ_formula      = formula,
    det_formula      = vd_det$det_formula,
    pos_formula      = vd_pos$det_formula,
    data             = data,
    y                = y,
    y_pos            = dots$y_pos,
    positive         = family$params$positive,
    species          = dots$species,
    det_visit_formula = vd_det$det_visit_formula,
    det_visit_data    = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula,
    pos_visit_data    = vd_pos$visits
  )

  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_occu_cover, fit_args)
}


.dispatch_ms_dyn_occu <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_dyn_occu() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_dyn_occu() requires `y` (a 4D array [n_sites x max_visits x ",
         "n_seasons x n_species] or a named list of 3D arrays).", call. = FALSE)
  }
  # Resolve structured terms: a shared areal field on the first-season occupancy
  # formula routes to the community dynamic-spatial fitter (stMsPGOcc /
  # svcTMsPGOcc; gcol33/tulpaObs#123). Bind the occupancy AND detection formulas
  # so a field on either arm is parsed (the FE part is what model.matrix sees);
  # `shared[2]` then identifies a detection-arm field to reject.
  bind <- .tobs_bind_formulas(list(psi1 = formula, p = detection), data)
  model <- .tobs_build_ms_dyn_occu(
    occ_formula = bind$fe$psi1, det_formula = bind$fe$p,
    col_formula = dots$colonization, ext_formula = dots$extinction,
    data = data, y = y, species = dots$species,
    structured_terms = bind$terms)
  structs <- .tobs_structures_from_model(model)

  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc) || !is.null(structs$latent)) {
    stop("ms_dyn_occu(): temporal / re / svc / latent terms are not wired for ",
         "the community dynamic occupancy family; a shared areal field ",
         "(icar()) on the first-season occupancy formula is the structured ",
         "term supported (gcol33/tulpaObs#123).", call. = FALSE)
  }

  if (!is.null(structs$spatial)) {
    if (!identical(engine, "nested_laplace")) {
      stop("a shared areal field on the ms_dyn_occu() occupancy formula needs ",
           "method = \"nested_laplace\" (drop the icar() term for the ",
           "non-spatial community dynamic fit).", call. = FALSE)
    }
    if (isTRUE(structs$spatial$shared[2L])) {
      stop("ms_dyn_occu() areal field sits on the first-season occupancy arm ",
           "only; a field on the detection formula is not supported.",
           call. = FALSE)
    }
    return(.tobs_fit_ms_dyn_occu_field(
      model, spatial = structs$spatial,
      max.iter   = control[["max.iter"]] %||% 200L,
      tol        = control[["tol"]] %||% 1e-4,
      sigma.beta = control[["sigma.beta"]] %||% 5,
      priors     = priors,
      max.outer  = control[["max.outer"]] %||% 25L,
      verbose    = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nested_laplace")) {
    stop("ms_dyn_occu(): method = \"nested_laplace\" needs a shared areal field ",
         "icar() on the first-season occupancy formula. For the non-spatial ",
         "community dynamic fit use method = \"laplace\".", call. = FALSE)
  }

  # Polya-Gamma Gibbs (method = "pg_gibbs", spOccupancy tMsPGOcc; tulpaObs#115,
  # #126): the community PG machinery (msPGOcc) + a 2-state HMM FFBS latent step,
  # giving a calibrated community-variance posterior (vs the attenuated
  # Laplace-EM). Constant transitions, site-level detection, no structured terms.
  if (identical(engine, "pg_gibbs")) {
    return(.tobs_fit_ms_dyn_occu_pg_gibbs(
      model, priors = priors,
      sigma.beta = control[["sigma.beta"]] %||% 2.5,
      n.iter   = as.integer(control[["n.iter"]]   %||% 3000L),
      n.warmup = as.integer(control[["n.warmup"]] %||% 1500L),
      n.chains = max(as.integer(control[["n.chains"]] %||% 2L), 2L),
      n.thin   = as.integer(control[["n.thin"]]   %||% 1L),
      seed     = as.integer(control[["seed"]]     %||% 1L),
      verbose  = isTRUE(control[["verbose"]])))
  }

  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_dyn_occu, fit_args)
}


.dispatch_ms_int_occu <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("ms_int_occu() requires a `detection` formula (or a list of formulas, ",
         "one per source).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("ms_int_occu() requires `y` (a list of per-source 3D arrays ",
         "[n_sites x J_d x n_species]).", call. = FALSE)
  }
  model <- .tobs_build_ms_int_occu(
    occ_formula = formula, det_formula = detection,
    data = data, y = y, species = dots$species, site_map = dots$site_map)
  fit_args <- c(list(model = model, priors = priors), control)
  do.call(.tobs_fit_ms_int_occu, fit_args)
}


