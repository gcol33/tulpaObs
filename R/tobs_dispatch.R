# ---------------------------------------------------------------------------
# Internal dispatchers — one per working family. Each wraps the internal
# model builder + fitter for that family's structure.
# ---------------------------------------------------------------------------

# Resolve an extra-arm formula from `...`, accepting the current bare
# process/symbol name and the deprecated `<name>_formula` spelling. The old
# name still works but emits a one-time deprecation warning pointing to the new
# one. Returns `default` when neither is supplied.
.tobs_arm_formula <- function(dots, new, old, default = ~1) {
  if (!is.null(dots[[new]])) return(dots[[new]])
  if (!is.null(dots[[old]])) {
    .Deprecated(new = new, old = old,
                msg = sprintf(
                  "The `%s` argument of tobs() is deprecated; use `%s` instead.",
                  old, new))
    return(dots[[old]])
  }
  default
}

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
  col_f <- .tobs_arm_formula(dots, "colonization", "col_formula", default = NULL)
  ext_f <- .tobs_arm_formula(dots, "extinction", "ext_formula", default = NULL)
  if (is.null(col_f)) {
    stop("dyn_occu() requires a `colonization = ~ ...` argument.", call. = FALSE)
  }
  if (is.null(ext_f)) {
    stop("dyn_occu() requires an `extinction = ~ ...` argument.", call. = FALSE)
  }
  model <- .tobs_build_model(
    occ_formula  = formula,
    det_formula  = detection,
    data         = data,
    y            = y,
    col_formula  = col_f,
    ext_formula  = ext_f
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
    return(.tobs_fit_ms_occu_spatial(
      model, spatial = structs$spatial,
      max.iter = control[["max.iter"]] %||% 100L,
      verbose  = isTRUE(control[["verbose"]])))
  }

  # Any other structured term (temporal / re / svc / latent) is not wired.
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc) || !is.null(structs$latent)) {
    stop("ms_occu() supports a shared areal field (icar()/bym2()/car_proper()) ",
         "on the occupancy formula under method = \"nested_laplace\"; temporal / ",
         "re / svc / latent terms are not wired.", call. = FALSE)
  }
  if (identical(engine, "nested_laplace")) {
    stop("method = \"nested_laplace\" for ms_occu() needs a shared areal field ",
         "(icar()/bym2()/car_proper()) on the occupancy formula. For the ",
         "non-spatial community fit use method = \"laplace\".", call. = FALSE)
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
  model <- .tobs_build_model(
    occ_formula = formula,
    data        = data,
    y           = y,
    jsdm        = TRUE,
    species     = dots$species
  )

  # A shared AREAL field (icar()/bym2()/car_proper()) on the occupancy formula
  # routes to the nested-Laplace JSDM fitter (the JSDM analogue of the
  # community-occupancy areal field; tulpaObs#76): one shared field on the latent
  # occupancy, atop the shared fixed effects and the per-species random intercept.
  # Continuous-mesh (spde()/gp()) terms fall through to the existing EM-Laplace
  # state field; only areal terms have a `$graph`.
  structs <- .tobs_structures_from_model(model)
  spatial_areal <- !is.null(structs$spatial) &&
    (structs$spatial$type %in% c("icar", "bym2", "car", "car_proper"))
  if (spatial_areal) {
    if (identical(engine, "nuts")) {
      stop("method = \"nuts\" for jsdm() is the non-spatial / continuous-mesh ",
           "sampler; a shared areal field (icar()/bym2()/car_proper()) uses ",
           "method = \"nested_laplace\".", call. = FALSE)
    }
    if (!identical(engine, "nested_laplace")) {
      stop("a shared areal field on the jsdm() formula needs ",
           "method = \"nested_laplace\". For the non-spatial JSDM use ",
           "method = \"laplace\".", call. = FALSE)
    }
    return(.tobs_fit_jsdm_spatial(
      model, spatial = structs$spatial,
      max.iter = control[["max.iter"]] %||% 100L,
      verbose  = isTRUE(control[["verbose"]])))
  }
  if (identical(engine, "nested_laplace")) {
    stop("method = \"nested_laplace\" for jsdm() needs a shared areal field ",
         "(icar()/bym2()/car_proper()) on the occupancy formula. For the ",
         "non-spatial JSDM use method = \"laplace\".", call. = FALSE)
  }

  do.call(.tobs_fit_model, c(
    list(model = model,
         method = .map_engine(engine, family = "jsdm"), priors = priors,
         approx = approx, correction = correction),
    control
  ))
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
    omega_formula = .tobs_arm_formula(dots, "omega", "omega_formula"),
    gamma_formula = .tobs_arm_formula(dots, "gamma", "gamma_formula"),
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
    fp_formula = .tobs_arm_formula(dots, "p10", "fp_formula"),
    b_formula = .tobs_arm_formula(dots, "certainty", "b_formula"))
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
      !is.null(structs$svc) || !is.null(structs$latent)) {
    stop("ms_abun() supports a spatial term (icar()/bym2()/car_proper()) on the ",
         "abundance formula; temporal / re / svc / latent terms are not yet ",
         "wired for the community N-mixture.", call. = FALSE)
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
  model <- .tobs_build_ms_dyn_occu(
    occ_formula = formula, det_formula = detection,
    col_formula = .tobs_arm_formula(dots, "colonization", "col_formula", NULL),
    ext_formula = .tobs_arm_formula(dots, "extinction", "ext_formula", NULL),
    data = data, y = y, species = dots$species)
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


