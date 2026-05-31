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
#'   The spatial fields also have a single-verb umbrella `spatial(...,
#'   model = ...)` that selects the field type by name, mirroring
#'   `temporal(time, type = ...)` and `INLA`'s `f(i, model = ...)`:
#'   `spatial(graph = adj, model = "bym2")` is `bym2(graph = adj)` and
#'   `spatial(lon, lat, model = "spde")` is `spde(lon, lat)`. `model` is one
#'   of `"icar"`, `"bym2"`, `"car"`, `"car_proper"`, `"gp"`, `"multiscale_gp"`,
#'   `"spde"`; per-model arguments pass through unchanged.
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
#'     multiple-imputation correction. The fixed-effect prior threads into the
#'     correction refits (gcol33/tulpa#27), so these use the same
#'     weakly-informative default prior as `"laplace"`; pass `priors = FALSE`
#'     for the unpenalised correction.
#'   * `"nested_laplace"` — multi-block nested Laplace (single-season
#'     occupancy and cover-hurdle joint).
#'   * `"nested_laplace_sla"` — nested Laplace with skew-corrected marginals.
#'   * `"nuts"` — HMC / NUTS sampler (every structure; reports Rhat / ESS).
#'   Not every method is available for every family (e.g. the cover hurdle has
#'   no `"nuts"` path; `"nested_laplace"` is occupancy- and cover-only). An
#'   unsupported method errors with the list of methods that family supports.
#' @param priors optional prior specification. For occupancy families fit
#'   with a Laplace method (`method = "laplace"`, `"laplace_sla"`,
#'   `"nested_laplace"`), pass a list or [occu_priors()] object to set
#'   weakly-informative quadratic priors on the fixed-effect coefficients
#'   (defaults pull the detection intercept toward `p = 0.5` and break the
#'   psi-p identifiability ridge at small `J`). Pass `priors = FALSE` to
#'   disable the default prior and recover the unpenalised MAP. The
#'   `"laplace_gibbs"` / `"laplace_mi"` routes apply the same default prior and
#'   thread it through the correction refits (gcol33/tulpa#27). For NUTS, this
#'   is forwarded to the underlying tulpa engine.
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
#'   * `n.seeds` — number of seed-offset refits to fit and LOO-stack into a
#'     `tobs_stack` ensemble (default 1, a single fit). Member `k` uses base
#'     seed `seed + k - 1`. Only meaningful for the stochastic routes
#'     (`"nuts"`, `"laplace_gibbs"`, `"laplace_mi"`); the deterministic Laplace
#'     methods reject it. Seed-variants are statistically identical, so their
#'     stacking weights come out roughly uniform (this is a Monte-Carlo
#'     robustness device) -- pass distinct fits to [tobs_stack()] for a genuine
#'     model average.
#'
#'   Laplace controls (`method = "laplace"` / `"laplace_sla"` /
#'   `"nested_laplace"`): `max.iter`, `tol`, `damping`, `sigma.beta`.
#'   * `re.aghq` — for a formula random effect under `method = "laplace"`, run
#'     the adaptive Gauss-Hermite debias of the variance components after the
#'     EM converges (default `TRUE`). Removes the Laplace small-cluster
#'     attenuation of `sigma` / the RE correlation for binary occupancy; set
#'     `FALSE` for the raw EM (Laplace, `nAGQ = 1`) fit.
#'   * `n.quad` — quadrature points per random-effect dimension for `re.aghq`
#'     (default 9). `n.quad = 1` is the plain Laplace (`nAGQ = 1`) marginal;
#'     higher values refine it toward the exact marginal.
#'   * `re.lkj` — LKJ shape (`eta`) regularizing a *correlated* random slope's
#'     correlation in the `re.aghq` refine (default 1.5). Pulls a
#'     weakly-identified RE correlation off the `+-1` boundary toward 0 without
#'     touching the marginal SDs; `re.lkj = 1` disables it (uniform). No effect
#'     on intercept / uncorrelated terms.
#'   * `inner_solver` — for a spatial community N-mixture (`ms_abun()` with an
#'     `icar()` / `bym2()` / `car_proper()` field on the abundance arm), the
#'     inner solver integrating the shared field given the community: `"em"`
#'     (default) the closed-form Laplace-EM M-step, or `"newton"` the exact-
#'     Newton shared-field solve alternated with a tulpa AGHQ community debias.
#'     Both integrate the field hyperparameter on the outer grid and return the
#'     same fit object; `"newton"` is Poisson- and areal-only, and markedly
#'     slower (an FD-gradient profile loop per grid node) -- an accuracy /
#'     validation alternative, not the production default.
#'   Stochastic-correction controls (`"laplace_gibbs"` / `"laplace_mi"`):
#'   `n.gibbs` / `n.imputations` (Rubin-pooled draw count) and `seed` (stored
#'   on `$seed`).
#'
#'   Control names are validated against the chosen `method`: passing a
#'   sampler control (e.g. `n.chains`) to a Laplace method, a Laplace control
#'   (e.g. `max.iter`) to `"nuts"`, `seed` to a deterministic route, or an
#'   unrecognized name raises an error rather than being silently ignored.
#' @param ... family-specific named arguments forwarded to the underlying
#'   engine builder.
#'
#' @return An object of class `c("tobs_fit", "<family>_fit", "tulpa_fit")`.
#'   When `control$n.seeds > 1`, a `tobs_stack` ensemble of the seed-offset
#'   refits is returned instead (see [tobs_stack()]).
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

  # Reject control options that the resolved method does not use, rather than
  # silently swallowing them via the splat into `.tobs_fit_model()`'s `...`.
  .tobs_validate_control(control, route, family)

  # Outer-grid progress (tulpaObs#25). Surface control$progress[.every/.throttle/
  # .file] to every nested-Laplace grid below -- tulpa's nested fitters and the
  # tulpaObs N-mixture spatial fitters read the scoped `tulpa.nl_progress` option.
  # Restored on exit so it never leaks past this fit.
  .op_nl_progress <- options(tulpa.nl_progress = .tobs_progress_opt(control))
  on.exit(options(.op_nl_progress), add = TRUE)

  # Gibbs / MI corrections thread the fixed-effect prior through their refits
  # (gcol33/tulpa#27), so the `"laplace_gibbs"` / `"laplace_mi"` routes use the
  # same weakly-informative default prior as `"laplace"`. Pass `priors = FALSE`
  # to recover the unpenalised correction.

  if (family$status == "planned") {
    .stop_planned_family(family)
  }

  # Backend coverage is enforced centrally: each working family declares the
  # methods it actually supports (`.tobs_family_methods`). Reject an unsupported
  # method with a pointer to the supported set, rather than silently downgrading
  # the engine (e.g. nested_laplace -> single-Laplace) and then mislabelling
  # `fit$method`.
  .tobs_validate_family_method(route$method, family)

  # n.seeds > 1: fit K members under offset RNG seeds and LOO-stack them (see
  # `tobs_stack()`). Control validation has already rejected `n.seeds` on the
  # deterministic Laplace routes, where seed-variants would be identical.
  n_seeds <- control[["n.seeds"]]
  if (!is.null(n_seeds) && as.integer(n_seeds) > 1L) {
    K           <- as.integer(n_seeds)
    base_seed   <- as.integer(control[["seed"]] %||% 42L)
    member_ctrl <- control
    member_ctrl[["n.seeds"]] <- NULL
    members <- lapply(seq_len(K), function(i) {
      ci <- member_ctrl
      ci[["seed"]] <- base_seed + i - 1L
      tobs(formula = formula, data = data, family = family,
           detection = detection, y = y, visits = visits,
           method = method, priors = priors, control = ci, ...)
    })
    names(members) <- paste0("seed", base_seed + seq_len(K) - 1L)
    return(tobs_stack(members))
  }
  control[["n.seeds"]] <- NULL   # orchestration knob, not a fitter argument

  dispatch <- switch(
    family$name,
    occu     = .dispatch_occu,
    dyn_occu = .dispatch_dyn_occu,
    ms_occu  = .dispatch_ms_occu,
    int_occu = .dispatch_int_occu,
    jsdm     = .dispatch_jsdm,
    abun     = .dispatch_abun,
    ms_abun  = .dispatch_ms_abun,
    cover    = .dispatch_cover,
    occu_cover = .dispatch_occu_cover,
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


# Normalize the outer-grid progress knobs from a `control` list into the scoped
# `tulpa.nl_progress` option value read by tulpa's nested-Laplace fitters and the
# tulpaObs N-mixture spatial fitters (tulpaObs#25 / gcol33/tulpa#45). Off by
# default; `progress = TRUE` turns on the flushed cell-k/n_grid + ETA reporter,
# `progress.file` adds a heartbeat file for detached runs.
.tobs_progress_opt <- function(control) {
  list(
    progress          = isTRUE(control$progress),
    progress_every    = as.integer(control$progress.every    %||% 0L),
    progress_throttle = as.numeric(control$progress.throttle %||% 2),
    progress_file     = as.character(control$progress.file    %||% "")
  )
}

# The four cpp-side progress arguments for the tulpaObs N-mixture spatial
# entries (which run their own outer grids, not tulpa's driver). Reads the same
# scoped `tulpa.nl_progress` option `tobs()` sets, so a single control surfaces
# to every spatial backend.
.tobs_nl_progress_cpp <- function() {
  p <- getOption("tulpa.nl_progress", NULL)
  if (!is.list(p)) p <- .tobs_progress_opt(list())
  list(progress          = isTRUE(p$progress),
       progress_every    = as.integer(p$progress_every),
       progress_throttle = as.numeric(p$progress_throttle),
       progress_file     = as.character(p$progress_file))
}

# Call an N-mixture spatial cpp entry with the four progress arguments appended
# from the scoped option. One injection point for every nmix spatial backend.
.cpp_nmix_progress <- function(.fn, ...) {
  do.call(.fn, c(list(...), .tobs_nl_progress_cpp()))
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
  route$method <- method   # resolved public name (auto -> concrete) for messages
  route
}

# ---------------------------------------------------------------------------
# Per-family backend coverage
#
# Single source of truth for which `method` each working family supports.
# `tobs()` validates the resolved (auto -> concrete) method against this set
# and errors with a pointer to the supported methods, rather than silently
# downgrading the engine (the old `.map_engine()` nested_laplace -> single-
# Laplace fall-back, which also mislabelled `fit$method`) or scattering the
# rejection across each family's dispatcher (the cover hurdle's bespoke
# `stop()`s). Gating mirrors the (engine, approx, correction) architecture: a
# method is listed iff its engine has a real execution path for the family.
#
#   * nested_laplace -- the nested-Laplace engine assembles a multi-block latent
#     prior (spatial / temporal / iid) and routes the state ("occ") M-step block
#     through `tulpa::tulpa_nested_laplace()`. Wired for single-season,
#     integrated, community, and dynamic occupancy (`.tobs_em_nested_laplace()`,
#     which shares the per-model-type callbacks with the Laplace path) and for
#     the cover hurdle's joint path (`tulpa_nested_laplace_joint()`). It also
#     supports INLA-style NA-response prediction (held-out sites), so the latent
#     field interpolates occupancy at unsurveyed sites.
#   * nested_laplace_sla -- the skew correction on the nested path is wired for
#     single-season occupancy and the cover hurdle only.
#   * laplace / laplace_sla / laplace_gibbs / laplace_mi -- run on tulpa's
#     EM+Laplace engine, which has callbacks for every occupancy family. The
#     cover hurdle is fit by a separate two-Laplace dispatcher with no EM
#     correction engine, so it offers laplace / laplace_sla only (no gibbs/mi).
#   * nuts -- the C++ sampler covers single / dynamic / community / integrated /
#     jsdm; the cover hurdle has no HMC likelihood yet.
#
# Planned families (status == "planned") have no entry and error earlier via
# `.stop_planned_family()`; the validator is a no-op for them.
.tobs_family_methods <- list(
  occu     = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nested_laplace_sla", "nuts"),
  dyn_occu = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nuts"),
  ms_occu  = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nuts"),
  int_occu = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nuts"),
  jsdm     = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi", "nuts"),
  # abun: non-spatial Poisson N-mixture (laplace) + areal-spatial offset
  # (nested_laplace: icar / bym2 / car_proper on the abundance arm). tulpa's
  # spatial fitters return the grid-integrated coefficient covariance, so the
  # spatial SEs are calibrated (law-of-total-covariance over the hyperparameter
  # grid). negbin and NUTS remain upstream-pending.
  abun     = c("laplace", "nested_laplace"),
  # ms_abun: community / multispecies N-mixture via the in-tree C++ Laplace-EM
  # (per-species coefficient RE with Gaussian community covariances). A shared
  # areal field (icar / bym2 / car_proper) on the abundance arm fits under
  # nested_laplace (gcol33/tulpaObs#12); Poisson or grid-integrated negbin size.
  ms_abun  = c("laplace", "nested_laplace"),
  cover    = c("laplace", "laplace_sla", "nested_laplace", "nested_laplace_sla"),
  # occu_cover: non-spatial Laplace via direct optim on the exact two-state
  # marginal (v1); nested-Laplace adds a cell-level ICAR field shared across
  # psi and cover arms with scaling alpha (v2, the mod.joint analogue with
  # `copy = "cell.occ"`). v2 currently reads bym2() as ICAR (rho fixed to 1);
  # free-rho BYM2 + outer-grid integration of (sigma, alpha) is v3.
  occu_cover = c("laplace", "nested_laplace")
)

# Validate a resolved public method name against the family's supported set.
# `method` is the concrete name ("auto" already resolved upstream). No-op for a
# family with no entry (planned families, which error before reaching here).
.tobs_validate_family_method <- function(method, family) {
  supported <- .tobs_family_methods[[family$name]]
  if (is.null(supported) || method %in% supported) return(invisible(NULL))
  stop(
    sprintf(
      "method = \"%s\" is not available for %s() (%s). Supported: %s.",
      method, family$name, family$class_long,
      paste0("\"", supported, "\"", collapse = ", ")
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Control-option validation
#
# `control` is splatted as named args onto `.tobs_fit_model()`, whose formals
# are the union of every knob across all methods (plus a trailing `...`). That
# means a control that does not apply to the chosen method is silently ignored
# (e.g. `n.chains` under `"laplace"`) and a typo (`niter`) vanishes into `...`.
# We validate names up front against a per-route allowlist so misapplied or
# misspelled controls error instead.
#
# Single source of truth: control names are grouped by capability, and each
# engine/correction route admits a set of groups. `sigma.beta` is shared by the
# Laplace and NUTS paths; `seed` and `n.seeds` by the stochastic-correction and
# NUTS paths (the deterministic Laplace routes reject `n.seeds` here, since
# seed-variant fits would be identical -- see the ensemble branch in tobs()).
# ---------------------------------------------------------------------------
.tobs_control_groups <- list(
  laplace_em = c("max.iter", "tol", "damping", "sigma.beta",
                 "re.aghq", "n.quad", "re.lkj", "optimizer", "hessian",
                 "inner_solver"),
  correction = c("n.gibbs", "n.imputations", "seed", "n.seeds"),
  sampler    = c("n.iter", "n.warmup", "n.thin", "n.chains", "n.threads",
                 "adapt.delta", "max.treedepth", "seed", "sigma.beta",
                 "sigma.re.scale", "n.seeds"),
  universal  = c("verbose",
                 "progress", "progress.every", "progress.throttle",
                 "progress.file")
)

# Capability groups admitted by a resolved (engine, correction) route.
.tobs_control_allow <- function(engine, correction) {
  switch(
    engine,
    laplace        = c("laplace_em", if (correction != "none") "correction"),
    nested_laplace = "laplace_em",
    nuts           = "sampler",
    character(0)
  )
}

# Public method names that accept a given control key (for "wrong method"
# hints). Derived from the route table + allowlist so it stays in sync.
.tobs_methods_for_control <- function(key) {
  in_group <- vapply(.tobs_control_groups, function(g) key %in% g, logical(1))
  groups   <- names(.tobs_control_groups)[in_group]
  if (!length(groups)) return(character(0))
  methods <- names(.tobs_method_table)
  keep <- vapply(methods, function(m) {
    r <- .tobs_method_table[[m]]
    any(.tobs_control_allow(r$engine, r$correction) %in% groups)
  }, logical(1))
  methods[keep]
}

# Validate `control` names against the resolved route. Errors on (a) a known
# control that the chosen method does not use, or (b) an unrecognized name
# (with a fuzzy "did you mean" suggestion). Collects all offenders into one
# message. No-op for a valid (or empty) control list.
.tobs_validate_control <- function(control, route, family = NULL) {
  if (length(control) == 0L) return(invisible(NULL))
  nms <- names(control)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("`control` must be a fully named list, e.g. ",
         "control = list(n.iter = 4000).", call. = FALSE)
  }

  # Family-specific dispatchers (e.g. the cover hurdle) declare extra control
  # names via family$control_keys; admit those alongside the engine controls.
  family_keys <- if (!is.null(family)) family$control_keys %||% character(0)
                 else character(0)

  allowed_groups <- c("universal",
                      .tobs_control_allow(route$engine, route$correction))
  allowed_keys <- unique(c(unlist(.tobs_control_groups[allowed_groups],
                                  use.names = FALSE),
                           family_keys))
  vocabulary   <- unique(c(unlist(.tobs_control_groups, use.names = FALSE),
                           family_keys))

  bad <- setdiff(nms, allowed_keys)
  if (!length(bad)) return(invisible(NULL))

  method <- route$method %||% route$engine
  msgs <- vapply(bad, function(key) {
    if (key %in% vocabulary) {
      uses <- .tobs_methods_for_control(key)
      sprintf("  - '%s' is not used by method = \"%s\"; it applies to %s.",
              key, method,
              paste0("method = \"", uses, "\"", collapse = " / "))
    } else {
      near <- agrep(key, vocabulary, value = TRUE, max.distance = 0.34)
      hint <- if (length(near))
        sprintf(" Did you mean %s?",
                paste0("'", near, "'", collapse = " / ")) else ""
      sprintf("  - '%s' is not a known control option.%s", key, hint)
    }
  }, character(1))

  stop("Invalid `control` option(s) for method = \"", method, "\":\n",
       paste(msgs, collapse = "\n"),
       "\nSee `?tobs` for the controls each method uses.", call. = FALSE)
}

.map_engine <- function(engine, family = NULL) {
  # Engine name translation between the tobs vocabulary and what the underlying
  # fitter currently understands. The nested-Laplace engine (`.tobs_fit_model()`
  # -> `.tobs_em_nested_laplace()`) is wired for single-season, integrated,
  # community, and dynamic occupancy; the per-family method registry
  # (`.tobs_family_methods`) rejects `nested_laplace` for every other family
  # before dispatch, so reaching here with an unsupported family is an internal
  # mis-wire rather than a user error to downgrade silently.
  if (engine == "nested_laplace") {
    if (family %in% c("occu", "int_occu", "ms_occu", "dyn_occu", "abun",
                       "occu_cover")) {
      return("nested_laplace")
    }
    stop(sprintf(
      "Internal error: nested_laplace reached .map_engine for family '%s'; the method registry should have rejected it.",
      family %||% "(unspecified)"), call. = FALSE)
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
    ms_abun  = "Phase 2 (multispecies N-mixture)",
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
    cat(sprintf("  default method : %s\n", fam$default_engine))
    cat("\n")
  }
  model <- x$model
  if (!is.null(model)) {
    if (model$model_type == "single" || model$model_type == "nmix") {
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
  if (!is.null(x$divergent) && isTRUE(sum(x$divergent) > 0)) {
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
        epsilon = "Mean extinction (intercept)",
        lambda  = "Mean abundance (intercept)"
      )
      if (!is.null(label)) {
        cat(sprintf("%s: %.3f\n", label, x$intercepts[[nm]]))
      }
    }
  }
  if (!is.null(x$nmix_dispersion)) {
    d <- x$nmix_dispersion
    if (isTRUE(is.finite(d$r_sd))) {
      cat(sprintf("NB dispersion (size r): %.3f (SE %.3f)\n", d$r, d$r_sd))
    } else {
      cat(sprintf("NB dispersion (size r): %.3f\n", d$r))
    }
  }
  invisible(x)
}
