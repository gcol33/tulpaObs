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
#' @param occurrence state-process (occupancy / latent-presence) formula for the
#'   cover hurdle ([occu_cover()] / [cover()]), reading symmetrically with
#'   `detection` and `positive`. The front-door name for `formula`; supply one
#'   of the two.
#' @param positive positive-arm (cover) formula for [occu_cover()], e.g.
#'   `~ time.sc + habitat`. The occurrence spatial field is carried onto this arm
#'   with a `copy()` selector, the INLA-style cross-arm edge written in the
#'   formula. The selector is a constructor, so no field name is needed in the
#'   common case:
#'
#'   * `copy(spatial(), alpha = grid(g))` copies the occurrence arm's spatial
#'     effect and marginalizes the coupling amplitude over the grid `g`;
#'     `alpha = <scalar>` fixes it.
#'   * `copy(spatial(), terms = list(intercept = grid(g0), time.sc = grid(g1)))`
#'     gives a per-component amplitude, keyed by the field's own block names (the
#'     intercept block and a `||`-declared trend column, or its alias `trend`);
#'     `terms =` must address every block.
#'   * `copy(spatial(cell_idx), ...)` disambiguates by grouping variable when the
#'     occurrence arm carries several spatial effects.
#'   * decoupling is structural: omit `copy()` so the field rides occupancy only
#'     (a block with no `copy()` is pinned at zero coupling), or write a
#'     `spatial()` term for the cover arm's own field.
#'
#'   Defaults to `detection` when unset (a per-visit cover design matching the
#'   detection design).
#' @param formula state-process formula, e.g. `~ elev + forest`. For occupancy
#'   this is the occupancy probability formula; for N-mixture the abundance
#'   formula; for the cover hurdle the latent-presence formula (also given as
#'   `occurrence`, which reads symmetrically with `detection` and `positive`).
#'
#'   Single-vector-response families (the cover hurdle, [cover()]) also accept
#'   the response on the left-hand side, `response ~ predictors`, in which case
#'   `y =` is omitted (e.g. `cover.flat ~ time + habitat`). The LHS is evaluated
#'   against `data` (then the calling environment), so it may be a bare column
#'   or an expression. Matrix / array / list response families ([occu()],
#'   [abun()], the `ms_*` families, ...) keep the one-sided form and supply the
#'   response via `y =`; a two-sided formula for those errors.
#'
#'   Structured effects are written as terms inside the formula, the way
#'   `lme4`, `mgcv`, and `INLA` do: spatial fields `icar(graph = adj)`,
#'   `bym2(graph = adj)`, `gp(lon, lat, prior_range = c(r0, alpha))`,
#'   `spde(lon, lat)`; random effects `re(group)`; temporal fields
#'   `temporal(time)`; spatially varying coefficients `spatial(~ 1 + w || cell,
#'   graph = adj)` over a graph or `spatial(lon, lat, model = "svc",
#'   coefficients = ..., prior_range = c(r0, alpha))` over coordinates;
#'   community latent factors
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
#'   `"spde"`, `"svc"`; per-model arguments pass through unchanged. The
#'   umbrella also reads a varying-coefficient bar
#'   (`spatial(~ 1 + w || cell, graph = adj)`), so both flavours of
#'   spatially-varying coefficient -- areal over a graph, continuous over
#'   coordinates (`model = "svc"`) -- are written through the same verb.
#'
#'   The continuous fields (`gp()`, `svc()`, `spde()`) require `prior_range =
#'   c(r0, alpha)`, a penalized-complexity prior on the spatial range encoding
#'   `P(range < r0) = alpha` (Fuglstad et al. 2019). The range is in the units
#'   of the coordinates -- the kernel is `exp(-d / range)` -- so choose `r0` as
#'   a distance below which the field's correlation would be surprisingly short:
#'   on unit-square coordinates `prior_range = c(0.1, 0.05)` reads as "a 5%
#'   chance the range is under 0.1". There is no default, deliberately. The
#'   range is weakly identified by the likelihood alone, so a default would be
#'   an invented prior doing real work on the posterior rather than a
#'   convenience.
#'
#'   Random effects also accept `lme4` bar syntax as shorthand for `re()`:
#'   `(1 | g)` is `re(g)`, `(x | g)` is a correlated random intercept and
#'   slope `re(g, type = "slope", covariate = x)`, and `(x || g)` drops the
#'   correlation (`correlated = FALSE`). Use `re()` directly for AR1/RW
#'   structures or other options.
#' @param data data frame of site-level covariates with `nrow(data) ==
#'   nrow(y)`, or a `tobs_data` frame (from [tobs_data()] / [tobs_format()])
#'   bundling the response, site covariates, and visit covariates. When a frame
#'   is passed, `y` and `visits` are taken from it and must not also be supplied.
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
#'
#'   For a single-vector-response family ([cover()]) the response may instead
#'   be written on the `formula` left-hand side (`response ~ predictors`), in
#'   which case `y =` is omitted. Supplying the response both on the LHS and via
#'   `y =` errors. Matrix / array / list response families take the response
#'   via `y =` only.
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
#'     correction refits, so these use the same
#'     weakly-informative default prior as `"laplace"`; pass `priors = FALSE`
#'     for the unpenalised correction.
#'   * `"pg_gibbs"` — a Polya-Gamma Gibbs sampler over the exact single-season
#'     occupancy posterior (the spOccupancy `PGOcc` engine). A real MCMC chain
#'     (with `Rhat` / `ESS` diagnostics), distinct from `"laplace_gibbs"` (a
#'     stochastic-EM variance correction). Sampler controls (`n.iter`,
#'     `n.warmup`, `n.chains`, `n.thin`, `seed`, `sigma.beta`). v1: single-season
#'     `occu()`, site-level detection, no structured terms.
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
#'   thread it through the correction refits. For NUTS, this
#'   is forwarded to the underlying tulpa engine.
#' @param control list of low-level engine controls. Names follow the
#'   dotted-separator convention. Sampler controls (`method = "nuts"`):
#'   * `n.iter` — post-warmup sampling iterations kept per chain (default 2000);
#'     the total run per chain is `n.iter + n.warmup`.
#'   * `n.warmup` — warmup / adaptation iterations per chain, discarded
#'     (default 1000).
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
#'   * `integration` — how the in-package spatial / community nested-Laplace
#'     fitters integrate the outer field hyperparameters (`tau`, `rho`, `sigma`,
#'     `range`): `"grid"` (default) a fixed tensor grid, or `"ccd"` a mode-centred
#'     central-composite design placed at the marginal-likelihood mode and scaled
#'     by the outer posterior covariance, with the outer PSIS Pareto-k reported on
#'     `fit$spatial_pareto_k`. CCD declines to the grid when the outer curvature
#'     is ill-conditioned (a weakly-identified axis) or for a single positive
#'     hyperparameter; `fit$spatial_integration` records which ran. Each outer
#'     node is a full inner solve, so `"ccd"` adds a mode-find without a node
#'     saving on these coarse grids -- it is opt-in, most useful when a
#'     multi-axis hyperparameter posterior is well identified.
#'   * `diagnose.k` — for the joint-coupled spatial families (`occu_cover()`,
#'     `occu()` spatial, `occu_multiscale_cover()`), whether to score the outer
#'     hyperparameter Gaussian summary with an importance-sampling Pareto-k.
#'     Defaults `FALSE` (it re-solves the inner Laplace on the full field
#'     `k.samples` times; gcol33/tulpaObs#101). When `TRUE`, the fit carries
#'     `pareto_k` (`< 0.7` = reliable summary), `pareto_k_is_ess` (the
#'     importance-sampling ESS on the PSIS-smoothed weights; `/ k.samples` is the
#'     relative IS efficiency), and
#'     `pareto_k_proposal_source` at the top level (and in [glance()]).
#'     `pareto_k_proposal_source` (gcol33/tulpa#116) is `"mode_hessian"` when the
#'     importance proposal came from the Laplace curvature at the hyperparameter
#'     mode -- curvature-backed, so the k-hat stays trustworthy even when a sharp
#'     posterior collapses the integration grid -- or `"grid_moment"` when it came
#'     from the grid-weighted node covariance, the regime to watch as the grid
#'     concentrates.
#'   Stochastic-correction controls (`"laplace_gibbs"` / `"laplace_mi"`):
#'   `n.gibbs` / `n.imputations` (Rubin-pooled draw count) and `seed` (stored
#'   on `$seed`).
#'
#'   Control names are validated against the chosen `method`: passing a
#'   sampler control (e.g. `n.chains`) to a Laplace method, a Laplace control
#'   (e.g. `max.iter`) to `"nuts"`, `seed` to a deterministic route, or an
#'   unrecognized name raises an error rather than being silently ignored.
#' @param by optional name of a species column for a per-species batched fit
#'   (`occu_cover()` / `cover()` only). When supplied, `data` is a long,
#'   plot-level frame: one row per site-visit (`occu_cover()`) or per site
#'   (`cover()`), with `by` giving the species. `tobs()` splits `data` by that
#'   column, builds each species' response onto one shared site x visit grid
#'   (via [tobs_data()]), fits the B species independently, and returns a
#'   `tobs_batch`. Per-species fits are statistically independent (no pooling) --
#'   identical to fitting each species with a separate single-species `tobs()`
#'   call; for a community model that borrows strength across species use
#'   [ms_occu_cover()]. The long -> response pivot needs the column names: pass
#'   `site = ` and `response = ` (the detection 0/1 column for `occu_cover()`,
#'   the cover column for `cover()`), plus, for `occu_cover()`, `visit = ` (the
#'   replicate column) and `y_pos = ` (the cover column) and any visit-level
#'   `det.covs = `. Site-level covariates (those the `formula` / `detection`
#'   reference at the cell level) are read as the first row per site, the way
#'   [tobs_data()]`(occ.covs = )` does.
#'
#'   The same long-frame contract drives a SINGLE `occu_cover()` fit when `by`
#'   is omitted: pass `site = `, `visit = `, `response = ` (the 0/1 detection
#'   column) and `y_pos = ` (the cover column), plus any visit-level
#'   `det.covs = `, with a long, plot-level `data`, and `tobs()` builds the
#'   paired occurrence / cover arms and the site-level design for you -- the
#'   by= batch reduced to one species, so you no longer hand-build the two
#'   responses and align them. The arms default to the compact (ragged) layout
#'   on the nested-Laplace route (no per-site visit cap); set
#'   `control$compact = FALSE` for the dense grid. See [occu_cover_inputs()] to
#'   build and inspect the arms without fitting.
#'
#'   The long-frame `response =` here names a COLUMN of `data` (the pivot key
#'   holding the detection / cover values). It is a separate argument from the
#'   positive-part distribution the family constructor takes
#'   (`cover(response = "beta")`, `occu_cover(response = "lognormal")`), which
#'   selects the cover-arm likelihood rather than a data column. Both may appear
#'   in one call -- `tobs(long, family = cover(response = "beta"), site = ,
#'   response = "cover.flat", ...)` -- and do not interact.
#' @param ... family-specific named arguments forwarded to the underlying
#'   engine builder.
#'
#' @return An object of class `c("tobs_fit", "<family>_fit", "tulpa_fit")`.
#'   When `control$n.seeds > 1`, a `tobs_stack` ensemble of the seed-offset
#'   refits is returned instead (see [tobs_stack()]). When `by` is supplied, a
#'   `tobs_batch` of per-species fits is returned (see [tobs_get()]).
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
                 occurrence = NULL,
                 detection  = NULL,
                 positive   = NULL,
                 y          = NULL,
                 visits     = NULL,
                 method     = c("auto", "laplace", "laplace_sla",
                                "laplace_gibbs", "laplace_mi", "pg_gibbs",
                                "nested_laplace", "nested_laplace_sla", "nuts"),
                 priors     = NULL,
                 control    = list(),
                 by         = NULL,
                 ...) {

  method <- match.arg(method)

  # A pre-built `tobs_data` frame stands in for the (data, y, visits) triple:
  # unpack it into the same arguments raw inputs use, so a frame routes through
  # one validated entry rather than a parallel pipeline. The frame's site-level
  # `occ.covs` becomes `data`, its response `y`, and its `det.covs` the `visits`.
  if (inherits(data, "tobs_data")) {
    unpacked <- .tobs_unpack_frame(data, y = y, visits = visits)
    data   <- unpacked$data
    y      <- unpacked$y
    visits <- unpacked$visits
  }

  # `occurrence` is the state-process formula's name for the occu_cover() /
  # cover() hurdle, reading symmetrically with `detection` and `positive`. It is
  # the front-door name for `formula`; give one of the two.
  if (!is.null(occurrence)) {
    if (!missing(formula) && !is.null(formula)) {
      stop("Give the state-process formula as `occurrence`, not both ",
           "`occurrence` and `formula`.", call. = FALSE)
    }
    formula <- occurrence
  }
  if (missing(formula) || is.null(formula)) {
    # cover() per-arm mode: `presence =` / `positive =` formulas stand in for a
    # shared state formula (arm = formula). The presence formula seeds the
    # pipeline; the cover encoder reads both per-arm formulas.
    .tobs_dots <- list(...)
    if (inherits(family, "tobs_family") && identical(family$name, "cover") &&
        !is.null(.tobs_dots$presence)) {
      formula <- .tobs_dots$presence
    } else {
      stop("A state-process formula is required (`occurrence =` for the cover ",
           "hurdle, `formula =` otherwise).", call. = FALSE)
    }
  }

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

  # Per-species batched fit (`by = "<species_col>"`). `data` is long / plot-level
  # here; split it by species, build each species' response onto one shared
  # site x visit grid by REUSING the single-species long -> matrix construction
  # (tobs_data()), then route the per-species responses through the same batched-
  # independent driver as the hand-built multi-response `y` list below. Returns a
  # `tobs_batch`. Scoped to occu_cover() / cover() (the families with a per-plot
  # response that a species column splits); errors for any other family.
  # `positive` is a formal of tobs() but the per-species / batch / ensemble
  # sub-fits below thread the cover positive-arm formula through their `dots`
  # (the old `...` spelling). Fold it back in so those paths still receive it.
  fwd_dots <- list(...)
  if (!is.null(positive)) fwd_dots$positive <- positive

  if (!is.null(by)) {
    return(.tobs_fit_by_species(
      formula = formula, data = data, family = family, detection = detection,
      visits = visits, method = method, priors = priors, control = control,
      by = by, dots = fwd_dots))
  }

  # Single occu_cover() fit from a long / plot-level frame (gcol33/tulpaObs#107):
  # the same long-frame contract the by= batch path accepts, minus the species
  # split. The signal is `response = "<col>"` (a column name) on an occu_cover()
  # family handed a data-frame `data`; build the paired occurrence / cover arms
  # and the shared site / visit design with the SAME builder the by= loop uses,
  # then fall through to the normal single-fit dispatch -- so a user no longer
  # hand-rolls tobs_data() twice plus the alignment check. `compact` defaults on
  # for the nested-Laplace route (the joint engine reads the ragged arms with no
  # per-site visit cap) and is overridable via control$compact.
  if (is.null(by) && identical(family$name, "occu_cover") &&
      is.data.frame(data) && !is.null(fwd_dots$response)) {
    if (!is.null(y)) {
      stop("occu_cover(): supply the long-frame `response = ` OR a pre-built ",
           "`y = `, not both.", call. = FALSE)
    }
    eng     <- .tobs_resolve_method(method, family)$engine
    compact <- control[["compact"]] %||% identical(eng, "nested_laplace")
    control[["compact"]] <- NULL        # long-frame build knob, not a fitter control
    pos_type <- .occu_cover_pos_type(family$params$positive)
    arms <- .occu_cover_arms_from_long_call(data, visits, fwd_dots, compact,
                                            pos_type = pos_type)
    data   <- arms$site_data
    y      <- arms$y
    visits <- arms$visits
    fwd_dots$y_pos    <- arms$y_pos
    fwd_dots$response <- NULL
    fwd_dots$site     <- NULL
    fwd_dots$det.covs <- NULL
    fwd_dots$occ.covs <- NULL
    fwd_dots$coords   <- NULL
  }

  # Response on the top formula LHS (gcol33/tulpaObs#66). A single-vector-response
  # family (cover hurdle; family$response == "vector") may carry its response on
  # the formula LHS -- `cover.flat ~ predictors` -- and drop `y =`. Resolve the
  # LHS to `y` and strip `formula` to one-sided here, so the rest of the pipeline
  # (dispatchers, encoders, structured-term parser) sees the unchanged one-sided
  # interface. Matrix / array / list families take their response via `y =` only;
  # a two-sided formula for those errors with a pointer to `?tobs`.
  resolved <- .tobs_resolve_response_lhs(formula, y, family, data)
  formula  <- resolved$formula
  y        <- resolved$y

  # Batched multi-response (gcol33/tulpa#66): occu_cover with `y` a list of >= 2
  # response matrices (or a 3D array [n_sites x max_visits x B]) fits B species,
  # each with the per-species model, and returns a `tobs_batch`. Intercept here
  # so each species replays the full single-species tobs() pipeline below --
  # making every per-species fit byte-identical to an independent fit (the
  # validation oracle for the fused block-diagonal backend; see
  # R/occu_cover_batch.R).
  if (identical(family$name, "occu_cover")) {
    B <- .tobs_multiresponse_n(y)
    if (!is.null(B) && B >= 2L) {
      return(.tobs_fit_occu_cover_batch(
        tobs_args = list(formula = formula, data = data, family = family,
                         detection = detection, visits = visits,
                         method = method, priors = priors, control = control,
                         dots = fwd_dots),
        y = y, B = B
      ))
    }
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
           detection = detection, positive = positive, y = y, visits = visits,
           method = method, priors = priors, control = ci, ...)
    })
    names(members) <- paste0("seed", base_seed + seq_len(K) - 1L)
    return(tobs_stack(members))
  }
  control[["n.seeds"]] <- NULL   # orchestration knob, not a fitter argument

  # Canonical response totals (sites x visits, plus sources for the integrated
  # families), reported once before dispatch under control$verbose and stored on
  # the returned fit. n_sites is anchored to the site-level data the design binds
  # to. Computed only on the single-fit path; the n.seeds / by / multi-response
  # branches above each return their own container, and every leaf fit reports
  # itself when it reaches here.
  dims <- .tobs_input_dims(y)
  if (is.data.frame(data)) dims$n_sites <- nrow(data)
  if (isTRUE(control[["verbose"]])) {
    message(.tobs_input_message(family$name, dims$n_sites,
                                dims$max_visits, dims$n_sources))
  }

  dispatch <- switch(
    family$name,
    occu     = .dispatch_occu,
    dyn_occu = .dispatch_dyn_occu,
    ms_occu  = .dispatch_ms_occu,
    int_occu = .dispatch_int_occu,
    jsdm     = .dispatch_jsdm,
    count    = .dispatch_count,
    ms_count = .dispatch_ms_count,
    abun     = .dispatch_abun,
    royle_nichols = .dispatch_royle_nichols,
    occu_ttd        = .dispatch_occu_ttd,
    occu_multi      = .dispatch_occu_multi,
    double_observer = .dispatch_double_observer,
    gdistremoval    = .dispatch_gdistremoval,
    distsamp_open   = .dispatch_distsamp_open,
    dyn_int_occu    = .dispatch_dyn_int_occu,
    t_occu          = .dispatch_t_occu,
    ms_abun  = .dispatch_ms_abun,
    removal  = .dispatch_removal,
    distance = .dispatch_distance,
    ms_distance = .dispatch_ms_distance,
    fp_occu  = .dispatch_fp_occu,
    dyn_abun = .dispatch_dyn_abun,
    cover    = .dispatch_cover,
    occu_cover = .dispatch_occu_cover,
    occu_multiscale_cover = .dispatch_occu_multiscale_cover,
    ms_occu_cover = .dispatch_ms_occu_cover,
    occu_categorical = .dispatch_occu_categorical,
    ms_dyn_occu = .dispatch_ms_dyn_occu,
    ms_int_occu = .dispatch_ms_int_occu,
    stop(sprintf(
      "Internal error: family '%s' has status 'working' but no dispatcher.",
      family$name
    ), call. = FALSE)
  )

  # `positive` (the cover-hurdle positive-arm formula) is a formal argument that
  # the cover dispatchers consume as `dots$positive`; `fwd_dots` folded it back
  # into the splat above.
  fit <- do.call(dispatch, c(
    list(formula    = formula,
         data       = data,
         family     = family,
         detection  = detection,
         y          = y,
         visits     = visits,
         engine     = engine,
         approx     = approx,
         correction = route$correction,
         priors     = priors,
         control    = control),
    fwd_dots
  ))

  if (!inherits(fit, "tobs_fit")) {
    class(fit) <- c("tobs_fit", class(fit))
  }
  attr(fit, "tobs_family") <- family
  # Record the resolved public route for provenance / reproducibility.
  fit$method <- method
  fit$dims   <- dims
  fit
}


# Normalize the outer-grid progress knobs from a `control` list into the scoped
# `tulpa.nl_progress` option value read by tulpa's nested-Laplace fitters and the
# tulpaObs N-mixture spatial fitters (tulpaObs#25 / gcol33/tulpa#45). The
# flushed cell-k/n_grid + ETA reporter is ON by default (tulpaObs#43), matching
# the cover()/occu_cover() hurdle paths; set control$progress = FALSE to silence
# it. `progress.file` adds a heartbeat file for detached runs, written whenever
# it is non-empty regardless of the console bar.
.tobs_progress_opt <- function(control) {
  list(
    # `[[` (exact) not `$`: `control$progress` prefix-matches `progress.file`,
    # so a fit that sets only progress.file would otherwise read the file path
    # string as the console flag.
    progress          = control[["progress"]] %||% TRUE,
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


