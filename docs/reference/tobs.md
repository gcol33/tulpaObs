# Fit a hierarchical latent-state observation model

Unified entry point for occupancy, abundance, distance, removal, and
related models that share the latent-state-plus-imperfect-detection
generative template. The specific model is chosen via the `family`
argument; the engine (Laplace / nested Laplace / NUTS) via `engine`.

## Usage

``` r
tobs(
  formula,
  data,
  family,
  occurrence = NULL,
  detection = NULL,
  positive = NULL,
  y = NULL,
  visits = NULL,
  method = c("auto", "laplace", "laplace_sla", "laplace_gibbs", "laplace_mi", "pg_gibbs",
    "nested_laplace", "nested_laplace_sla", "nuts"),
  priors = NULL,
  control = list(),
  by = NULL,
  ...
)
```

## Arguments

- formula:

  state-process formula, e.g. `~ elev + forest`. For occupancy this is
  the occupancy probability formula; for N-mixture the abundance
  formula; for the cover hurdle the latent-presence formula (also given
  as `occurrence`, which reads symmetrically with `detection` and
  `positive`).

  Single-vector-response families (the cover hurdle,
  [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md))
  also accept the response on the left-hand side,
  `response ~ predictors`, in which case `y =` is omitted (e.g.
  `cover.flat ~ time + habitat`). The LHS is evaluated against `data`
  (then the calling environment), so it may be a bare column or an
  expression. Matrix / array / list response families
  ([`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
  [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md), the
  `ms_*` families, ...) keep the one-sided form and supply the response
  via `y =`; a two-sided formula for those errors.

  Structured effects are written as terms inside the formula, the way
  `lme4`, `mgcv`, and `INLA` do: spatial fields `icar(graph = adj)`,
  `bym2(graph = adj)`, `gp(lon, lat, prior_range = c(r0, alpha))`,
  `spde(lon, lat)`; random effects `re(group)`; temporal fields
  `temporal(time)`; spatially varying coefficients
  `spatial(~ 1 + w || cell, graph = adj)` over a graph or
  `spatial(lon, lat, model = "svc", coefficients = ..., prior_range = c(r0, alpha))`
  over coordinates; community latent factors `latent(k)`. A term enters
  whichever linear predictor it is written in (occupancy `formula` or
  `detection`). To share one realization across both predictors, tag the
  term with `id = "u"` and write `share("u")` in the other formula.

  The spatial fields also have a single-verb umbrella
  `spatial(..., model = ...)` that selects the field type by name,
  mirroring `temporal(time, type = ...)` and `INLA`'s
  `f(i, model = ...)`: `spatial(graph = adj, model = "bym2")` is
  `bym2(graph = adj)` and `spatial(lon, lat, model = "spde")` is
  `spde(lon, lat)`. `model` is one of `"icar"`, `"bym2"`, `"car"`,
  `"car_proper"`, `"gp"`, `"multiscale_gp"`, `"spde"`, `"svc"`;
  per-model arguments pass through unchanged. The umbrella also reads a
  varying-coefficient bar (`spatial(~ 1 + w || cell, graph = adj)`), so
  both flavours of spatially-varying coefficient – areal over a graph,
  continuous over coordinates (`model = "svc"`) – are written through
  the same verb.

  The continuous fields
  ([`gp()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
  [`svc()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
  [`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md))
  require `prior_range = c(r0, alpha)`, a penalized-complexity prior on
  the spatial range encoding `P(range < r0) = alpha` (Fuglstad et al.
  2019). The range is in the units of the coordinates – the kernel is
  `exp(-d / range)` – so choose `r0` as a distance below which the
  field's correlation would be surprisingly short: on unit-square
  coordinates `prior_range = c(0.1, 0.05)` reads as "a 5% chance the
  range is under 0.1". There is no default, deliberately. The range is
  weakly identified by the likelihood alone, so a default would be an
  invented prior doing real work on the posterior rather than a
  convenience.

  Random effects also accept `lme4` bar syntax as shorthand for
  [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md):
  `(1 | g)` is `re(g)`, `(x | g)` is a correlated random intercept and
  slope `re(g, type = "slope", covariate = x)`, and `(x || g)` drops the
  correlation (`correlated = FALSE`). Use
  [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  directly for AR1/RW structures or other options.

- data:

  data frame of site-level covariates with `nrow(data) == nrow(y)`, or a
  `tobs_data` frame (from
  [`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)
  /
  [`tobs_format()`](https://gillescolling.com/tulpaObs/reference/tobs_format.md))
  bundling the response, site covariates, and visit covariates. When a
  frame is passed, `y` and `visits` are taken from it and must not also
  be supplied.

- family:

  a `tobs_family` object (see
  [`obs_family()`](https://gillescolling.com/tulpaObs/reference/obs_family.md)
  and the concrete constructors
  [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
  [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
  [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md),
  ...).

- occurrence:

  state-process (occupancy / latent-presence) formula for the cover
  hurdle
  ([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  / [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)),
  reading symmetrically with `detection` and `positive`. The front-door
  name for `formula`; supply one of the two.

- detection:

  detection-process formula, e.g. `~ observer + effort`.
  Family-dependent: required for
  [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) and
  [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
  ignored for
  [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md) and
  (currently)
  [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md).

- positive:

  positive-arm (cover) formula for
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
  e.g. `~ time.sc + habitat`. The occurrence spatial field is carried
  onto this arm with a `share()` selector, the INLA-style cross-arm edge
  written in the formula. The selector is a constructor, so no field
  name is needed in the common case:

  - `share(spatial(), alpha = grid(g))` copies the occurrence arm's
    spatial effect and marginalizes the coupling amplitude over the
    nodes `g`; `alpha = grid(n = 9)` marginalizes it over the engine's
    own axis read at `n` slab nodes instead, which sharpens the axis
    without restating the prior structure it carries; `alpha = <scalar>`
    fixes it.

  - `share(spatial(), terms = list(intercept = grid(g0), time.sc = grid(g1)))`
    gives a per-component amplitude, keyed by the field's own block
    names (the intercept block and a `||`-declared trend column, or its
    alias `trend`); `terms =` must address every block.

  - `share(spatial(), prior = list("pc.prec", c(4, 0.01)))` regularizes
    the copy coefficient itself (`prior_alpha` in
    [`tulpa::tulpa_nested_laplace_joint()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace_joint.html)).
    One reaches the engine per fit.

  - `share(spatial(cell_idx), ...)` disambiguates by grouping variable
    when the occurrence arm carries several spatial effects.

  - decoupling is structural: omit `share()` so the field rides
    occupancy only (a block with no `share()` is pinned at zero
    coupling), or write a
    [`spatial()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    term for the cover arm's own field.

  Defaults to `detection` when unset (a per-visit cover design matching
  the detection design).

- y:

  response. Shape depends on family:

  - [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) – N
    x J detection-history matrix.

  - [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) – N
    x J integer count matrix.

  - `ms_*` – S x N x J array.

  - [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md) –
    length-N vector of cover proportions in \[0, 1\].

  For a single-vector-response family
  ([`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md))
  the response may instead be written on the `formula` left-hand side
  (`response ~ predictors`), in which case `y =` is omitted. Supplying
  the response both on the LHS and via `y =` errors. Matrix / array /
  list response families take the response via `y =` only.

- visits:

  optional visit-level detection covariates. Accepts either:

  - a named list of `[n_sites, max_visits]` matrices (the shape returned
    by
    [`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)
    in `det.covs`) – flattened internally to a long data frame in
    site-major order;

  - a data frame with `nrow(y) * ncol(y)` rows in site-major order.

  When `visits` is provided without a `"formula"` attribute, the
  `detection` argument is interpreted as the visit-level detection
  formula and the site-level detection design matrix is an intercept
  only. To split visit-level and site-level detection covariates (e.g.
  visit-level effort plus site-level observer category), pass a long
  data frame with `attr(visits, "formula") = ~ effort` and use
  `detection = ~ observer` for the site-level terms.

  For [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
  a random-effect bar in the visit-level formula is a random effect over
  visit rows: `detection = ~ effort + (1 | observer)` with `observer` a
  column of `visits` gives each observer its own shift of the detection
  logit on the visits it made. It is fitted with `method = "laplace"` or
  `"nuts"`, and other families refuse it.

- method:

  inference route, naming a fully-specified path rather than a pair of
  orthogonal knobs:

  - `"auto"` – the family's default route (see `default_engine`).

  - `"laplace"` – EM + Laplace with Gaussian marginals (fast default).

  - `"laplace_sla"` – Laplace with skew-corrected (simplified-Laplace)
    marginals.

  - `"laplace_gibbs"` / `"laplace_mi"` – Laplace with a post-EM Gibbs /
    multiple-imputation correction. The fixed-effect prior threads into
    the correction refits, so these use the same weakly-informative
    default prior as `"laplace"`; pass `priors = FALSE` for the
    unpenalised correction.

  - `"pg_gibbs"` – a Polya-Gamma Gibbs sampler over the exact
    single-season occupancy posterior (the spOccupancy `PGOcc` engine).
    A real MCMC chain (with `Rhat` / `ESS` diagnostics), distinct from
    `"laplace_gibbs"` (a stochastic-EM variance correction). Sampler
    controls (`n.iter`, `n.warmup`, `n.chains`, `n.thin`, `seed`,
    `sigma.beta`). v1: single-season
    [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
    site-level detection, no structured terms.

  - `"nested_laplace"` – multi-block nested Laplace (single-season
    occupancy and cover-hurdle joint).

  - `"nested_laplace_sla"` – nested Laplace with skew-corrected
    marginals.

  - `"nuts"` – HMC / NUTS sampler (every structure; reports Rhat / ESS).
    Not every method is available for every family (e.g. the cover
    hurdle has no `"nuts"` path; `"nested_laplace"` is occupancy- and
    cover-only). An unsupported method errors with the list of methods
    that family supports.

- priors:

  optional prior specification. For occupancy families fit with a
  Laplace method (`method = "laplace"`, `"laplace_sla"`,
  `"nested_laplace"`), pass a list or
  [`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md)
  object to set weakly-informative quadratic priors on the fixed-effect
  coefficients (defaults pull the detection intercept toward `p = 0.5`
  and break the psi-p identifiability ridge at small `J`). Pass
  `priors = FALSE` to disable the default prior and recover the
  unpenalised MAP. The `"laplace_gibbs"` / `"laplace_mi"` routes apply
  the same default prior and thread it through the correction refits.
  For NUTS, this is forwarded to the underlying tulpa engine.

- control:

  list of low-level engine controls. Names follow the dotted-separator
  convention. Every default below is resolved from one table
  (`.TOBS_ENGINE_DEFAULTS`), so it is a property of the fitting engine
  rather than of the family. Where a default differs between the
  single-species and community entries, both are given.

  Sampler controls (`method = "nuts"`):

  - `n.iter` – post-warmup sampling iterations kept per chain (default
    1000); the total run per chain is `n.iter + n.warmup`.

  - `n.warmup` – warmup / adaptation iterations per chain, discarded
    (default 1000).

  - `n.thin` – keep every `n.thin`-th post-warmup draw (default 1). The
    kept draws and the per-iteration diagnostics (`divergent`,
    `accept_prob`, `treedepth`) are thinned by the same stride.

  - `n.chains` – number of chains, run with offset seeds and pooled
    (default 1). Split-Rhat / bulk / tail ESS are reported on
    `$convergence`.

  - `n.threads` – chains to run in parallel (default 1, sequential).
    Values `> 1` use a PSOCK cluster and require tulpaObs to be
    installed. This is chain parallelism; it does not change the thread
    count inside a single gradient evaluation.

  - `n.threads.grad` – OpenMP threads inside ONE gradient evaluation of
    a community NUTS target
    ([`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
    [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
    / [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md),
    [`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md),
    [`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md)),
    whose per-species loop is parallel. Default 0 leaves the count to
    OpenMP. The per-species reduction is serial and order-fixed, so the
    gradient is the same at any count.

  - `adapt.delta` – target acceptance probability (default 0.8 on the
    single-species families, 0.9 on the community samplers).

  - `max.treedepth` – NUTS maximum tree depth (default 10).

  - `seed` – base RNG seed; chain `c` uses `seed + c - 1` (default 42 on
    the single-species families, 1 on the community samplers). The
    resolved per-chain seeds are stored on `$seeds`.

  - `sigma.beta` – prior SD on the coefficients (default 10 on the
    single-species families and on the log-link community samplers
    [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
    / [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md) /
    [`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md),
    whose unit change is multiplicative; 5 on the logit-link community
    samplers
    [`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
    [`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
    [`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md),
    [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)).

  - `sigma.logr` – prior SD on the community-mean log-dispersion
    `mu_log_r` (default 1.5), on the negative-binomial samplers that
    carry one
    ([`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md),
    [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md),
    [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md)).
    At the default this is an informative prior on the dispersion scale,
    so raise it to compare the sampler against a maximum-likelihood
    target. Ignored where the family or mixture has no log-dispersion
    arm.

  Sampler controls (`method = "pg_gibbs"`). A conjugate sweep is far
  cheaper than a NUTS trajectory, so the chain is longer and two chains
  run by default. **`n.iter` counts differently here**: it is the TOTAL
  number of sweeps and warmup comes out of it, so the chain keeps
  `n.iter - n.warmup` draws, where under `"nuts"` it is the kept count
  and the run is `n.iter + n.warmup` long.

  - `n.iter` – total sweeps per chain (default 3000).

  - `n.warmup` – sweeps discarded from the front (default 1500), leaving
    1500 kept draws.

  - `n.thin` – keep every `n.thin`-th post-warmup sweep (default 1).

  - `n.chains` – number of chains (default 2, so split-Rhat is available
    without a second call).

  - `seed` – base RNG seed (default 1).

  - `sigma.beta` – coefficient prior SD (default 2.5; tighter than the
    NUTS one because a conjugate update has no step-size adaptation to
    absorb a wide prior). There is deliberately no `adapt.delta` /
    `max.treedepth`: those are HMC knobs.

  - `n.seeds` – number of seed-offset refits to fit and LOO-stack into a
    `tobs_stack` ensemble (default 1, a single fit). Member `k` uses
    base seed `seed + k - 1`. Only meaningful for the stochastic routes
    (`"nuts"`, `"laplace_gibbs"`, `"laplace_mi"`); the deterministic
    Laplace methods reject it. Seed-variants are statistically
    identical, so their stacking weights come out roughly uniform (this
    is a Monte-Carlo robustness device) – pass distinct fits to
    [`tobs_stack()`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)
    for a genuine model average.

  Laplace controls (`method = "laplace"` / `"laplace_sla"` /
  `"nested_laplace"`): `max.iter`, `tol`, `damping`, `sigma.beta`.

  - `logr.sigma.prior` – Penalized-Complexity prior `c(U, alpha)`
    (`P(sigma_log_r > U) = alpha`) on the per-species log-dispersion SD
    under `ms_abun(mixture = "negbin" / "zinb")`. Default `NULL`, pure
    maximum likelihood. `sigma_log_r` is one scalar variance over
    species and at few species can settle near its lower boundary, which
    shifts `mu_log_r` and narrows its interval together (see
    [`?ms_abun`](https://gillescolling.com/tulpaObs/reference/ms_abun.md));
    the prior adds curvature there. `omega.sigma.prior` is the same knob
    for the structural-zero variance and does default to `c(1, 0.05)`;
    when both are set they must be equal, since one Penalized-Complexity
    prior is applied across every regularized block.

  - `re.aghq` – for a formula random effect under `method = "laplace"`,
    run the adaptive Gauss-Hermite debias of the variance components
    after the EM converges (default `TRUE`). Removes the Laplace
    small-cluster attenuation of `sigma` / the RE correlation for binary
    occupancy; set `FALSE` for the raw EM (Laplace, `nAGQ = 1`) fit.

  - `n.quad` – quadrature points. One name across several routes, each
    integrating a different marginal over a different latent dimension,
    so the default is per route rather than one number. `n.quad = 1` is
    always the plain Laplace (`nAGQ = 1`) marginal; higher values refine
    it toward the exact one. Formula random effect under
    `method = "laplace"`: 9. Adaptive Gauss-Hermite over the exact
    per-group marginal, debiasing the Laplace small-cluster attenuation
    of a binary occupancy variance component. Binary data carries little
    information per group, so the refine wants many nodes. Community
    N-mixture
    ([`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)):
    1, i.e. the EM default. Each species' count marginal is already
    informative, so the AGHQ refine barely moves the community
    covariances; it is opt-in via `optimizer = "joint_fd"` for the
    sparse / rare-species regime. `n.quad.scalar` (default 3, and
    floored there) is the trailing per-species log-dispersion /
    structural-zero-logit coordinate, integrated separately. Community
    joint occupancy-cover
    ([`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)): 5.
    Tensor AGHQ over the joint per-species RE vector, so the node count
    is raised to a power of the RE dimension and stays small. Latent
    cover-per-unit (`cover_aggregate = "latent"`): 15 for a beta cover
    arm, 1 for lognormal, whose per-unit marginal is closed form and
    needs no quadrature at all. Community
    [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    factors: 5, the Gauss-Hermite nodes the joint site marginal
    integrates the factor scores on. Both the loading magnitude and the
    score-matched offset are insensitive to it (argmax stable to \< 0.4%
    against 21 nodes).

  - `max.outer` – for a community family whose
    [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    factors or shared areal field are fit by block coordinate ascent,
    the cap on the outer alternation between the community EM and the
    field / factor update. A field block reaches `tol` and stops early,
    so its default 25 is only a cap; a factor block does not, and each
    family sets its own budget from a measured bias curve (150 on
    [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
    / [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md) /
    [`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
    25 elsewhere).

  - `factor.starts` – candidate starting directions the first factor
    pass selects over, on the joint marginal. Each costs a full loading
    EM against that family's oracle, so the default is per family: 1 on
    [`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
    (measured to buy nothing there against a 2.0-2.3x cost), 8
    elsewhere. Accepted only by the families fit this way; the resolved
    value, and `max.outer` / `n.quad` alongside it, is reported as
    `fit$latent_control`.

  - `re.lkj` – LKJ shape (`eta`) regularizing a *correlated* random
    slope's correlation in the `re.aghq` refine (default 1.5). Pulls a
    weakly-identified RE correlation off the `+-1` boundary toward 0
    without touching the marginal SDs; `re.lkj = 1` disables it
    (uniform). No effect on intercept / uncorrelated terms.

  - `sd.load` – prior SD on a spatial-factor loading in the community
    occupancy-cover fit (default 1). The auto-rank ladder selects `K` by
    marginal evidence under this prior, so the selection fit and the
    final fit necessarily read the same value.

  - `inner.solver` – for a spatial community N-mixture
    ([`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
    with an
    [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    /
    [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    /
    [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
    field on the abundance arm), the inner solver integrating the shared
    field given the community: `"em"` (default) the closed-form
    Laplace-EM M-step, or `"newton"` the exact- Newton shared-field
    solve alternated with a tulpa AGHQ community debias. Both integrate
    the field hyperparameter on the outer grid and return the same fit
    object; `"newton"` is Poisson- and areal-only, and markedly slower
    (an FD-gradient profile loop per grid node) – an accuracy /
    validation alternative, not the production default.

  - `integration` – how the in-package spatial / community
    nested-Laplace fitters integrate the outer field hyperparameters
    (`tau`, `rho`, `sigma`, `range`): `"grid"` (default) a fixed tensor
    grid, or `"ccd"` a mode-centred central-composite design placed at
    the marginal-likelihood mode and scaled by the outer posterior
    covariance, with the outer PSIS Pareto-k reported on
    `fit$spatial_pareto_k`. CCD declines to the grid when the outer
    curvature is ill-conditioned (a weakly-identified axis) or for a
    single positive hyperparameter; `fit$spatial_integration` records
    which ran. Each outer node is a full inner solve, so `"ccd"` adds a
    mode-find without a node saving on these coarse grids – it is
    opt-in, most useful when a multi-axis hyperparameter posterior is
    well identified.

  - `adaptive.grid.cutoff`, `adaptive.grid.stride`,
    `adaptive.grid.max.frac`, `adaptive.grid.min.cells` – tuning for
    `integration = "grid_adaptive"` on the joint-coupled spatial
    families
    ([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
    [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md),
    [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md)
    spatial,
    [`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md)),
    which the engine hosts alongside `"grid"` and `"ccd"`. That
    integrator evaluates a strict subset of the same tensor lattice – it
    floods outward from the posterior mode and keeps the cells within
    `adaptive.grid.cutoff` log-density of the peak – and declines back
    to the dense tensor when the kept region would rival it, so it
    trades inner solves for nothing but resolution of the hyperparameter
    tails. `adaptive.grid.cutoff` is the keep / expand radius (larger
    keeps more cells, closer to the dense tensor),
    `adaptive.grid.stride` the coarse-seed subsample stride per axis,
    `adaptive.grid.max.frac` the kept-fraction ceiling past which it
    declines, and `adaptive.grid.min.cells` the smallest dense tensor
    worth the machinery. Unset leaves each at the engine default. These
    are a DIFFERENT mechanism from `adaptive.grid` /
    `adaptive.grid.edge.thresh` / `adaptive.grid.max.passes`, which
    refine an already-integrated grid; the two compose and neither
    implies the other. The saving is largest where a fit carries several
    outer axes – two spatial arms (a shared field plus an arm-specific
    one) put three on the grid, which `integration = "auto"` resolves to
    the dense tensor.

  - `diagnose.k` – for the joint-coupled spatial families
    ([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
    [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md)
    spatial,
    [`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md)),
    whether to score the outer hyperparameter Gaussian summary with an
    importance-sampling Pareto-k. Defaults `FALSE` (it re-solves the
    inner Laplace on the full field `k.samples` times). When `TRUE`, the
    fit carries `pareto_k` (`< 0.7` = reliable summary),
    `pareto_k_is_ess` (the importance-sampling ESS on the PSIS-smoothed
    weights; `/ k.samples` is the relative IS efficiency), and
    `pareto_k_proposal_source` at the top level (and in
    [`glance()`](https://generics.r-lib.org/reference/glance.html)).
    `pareto_k_proposal_source` is `"mode_hessian"` when the importance
    proposal came from the Laplace curvature at the hyperparameter mode
    – curvature-backed, so the k-hat stays trustworthy even when a sharp
    posterior collapses the integration grid – or `"grid_moment"` when
    it came from the grid-weighted node covariance, the regime to watch
    as the grid concentrates. Stochastic-correction controls
    (`"laplace_gibbs"` / `"laplace_mi"`): `n.gibbs` / `n.imputations`
    (Rubin-pooled draw count) and `seed` (stored on `$seed`).

  Progress controls (every method): `progress` toggles the console
  iteration / grid bar (default `TRUE` in an interactive session,
  `FALSE` otherwise), `progress.every` and `progress.throttle` set its
  emit cadence, and `progress.file` names a heartbeat file rewritten
  with `"<done> <total> <elapsed_s> <eta_s>"`. The file is the channel
  that survives a detached run, where a console flush does not, and it
  is written whenever it is set regardless of `progress`. The
  environment variable `TULPAOBS_PROGRESS` sets the console default for
  the whole session: `0` (also `false` / `no` / `off`) turns it off, any
  other value turns it on – for a batch or CI run whose caller cannot
  pass `control` to each individual fit; an explicit `control$progress`
  still wins.

  Control names are validated against the chosen `method`: passing a
  sampler control (e.g. `n.chains`) to a Laplace method, a Laplace
  control (e.g. `max.iter`) to `"nuts"`, `seed` to a deterministic
  route, or an unrecognized name raises an error rather than being
  silently ignored.

- by:

  optional name of a species column for a per-species batched fit
  ([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  / [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
  only). When supplied, `data` is a long, plot-level frame: one row per
  site-visit
  ([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md))
  or per site
  ([`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)),
  with `by` giving the species. `tobs()` splits `data` by that column,
  builds each species' response onto one shared site x visit grid (via
  [`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)),
  fits the B species independently, and returns a `tobs_batch`.
  Per-species fits are statistically independent (no pooling) –
  identical to fitting each species with a separate single-species
  `tobs()` call; for a community model that borrows strength across
  species use
  [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md).
  The long -\> response pivot needs the column names: pass `site = ` and
  `response = ` (the detection 0/1 column for
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
  the cover column for
  [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)),
  plus, for
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
  `visit = ` (the replicate column) and `y_pos = ` (the cover column)
  and any visit-level `det.covs = `. Site-level covariates (those the
  `formula` / `detection` reference at the cell level) are read as the
  first row per site, the way
  [`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)`(occ.covs = )`
  does.

  The same long-frame contract drives a SINGLE
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  fit when `by` is omitted: pass `site = `, `visit = `, `response = `
  (the 0/1 detection column) and `y_pos = ` (the cover column), plus any
  visit-level `det.covs = `, with a long, plot-level `data`, and
  `tobs()` builds the paired occurrence / cover arms and the site-level
  design for you – the by= batch reduced to one species, so you no
  longer hand-build the two responses and align them. The arms default
  to the compact (ragged) layout on the nested-Laplace route (no
  per-site visit cap); set `control$compact = FALSE` for the dense grid.
  See
  [`occu_cover_inputs()`](https://gillescolling.com/tulpaObs/reference/occu_cover_inputs.md)
  to build and inspect the arms without fitting.

  The long-frame `response =` here names a COLUMN of `data` (the pivot
  key holding the detection / cover values). It is a separate argument
  from the positive-part distribution the family constructor takes
  (`cover(response = "beta")`, `occu_cover(response = "lognormal")`),
  which selects the cover-arm likelihood rather than a data column. Both
  may appear in one call –
  `tobs(long, family = cover(response = "beta"), site = , response = "cover.flat", ...)`
  – and do not interact.

- ...:

  family-specific named arguments forwarded to the underlying engine
  builder.

## Value

An object of class `c("tobs_fit", "<family>_fit", "tulpa_fit")`. When
`control$n.seeds > 1`, a `tobs_stack` ensemble of the seed-offset refits
is returned instead (see
[`tobs_stack()`](https://gillescolling.com/tulpaObs/reference/tobs_stack.md)).
When `by` is supplied, a `tobs_batch` of per-species fits is returned
(see
[`tobs_get()`](https://gillescolling.com/tulpaObs/reference/tobs_get.md)).

## Examples

``` r
# \donttest{
# Single-season occupancy
sim <- simulate_occu(N = 100, J = 3, n_occ_covs = 1, n_det_covs = 1,
                     seed = 1)
fit <- tobs(
  formula   = ~ occ_cov1,
  data      = sim$data,
  family    = occu(),
  detection = ~ det_cov1,
  y         = sim$y,
  method    = "laplace",
  control   = list(verbose = FALSE, progress = FALSE)
)
summary(fit)
# }
```
