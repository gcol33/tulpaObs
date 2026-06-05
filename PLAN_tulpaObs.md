# tulpaObs — Roadmap

`tulpaOcc` has been renamed to **tulpaObs**: a unified framework for
hierarchical latent-state observation models on the `tulpa` backend. The
rename (Phase 5 in the rollout below) was completed ahead of schedule to
keep the rename diff mechanical before further families ship.

This document captures the target architecture, the family roster, the
engine matrix, and the phased rollout. Nothing in this file ships in v0.1
beyond the new `tulpa_obs()` dispatcher and the family-object stubs.

---

## 1. Conceptual identity

tulpaObs models share four assumptions:

1. A latent ecological state exists (`z`, `N`, distance-attenuated density, ...).
2. Observations are imperfect.
3. Detection is probabilistic.
4. Replication identifies detectability.

What lives outside tulpaObs:

- Pure GLMM / GLM regression without a latent observation process — that is
  `tulpaglmm`.
- Process-only spatial smoothing without detection — that is `tulpa` directly.
- Compositional / Dirichlet community-cover models — separate scope.

What lives inside tulpaObs (in tension, by user request):

- `cover_hurdle()` family. Strictly, vegetation cover does not share the
  replicate-detection assumption (single survey per plot, the cover scale
  *is* the observation). It is included as a family because (a) the
  hurdle-on-cover model fits the same latent-state-plus-link generative
  template, with the latent state being "did the species pass the detection
  threshold" and the positive part being conditional cover, and (b) it is a
  concrete user need (MOTIVATE resurvey workflow). The framework remains
  honest: `cover_hurdle()` declares `replicates = "single"` and the family
  documentation flags the conceptual difference.

---

## 2. Architecture

```
+-----------------------------------------------------------------+
| tulpaObs — latent-state observation models                      |
|                                                                 |
|   tulpa_obs(                                                    |
|     formula,           # state-process formula                  |
|     data,                                                       |
|     family = occ(),    # latent-state family                    |
|     detection = ~ ...,  # detection-process formula              |
|     spatial = ...,      # spatial spec from tulpa / tulpaMesh   |
|     temporal = ...,     # temporal spec from tulpa              |
|     engine = "auto"     # "laplace", "nested_laplace", "nuts"   |
|   )                                                             |
+-----------------------------------------------------------------+
              |                          |                  |
              v                          v                  v
+----------------------+    +-------------------+   +----------------+
| Family object        |    | Engine selection  |   | Output         |
|  - latent type       |    |  Laplace (EM)     |   | tulpa_obs_fit  |
|  - obs likelihood    |    |  Nested Laplace   |   |  inherits      |
|  - replicate rule    |    |  NUTS             |   |  tulpa_fit     |
|  - encode/decode     |    |  (auto)           |   |                |
+----------------------+    +-------------------+   +----------------+
              |
              v
+-----------------------------------------------------------------+
| tulpa engine (latent priors, EM+Laplace, NUTS, nested Laplace)  |
+-----------------------------------------------------------------+
```

**Package boundaries**

- `tulpa` — engine. Latent priors (SPDE, NNGP, BYM2, AR1, ...), Laplace,
  EM+Laplace, MI/Gibbs correction, NUTS, nested Laplace, generic S3.
- `tulpaMesh` — mesh + SPDE precision builders.
- `tulpaglmm` — GLMM with the same engines, no latent observation layer.
- `tulpaObs` (this package, formerly `tulpaOcc`) — families, family-
  specific E-step weights, family-specific encode/decode, prediction,
  family-specific diagnostics.

---

## 3. Family roster

| family               | latent state          | observation likelihood        | replicates | spatial / temporal | engines    | status |
|----------------------|-----------------------|-------------------------------|-----------|--------------------|------------|--------|
| `occ()`              | Bernoulli z           | Binomial(p) per visit         | required  | all                | L, NL, NUTS| working (single-season) |
| `dynamic_occ()`      | Bernoulli z_t, HMM    | Binomial(p_{i,t,j})           | required  | all                | L, NUTS    | working |
| `multispecies_occ()` | z_{s,i}               | Binomial(p_{s,i,j})           | required  | all                | L, NUTS    | working (community RE)  |
| `integrated_occ()`   | shared z              | multi-source likelihoods      | required  | all                | L, NUTS    | working |
| `jsdm()`             | latent factor         | multivariate Bernoulli/Probit | no        | spatial            | NUTS       | working |
| `abun()` / `nmixture`| Poisson / NB N        | Binomial(N, p) per visit      | required  | all                | L, NL, NUTS| working (Poisson + negbin; non-spatial L + areal-spatial NL + non-spatial NUTS via the in-tree FullGradFn, tulpaObs#41) |
| `multispecies_nmix()`| Poisson N_{s,i}       | Binomial(N, p)                | required  | all                | L, NL, NUTS| working (`ms_abun()` on L via C++ community Laplace-EM; non-spatial Pois + negbin, areal-spatial Pois/NB via nested-Laplace, tulpaObs#12; community NUTS pending) |
| `dyn_abun()`         | Dail-Madsen N_t       | Binomial(N_t, p)              | required  | non-spatial        | L, NUTS    | working (exact HMM forward marginal, forward-mode-diff analytic gradients; Poisson init + constant recruitment; analytic-grad BFGS + in-tree NUTS; tulpaObs#37) |
| `distance()`         | Poisson / NB N        | hazard / half-normal binned   | replaced by distance bins | non-spatial | L, NUTS  | working (closed-form multinomial-over-N marginal on L + in-tree NUTS; half-normal / hazard-rate, line / point, Poisson + negbin; tulpaObs#38) |
| `removal()`          | Poisson / NB N        | sequential removal            | required  | non-spatial        | L, NUTS    | working (closed-form depleting-binomial marginal on L + in-tree NUTS; Poisson + negbin; tulpaObs#39) |
| `fp_occu()`          | z + classification    | confirmed / ambiguous (0/1/2) | required  | non-spatial        | L, NUTS    | working (Miller 2011 multistate; closed-form two-state marginal, analytic-grad BFGS + in-tree NUTS; tulpaObs#40) |
| `cover_hurdle()`     | latent presence + mu  | Binomial(occur) + Beta or LN  | single    | CAR, BYM2, SPDE    | L, NL      | working (lognormal + beta on L; lognormal + beta on NL via shared spatial field) |

Engines: **L** = single Laplace via tulpa, **NL** = nested Laplace via
tulpa, **NUTS** = HMC via tulpa.

---

## 4. API — `tulpa_obs()`

The new entry point. Existing entry points (`occu()`, `occu_fit()`) keep
working through the migration period and continue to be exported.

```r
fit <- tulpa_obs(
  formula     = ~ elev + forest,
  data        = sites,
  family      = occ(),
  detection   = ~ observer + effort,
  visit_data  = visits,
  y           = y_matrix,
  spatial     = tulpa::spatial_bym2(adj),
  temporal    = tulpa::temporal_ar1(year),
  engine      = "nested_laplace",
  priors      = tulpa::tulpa_priors(beta = "normal(0, 2.5)"),
  control     = list(n_threads = 4)
)
```

Family-specific kwargs (e.g. `K_max` for N-mixture, `positive` for
`cover_hurdle`) are passed via the family object:

```r
fit <- tulpa_obs(
  formula  = ~ elev,
  data     = sites,
  family   = nmixture(K_max = 80),
  detection = ~ effort,
  y        = counts,
  engine   = "nested_laplace"
)
```

```r
fit <- tulpa_obs(
  formula = ~ year_centred + habitat,
  data    = plots,
  family  = cover_hurdle(positive = "beta"),
  y       = cover,
  spatial = tulpa::spatial_bym2(adj),
  engine  = "nested_laplace"
)
```

**Migration**: `occu()` and `occu_fit()` remain the canonical single-season
occupancy entry until Phase 1 is complete. `tulpa_obs(..., family = occ())`
dispatches to `occu()` internally during the transition.

---

## 5. Engine matrix

| family            | Laplace (EM) | nested Laplace | NUTS |
|-------------------|--------------|----------------|------|
| occ               | working      | partial (needs hyperparam grid for spatial) | working |
| dynamic_occ       | working      | not yet        | working |
| multispecies_occ  | working      | not yet        | working |
| nmixture          | needs port from INLAabun | needs likelihood added to tulpa nested-Laplace registry | planned |
| cover_hurdle      | n/a (no E-step) | shipped (BYM2 / ICAR / CAR_proper, lognormal *and* beta positive; phi profiled via `tulpa_laplace_beta()` pre-fit, full posterior integration over phi is Phase 3) | planned |
| distance          | closed-form multinomial-over-N marginal (no E-step) | non-spatial shipped (tulpaObs#38); spatial offset not yet | working |
| dyn_abun          | exact HMM forward marginal, forward-mode-diff gradients (no E-step) | non-spatial shipped (tulpaObs#37); spatial / negbin / season-varying not yet | working |
| fp_occu           | closed-form two-state marginal (no E-step) | non-spatial shipped (tulpaObs#40); spatial not yet | working |

---

## 6. Nested Laplace gaps in `tulpa`

What `tulpa_nested_laplace()` already covers (registry in `nested_laplace.R`):

- Areal: `icar`, `bym2`, `car_proper`.
- Continuous: `nngp`, `hsgp`.
- Temporal: `rw1`, `rw2`, `ar1`.
- SPDE: `cpp_nested_laplace_spde()`.
- Likelihoods: binomial, Poisson, negbin-2.
- Single latent prior block per call.

What's missing for full INLA parity:

- **Joint multi-likelihood models** — e.g. binomial + beta in one fit with
  shared latent fields (the hurdle pattern, INLA `family = c(...)`).
  *Phase 1c shipped* `tulpa::tulpa_nested_laplace_joint()` for areal
  spatial priors (BYM2, ICAR, CAR_proper) with binomial / gaussian / beta
  arms and INLA-style `copy=` scaling on a designated arm. NNGP/HSGP/
  RW1/RW2/AR1 joint variants land under Phase 3.
- **Beta and lognormal likelihoods**.
- **Multiple latent prior blocks** in one call — e.g. spatial BYM2 + temporal
  AR1 + IID habitat. Currently single-block; extend the grid construction
  and inner Laplace to handle independent product blocks.
- **PC priors** as defaults across all priors.
- **Posterior linear predictor at held-out rows** — `INLA::inla()` returns
  fitted-value summaries at NA-response rows. Add a `predict_at` argument.
- **Skewness correction** for non-Gaussian likelihood marginals (the
  "simplified Laplace" of INLA).

These are scheduled under Phase 3 below.

---

## 7. Phased rollout

**Phase 0 — Architecture scaffold (this commit)**

- `tulpa_obs()` dispatcher introduced as a thin wrapper around `occu()`
  for `family = occ()`. All other families raise informative errors.
- Family-object infrastructure: `obs_family()` constructor + concrete
  family functions (`occ()`, `nmixture()`, `cover_hurdle()`, ...).
- No new functionality, no rename, no breaking changes.

**Phase 1 — Cover hurdle**

- *(Phase 1a — done)* Lognormal-positive variant shipped via two independent
  `tulpa::tulpa_laplace()` calls (binomial on occurrence, Gaussian on
  `log(cover)` over positive sites). Sigma estimated post-hoc as residual SE;
  Gaussian-arm Hessian SEs are scaled by `sigma_pos^2`. Files:
  `R/family_cover_hurdle.R`, `R/sim_cover_hurdle.R`,
  `tests/testthat/test-cover-hurdle-lognormal.R`, `vignettes/cover-hurdle.Rmd`.
- *(Phase 1b — done)* Beta likelihood in tulpa via `tulpa_laplace_beta()`.
- *(Phase 1c — shipped)* Joint multi-likelihood
  `tulpa_nested_laplace_joint()` with shared latent fields. Areal
  backends: BYM2, ICAR, CAR_proper. Likelihoods available in the kernel:
  binomial, gaussian, poisson, neg_binomial_2, beta (and the rest of
  `grad_hess_for_family`). cover_hurdle wired through with
  `engine = "nested_laplace"` for the lognormal-positive variant on any
  of the three areal spatial specs. Remaining: beta-positive variant of
  cover_hurdle, NNGP/HSGP/RW1/RW2/AR1 joint backends — Phase 3.
- *(Phase 1d — shipped)* Beta-positive variant of `cover_hurdle()` on the
  joint nested-Laplace engine. Precision `phi` is **profiled**: pre-fit via
  `tulpa::tulpa_laplace_beta()` on the positive subset (no spatial), then
  held fixed while the joint engine integrates over the spatial
  hyperparameters. Mirrors the `sigma_pos` handling for lognormal. Full
  posterior integration over `phi` is Phase 3. MOTIVATE-style
  within/between (Mundlak) decomposition helper shipped as
  `tulpaObs::within_between(data, group, vars)` — adds `<var>_btw`
  (per-group mean) and `<var>_wtn` (deviation from per-group mean) so
  longitudinal formulas like `~ year_btw + year_wtn` separate cross-plot
  baseline heterogeneity from within-plot temporal trend.
- *(Phase 1e — shipped)* Reproducible reduction of the MOTIVATE example as
  a vignette. `vignettes/cover-hurdle-motivate.Rmd` exercises the joint
  nested-Laplace cover-hurdle on a synthetic plot × year panel (chain
  BYM2, 25 plots, 4 visits each) for both `cover("beta")` and
  `cover("lognormal")` positive parts, using `within_between()` to split
  raw `year` into the cross-plot baseline and the within-plot deviation.
  Both within-plot slopes recover to truth at z > 2.4 on n_positive ≈ 37.

**Phase 2 — N-mixture**

- *(Phase 2a — shipped)* Single-species Poisson `abun()` on the **closed-form
  marginal** Laplace path. tulpa grew a direct N-mixture engine
  (`tulpa_nmix_laplace()`, analytical gradients + observed-Fisher curvature,
  exact sum over N to `K_max`) — strictly better than EM-around-INLA, so
  tulpaObs wires the family rather than porting the EM. `R/abun.R`:
  `.tobs_build_abun()` (`model_type = "nmix"`), `.tobs_fit_nmix()` ->
  `tulpa_nmix_laplace()`, `build_nmix_fit()` (joint lambda/p vcov, MVN draws),
  nmix S3 (`fitted`/`predict`/`simulate`/`residuals`/`nobs`), `simulate_abun()`.
  `tests/testthat/test-abun.R`: point recovery + 95% CI coverage across 30 seeds.
- *(Phase 2b — shipped)* Areal-spatial Poisson `abun()` via
  `method = "nested_laplace"` (an `icar()` / `bym2()` / `car_proper()` term on
  the abundance formula). Required a tulpa engine addition: the spatial
  N-mixture kernels (`nmix_spatial.cpp`, `nmix_spatial_bym2.cpp`) now return
  per-grid `cov_blocks` — the beta-block of each grid mode's joint `H^{-1}`,
  computed under the sum-to-zero constraint for the rank-deficient intrinsic
  fields (penalty-method, so the intercept variance is the constrained one, not
  the flat (intercept, field-mean) confounding of the improper prior). The R
  wrappers assemble the grid-integrated coefficient covariance via the law of
  total covariance (`.nmix_grid_vcov()`), returned as `vcov`. tulpaObs ungates
  the path (`.tobs_fit_nmix_spatial()` -> `tulpa_nmix_laplace_{icar,bym2,
  car_proper}`, one spatial unit per site) and reads the returned `vcov`.
  Calibrated: slope 95% CIs cover at nominal rate over seeds; intercept SE is
  finite/sane (`test-abun.R`).
- *(Phase 2c — shipped)* Negative-binomial abundance mixture
  (`abun(mixture = "negbin")`). tulpa's N-mixture kernel grew a parameterized
  marginal (`mixture = "NB"`, analytic dispersion score, NB size integrated as
  an extra outer grid axis on the spatial fitters); tulpaObs maps the family
  `mixture` to the kernel code and threads it through both the non-spatial and
  areal-spatial paths. Non-spatial: `log_r` jointly estimated (trailing `vcov`
  coordinate, with SE). Spatial: `r` grid-integrated (`r_mean` / `r_sd`).
  `simulate*()` draw NB. Matches `unmarked::pcount(mixture = "NB")`; recovery /
  coverage in `test-abun.R`.
- *(Phase 2d — shipped)* Community / multispecies N-mixture (`ms_abun()`, the
  spAbundance `msNMix` model). Per-species abundance and detection coefficients
  are random effects with Gaussian community hyperpriors
  (`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`, `beta_p_s ~ N(mu_p, Sigma_p)`).
  Required a tulpa engine addition: a C++ community N-mixture Laplace-EM
  (`tulpa_nmix_laplace_re`, `src/nmix_re.cpp`) reusing the per-site marginal
  kernel -- per-species coefficient mode-finding (complete-data Fisher), a
  closed-form covariance M-step, and fixed-effect SEs from the marginal
  observed-information Schur complement. The kernel grew a cached path
  (`compute_nmix_site_cached`, eta-independent `lgamma` precompute) that the
  single-shot Poisson kernel also delegates to (byte-identical; nmix regression
  suite passes). tulpaObs wires the family (`R/ms_abun.R`,
  `.tobs_fit_ms_nmix()` -> `tulpa_nmix_laplace_re`). Poisson only; recovery /
  20-seed coverage / per-species recovery / S3 in `test-ms-abun.R`. The R-only
  alternative (`tulpa_nmix_site_marginal()` -> `tulpa_re_aghq(make_group=)`) is
  correct but slow (R-interpreter per-group Newton); the C++ path is production.
- *(Phase 2e — shipped, gcol33/tulpaObs#13)* Grouped random effects on the
  single-species N-mixture. `abun()` accepts `(1 | g)` / slope / correlated
  formula RE on either arm (site-level grouping; observer-per-site, station,
  site cluster). `R/nmix_re_aghq.R::.tobs_nmix_re_aghq()` is a thin `make_site`
  callback over `nmix_site_marginal()` driven by `tulpa::tulpa_re_aghq()` --
  `n_quad = 1` is the joint Laplace, `n_quad > 1` is the AGHQ debias of the
  small-cluster sigma attenuation. Poisson and NB (`log_r` jointly estimated
  as the trailing theta coordinate). Gated: RE + spatial, RE + visit-level
  detection covariates, RE shared across both arms. Cost is tens of seconds
  per fit (R-closure marginal eval); a native per-group oracle along the
  `NMixCommunityOracle` pattern would close that gap.
- *(Phase 2f — shipped, tulpaObs#12)* Areal-spatial community N-mixture
  (`ms_abun()` + `icar()` / `bym2()` / `car_proper()`, the spAbundance
  `sfMsNMix` model). A shared areal field on `log lambda` with per-species RE,
  fitted by an in-tree nested Laplace-EM (`nmix_community_spatial.cpp`); Poisson
  and NB (`r` grid-integrated). Recovery / coverage / S3 / interop smoke in
  `test-ms-abun-spatial.R`.
- *(Phase 2g — shipped, tulpaObs#41)* NUTS for the single-species N-mixture
  (`abun()`, `method = "nuts"`; Poisson and negbin). The family likelihood is a
  tulpaObs `FullGradFn` (`src/abun_nuts.cpp`) over the closed-form per-site
  marginal, driving tulpa's generic NUTS engine -- no upstream tulpa change, the
  same boundary as the spatial-factor community occu_cover NUTS (gcol33/tulpa#67).
  Warm-started at the Laplace mode with a diagonal Laplace metric; the draws give
  calibrated (non-Gaussian) intervals and the per-site WAIC / LOO the Gaussian
  Laplace draws cannot. R target FD-verified, C++ port byte-exact vs it,
  front-door recovery + WAIC + cross-check tests in `test-abun.R`. Community
  non-spatial NB (`ms_abun(mixture = "negbin")`, per-species `log_r_s` RE) is also
  in place (verified: `r` recovers, community means unbiased; `test-ms-abun.R`).
- **Pending** (not bugs): NUTS for the community / areal-spatial N-mixture
  (`ms_abun()`) -- the same `FullGradFn` approach extended with the per-species RE
  / shared-field blocks.
- spAbundance / unmarked benchmarks.

**Phase 3 — Nested-Laplace extensions in `tulpa`**

- Beta + lognormal families.
- Joint multi-likelihood fits with shared latent fields.
- Multi-block latent priors.
- PC priors as defaults.
- `predict_at` row API.

**Phase 4 — Other observation processes**

- *(removal shipped, tulpaObs#39)* `removal()` sequential-depletion sampling.
  Latent `N ~ Poisson(lambda)` / negbin observed through `K` ordered removal
  passes; the depleting-binomial product (pass `k` sees `N - sum_{l<k} y_l`
  trials) equals the multinomial-removal likelihood and the latent `N` is summed
  out in closed form (sum to `K_max`). Shares the count-marginal Laplace driver
  (`src/marginal_count_laplace.h`) and NUTS machinery
  (`src/marginal_count_nuts.h`) with `abun()`; the per-site math is
  `src/removal_kernel.h` (built on the shared `accumulate_count_moments` /
  `fill_nb_dispersion` extracted from `src/nmix_kernel.h`). `R/removal.R` (family
  wiring, S3, `simulate_removal()`), `R/removal_nuts.R`. Non-spatial Poisson + NB
  on `laplace`, NUTS over the same marginal; recovery / 95% coverage / NB
  dispersion / NUTS recovery + WAIC + Poisson-equivalence anchor in
  `test-removal.R`. Spatial / RE removal not yet wired.
- *(distance shipped, tulpaObs#38)* `distance()` binned distance sampling.
  Latent `N ~ Poisson(lambda)` / negbin in a covered region; the per-bin detected
  counts are multinomial over `(bin 1, ..., bin B, undetected)` with cell
  probabilities `pi_b = integral_bin g(x; sigma) f(x) dx` (half-normal /
  hazard-rate key, line / point transect distance density `f`), integrated by
  Gauss-Legendre quadrature, and the latent `N` is summed out in closed form (sum
  to `K_max`). Reuses the shared count-marginal core (`accumulate_count_moments` /
  `fill_nb_dispersion` from `src/nmix_kernel.h`) for the abundance / NB arm; the
  detection arm (site-level `log sigma`, optional scalar hazard shape, bin
  integrals + first/second eta-derivatives) is `src/distance_quad.h` /
  `src/distance_kernel.h` with its own Laplace driver (`src/distance_laplace.cpp`)
  and NUTS target (`src/distance_nuts.cpp`) over the shared tulpa-NUTS engine
  driver (`src/nuts_engine.h`). `R/distance.R` (family wiring, S3,
  `simulate_distance()`), `R/distance_nuts.R`. Non-spatial Poisson + NB on
  `laplace`, NUTS over the same marginal; recovery / 95% coverage / hazard-shape /
  NB dispersion / point-transect / NUTS recovery + WAIC + a Poisson-thinning
  closed-form anchor + an analytic-vs-FD observed-information check in
  `test-distance.R`. Spatial / RE distance not yet wired.
- *(fp_occu shipped, tulpaObs#40)* `fp_occu()` multistate false-positive
  occupancy (Miller et al. 2011 confirmed-detection design). Detection states
  `y in {0,1,2}` (none / ambiguous / certain); certain detections only at
  occupied sites identify the model. Latent occupancy `z` summed out in closed
  form (two states); four site-level logit arms (psi, true detection p11,
  false-positive p10, certain-classification b). Per-site math
  `src/fp_occu_kernel.h` (analytic gradient, no second derivatives needed for the
  vcov -- the Laplace fit is analytic-gradient BFGS with an observed-information
  covariance from the FD-Jacobian of the analytic gradient at the mode), NUTS
  target `src/fp_occu_nuts.cpp` over the shared `src/nuts_engine.h`. `R/fp_occu.R`
  (family wiring, S3, `simulate_fp_occu()`), `R/fp_occu_nuts.R`. Recovery / 95%
  coverage / fp-arm covariate / NUTS recovery + WAIC + a direct-marginal anchor
  in `test-fp_occu.R`. Spatial / RE not yet wired.
- *(dyn_abun shipped, tulpaObs#37)* `dyn_abun()` Dail-Madsen open-population
  N-mixture. `N_1 ~ Poisson(lambda)`; `N_t = Binomial(N_{t-1}, omega) +
  Poisson(gamma)`; `Binomial(N_t, p)` observation. The latent abundance sequence
  is summed out by an exact HMM forward recursion over states `0..K_max` (not
  closed form); analytic gradients by forward-mode differentiation of the scaled
  forward algorithm (`src/dyn_abun_kernel.h`). Analytic-grad BFGS Laplace with an
  observed-information vcov (FD-Jacobian of the analytic gradient) and NUTS over
  the same marginal (`src/dyn_abun_nuts.cpp` via shared `src/nuts_engine.h`).
  `R/dyn_abun.R` (family wiring, S3, `simulate_dyn_abun()`), `R/dyn_abun_nuts.R`.
  Recovery / 95% coverage / NUTS recovery + WAIC + a forward-recursion anchor
  (C++ vs an independent R forward, exact to 1e-9) in `test-dyn_abun.R`. Poisson
  initial + constant recruitment this round; negbin / season-varying dynamics /
  spatial / RE not yet wired.
- All five filed observation-family issues (tulpaObs#37/#38/#39/#40) shipped;
  remaining open work is spatial / RE / negbin extensions and the threaded
  nested-Laplace memory issue (tulpaObs#42).

**Phase 5 — Package rename (DONE ahead of schedule, before Phase 1b)**

- `tulpaOcc` → `tulpaObs`: Package field, `useDynLib`, all `tulpaOcc_*`
  S3 class names → `tulpaObs_*`, C++ `namespace tulpaOcc` → `namespace
  tulpaObs`, DLL symbol via `Rcpp::compileAttributes()`. Single mechanical
  commit. Local folder renamed `tulpaOcc/` → `tulpaObs/`; GitHub repo
  renamed `gcol33/tulpaOcc` → `gcol33/tulpaObs` on 2026-05-15 (GitHub
  301-redirects the old URL, so prior install_github clones keep working).
- **Deferred**: deprecation shim package `tulpaOcc` that re-exports the
  public API. Not needed yet because INLAabun and other downstream
  consumers don't import `tulpaOcc::` directly. File when the first
  external consumer breaks.

---

## 8. Non-goals

- Generic compositional / Dirichlet vegetation modelling.
- Bespoke functional-data ecology.
- General-purpose GLMM API — that belongs in `tulpaglmm`.
- Replacement for `unmarked` interactive plotting and study-design tooling.
- Pure MCMC-first design — `tulpa` provides NUTS, but the default path is
  Laplace / nested Laplace.

---

## 9. Open questions

1. Should `cover_hurdle()` live in tulpaObs or in a sibling package
   `tulpaCover` that depends on `tulpa` + `tulpaMesh`? Decision: in
   tulpaObs for now, behind a clear family identity flag. Reconsider if
   the family grows compositional-data-specific machinery.
2. Should `tulpa_obs()` accept `engine = "auto"` and select per-family
   defaults? Tentative: yes, with explicit override.
3. Should the rename happen before or after Phase 1? **Resolved 2026-05-13**:
   done after Phase 1a so the rename commit is mechanical and isolated
   from Phase 1b feature work.
