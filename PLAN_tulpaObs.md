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
| `abun()` / `nmixture`| Poisson N             | Binomial(N, p) per visit      | required  | all                | L, NL, NUTS| working (Poisson, non-spatial L). negbin / areal-spatial / NUTS pending upstream tulpa |
| `multispecies_nmix()`| Poisson N_{s,i}       | Binomial(N, p)                | required  | all                | L, NL, NUTS| planned (Phase 2) |
| `dynamic_nmix()`     | Dail-Madsen N_t       | Binomial(N_t, p)              | required  | all                | NUTS       | planned (Phase 3) |
| `distance()`         | density               | hazard / half-normal binned   | replaced by distance bins | all | L, NUTS  | planned (Phase 4) |
| `removal()`          | N                     | sequential removal            | required  | spatial            | L, NUTS    | planned (Phase 4) |
| `false_positive()`   | z + classification    | confirmed / ambiguous         | required  | all                | NUTS       | planned (Phase 4) |
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
| distance          | needs E-step over distance bins | not yet | planned |

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
- **Pending upstream tulpa** (deferred, not bugs):
  - **negbin marginal likelihood** — `abun(mixture = "negbin")` errors until the
    NB marginal sum lands in tulpa's N-mixture kernel.
  - **N-mixture NUTS** — no HMC likelihood for N-mixture in tulpa yet.
- `multispecies_nmix()` (`ms_abun`) and the closed multi-season via season RE /
  AR1: after the negbin gap closes.
- spAbundance / unmarked benchmarks.

**Phase 3 — Nested-Laplace extensions in `tulpa`**

- Beta + lognormal families.
- Joint multi-likelihood fits with shared latent fields.
- Multi-block latent priors.
- PC priors as defaults.
- `predict_at` row API.

**Phase 4 — Other observation processes**

- `distance()`, `removal()`, `false_positive()`, `dynamic_nmix()`.

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
