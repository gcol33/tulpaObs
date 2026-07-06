# tulpaObs API

Hierarchical latent-state observation models (occupancy / abundance / cover)
fit through the [`tulpa`](https://github.com/gcol33/tulpa) inference engine.

The public surface is small and uniform: **one fitter** (`tobs()`), **one
family object per model type**, and **structured effects written inside the
formula** (lme4 / mgcv / INLA style). Everything else is S3 methods, data
helpers, simulators, and diagnostics.

```r
fit <- tobs(
  formula   = ~ elev + forest,   # state process (psi / abundance / latent cover)
  data      = sites,
  family    = occu(),            # model type
  detection = ~ effort,          # detection process
  y         = y_matrix           # response (shape depends on family)
)
```

---

## 1. The fitter: `tobs()`

```r
tobs(formula, data, family,
     detection = NULL, y = NULL, visits = NULL,
     method = c("auto", "laplace", "laplace_sla", "laplace_gibbs",
                "laplace_mi", "nested_laplace", "nested_laplace_sla", "nuts"),
     priors = NULL, control = list(), ...)
```

| Argument    | Meaning |
|-------------|---------|
| `formula`   | State-process formula. Structured terms (spatial / RE / temporal / SVC / latent) go *inside* it — see §3. |
| `data`      | Site-level covariate data frame, `nrow(data) == nrow(y)`. |
| `family`    | A `tobs_family` object (`occu()`, `abun()`, ...). Required. |
| `detection` | Detection-process formula. Required for `occu()`/`abun()`; ignored by `jsdm()` and `cover()`. |
| `y`         | Response; shape depends on family (see §2). |
| `visits`    | Optional visit-level detection covariates (named list of `[n_sites, max_visits]` matrices, or a long data frame in site-major order). |
| `method`    | Inference route — one fully-specified path, not orthogonal knobs (see §4). |
| `priors`    | `occu_priors()` / list for the Laplace methods; `FALSE` disables the default penalty; forwarded to tulpa under NUTS. The `laplace_gibbs`/`laplace_mi` routes disable it automatically. |
| `control`   | Low-level engine controls (see §5). |
| `...`       | Family-specific args, e.g. `species=` (`ms_occu`), `col_formula=`/`ext_formula=` (`dyn_occu`). |

The resolved route is recorded on `fit$method` for provenance.

Returns an object of class `c("tobs_fit", "<family>_fit", "tulpa_fit")`.

### Family-specific `...` arguments

| Family       | Extra required `...` |
|--------------|----------------------|
| `dyn_occu()` | `col_formula = ~ ...`, `ext_formula = ~ ...` (colonisation / extinction) |
| `ms_occu()`  | `species = <id>` |
| `jsdm()`     | `species = <id>` (optional) |

---

## 2. Family constructors

A family object carries the latent-state type, observation likelihood,
replicate requirement, default engine, and status. Families do not fit —
`tobs()` reads the object and dispatches.

All families below carry `status = "working"` and are recovery-tested.

### Occupancy

| Constructor       | Model                          | Response `y`                       | Default engine |
|-------------------|--------------------------------|------------------------------------|----------------|
| `occu()`          | Single-season occupancy        | `N x J` detection matrix (0/1/NA)  | laplace |
| `dyn_occu()`      | Dynamic (multi-season) HMM     | `N x J x T` array                  | laplace |
| `ms_occu()`       | Multispecies / community       | `S x N x J` array                  | laplace |
| `ms_dyn_occu()`   | Dynamic multispecies           | `S x N x J x T` array              | laplace |
| `int_occu()`      | Integrated multi-source        | list of `N_s x J_s` matrices       | laplace |
| `ms_int_occu()`   | Integrated multispecies        | list of `S x N_s x J_s` arrays     | laplace |
| `jsdm()`          | Joint species distribution     | `N x S` presence matrix (no detection) | nuts |
| `fp_occu()`       | False-positive occupancy (Miller 2011) | `N x J` in `{0,1,2}`      | laplace |
| `occu_categorical(classes)` | Multi-state categorical detection | `N x J` class matrix    | laplace |

### Abundance

| Constructor                     | Model                        | Response `y`                | Default engine |
|---------------------------------|------------------------------|-----------------------------|----------------|
| `abun(K_max, mixture)`          | N-mixture (Royle 2004)       | `N x J` counts              | laplace |
| `ms_abun(K_max, mixture)`       | Community N-mixture (`msNMix`) | `S x N x J` counts        | laplace |
| `dyn_abun(K_max, mixture)`      | Open N-mixture (Dail-Madsen) | `N x J x T` counts          | laplace |
| `distance(key, transect, cutpoints, ...)` | Binned distance sampling | `N x B` bin counts   | laplace |
| `removal(K_max, mixture)`       | Sequential removal           | `N x K` removal-pass counts | laplace |

`mixture = c("poisson", "negbin")` on the count families; `distance(key = c("halfnorm", "hazard"))`.

### Cover

| Constructor                | Model                              | Response `y`                | Default engine |
|----------------------------|------------------------------------|-----------------------------|----------------|
| `cover(response, breaks)`  | Vegetation cover hurdle            | length-`N` cover in `[0, 1]`| laplace |
| `occu_cover(response)`     | Joint occupancy + cover hurdle     | detection matrix + cover    | laplace |
| `ms_occu_cover(response)`  | Community joint occu + cover       | per-species                 | laplace |
| `occu_multiscale_cover(response)` | Three-level (cell/plot/visit) + cover | plot detection + cover | nested_laplace |

`cover(response = c("beta", "beta_oi", "lognormal", "lognormal_trunc", "ordinal"))`
chooses the positive-part likelihood; `"ordinal"` also needs `breaks =` (interior
Braun-Blanquet class boundaries). The joint cover families take
`response = c("beta", "lognormal")`.

`obs_family()` is the low-level constructor behind all of these (exported but
`@keywords internal`; use the concrete constructors).

---

## 3. Structured terms (inside the formula)

Spatial fields, random effects, temporal structure, spatially varying
coefficients, and latent factors are **terms written in the process formula**,
not arguments to `tobs()`. The parser resolves bare symbols against `data`
columns exactly like `s()` does inside `gam()`.

A term enters whichever linear predictor it appears in (`formula` vs
`detection`). To share one realization across both predictors, tag it with
`id = "u"` in one formula and write `copy("u")` in the other.

```r
tobs(~ elev + bym2(graph = adj) + re(region),
     detection = ~ effort + re(observer),
     data = sites, family = occu(), y = y)
```

### Registry

| Term                              | Effect |
|-----------------------------------|--------|
| `icar(graph)`                     | Intrinsic CAR over an adjacency matrix |
| `bym2(graph, scale_factor)`       | BYM2 (ICAR + IID) reparameterization |
| `car(graph)` / `car_proper(graph)`| (Improper / proper) CAR areal field |
| `gp(lon, lat, ...)`               | NNGP-approximated Gaussian process |
| `multiscale_gp(lon, lat, ...)`    | Two-scale (local + regional) NNGP |
| `spde(lon, lat, ...)`             | Continuous Matern field via triangular mesh |
| `re(group, ...)`                  | Grouped random effect (intercept / slope / iid) |
| `temporal(time, ...)`             | AR1 / RW1 / RW2 / IID temporal field |
| `svc(lon, lat, indices)`          | Spatially varying coefficients on design columns |
| `latent(n_factors, ...)`          | Latent factors for community models |
| `copy("id")`                      | Share a named term's realization across processes |

Common term options:

- `re(group, type = c("intercept","slope","iid"), covariate = NULL, model = "iid", correlated = TRUE, intercept = TRUE)` — `covariate` may be a column name, several names, a `cbind()` matrix, or a one-sided formula (multi-slope block); `model` is `iid`/`ar1`/`rw1`/`rw2`.
- `temporal(time, type = c("ar1","rw1","rw2","iid"), group = NULL, cyclic = FALSE)`.
- `gp(..., cov = "exponential", nu = 1.5, nn = 15)` — `cov` is `exponential`/`matern`/`gaussian`/`spherical`.

Areal terms (`icar`/`bym2`/`car`/`car_proper`) accept `group_var =` to map
observations to graph nodes when the graph is over regions rather than rows.

### lme4 bar syntax (sugar over `re()`)

Bar terms are rewritten to `re()` calls on the AST before `terms()` runs, so
there is one parser and one term type:

| Bar form            | Equivalent             | Engine support |
|---------------------|------------------------|----------------|
| `(1 \| g)`          | `re(g)` (iid intercept)| Laplace + NUTS |
| `(x \|\| g)`        | uncorrelated slope     | Laplace + NUTS |
| `(0 + x \| g)`      | slope only, no intercept | Laplace + NUTS |
| `(1 + x \|\| g)`    | uncorrelated int+slope | Laplace + NUTS |
| `(x \| g)`          | correlated int + slope | Laplace + NUTS |
| `(1 + x + z \| g)`  | multi-slope, correlated| Laplace + NUTS (RE dim <= 3 for the AGHQ debias) |
| `(1 \| g:h)`, `(1 \| g/h)` | crossed / nested | one `re()` per implied factor |

Each form fits on EITHER the occupancy or the detection predictor (its position
in the formula sets the arm: `detection = ~ (1 | g)` puts the RE on detection).
The Laplace methods fit iid intercept RE, *uncorrelated* slopes, and *correlated*
slopes (a full RE covariance) on a single-season model. A single RE shared across
BOTH predictors, RE + visit-level detection, RE + spatial, and RE on non-single
families error from `.validate_re_laplace()` with a pointer to `method = "nuts"`
rather than being dropped. The raw EM variance components (sigma, correlation)
carry the Laplace small-cluster bias for binary data (the glmer nAGQ=1 regime,
not Breslow-Clayton PQL); by default (`re.aghq = TRUE`, `n.quad = 9`) they are
debiased after the EM by an adaptive Gauss-Hermite refinement on the exact
per-group marginal (the nAGQ > 1 fix), branching on the arm (an occupancy-arm RE
moves psi, a detection-arm RE moves p). It cuts the occupancy per-group-n = 8
sigma bias from ~18% to ~4% (matching NUTS), and the detection-arm bias from
~70% (only occupied sites inform p) to ~1%. A default LKJ (`re.lkj = 1.5`)
penalty regularizes a weakly-identified RE correlation off the `+-1` boundary
(without touching the marginal SDs); NUTS is available for a full posterior
treatment of the correlation.

---

## 4. Methods

`method` names a fully-specified route. Each maps internally to an orthogonal
`(engine, approx, correction)` triple — collapsed into one name so no argument
is silently ignored (e.g. there is no NUTS-with-SLA-marginal combination).

| `method`               | Engine | Marginal | Correction | Notes |
|------------------------|--------|----------|------------|-------|
| `"laplace"`            | Laplace | Gaussian | none | Default. EM + tulpa Laplace. Fast. |
| `"laplace_sla"`        | Laplace | skew (SLA) | none | Skew-corrected marginals; `summary()` gains a `skew` column. |
| `"laplace_gibbs"`      | Laplace | Gaussian | Gibbs | Post-EM Gibbs chain, Rubin-pooled. Runs unpenalised (prior auto-disabled, gcol33/tulpa#27); seedable + reproducible. |
| `"laplace_mi"`         | Laplace | Gaussian | MI | Post-EM multiple imputation, Rubin-pooled. Same prior/seed behaviour. |
| `"nested_laplace"`     | nested Laplace | Gaussian | none | Multi-block (non-conjugate hyperpriors, areal spatial fields, cover-hurdle joint). Supported by most families (`occu`, `dyn_occu`, `int_occu`, `ms_occu`, `jsdm`, `abun`, `ms_abun`, `removal`, `distance`, `fp_occu`, `dyn_abun`, `cover`, `occu_cover`, `ms_occu_cover`, `occu_multiscale_cover`). The per-family support set is `.tobs_family_methods`; an unsupported `method` errors with a pointer rather than silently downgrading. |
| `"nested_laplace_sla"` | nested Laplace | skew (SLA) | none | Nested Laplace with skew-corrected marginals. |
| `"nuts"`               | NUTS | — | — | HMC via tulpa. Fits every structure, incl. correlated slopes and stacked spatial+RE. Reports split-Rhat / bulk / tail ESS on `$convergence`. |
| `"auto"`               | — | — | — | Resolves to the family's `default_engine`. |

`fit$method` records the resolved route.

---

## 5. `control` list

Sampler controls (`method = "nuts"`):

| Name             | Default | Meaning |
|------------------|---------|---------|
| `n.iter`         | 2000    | Total iterations per chain (incl. warmup) |
| `n.warmup`       | 1000    | Warmup / adaptation iterations |
| `n.thin`         | 1       | Keep every `n.thin`-th post-warmup draw |
| `n.chains`       | 1       | Chains (offset seeds, pooled). Resolved seeds stored on `$seeds`. |
| `n.threads`      | 1       | Chains run in parallel (PSOCK cluster if `> 1`) |
| `adapt.delta`    | 0.8     | Target acceptance probability |
| `max.treedepth`  | 10      | NUTS maximum tree depth |
| `seed`           | 42      | Base RNG seed; chain `c` uses `seed + c - 1` |
| `n.seeds`        | 1       | Fit K seed-offset refits and LOO-stack them into a `tobs_stack` (see §9). Member `k` uses base seed `seed + k - 1`. |

Laplace controls (`method = "laplace"` / `"laplace_sla"` / `"nested_laplace"` /
`"nested_laplace_sla"`): `max.iter`, `tol`, `damping`, `sigma.beta`,
`sigma.re.scale`. Stochastic-correction controls (`"laplace_gibbs"` /
`"laplace_mi"`): `n.gibbs` / `n.imputations` (Rubin-pooled draw count) and
`seed` (stored on `$seed`).

`n.seeds > 1` is only meaningful for the stochastic routes (`"nuts"`,
`"laplace_gibbs"`, `"laplace_mi"`) — seed-variants would be identical under the
deterministic Laplace methods, which reject it. The seed-offset members are
statistically identical, so this is a Monte-Carlo robustness device with
roughly-uniform stacking weights; pass distinct fits to `tobs_stack()` for a
genuine model average.

---

## 6. Priors

`occu_priors()` builds the weakly-informative quadratic prior used by the
occupancy Laplace path (breaks the psi-p identifiability ridge at small `J`):

```r
occu_priors(p_intercept        = list(mean = 0, sd = 1.5),
            p_slope            = list(mean = 0, sd = 2.5),
            beta_occ_intercept = list(mean = 0, sd = 2),
            beta_occ_slope     = list(mean = 0, sd = 5))
```

Set any `sd = Inf` to disable that component. Pass `priors = FALSE` to `tobs()`
to recover the unpenalised MAP. For NUTS, the occupancy path does not apply a
user `priors` argument — the sampler uses the weakly-informative
`N(0, sigma.beta)` coefficient prior, and a supplied `priors` now raises a
warning rather than being silently ignored.

`cover_priors()` is the opt-in fixed-effect prior for the cover hurdle, with one
intercept/slope bucket per arm:

```r
cover_priors(occ_intercept = list(mean = 0, sd = 2),
             occ_slope     = list(mean = 0, sd = 2.5),
             pos_intercept = list(mean = 0, sd = 3),
             pos_slope     = list(mean = 0, sd = 2.5))
```

Cover has no psi-p-style ridge, so priors are off by default (a cover fit with
`priors = NULL` is unpenalised); the main use is taming perfect separation in
the occurrence arm at small `N`. Both fit paths thread the prior. On the
separate-Laplace path (`method = "laplace"`/`"laplace_sla"`, no spatial term)
both arms are penalised — occurrence and lognormal-positive through
`tulpa_laplace()`, beta-positive through `tulpa_laplace_beta()`. On the joint
nested-Laplace / spatial path the per-arm `cover_priors()` buckets are passed to
`tulpa_nested_laplace_joint()` as per-response `beta_prior_mean` /
`beta_prior_prec` (natural scale, applied to the autoscaled O(1) design), with
`priors = NULL`/`FALSE`/`"none"` leaving both arms unpenalised (gcol33/tulpaObs#54).
The `occu_cover()` joint path threads the prior the same way. Each prior
constructor (`occu_priors`, `cover_priors`) carries the natural parameters of
its family group — there is no generic prior object.

---

## 7. Data helpers

| Function                                          | Purpose |
|---------------------------------------------------|---------|
| `tobs_format(y, occ.covs, det.covs, coords, species)` | Build a `tobs_data` object from matrices/lists |
| `tobs_data(df, y, site, visit, occ.covs, det.covs, coords)` | Convert long format (one row per site-visit) |
| `tobs_format_ms(y, occ.covs, det.covs, coords, species_names)` | Multi-species (3D array or list of matrices) |

`summary()` / `plot()` / `print()` methods on `tobs_data` report naive
occupancy/detection, per-visit rates, completeness, and (with coordinates) a
detection map.

### Bundled example datasets

Three synthetic datasets (fixed-seed generative models with a known `truth`,
built by `data-raw/make_datasets.R`) ship for load-and-run examples:

| Dataset          | Family it feeds       | Shape |
|------------------|-----------------------|-------|
| `peatland_occu`  | `occu()`              | list: `y` (120 x 4), `occ.covs`, `det.covs`, `coords`, `truth` |
| `foray_counts`   | `abun()`              | list: `y` (100 x 3 counts), `occ.covs`, `det.covs`, `truth` |
| `meadow_cover`   | `cover()` / `within_between()` | data frame, 150 plots x (plot, year, year_c, moisture, grazing, cover) |

Load with `data(peatland_occu)`; see `?peatland_occu` for the generative model.

---

## 8. Simulators

All return a list with `y`, `data`, and `truth` (component truths for
parameter-recovery tests).

| Function                  | Generates |
|---------------------------|-----------|
| `simulate_occu()`         | Single-season occupancy |
| `simulate_ms_occu()`      | Multispecies occupancy (3D array) |
| `simulate_dyn_occu()`     | Dynamic occupancy (colonisation / extinction) |
| `simulate_int_occu()`     | Integrated multi-source (list of matrices + `site_maps`) |
| `simulate_dyn_ms_occu()`  | Dynamic multispecies (4D array) |
| `simulate_int_ms_occu()`  | Integrated multispecies (list of 3D arrays) |
| `simulate_cover()`        | Cover hurdle (optional exponential-kernel spatial field) |
| `simulate_cover_joint()`  | Cover hurdle with a **shared demeaned BYM2 field** (matches the joint nested-Laplace parameterisation) |

---

## 9. S3 methods on `tobs_fit`

Generic methods (`coef`, `confint`, `vcov`, `logLik`, `tidy`, `glance`) are
inherited from `tulpa::tulpa_fit`. `logLik(fit)` carries `df` (the fixed-effect
count) and `nobs` attributes, so `AIC(fit)` / `BIC(fit)` resolve through the
`stats` defaults. Note the two caveats: `df` counts only fixed effects (variance
components and spatial-field hyperparameters are not penalised), and under
`method = "nuts"` the log-likelihood is the mean posterior log-density (prior
included), not a maximised likelihood — so for hierarchical model comparison
prefer `tobs_waic()` / `tobs_dic()` / `tobs_cpo()`. tulpaObs overrides or adds:

| Method                          | Notes |
|---------------------------------|-------|
| `print()`                       | Family, dims, convergence, intercepts |
| `summary()`                     | Adds `skew` (SLA) and `rhat`/`ess_bulk`/`ess_tail` (NUTS) columns |
| `coef()`                        | Per-process list; appends visit-level `p_visit_*` coefficients |
| `ranef()`                       | Per-group BLUPs (`group`, `level`, `term`, `estimate`, `std.error`) |
| `fitted()`                      | `list(psi, p, z)` at posterior mean |
| `residuals(type=)`              | `deviance` (default) / `pearson` / `response`; `list(occ, det)` |
| `predict()`                     | In-sample / design-matrix (`X.0=`) / terms-based (`terms=`) |
| `simulate(nsim, seed)`          | Posterior replicate datasets (single-season) |
| `update(...)`                   | Refit with overridden controls; structured terms travel with the model |
| `nobs()`                        | Non-missing detection-history count |
| `converged()` / `convergence()` | Convergence flag / full convergence diagnostics (split-Rhat, ESS under NUTS) |
| `$`                             | spOccupancy-compatible accessors: `beta.samples`, `alpha.samples`, `psi.samples`, `p.samples`, `z.samples`, `run.time` |

`tobs_marginal_effect(object, covariate, process, n_points)` and
`tobs_richness(object)` (community models) are standalone exported helpers;
`predict(..., terms=)` returns a `tobs_prediction` with a `plot()` method.

`tobs_stack(..., method = c("stacking", "pseudobma"))` combines two or more
`tobs_fit` objects (fit to the same observation set, otherwise free to differ
in covariates / priors / `method`) into a LOO-weighted ensemble. Returns a
`tobs_stack` with `weights`, `fits`, `loo`, `comparison`, and `method`;
`predict()` / `fitted()` return the weight-combined predictive. Works across
every family. `tobs(..., control = list(n.seeds = K))` produces one of these
directly from seed-offset refits (§5).

---

## 10. Diagnostics

| Function                       | Returns |
|--------------------------------|---------|
| `tobs_waic()`                  | `waic`, `elpd`, `p_waic`, `lppd` |
| `tobs_dic()`                   | DIC, effective parameter count `p_D` |
| `tobs_cpo()`                   | Conditional predictive ordinate / LOO log-score (`loo.unit = "obs"`/`"cell"`) |
| `tobs_ppc(fit.stat, n.samples)`| Posterior predictive check + Bayesian p-value |
| `tobs_pit_residuals(n.samples)`| PIT residual vector |
| `tobs_test_uniformity(pit)`    | KS test of PIT residuals against uniform |
| `tobs_test_dispersion()`       | Observed vs simulated variance ratio + p |
| `tobs_test_zero_inflation()`   | Observed vs simulated all-zero sites + p |
| `tobs_test_outliers()`         | Count of sites outside simulated 95% envelope + p |
| `tobs_check(coords, n.samples)`| Roll-up report (sampler, WAIC, PPC, zero-inflation, Moran's I) |
| `tobs_check_id(model, fit)`    | Pre/post-fit identifiability checks (confounding, low detection, sparse data) |

Generic spatial/temporal diagnostics (`moranI`, `durbinWatson`, `variogram`,
`compare_models`) are inherited from tulpa. Many tobs diagnostics currently
support single-season models only.

---

## 11. Spatial prediction & utilities

- `tobs_predict_spatial(object, newcoords, newocc.covs, quantiles)` — occupancy
  at new coordinates, interpolating the fitted spatial field via IDW (k=5).
  Requires a fit with a spatial component.
- `within_between(data, group, vars, suffix, na.rm)` — Mundlak within/between
  decomposition: splits each `var` into a per-group mean (`<var>_btw`) and the
  within-group deviation (`<var>_wtn`), so `~ x_btw + x_wtn` separates
  cross-group from within-group association (MOTIVATE resurvey workflow).

---

## 12. Batch fits (multi-response)

Passing a multi-response `y` (e.g. a per-species list) to `tobs()` with a cover
family fits each response independently on the shared design and returns a
`tobs_batch` — a container of per-response `tobs_fit` objects, each byte-identical
to a separate single-response call.

| Function / method       | Purpose |
|-------------------------|---------|
| `occu_cover_inputs(data, site, visit, response, y_pos, ...)` | Build the ragged / dense occupancy + cover arms from a long (one row per site-visit) plot-level data frame, ready to hand to `tobs(family = occu_cover())`. |
| `tobs_get(x, species)`  | Extract one response's `tobs_fit` from a `tobs_batch` (by label or index). |
| `coef(<tobs_batch>)`    | Named list of per-response coefficient vectors. |
| `print(<tobs_batch>)`   | Species count, family, and where the per-response fits live (`$fits`). |

`tobs_associations(object)` returns the residual species-association matrix from
a JSDM / community-factor fit (the spatial-factor `ms_occu_cover()` loadings).
`occu_aggregation_scan()` sweeps cover-arm aggregation choices for a joint fit.

---

## Quick reference: full export list

**Fitter** `tobs` · **Families** `occu` `dyn_occu` `ms_occu` `ms_dyn_occu`
`int_occu` `ms_int_occu` `jsdm` `fp_occu` `occu_categorical` `abun` `ms_abun`
`dyn_abun` `distance` `removal` `cover` `occu_cover` `ms_occu_cover`
`occu_multiscale_cover` `obs_family` · **Priors** `occu_priors` `cover_priors` ·
**Data** `tobs_format` `tobs_data` `tobs_format_ms` `occu_cover_inputs` `tobs_get`
· **Simulators** `simulate_occu` `simulate_ms_occu` `simulate_dyn_occu`
`simulate_dyn_ms_occu` `simulate_int_occu` `simulate_int_ms_occu`
`simulate_occu_categorical` `simulate_abun` `simulate_ms_abun`
`simulate_dyn_abun` `simulate_distance` `simulate_removal` `simulate_fp_occu`
`simulate_cover` `simulate_cover_joint` `simulate_occu_cover`
`simulate_ms_occu_cover` `simulate_ms_occu_cover_spatial`
`simulate_occu_multiscale_cover` · **Diagnostics** `tobs_waic` `tobs_cpo`
`tobs_dic` `tobs_ppc` `tobs_pit_residuals` `tobs_test_uniformity`
`tobs_test_dispersion` `tobs_test_zero_inflation` `tobs_test_outliers`
`tobs_check` `tobs_check_id` · **Prediction / effects** `tobs_predict_spatial`
`tobs_marginal_effect` `tobs_richness` `tobs_associations`
`occu_aggregation_scan` `within_between` · **Ensembles** `tobs_stack` ·
**Convergence** `converged` `convergence` · **Generic re-exports** `ranef`
`tidy` `glance`

> Structured terms (`icar`, `bym2`, `car`, `car_proper`, `gp`, `multiscale_gp`,
> `spde`, `re`, `temporal`, `svc`, `latent`, `copy`) are **not exported** — they
> are only meaningful inside a `tobs()` formula and resolved through the
> internal registry.
