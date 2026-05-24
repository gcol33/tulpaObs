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

### Working

| Constructor       | Model                          | Response `y`                       | Default engine |
|-------------------|--------------------------------|------------------------------------|----------------|
| `occu()`          | Single-season occupancy        | `N x J` detection matrix (0/1/NA)  | laplace |
| `dyn_occu()`      | Dynamic (multi-season) HMM     | `N x J x T` array                  | laplace |
| `ms_occu()`       | Multispecies / community       | `S x N x J` array                  | laplace |
| `int_occu()`      | Integrated multi-source        | list of `N_s x J_s` matrices       | laplace |
| `jsdm()`          | Joint species distribution     | `N x S` presence matrix (no detection) | nuts |
| `cover(positive)` | Vegetation cover hurdle        | length-`N` cover in `[0, 1]`       | laplace |

`cover(positive = c("beta", "lognormal"))` chooses the positive-part likelihood.

### Planned (error informatively until implemented)

`abun(K_max, family)`, `ms_abun(...)`, `dyn_abun(...)` (Dail-Madsen),
`distance(key)`, `removal()`, `fp_occu()` (false-positive occupancy).

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
| `(x \| g)`          | correlated int + slope | **NUTS only**  |
| `(1 + x + z \| g)`  | multi-slope, correlated| **NUTS only**  |
| `(1 \| g:h)`, `(1 \| g/h)` | crossed / nested | one `re()` per implied factor |

The Laplace methods fit iid intercept RE and *uncorrelated* slopes on the
occupancy predictor of a single-season model. Correlated slopes, RE on
detection, RE + spatial, and RE on non-single families error from
`.validate_re_laplace()` with a pointer to `method = "nuts"` rather than being
dropped (correlated slopes are blocked upstream by tulpa's diagonal RE
precision, gcol33/tulpa#28). Deterministic Laplace variance estimates for
binary occupancy carry the usual small-cluster (PQL) bias; NUTS is the
calibrated route.

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
| `"nested_laplace"`     | nested Laplace | Gaussian | none | Multi-block (non-conjugate hyperpriors / cover-hurdle joint). Wired for single-season `occu()` and `cover()`; other families fall back with a message. |
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

Laplace controls (`method = "laplace"` / `"laplace_sla"` / `"nested_laplace"` /
`"nested_laplace_sla"`): `max.iter`, `tol`, `damping`, `sigma.beta`,
`sigma.re.scale`. Stochastic-correction controls (`"laplace_gibbs"` /
`"laplace_mi"`): `n.gibbs` / `n.imputations` (Rubin-pooled draw count) and
`seed` (stored on `$seed`).

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
to recover the unpenalised MAP. `tobs_priors()` is the generic prior container
(`beta.normal`, `alpha.normal`, `sigma.sq.psi`, `sigma.sq.p`).

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
inherited from `tulpa::tulpa_fit`. tulpaObs overrides or adds:

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
| `$`                             | spOccupancy-compatible accessors: `beta.samples`, `alpha.samples`, `psi.samples`, `p.samples`, `z.samples`, `run.time` |

`tobs_marginal_effect(object, covariate, process, n_points)` and
`tobs_richness(object)` (community models) are standalone exported helpers;
`predict(..., terms=)` returns a `tobs_prediction` with a `plot()` method.

---

## 10. Diagnostics

| Function                       | Returns |
|--------------------------------|---------|
| `tobs_waic()`                  | `waic`, `elpd`, `p_waic`, `lppd` |
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

## Quick reference: full export list

**Fitter** `tobs` · **Families** `occu` `dyn_occu` `ms_occu` `int_occu` `jsdm`
`cover` `abun` `ms_abun` `dyn_abun` `distance` `removal` `fp_occu` `obs_family`
· **Priors** `occu_priors` `tobs_priors` · **Data** `tobs_format` `tobs_data`
`tobs_format_ms` · **Simulators** `simulate_occu` `simulate_ms_occu`
`simulate_dyn_occu` `simulate_int_occu` `simulate_dyn_ms_occu`
`simulate_int_ms_occu` `simulate_cover` `simulate_cover_joint` · **Diagnostics**
`tobs_waic` `tobs_ppc` `tobs_pit_residuals` `tobs_test_uniformity`
`tobs_test_dispersion` `tobs_test_zero_inflation` `tobs_test_outliers`
`tobs_check` `tobs_check_id` · **Prediction / effects** `tobs_predict_spatial`
`tobs_marginal_effect` `tobs_richness` `within_between` · **Generic re-export**
`ranef`

> Structured terms (`icar`, `bym2`, `car`, `car_proper`, `gp`, `multiscale_gp`,
> `spde`, `re`, `temporal`, `svc`, `latent`, `copy`) are **not exported** — they
> are only meaningful inside a `tobs()` formula and resolved through the
> internal registry.
