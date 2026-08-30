# tulpaObs

*false zeros from imperfect detection*

[![Lifecycle: stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/gcol33/tulpaObs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/tulpaObs/actions/workflows/R-CMD-check.yaml)
[![smoke](https://github.com/gcol33/tulpaObs/actions/workflows/smoke.yaml/badge.svg)](https://github.com/gcol33/tulpaObs/actions/workflows/smoke.yaml)

**Hierarchical observation models on the [tulpa](https://github.com/gcol33/tulpa) inference engine: occupancy, N-mixture abundance, distance, removal, cover.**

A site reads zero because the species is absent, or because you visited on the wrong day.
`tulpaObs` fits the models that hold those two apart: a latent ecological state under an
imperfect detection process, with the state integrated out exactly so occupancy and
detection come back as separate estimates. One fitter, one family object per model, and
the spatial fields, temporal structure, random effects, and latent factors written inside
the formula.

```r
library(tulpaObs)

fit <- tobs(
  ~ elevation + spatial(graph = adj), # occupancy: fixed effects + a spatial field
  data      = camera_data,
  family    = occu(),                # MacKenzie et al. (2002) single-season
  detection = ~ wind + time,         # the detection process
  y         = y                      # N x J detection-history matrix
)

summary(fit)
predict(fit, newdata = grid)         # marginalized per-site psi
```

## One fitter, every family

`tobs()` is the only fitter; the family object chooses the model, and the response shape
follows from it.

```r
tobs(~ forest, data = sites, family = abun(K_max = 50), detection = ~ effort, y = counts)
tobs(~ moisture, data = plots, family = cover(response = "beta"), y = pct_cover)
tobs(~ elev, data = sites, family = ms_occu(), detection = ~ effort, y = y_array,
     species = species_id)
```

| Constructor               | Model                                          | Reference                  |
|---------------------------|------------------------------------------------|----------------------------|
| `occu()`                  | Single-season occupancy                        | MacKenzie et al. (2002)    |
| `dyn_occu()`              | Dynamic occupancy (colonization / extinction)  | MacKenzie et al. (2003)    |
| `t_occu()`                | Multi-season occupancy, AR1 year effect        | spOccupancy `tPGOcc`       |
| `int_occu()`              | Integrated multi-source occupancy              | Miller et al. (2019)       |
| `dyn_int_occu()`          | Multi-season integrated occupancy              | spOccupancy `tIntPGOcc`    |
| `fp_occu()`               | False-positive occupancy                       | Miller et al. (2011)       |
| `royle_nichols()`         | Abundance-induced detection heterogeneity      | Royle & Nichols (2003)     |
| `occu_ttd()`              | Time-to-detection occupancy                    | Garrard et al. (2008)      |
| `occu_multi()`            | Multi-species co-occurrence occupancy          | Rota et al. (2016)         |
| `jsdm()`                  | Joint species distribution (no detection)      |                            |
| `abun()`                  | N-mixture abundance (Poisson / negbin)         | Royle (2004)               |
| `dyn_abun()`              | Open (Dail-Madsen) N-mixture                   | Dail & Madsen (2011)       |
| `count()`                 | Count / relative-abundance GLMM (no detection) | spAbundance `abund`        |
| `ms_count()`              | Community relative-abundance GLMM (no detection) | spAbundance `msAbund`    |
| `distance()`              | Binned distance sampling                       | Buckland et al. (2001)     |
| `distsamp_open()`         | Open-population distance sampling              | Sollmann et al. (2015)     |
| `gdistremoval()`          | Joint distance + removal                       | Amundson et al. (2014)     |
| `removal()`               | Removal sampling                               | Dorazio et al. (2005)      |
| `double_observer()`       | Double-observer abundance                      | Nichols et al. (2000)      |
| `cover()`                 | Cover hurdle (beta / lognormal / ordinal)      |                            |
| `occu_categorical()`      | Presence + nominal-class hurdle                |                            |
| `occu_cover()`            | Joint occupancy + cover hurdle                 |                            |
| `occu_multiscale_cover()` | Three-level (cell / plot / visit) + cover      | Nichols et al. (2008)      |
| `ms_occu()`               | Community occupancy                            | spOccupancy `msPGOcc`      |
| `ms_dyn_occu()`           | Community dynamic occupancy                    |                            |
| `ms_int_occu()`           | Community integrated occupancy                 |                            |
| `ms_abun()`               | Community N-mixture                            | spAbundance `msNMix`       |
| `ms_distance()`           | Community distance sampling                    | Sollmann et al. (2016)     |
| `ms_occu_cover()`         | Community joint occupancy + cover              |                            |

Every family carries the same structured-effect vocabulary, so any of them composes with
the engine's spatial, temporal, and latent structure.

## Structured effects inside the formula

Structured terms are written *in* the formula, the way `lme4`, `mgcv`, and `INLA` do. The
parser resolves bare symbols against `data` columns, so a spatial field is one term in the
formula:

```r
~ elevation +
  spatial(graph = adj) +               # areal: icar / bym2 / car / car_proper
  spatial(lon, lat, model = "spde") +  # continuous: spde / gp / multiscale_gp
  temporal(year, type = "ar1") +       # temporal: ar1 / rw1 / rw2 / iid
  (1 | observer)                       # random intercept (lme4 bar syntax)
```

`spatial()` and `temporal()` are the two front doors: `model =` picks the spatial
model (`"icar"` by default), `type =` the temporal one.

A coefficient that varies over space is a spatial field, so it is written under
`spatial()` and `model =` says whether space is a graph or a set of coordinates.
The areal form is a bar over graph nodes and has the broad family coverage:

```r
~ spatial(~ 1 + year || cell, graph = adj)   # areal SVC, method = "nested_laplace"
```

The continuous form names the varying coefficients directly and asks for a range
prior (`P(range < r0) = alpha`):

```r
~ elevation + spatial(lon, lat, model = "svc",
                      coefficients = "elevation", prior_range = c(50, 0.05))
```

The continuous form fits on single-season `occu()` under `laplace` /
`nested_laplace` / `nuts`, and on `removal()`, `distance()`, `fp_occu()` and
`dyn_abun()` under `laplace` / `nested_laplace`; the surfaces come back as
`fit$svc_field`. Elsewhere the term errors and points at the areal bar.
`svc(lon, lat, coefficients = )` is the direct constructor, as `icar()` is for
`spatial(model = "icar")`.

A term enters whichever process it is written in, so `detection = ~ (1 | observer)` puts
the random effect on detection. `copy("id")` shares one realization across both processes.
`?tobs_terms` documents every term and its arguments; `?svc` and `?icar` land there too.
Bar syntax follows `lme4` throughout: `(1 | site)`, `(x | site)` for a correlated random
slope, `(x || site)` uncorrelated, `(1 | g:h)` and `(1 | g/h)` for crossed and nested
factors.

Advanced usage: each spatial model also has a direct constructor, and `spatial()` forwards
its arguments to it, so the two forms build the same term.

```r
icar(graph = adj)                     # spatial(graph = adj)
bym2(graph = adj)                     # spatial(graph = adj, model = "bym2")
spde(lon, lat)                        # spatial(lon, lat, model = "spde")
gp(lon, lat, cov = "matern",          # spatial(lon, lat, model = "gp", cov = "matern",
   prior_range = c(0.1, 0.05))        #         prior_range = c(0.1, 0.05))
```

Model-specific arguments are checked against the model named in `model =`, so an argument
that model does not take is an error. `gp()` requires a PC prior on the range:
`prior_range = c(r0, alpha)` encodes `P(range < r0) = alpha` in the units of the
coordinates.

## The latent state integrates out

An occupancy site contributes a mixture over its latent state: occupied with probability
`psi` and detected on a Bernoulli schedule, or unoccupied and silent. `tulpaObs` evaluates
that marginalized likelihood directly, and the same holds for the N-mixture abundance `N`,
so `psi` and `p` come back as separate processes. Fitting the detection histories as one
stacked binomial identifies only their product `psi * p`; the
[occupancy-vs-INLA vignette](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy-vs-inla.Rmd)
runs that head-to-head on simulated data with known truth.

## Choose the inference route

`method` names one fully specified route. Laplace is fast and the default; the stochastic
corrections debias it where it is biased; NUTS gives the full posterior and fits every
structure.

```r
tobs(..., method = "laplace")         # EM + Laplace, deterministic
tobs(..., method = "laplace_sla")     # skew-corrected marginals, adds a `skew` column
tobs(..., method = "laplace_gibbs")   # post-EM Gibbs chain, Rubin-pooled
tobs(..., method = "nested_laplace")  # multi-block: areal fields, joint cover hurdle
tobs(..., method = "nuts")            # HMC, reports split-Rhat and bulk / tail ESS
```

Occupancy families also take `method = "pg_gibbs"`, an exact Polya-Gamma Gibbs sampler that
recovers the community variance components the Laplace EM attenuates. `fit$method` records
the resolved route, and an unsupported combination errors with a pointer to the methods
that family accepts.

Random-effect variance components carry the Laplace small-cluster bias for binary data (the
`glmer` `nAGQ = 1` regime). An adaptive Gauss-Hermite refinement on the exact per-group
marginal runs by default (`re.aghq`), cutting the occupancy per-group-`n = 8` sigma bias
from ~18% to ~4% and the detection-arm bias from ~70% to ~1%, both against NUTS.

## The data shape you already have

Detection histories arrive as matrices, as long tables with one row per site-visit, or as
species arrays. Three helpers build the same `tobs_data` object from any of them:

```r
tobs_format(y, occ.covs = sites, det.covs = list(wind = wind_mat), coords = xy)
tobs_data(visits_df, y = "detected", site = "site_id", visit = "visit")
tobs_format_ms(y_array, occ.covs = sites, species_names = spp)

occu_cover_inputs(plot_df, site = "plot", visit = "visit",
                  response = "present", y_pos = "cover")
```

`summary()`, `plot()`, and `print()` on a `tobs_data` report naive occupancy and detection,
per-visit rates, completeness, and a detection map when coordinates are present. Three
fixed-seed example datasets ship with a known `truth` attached: `peatland_occu`,
`foray_counts`, and `meadow_cover`.

## Simulate, fit, recover

Every family has a matching simulator that returns `y`, `data`, and the `truth` that
generated them, so a model can be checked against known parameters before it meets real
data. The package's own recovery tests are built on these:

```r
sim <- simulate_occu(N = 200, J = 4,
                     beta_occ = c(-0.5, 1.2, 0.4),
                     beta_det = c(0.2, -0.8),
                     seed = 1)

fit <- tobs(~ occ_cov1 + occ_cov2, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y)

coef(fit)
sim$truth
```

`simulate_abun()`, `simulate_cover()`, `simulate_dyn_occu()`, `simulate_ms_occu()`,
`simulate_distance()`, and the rest cover the remaining families, including a
`simulate_cover_joint()` that shares a demeaned BYM2 field across both hurdle arms.

## The usual R surface

A fit is an ordinary R model object, whichever route produced it:

```r
coef(fit); confint(fit); vcov(fit); logLik(fit)
fitted(fit); residuals(fit); predict(fit); simulate(fit)
ranef(fit); tidy(fit); glance(fit); update(fit)
converged(fit); convergence(fit)
```

`fitted()` returns `list(psi, p, z)` at the posterior mean, `ranef()` per-group BLUPs with
standard errors, and `predict()` works in-sample, on a design matrix (`X.0=`), or on terms
(`terms=`, which returns a `tobs_prediction` with a `plot()` method). Spatial fits
interpolate to new coordinates through `tobs_predict_spatial()`. `$` accessors
(`beta.samples`, `psi.samples`, `z.samples`, `run.time`) follow spOccupancy naming.

## Diagnostics and identifiability

```r
check_model(fit, coords = xy)    # roll-up: sampler, WAIC, PPC, zero-inflation, Moran's I
tobs_check_id(fit$model)         # confounding, low detection, sparse data

waic(fit); dic(fit); cpo(fit)
ppc(fit, fit.stat = "chi-squared")
test_dispersion(fit); test_zero_inflation(fit); test_outliers(fit)
test_uniformity(pit_residuals(fit))
```

`waic()` and `loo()` are the loo generics, so `loo::loo_compare()` reads a `tobs_fit`
directly.

Moran's I, Durbin-Watson, variograms, and `spatial_range` / `temporal_corr` for the latent
structure come through from tulpa. `tobs_stack()` LOO-weights two or more fits on the same
observations into an ensemble whose `predict()` returns the combined predictive.

## Longitudinal and community extras

`within_between()` splits a covariate into its per-group mean and within-group deviation, so
`~ x_btw + x_wtn` separates the cross-site gradient from change inside a site across a
resurvey:

```r
d <- within_between(plots, group = "plot", vars = c("temperature", "grazing"))
tobs(~ temperature_btw + temperature_wtn, data = d, family = cover(), y = d$cover)
```

For community fits, `tobs_richness()` returns posterior species richness and
`tobs_associations()` the residual species-association matrix from the latent factors.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/tulpaObs")            # latest from GitHub
pak::pak("gcol33/tulpaObs@v0.1.1")     # a specific tagged release
```

Tagged releases are listed at <https://github.com/gcol33/tulpaObs/releases>.

pak resolves the dependency tree, pulling `tulpa` and `tulpaMesh` from GitHub (declared in
`Remotes:`). A C++17 toolchain is needed (Rtools on Windows, Xcode CLI tools on macOS,
`r-base-dev` on Linux); both `tulpa` and `tulpaObs` compile their backends on first install.

## Documentation

Getting started

- [Quickstart](https://github.com/gcol33/tulpaObs/blob/main/vignettes/quickstart.Rmd)
- [Data formatting](https://github.com/gcol33/tulpaObs/blob/main/vignettes/data-formatting.Rmd)
- [Occupancy](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy.Rmd)
- [Occupancy vs INLA](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy-vs-inla.Rmd)

Model families

- [Abundance (N-mixture)](https://github.com/gcol33/tulpaObs/blob/main/vignettes/abundance.Rmd)
- [Dynamic occupancy](https://github.com/gcol33/tulpaObs/blob/main/vignettes/dynamic-occupancy.Rmd)
- [Integrated occupancy](https://github.com/gcol33/tulpaObs/blob/main/vignettes/integrated-occupancy.Rmd)
- [Community models](https://github.com/gcol33/tulpaObs/blob/main/vignettes/community-models.Rmd)
- [Cover hurdle](https://github.com/gcol33/tulpaObs/blob/main/vignettes/cover-hurdle.Rmd)
- [Cover hurdle, motivated](https://github.com/gcol33/tulpaObs/blob/main/vignettes/cover-hurdle-motivate.Rmd)
- [Cover hurdle vs INLA](https://github.com/gcol33/tulpaObs/blob/main/vignettes/cover-hurdle-vs-inla.Rmd)

Structure

- [Spatial occupancy](https://github.com/gcol33/tulpaObs/blob/main/vignettes/spatial-occupancy.Rmd)
- [Spatial occupancy with SPDE](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy-spatial-spde.Rmd)
- [Random effects](https://github.com/gcol33/tulpaObs/blob/main/vignettes/random-effects.Rmd)
- [Diagnostics](https://github.com/gcol33/tulpaObs/blob/main/vignettes/diagnostics.Rmd)

[API.md](https://github.com/gcol33/tulpaObs/blob/main/API.md) is the full reference: every
argument of `tobs()`, the term registry, the method table, the `control` list, and the
complete export list.

## Status

The public API (`tobs()`, the family constructors, and the `tobs_*` S3 classes) is
stable. Internal entry points (`.tobs_build_model()`, `.tobs_fit_model()`,
`.tobs_laplace()`) and the C++ surface follow tulpa's ABI version and may change between
releases.

## Support

> "Software is like sex: it's better when it's free."
>
> Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools
should be free and open. I started these projects for my own work and figured others
might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks.
It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (see the LICENSE file)

## Citation

```bibtex
@software{tulpaObs,
  author = {Colling, Gilles},
  title  = {tulpaObs: Hierarchical Latent-State Observation Models via tulpa},
  year   = {2026},
  url    = {https://github.com/gcol33/tulpaObs}
}
```
