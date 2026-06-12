# tulpaObs

> Small exact engines for scientific computing in R.

*false zeros from imperfect detection*

[![Lifecycle: experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Latent-state observation models (occupancy, N-mixture abundance, distance, removal, cover) on the [tulpa](https://github.com/gcol33/tulpa) engine.**

A site can read zero because the species is absent, or because you were there on the wrong day. `tulpaObs` fits the hierarchical models that separate the two: a latent ecological state (occupied, abundant, present-and-covering) under an imperfect detection process. You write one model menu's worth of likelihoods, and they compose with the spatial fields, temporal structure, random effects, and latent factors that the tulpa engine already carries.

```r
library(tulpaObs)

fit <- tobs(
  ~ elevation + icar(graph = adj),   # occupancy: fixed effects + a spatial field
  data      = camera_data,
  family    = occu(),                # MacKenzie et al. (2002) single-season
  detection = ~ wind + time,         # the detection process
  y         = y                      # N x J detection-history matrix
)

summary(fit)
predict(fit, newdata = grid)         # marginalized per-site psi, held-out sites interpolated
```

## The likelihood INLA can't write row-wise

A single-season occupancy site contributes a *mixture* over its latent state: occupied with probability `psi` and detected on a Bernoulli schedule, or unoccupied and silent. That mixture is not a row-wise GLM likelihood, so `INLA::inla()` cannot fit it without writing a custom `inla.rgeneric`. The common workaround stacks the detection histories into one binomial, which only ever identifies the product `psi * p` and confounds occupancy change with detection change.

`tulpaObs` implements the marginalized likelihood directly. The latent state integrates out exactly (the occupancy `z`, the N-mixture abundance `N`), so the fit recovers `psi` and `p` as separate processes. The `occupancy-vs-inla` vignette runs the head-to-head against the binomial workaround on simulated data with known truth.

## Two inference backends, one model

Laplace by default, with Gibbs / multiple-imputation debiasing where the approximation is biased; HMC/NUTS on demand for a full posterior. Both backends share the same model, the same data, and the same S3 surface:

```r
coef(fit); confint(fit); vcov(fit); logLik(fit)
fitted(fit); residuals(fit); predict(fit); simulate(fit)
ranef(fit); tidy(fit); glance(fit); update(fit)
```

The random-effect variance components carry the Laplace small-cluster bias for binary data (the `glmer` `nAGQ = 1` regime); an adaptive Gauss-Hermite refinement (`re.aghq`, on by default) corrects it back toward the NUTS estimate, and NUTS remains the calibrated route for the full covariance.

## The model menu

`tobs()` is the single dispatcher; the family object chooses the model.

| Constructor    | Model                                          | Reference                  |
|----------------|------------------------------------------------|----------------------------|
| `occu()`       | Single-season occupancy                        | MacKenzie et al. (2002)    |
| `dyn_occu()`   | Dynamic occupancy (colonization / extinction)  |                            |
| `ms_occu()`    | Multi-species / community occupancy            |                            |
| `int_occu()`   | Integrated multi-source occupancy              |                            |
| `jsdm()`       | Joint species distribution                     |                            |
| `abun()`       | N-mixture abundance (Poisson / negbin)         | Royle (2004)               |
| `ms_abun()`    | Community N-mixture                            | spAbundance `msNMix`       |
| `dyn_abun()`   | Dynamic N-mixture                              |                            |
| `distance()`   | Distance sampling                              |                            |
| `removal()`    | Removal sampling                               |                            |
| `fp_occu()`    | False-positive occupancy                       |                            |
| `cover()`      | Hurdle-Beta / lognormal cover                  |                            |

`unmarked`, `spOccupancy`, and `spAbundance` cover these models with a fixed per-model menu of structure. Here the same likelihoods sit on the tulpa engine, so any of them composes with the engine's latent structure.

## Structured effects inside the formula

Structured terms are written *in* the formula, the way `lme4`, `mgcv`, and `INLA` do:

```r
~ elevation +
  icar(graph = adj) +          # spatial field: icar / bym2 / gp / spde
  temporal(year, type = "ar1") +  # temporal: rw1 / rw2 / ar1
  svc(slope, graph = adj) +    # spatially varying coefficient
  (1 | observer)               # random intercept (lme4 bar syntax)
```

A term enters whichever process it is written in (occupancy or detection); `copy("id")` shares one realization across both. Bar syntax follows `lme4`: `(1 | site)` is `re(site)`, `(x | site)` a correlated random slope, `(x || site)` an uncorrelated one. `spatial(..., model = ...)` is a single-verb umbrella over the areal and continuous spatial terms.

## Diagnostics

WAIC, posterior predictive checks, PIT residuals, over-dispersion and zero-inflation tests, Moran's I, Durbin-Watson, variograms, and `spatialRange` / `temporalCorr` for the latent structure.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/tulpaObs@v0.0.1")     # tagged release
pak::pak("gcol33/tulpaObs")            # development version
```

pak resolves the dependency tree, pulling `tulpa` and `tulpaMesh` from GitHub (declared in `Remotes:`). A C++17 toolchain is needed (Rtools on Windows, Xcode CLI tools on macOS, `r-base-dev` on Linux); both `tulpa` and `tulpaObs` compile their backends on first install.

## Documentation

- [Occupancy vs INLA](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy-vs-inla.Rmd)
- [Cover hurdle vs INLA](https://github.com/gcol33/tulpaObs/blob/main/vignettes/cover-hurdle-vs-inla.Rmd)
- [Spatial occupancy with SPDE](https://github.com/gcol33/tulpaObs/blob/main/vignettes/occupancy-spatial-spde.Rmd)
- [Cover hurdle, motivated](https://github.com/gcol33/tulpaObs/blob/main/vignettes/cover-hurdle-motivate.Rmd)

## Status

Experimental. The public API (`tobs()` plus the family constructors and `tobs_*` S3 classes) is stable. Internal entry points (`.tobs_build_model()`, `.tobs_fit_model()`, `.tobs_laplace()`) and the C++ surface follow tulpa's ABI version and may change between releases.

## License

MIT (see the LICENSE file)
