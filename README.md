# tulpaObs

[![Lifecycle: experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Hierarchical Latent-State Observation Models via tulpa**

`tulpaObs` fits Bayesian occupancy, abundance, distance, removal, and
cover models on top of the [tulpa](https://github.com/gcol33/tulpa)
engine. Models compose freely with spatial fields (ICAR, BYM2, NNGP),
temporal structure (RW1/RW2, AR1), random effects, spatially varying
coefficients, and latent factors.

## Quick start

```r
library(tulpaObs)

fit <- tobs(
  occu(detection ~ wind + time, occupancy ~ elevation + tobs_icar(adj)),
  data = camera_data
)

summary(fit)
predict(fit, newdata = grid)
```

`tobs()` is the single dispatcher; the family is chosen by the constructor
passed in:

| Constructor      | Model                                          |
|------------------|------------------------------------------------|
| `occu()`         | Single-season occupancy                        |
| `dyn_occu()`     | Dynamic occupancy (colonization / extinction)  |
| `ms_occu()`      | Multi-species / community occupancy            |
| `int_occu()`     | Integrated multi-source occupancy              |
| `jsdm()`         | Joint species distribution                     |
| `abun()`         | N-mixture abundance                            |
| `ms_abun()`      | Multi-species abundance                        |
| `dyn_abun()`     | Dynamic N-mixture                              |
| `distance()`     | Distance sampling                              |
| `removal()`      | Removal sampling                               |
| `fp_occu()`      | False-positive occupancy                       |
| `cover()`        | Hurdle-Beta cover                              |

## Inference

Laplace by default (with Gibbs / MI debiasing where required), HMC/NUTS
on demand. Both backends share the same model, data, and S3 surface
(`coef`, `confint`, `vcov`, `logLik`, `fitted`, `residuals`, `predict`,
`simulate`, `tidy`, `glance`, `ranef`, `update`, `summary`).

Diagnostics include WAIC, posterior predictive checks, PIT residuals,
dispersion and zero-inflation tests, Moran's I, Durbin–Watson, variograms,
and `spatialRange` / `temporalCorr` for the latent structure.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gcol33/tulpaMesh")
remotes::install_github("gcol33/tulpa")
remotes::install_github("gcol33/tulpaObs")
```

Requires a C++17 toolchain.

## Status

Experimental. The public API (`tobs()` + family constructors, `tobs_*`
S3 classes) is stable. Internal entry points (`.tobs_build_model()`,
`.tobs_fit_model()`, `.tobs_laplace()`) and the C++ surface follow
`tulpa`'s ABI version and may change between releases.
