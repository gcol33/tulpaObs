
# tulpaObs

*false zeros from imperfect detection*

<!-- badges: start -->

[![R-CMD-check](https://github.com/gcol33/tulpaObs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/tulpaObs/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/gcol33/tulpaObs/graph/badge.svg)](https://app.codecov.io/gh/gcol33/tulpaObs)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

tulpaObs fits Bayesian models for species surveys with imperfect
detection, such as site occupancy from repeated visits, N-mixture
abundance from replicated counts, distance and removal sampling, and
vegetation cover.

A site can read zero because the species is absent, or because it was
present and missed on every visit. Each model writes the ecological
state (occupied or not, how many individuals) and the detection process
as separate levels and integrates the latent state out of the
likelihood, so occupancy and detection come back as separate estimates.
Each model supplies its observation likelihood to the
[tulpa](https://CRAN.R-project.org/package=tulpa) engine, which fits it
by Laplace approximation, by nested Laplace integration over the
hyperparameters, or by the No-U-Turn sampler (NUTS). You choose the
model through a family object, write spatial, temporal and random
effects inside the formula, and pick the inference method for each fit.

tulpaObs supports:

- **Occupancy, abundance and cover models.** Single-season, dynamic,
  integrated, false-positive, time-to-detection and community occupancy;
  N-mixture, distance, removal and double-observer abundance; and a
  hurdle model for vegetation cover. One fitter, `tobs()`, takes a
  family constructor per model ([families in the
  reference](https://gillescolling.com/tulpaObs/reference/index.html)).
- **Structured effects inside the formula.** Areal and continuous
  spatial fields, temporal processes and lme4-style random effects enter
  the occupancy or the detection formula as terms.
- **Several inference routes.** Laplace approximation is the default;
  nested Laplace handles areal fields, NUTS gives the full posterior,
  and a Polya-Gamma Gibbs sampler is available for single-season and
  community models.
- **Simulation with known truth.** Every family has a simulator that
  returns the data with the parameters that generated them, so a model
  can be checked against known values before it meets field data.

The package’s recovery tests are built on these simulators. [Occupancy:
tulpaObs vs INLA
aggregated-binomial](https://gillescolling.com/tulpaObs/articles/occupancy-vs-inla.html)
and [Cover hurdle: tulpaObs vs INLA on simulated cover
data](https://gillescolling.com/tulpaObs/articles/cover-hurdle-vs-inla.html)
compare fits with INLA on simulated data with known truth.

## Installation

tulpaObs compiles C++17 code. Install the build tools for your operating
system first:

| Operating system | Build tools |
|----|----|
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version |
| macOS | Run `xcode-select --install` in Terminal |
| Debian / Ubuntu | Install `r-base-dev` |

Check that R can find the compiler:

``` r
pkgbuild::has_build_tools(debug = TRUE)
```

Then install tulpaObs:

``` r
# install.packages("pak")
pak::pak("gcol33/tulpaObs")
```

pak installs tulpa and
[tulpaMesh](https://CRAN.R-project.org/package=tulpaMesh) from CRAN
along with it. The build links against the BLAS and LAPACK libraries
that R itself uses, so no other system library is needed. Tagged
releases are listed on the [releases
page](https://github.com/gcol33/tulpaObs/releases).

### Optional: multithreading

OpenMP lets the community NUTS samplers (`ms_occu()`, `ms_count()`,
`jsdm()`, `ms_abun()`, `ms_dyn_occu()`) evaluate their gradient on
several CPU threads, set through `n.threads.grad` in the `control` list
of `tobs()`. The per-species reduction runs in a fixed order, so the
gradient is the same at any thread count.

Rtools on Windows and the default compiler on Linux include OpenMP.
Without it, the package installs and runs on one thread.

#### macOS setup

Apple’s compiler needs an additional OpenMP runtime. Install it in
Terminal:

``` sh
brew install libomp
```

Before installing tulpaObs, add the settings for your Mac to
`~/.R/Makevars`, R’s personal compiler-settings file. Create the
directory and file if needed.

**Apple Silicon:**

``` makefile
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
```

**Intel:**

``` makefile
CPPFLAGS += -I/usr/local/opt/libomp/include -Xclang -fopenmp
LDFLAGS += -L/usr/local/opt/libomp/lib -lomp
```

## Usage

1.  Fit a model with `tobs()`: the occupancy formula, the data, a family
    such as `occu()`, the detection formula and the detection histories.
2.  Read the estimates with `coef()` or `summary()`.
3.  Predict with `predict()` and check the fit with `check_model()`.

`peatland_occu` is a synthetic survey shipped with the package: 120
sites visited four times, with elevation and wetness per site and survey
effort per visit.

``` r
library(tulpaObs)

occ <- peatland_occu
fit <- tobs(
  ~ elevation + wetness,
  data = occ$occ.covs,
  family = occu(),
  detection = ~ effort,
  y = occ$y,
  visits = occ$det.covs
)
coef(fit)
#> psi_(Intercept)   psi_elevation     psi_wetness   p_(Intercept)  p_visit_effort 
#>     0.242198812    -0.505199672     0.691204862    -0.001590117     0.853642172
```

The `psi_` coefficients belong to the occupancy formula and the `p_`
coefficients to the detection formula, both on the logit scale.

The detection histories `y` are a matrix with one row per site and one
column per visit (1 detected, 0 not detected, `NA` missing). Site
covariates are a data frame with one row per site, and covariates that
change between visits are matrices of the same shape as `y`, passed
through `visits`. `tobs_format()` and `tobs_data()` build the same input
from matrices or from a long table with one row per site-visit; [Shaping
data for
tobs()](https://gillescolling.com/tulpaObs/articles/data-formatting.html)
covers both.

## Structured effects

Spatial fields, temporal processes and random effects are terms in the
formula, and each enters the process whose formula it is written in:

``` r
~ elevation +
  spatial(graph = adj) +
  temporal(year, type = "ar1") +
  (1 | observer)
```

`spatial()` takes an adjacency matrix for areal models (ICAR by default,
BYM2, CAR) or coordinates for continuous fields (`model = "spde"` or
`"gp"`). Random effects follow lme4’s bar syntax, including correlated
slopes, crossed and nested grouping. `?tobs_terms` documents every term.
[Areal spatial occupancy with ICAR, BYM2, and
CAR](https://gillescolling.com/tulpaObs/articles/spatial-occupancy.html),
[Occupancy with a continuous SPDE spatial
field](https://gillescolling.com/tulpaObs/articles/occupancy-spatial-spde.html)
and [Random effects with lme4 bar
syntax](https://gillescolling.com/tulpaObs/articles/random-effects.html)
work through each.

## Inference methods

`method` in `tobs()` selects the route: `"laplace"` (the default),
`"nested_laplace"` for areal fields and the joint cover hurdle, `"nuts"`
for the full posterior with split-Rhat and effective sample sizes, and
`"pg_gibbs"`, an exact Polya-Gamma Gibbs sampler, on single-season
occupancy, `t_occu()` and the community families. An unsupported
combination of family and method stops with the methods that family
accepts.

Variance components of binary random effects are biased under the
Laplace approximation when groups are small. An adaptive Gauss-Hermite
pass corrects them by default; with eight observations per group it
lowers the measured bias in the random-effect standard deviation from
about 18 % to about 4 %, matching a NUTS fit ([Random effects with lme4
bar
syntax](https://gillescolling.com/tulpaObs/articles/random-effects.html)).

## Checking and comparing fits

`check_model()` collects the sampler diagnostics, WAIC, a posterior
predictive check, a zero-inflation test and Moran’s I in one report.
`waic()` and `loo()` are the generics of the
[loo](https://mc-stan.org/loo/) package, so `loo::loo_compare()`
compares tulpaObs fits directly. [Diagnostics: checking, comparing, and
residuals](https://gillescolling.com/tulpaObs/articles/diagnostics.html)
covers the full workflow.

## Documentation

Start with [Getting started with
tulpaObs](https://gillescolling.com/tulpaObs/articles/quickstart.html):
one path from simulated detection histories to a fitted, compared and
predicted occupancy model.

- [Shaping data for
  tobs()](https://gillescolling.com/tulpaObs/articles/data-formatting.html):
  detection histories and covariates from matrices or long tables.
- [Single-season
  occupancy](https://gillescolling.com/tulpaObs/articles/occupancy.html):
  the base model end to end, including when occupancy and detection
  cannot be separated.
- [Dynamic (multi-season)
  occupancy](https://gillescolling.com/tulpaObs/articles/dynamic-occupancy.html):
  colonisation and extinction across seasons.
- [Integrated occupancy: combining data
  sources](https://gillescolling.com/tulpaObs/articles/integrated-occupancy.html):
  several detection processes sharing one occupancy state.
- [N-mixture abundance with
  tulpaObs](https://gillescolling.com/tulpaObs/articles/abundance.html):
  abundance from replicated counts.
- [Community models in
  tulpaObs](https://gillescolling.com/tulpaObs/articles/community-models.html):
  many species with shared community priors.
- [Cover hurdle (independent
  arms)](https://gillescolling.com/tulpaObs/articles/cover-hurdle.html)
  and [Cover hurdle on resurvey panels (MOTIVATE
  reduction)](https://gillescolling.com/tulpaObs/articles/cover-hurdle-motivate.html):
  vegetation cover with zeros.
- [Areal spatial occupancy with ICAR, BYM2, and
  CAR](https://gillescolling.com/tulpaObs/articles/spatial-occupancy.html)
  and [Occupancy with a continuous SPDE spatial
  field](https://gillescolling.com/tulpaObs/articles/occupancy-spatial-spde.html):
  spatial fields.
- [Random effects with lme4 bar
  syntax](https://gillescolling.com/tulpaObs/articles/random-effects.html):
  grouped effects and the small-group correction.
- [Diagnostics: checking, comparing, and
  residuals](https://gillescolling.com/tulpaObs/articles/diagnostics.html):
  checking and comparing fits.
- [Occupancy: tulpaObs vs INLA
  aggregated-binomial](https://gillescolling.com/tulpaObs/articles/occupancy-vs-inla.html)
  and [Cover hurdle: tulpaObs vs INLA on simulated cover
  data](https://gillescolling.com/tulpaObs/articles/cover-hurdle-vs-inla.html):
  comparisons with INLA.
- [Full
  reference](https://gillescolling.com/tulpaObs/reference/index.html):
  every function, grouped by job.

## Getting help

Questions and bug reports go to the [GitHub
issues](https://github.com/gcol33/tulpaObs/issues). A bug report is
easiest to act on with a small reproducible example, for instance made
with [reprex](https://reprex.tidyverse.org/).

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## Acknowledgements

Inference, the latent spatial and temporal structures and the generic
diagnostics come from [tulpa](https://CRAN.R-project.org/package=tulpa);
spatial meshes come from
[tulpaMesh](https://CRAN.R-project.org/package=tulpaMesh).

## License

tulpaObs is licensed under MIT.

## Citation

To cite the tulpaObs package in publications, use:

``` text
Colling G (2026). _tulpaObs: Bayesian Occupancy, Abundance and
Detection Models via 'tulpa'_. R package version 0.6.2,
<https://gillescolling.com/tulpaObs/>.
```

BibTeX:

``` bibtex
@Manual{colling2026tulpaobs,
  title = {{tulpaObs}: Bayesian Occupancy, Abundance and Detection Models via 'tulpa'},
  author = {Gilles Colling},
  year = {2026},
  note = {R package version 0.6.2},
  url = {https://gillescolling.com/tulpaObs/},
}
```

Download the citation as
[BibTeX](https://gillescolling.com/tulpaObs/tulpaObs.bib) or
[RIS](https://gillescolling.com/tulpaObs/tulpaObs.ris).
