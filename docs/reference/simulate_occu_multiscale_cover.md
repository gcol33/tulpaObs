# Simulate a three-level occupancy + cover hurdle data set

Draws cell occupancy `z_c ~ Bernoulli(psi_c)`, plot availability
`a_cj | z_c=1 ~ Bernoulli(theta_cj)`, visit detection
`y_cjv | a_cj=1 ~ Bernoulli(p_cjv)`, and the cover hurdle
`cover_cjv | y_cjv=1 ~ f_pos`, with a shared areal ICAR field on the
occupancy (`sigma`) and cover (`alpha * sigma`) arms. The
data-generating process for
[`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md)
recovery tests.

## Usage

``` r
simulate_occu_multiscale_cover(
  n.cells = 60L,
  plots.per.cell = 4L,
  visits.per.plot = 2L,
  beta.psi = c(0.4, 0.6),
  beta.theta = c(0.2, 0.5),
  beta.p = c(0, 0.5),
  beta.pos = c(log(0.1), -0.4),
  positive = c("lognormal", "beta", "gaussian"),
  phi = 0.35,
  adj = NULL,
  sigma = 0.7,
  alpha = 1,
  trend = FALSE,
  sigma.trend = 0.7,
  alpha.trend = 1,
  seed = NULL
)
```

## Arguments

- n.cells:

  number of areal cells (graph nodes).

- plots.per.cell:

  plots (the availability units) per cell.

- visits.per.plot:

  detection replicate visits per plot.

- beta.psi, beta.theta, beta.p, beta.pos:

  length-2 `c(intercept, slope)` coefficients for the four arms (one
  covariate each).

- positive:

  `"lognormal"` or `"beta"` cover arm.

- phi:

  cover dispersion (lognormal log-scale SD, or beta precision).

- adj:

  `n.cells x n.cells` adjacency; default a 1-D chain.

- sigma:

  areal-field marginal SD on the occupancy arm.

- alpha:

  cover-arm field scaling (`alpha * sigma` on the cover arm).

- trend:

  add a second, spatially-varying-coefficient (SVC) areal field weighted
  per cell by a cell-level covariate `tcov` (added to `data`):
  `tcov_c * sigma.trend * f_trend[c]` on occupancy and
  `tcov_c * alpha.trend * sigma.trend * f_trend[c]` on cover. Default
  `FALSE`.

- sigma.trend, alpha.trend:

  trend-field marginal SD and its cover-arm scaling (used only when
  `trend = TRUE`).

- seed:

  optional RNG seed.

## Value

A list with `y`, `y_pos` (`[n_plots x visits_per_plot]`), the plot-level
`data` (cell id `cell`, covariates), `adj`, and `truth`.

## Examples

``` r
sim <- simulate_occu_multiscale_cover(n.cells = 10, plots.per.cell = 2,
                                      visits.per.plot = 2, seed = 1)
dim(sim$y)
```
