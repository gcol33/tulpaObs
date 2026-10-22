# Simulate a joint distance + removal data set

Draws from the
[`gdistremoval()`](https://gillescolling.com/tulpaObs/reference/gdistremoval.md)
model: site abundance `N ~ Poisson(lambda)`, the detected birds
cross-classified by a distance band (half-normal key, scale `sigma`) and
by a removal period (per-period capture `r`).

## Usage

``` r
simulate_gdistremoval(
  N = 200,
  cutpoints = c(0, 10, 20, 30, 40),
  n.periods = 4L,
  transect = "line",
  n.abund.covs = 1,
  n.det.covs = 1,
  n.rem.covs = 1,
  beta.lambda = NULL,
  beta.sigma = NULL,
  beta.r = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- cutpoints:

  Distance-bin edges `0 = c_0 < ... < c_B`.

- n.periods:

  Number of removal periods.

- transect:

  `"line"` (default) or `"point"`.

- n.abund.covs, n.det.covs, n.rem.covs:

  Number of abundance / distance / removal covariates.

- beta.lambda, beta.sigma, beta.r:

  Coefficients on the log-abundance, log-scale, and logit-removal arms
  (`c(intercept, slopes...)`). Defaults give moderate abundance /
  detection.

- seed:

  Optional random seed.

## Value

A list with `y` (distance-band counts), `y_rem` (removal-period counts),
`data`, and `truth`.

## Examples

``` r
sim <- simulate_gdistremoval(N = 50, seed = 1)
dim(sim$y)
dim(sim$y_rem)
```
