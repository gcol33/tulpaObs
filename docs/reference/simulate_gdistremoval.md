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
  n_periods = 4L,
  transect = "line",
  n_abund_covs = 1,
  n_det_covs = 1,
  n_rem_covs = 1,
  beta_lambda = NULL,
  beta_sigma = NULL,
  beta_r = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- cutpoints:

  Distance-bin edges `0 = c_0 < ... < c_B`.

- n_periods:

  Number of removal periods.

- transect:

  `"line"` (default) or `"point"`.

- n_abund_covs, n_det_covs, n_rem_covs:

  Number of abundance / distance / removal covariates.

- beta_lambda, beta_sigma, beta_r:

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
