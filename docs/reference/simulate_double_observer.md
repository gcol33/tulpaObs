# Simulate a double-observer abundance data set

Draws from the
[`double_observer()`](https://gillescolling.com/tulpaObs/reference/double_observer.md)
model: site abundance `N ~ Poisson(lambda)` observed by two observers
with detection `p1` / `p2`. For `type = "independent"` each individual
is recorded by observer 1 only, observer 2 only, or both (three cell
counts). For `type = "dependent"` a primary observer records what it
detects and a secondary observer records only the primary's misses (two
cell counts, primary-detected and secondary-only); the primary observer
alternates across sites (role-swapping) so `p1` and `p2` are
identifiable, and the per-site primary indicator is returned as
`primary`.

## Usage

``` r
simulate_double_observer(
  N = 200,
  type = c("independent", "dependent"),
  n_abund_covs = 1,
  n_det_covs = 1,
  beta_lambda = NULL,
  beta_p1 = NULL,
  beta_p2 = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- type:

  `"independent"` (default) or `"dependent"` (role-swapping).

- n_abund_covs, n_det_covs:

  Number of abundance / detection covariates.

- beta_lambda:

  Log-abundance coefficients `c(intercept, slopes...)`. Default
  `c(log(8), runif(n_abund_covs, -0.5, 0.5))`.

- beta_p1, beta_p2:

  Per-observer detection coefficients (logit). Default moderate
  detection.

- seed:

  Optional random seed.

## Value

A list with `y` (`N x 3` cell counts for `"independent"`, `N x 2` for
`"dependent"`), `data`, `primary` (the per-site primary observer, for
`"dependent"`), and `truth`.

## Examples

``` r
sim <- simulate_double_observer(N = 50, seed = 1)
head(sim$y)
```
