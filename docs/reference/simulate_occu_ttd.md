# Simulate a time-to-detection occupancy data set

Draws from the
[`occu_ttd()`](https://gillescolling.com/tulpaObs/reference/occu_ttd.md)
model: site occupancy `z ~ Bernoulli(psi)` and, at occupied sites,
exponential time-to-detection at rate `lambda` over a survey of length
`Tmax` (a non-detection is censored at `Tmax`).

## Usage

``` r
simulate_occu_ttd(
  N = 200,
  J = 4,
  n_psi_covs = 1,
  n_rate_covs = 1,
  beta_psi = NULL,
  beta_rate = NULL,
  Tmax = 3,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- J:

  Number of replicate surveys per site (default 4).

- n_psi_covs, n_rate_covs:

  Number of occupancy / rate covariates.

- beta_psi:

  Occupancy coefficients `c(intercept, slopes...)` (logit). Default
  `c(qlogis(0.6), runif(n_psi_covs, -0.5, 0.5))`.

- beta_rate:

  Log detection-rate coefficients `c(intercept, slopes...)`. Default
  `c(log(0.6), runif(n_rate_covs, -0.4, 0.4))`.

- Tmax:

  Survey length (scalar). Default 3.

- seed:

  Optional random seed.

## Value

A list with `y` (N x J TTD matrix), `data`, `Tmax`, and `truth`.

## Examples

``` r
sim <- simulate_occu_ttd(N = 50, J = 3, seed = 1)
dim(sim$y)
```
