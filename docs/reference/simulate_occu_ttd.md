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
  n.psi.covs = 1,
  n.rate.covs = 1,
  beta.psi = NULL,
  beta.rate = NULL,
  Tmax = 3,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- J:

  Number of replicate surveys per site (default 4).

- n.psi.covs, n.rate.covs:

  Number of occupancy / rate covariates.

- beta.psi:

  Occupancy coefficients `c(intercept, slopes...)` (logit). Default
  `c(qlogis(0.6), runif(n.psi.covs, -0.5, 0.5))`.

- beta.rate:

  Log detection-rate coefficients `c(intercept, slopes...)`. Default
  `c(log(0.6), runif(n.rate.covs, -0.4, 0.4))`.

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
