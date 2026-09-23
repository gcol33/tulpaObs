# Simulate binned distance-sampling abundance data

Latent abundance `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, size)`)
with `log lambda_i = X_lambda beta_lambda`, observed through a
half-normal or hazard-rate detection function with
`log sigma_i = X_sigma beta_sigma`. Each individual lands in a distance
bin (or goes undetected) by the multinomial cell probabilities. Returns
an `N x n_bins` integer matrix of per-bin counts suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md).

## Usage

``` r
simulate_distance(
  N = 200,
  cutpoints = seq(0, 1, length.out = 6),
  key = c("halfnorm", "hazard"),
  transect = c("line", "point"),
  n_abund_covs = 1,
  n_sigma_covs = 1,
  beta_lambda = NULL,
  beta_sigma = NULL,
  shape = 3,
  mixture = c("poisson", "negbin"),
  size = 5,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- cutpoints:

  Distance-bin edges (length `n_bins + 1`). Default
  `seq(0, 1, length.out = 6)` (five bins out to 1).

- key:

  `"halfnorm"` (default) or `"hazard"`.

- transect:

  `"line"` (default) or `"point"`.

- n_abund_covs, n_sigma_covs:

  Number of abundance / detection covariates.

- beta_lambda:

  Abundance coefficients (log scale). Default
  `c(log(40), runif(n_abund_covs, -0.4, 0.4))`.

- beta_sigma:

  Detection-scale coefficients (log scale). Default
  `c(log(0.4), runif(n_sigma_covs, -0.3, 0.3))`.

- shape:

  Hazard-rate shape `b` (`key = "hazard"` only, default 3).

- mixture:

  `"poisson"` (default) or `"negbin"`.

- size:

  Negative-binomial size `r` (`mixture = "negbin"` only, default 5).

- seed:

  Optional random seed.

## Value

A list with `y` (N x n_bins count matrix), `data` (covariates),
`cutpoints`, and `truth` (coefficients, per-site `lambda` / `sigma`,
latent `N`, key/transect/mixture/shape/size).

## Examples

``` r
sim <- simulate_distance(N = 50, seed = 1)
head(sim$y)
```
