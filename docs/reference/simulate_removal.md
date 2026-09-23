# Simulate removal-sampling abundance data

Latent abundance `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, size)`)
with `log lambda_i = X_lambda beta_lambda`, observed through `K` ordered
removal passes: pass `k` removes
`Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the individuals still
present, with `logit p_{ik} = X_p beta_p` (site-level detection here).
The returned `y` is an `N x K` integer matrix of per-pass removals
suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md).

## Usage

``` r
simulate_removal(
  N = 100,
  K = 4,
  n_abund_covs = 2,
  n_det_covs = 1,
  beta_lambda = NULL,
  beta_p = NULL,
  mixture = c("poisson", "negbin"),
  size = 3,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- K:

  Number of removal passes (default 4).

- n_abund_covs:

  Number of abundance covariates (default 2).

- n_det_covs:

  Number of detection covariates (default 1).

- beta_lambda:

  Abundance coefficients on the log scale. Default
  `c(log(6), runif(n_abund_covs, -0.5, 0.5))`.

- beta_p:

  Detection coefficients on the logit scale. Default
  `c(0.4, runif(n_det_covs, -0.5, 0.5))`.

- mixture:

  `"poisson"` (default) or `"negbin"`.

- size:

  Negative-binomial size `r` (`mixture = "negbin"` only, default 3).

- seed:

  Optional random seed.

## Value

A list with `y` (N x K removal matrix), `data` (covariate data frame),
and `truth` (coefficients, per-site `lambda`, `p`, latent `N`,
mixture/size).

## Examples

``` r
sim <- simulate_removal(N = 50, K = 3, seed = 1)
dim(sim$y)
```
