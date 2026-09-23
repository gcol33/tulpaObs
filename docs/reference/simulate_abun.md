# Simulate Royle (2004) N-mixture abundance data

Latent abundance `N_i ~ Poisson(lambda_i)` (or, under
`mixture = "negbin"`, `N_i ~ NegBin(mean = lambda_i, size = size)` with
variance `lambda + lambda^2 / size`) with
`log lambda_i = X_lambda beta_lambda`, and replicate counts
`y_ij ~ Binomial(N_i, p_i)` with `logit p_i = X_p beta_p` (site-level
detection). The returned `y` is an `N x J` integer count matrix suitable
for [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
with [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md).

## Usage

``` r
simulate_abun(
  N = 100,
  J = 4,
  n_abund_covs = 2,
  n_det_covs = 1,
  beta_lambda = NULL,
  beta_p = NULL,
  mixture = c("poisson", "negbin", "zip", "zinb"),
  size = 2,
  omega = 0.3,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Number of replicate visits (default 4).

- n_abund_covs:

  Number of abundance covariates (default 2).

- n_det_covs:

  Number of detection covariates (default 1).

- beta_lambda:

  Abundance coefficients `c(intercept, slopes...)` on the log scale.
  Default `c(log(3), runif(n_abund_covs, -0.5, 0.5))`.

- beta_p:

  Detection coefficients `c(intercept, slopes...)` on the logit scale.
  Default `c(0, runif(n_det_covs, -0.5, 0.5))` (intercept 0 = p 0.5).

- mixture:

  Abundance mixing distribution: `"poisson"` (default), `"negbin"`
  (negative binomial, overdispersed), or their zero-inflated
  counterparts `"zip"` / `"zinb"` (a structural-zero share `omega` of
  sites have `N = 0` regardless of `lambda`).

- size:

  Negative-binomial size `r` (the `"negbin"` / `"zinb"` dispersion;
  variance `lambda + lambda^2 / r`). Smaller `r` means more
  overdispersion; as `r` grows large the model approaches Poisson.
  Default 2. Ignored under Poisson.

- omega:

  Structural-zero probability for `"zip"` / `"zinb"` (the share of sites
  with `N = 0` independent of `lambda`). Default 0.3. Ignored otherwise.

- seed:

  Optional random seed.

## Value

A list with `y` (N x J count matrix), `data` (covariate data frame), and
`truth` (the coefficients, per-site `lambda`, `p`, latent `N`, and the
`mixture` / `size` used).

## Examples

``` r
sim <- simulate_abun(N = 50, J = 3, seed = 1)
dim(sim$y)
```
