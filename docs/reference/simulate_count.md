# Simulate count / relative-abundance GLMM data

One observed value per site from a GLMM with no detection process,
matching the
[`count()`](https://gillescolling.com/tulpaObs/reference/count.md)
family: `log`-link Poisson / negative-binomial counts, an identity-link
Gaussian response, or a logit-link binomial response (`k` successes out
of `n` trials per site). The mean predictor is `X beta` with `beta` the
fixed-effect coefficients (intercept first).

## Usage

``` r
simulate_count(
  N = 200,
  beta = c(1, 0.5),
  response = c("poisson", "negbin", "gaussian", "binomial"),
  size = 2,
  sd = 1,
  trials = 10,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- beta:

  Fixed-effect coefficients, length `1 + n_covs` (intercept first).

- response:

  The response distribution: `"poisson"`, `"negbin"`, `"gaussian"`, or
  `"binomial"`.

- size:

  Negative-binomial size (dispersion); larger is closer to Poisson. Used
  only for `response = "negbin"`.

- sd:

  Gaussian residual SD. Used only for `response = "gaussian"`.

- trials:

  Per-site trial count for `response = "binomial"`: a scalar (recycled)
  or a length-`N` vector. Default 10 (a scalar 1 gives Bernoulli data,
  the `svcPGBinom` `trials = 1` setting).

- seed:

  Random seed.

## Value

A list with `y` (a length-`N` numeric response vector; success counts
for the binomial response), `data` (the site covariates), and `truth`
(the coefficients, response, dispersion, and the binomial `trials`).

## Examples

``` r
sim <- simulate_count(N = 100, beta = c(1, 0.5), seed = 1)
length(sim$y)
```
