# Simulate cover-hurdle data (lognormal positive part)

Generates synthetic data matching the `cover(response = "lognormal")`
generative model:

## Usage

``` r
simulate_cover(
  N = 200L,
  beta_occ = c(-0.5, 0.8),
  beta_pos = c(-1, 0.3),
  sigma_pos = 0.4,
  response = c("lognormal", "gaussian"),
  spatial_range = NULL,
  spatial_var = 1,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- beta_occ:

  Length-2 occurrence coefficients (intercept, slope on x).

- beta_pos:

  Length-2 log-cover coefficients (intercept, slope on x).

- sigma_pos:

  Lognormal residual standard deviation (default 0.4).

- response:

  Positive-arm likelihood: `"lognormal"` (default) or `"gaussian"`.

- spatial_range:

  Length scale of the exponential spatial field (in unit-square
  distance). `NULL` disables the spatial layer.

- spatial_var:

  Marginal variance of the spatial field (default 1).

- seed:

  Optional integer seed.

## Value

A list with:

- data:

  A data frame with `cover`, covariate `x`, and `lon`, `lat`.

- y:

  `data$cover` (length-N numeric vector, for passing to
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)).

- coords:

  An N x 2 matrix of coordinates.

- truth:

  A list with `beta_occ`, `beta_pos`, `sigma_pos`, `p`, `mu`, `occur`,
  `spatial_occ`, `spatial_pos`.

## Details

- eta_occ:

  `= X %*% beta_occ + spatial_occ`

- eta_pos:

  `= X %*% beta_pos + spatial_pos`

- p:

  = plogis(eta_occ)

- occur:

  ~ Bernoulli(p)

- log_cover:

  ~ Normal(eta_pos, sigma_pos^2) \| occur = 1

- cover:

  = exp(log_cover) when occur = 1, else 0

Both layers share the same design matrix by default (one continuous
covariate plus an intercept). When `spatial_range` is supplied, a simple
exponential-kernel Gaussian field on a random unit-square coordinate set
is added to both linear predictors; the two layers see independently
simulated draws of that field, which matches the two-independent-arm
fit. The joint engine's shared field is a different generative model.

## Examples

``` r
sim <- simulate_cover(N = 200, seed = 1)
head(sim$data)
```
