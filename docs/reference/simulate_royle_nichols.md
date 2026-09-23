# Simulate Royle-Nichols occupancy data

Latent abundance `N_i ~ Poisson(lambda_i)`,
`log lambda = X beta_lambda`, and per-visit detection
`y_ij ~ Bernoulli(1 - (1 - r_ij)^{N_i})` with per-individual detection
`logit r_ij = beta_r`. Detection is site-level by default; supplying
`beta_r_visit` draws a per-visit covariate `w` and makes
`logit r_ij = beta_r + beta_r_visit * w_ij` visit-varying. Matches the
[`royle_nichols()`](https://gillescolling.com/tulpaObs/reference/royle_nichols.md)
family; the response is an `N x J` 0/1 detection-history matrix.

## Usage

``` r
simulate_royle_nichols(
  N = 200,
  J = 5,
  beta_lambda = c(0.3, 0.5),
  beta_r = -0.8,
  beta_r_visit = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- J:

  Number of visits per site (default 5).

- beta_lambda:

  Log-abundance coefficients `c(intercept, slope_on_x)`.

- beta_r:

  Per-individual detection logit (a scalar intercept).

- beta_r_visit:

  Optional scalar slope on a drawn per-visit detection covariate `w`.
  When supplied, detection varies by visit and the returned list carries
  `visits = list(w = )` (an `N x J` matrix) for `tobs(..., visits =)`.

- seed:

  Random seed.

## Value

A list with `y` (`N x J` 0/1 matrix), `data`, `truth` (`beta_lambda`,
`beta_r`, `beta_r_visit`, realised abundance `N`, `lambda`, per-site or
per-visit `r`), and, when `beta_r_visit` is supplied, `visits`.

## Examples

``` r
sim <- simulate_royle_nichols(N = 100, J = 4, seed = 1)
dim(sim$y)
# Visit-varying detection:
sv <- simulate_royle_nichols(N = 100, J = 4, beta_r_visit = 0.8, seed = 1)
names(sv$visits)
```
