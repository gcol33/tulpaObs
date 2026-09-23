# Simulate a multi-species co-occurrence occupancy data set

Draws from the
[`occu_multi()`](https://gillescolling.com/tulpaObs/reference/occu_multi.md)
model: a joint occupancy state for `S` species from the log-linear
(first + second order) occupancy model, then per-species site-level
detections at the observed visit pattern.

## Usage

``` r
simulate_occu_multi(
  S = 2,
  N = 300,
  J = 4,
  n_state_covs = 1,
  beta_first = NULL,
  beta_second = NULL,
  beta_p = NULL,
  seed = NULL
)
```

## Arguments

- S:

  Number of species (default 2).

- N:

  Number of sites (default 300).

- J:

  Number of replicate visits (default 4).

- n_state_covs:

  Number of shared occupancy covariates (default 1).

- beta_first:

  Length-`S` list of first-order natural-parameter coefficients
  `c(intercept, slopes...)`. Default moderate occupancy.

- beta_second:

  Length-`choose(S, 2)` list of second-order (interaction) coefficients
  (in `combn(S, 2)` order). Default a single positive interaction for
  `S = 2`, else 0.

- beta_p:

  Length-`S` list of detection coefficients (logit). Default
  `c(qlogis(0.5), ...)`.

- seed:

  Optional random seed.

## Value

A list with `y` (a length-`S` list of `N x J` matrices), `data`,
`species`, and `truth`.

## Examples

``` r
sim <- simulate_occu_multi(S = 2, N = 50, J = 3, seed = 1)
length(sim$y)
```
