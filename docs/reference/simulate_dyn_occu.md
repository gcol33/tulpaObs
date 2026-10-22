# Simulate temporal (multi-season) occupancy data

Colonization and extinction are constant across the `n_seasons - 1`
transition intervals by default. Supplying `beta.gamma` and/or
`beta.epsilon` makes the corresponding rate SEASON-VARYING: a
per-`(site, interval)` covariate is drawn into an
`[N x (n_seasons - 1)]` matrix column (`gamma_cov` / `eps_cov`) of the
returned `data`, and the rate is `plogis(beta[1] + beta[2] * cov)`. Fit
these with `colonization = ~ gamma_cov` / `extinction = ~ eps_cov` (the
matrix column drives the interval-indexed rate).

## Usage

``` r
simulate_dyn_occu(
  N = 100,
  J = 4,
  n.seasons = 5,
  beta.occ = c(0.5),
  beta.det = c(0),
  gamma = 0.2,
  epsilon = 0.1,
  beta.gamma = NULL,
  beta.epsilon = NULL,
  beta.det.season = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Number of visits per season (default 4).

- n.seasons:

  Number of seasons (default 5).

- beta.occ:

  Initial occupancy coefficients.

- beta.det:

  Detection coefficients.

- gamma:

  Colonization probability (default 0.2); ignored when `beta.gamma` is
  given.

- epsilon:

  Extinction probability (default 0.1); ignored when `beta.epsilon` is
  given.

- beta.gamma:

  Optional `c(intercept, slope)` for a season-varying colonization logit
  driven by a drawn per-`(site, interval)` covariate.

- beta.epsilon:

  Optional `c(intercept, slope)` for a season-varying extinction logit
  driven by a drawn per-`(site, interval)` covariate.

- beta.det.season:

  Optional `c(intercept, slope)` for a season-varying detection logit
  driven by a drawn per-`(site, season)` covariate `det_cov` (a
  `[N x T]` matrix column of `data`); fit with `detection = ~ det_cov`.

- seed:

  Random seed.

## Value

A list with `y` (3D array), `data`, and `truth`.

## Examples

``` r
sim <- simulate_dyn_occu(N = 40, J = 3, n.seasons = 4, seed = 1)
dim(sim$y)
```
