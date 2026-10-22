# Simulate joint species distribution (presence/absence) data

Presence `y_{i,s} ~ Bernoulli(psi_{i,s})` with `logit psi = X beta_s`,
the per-species occupancy coefficients drawn from Gaussian community
hyperpriors. Matches the
[`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md) family,
which observes presence directly (no detection process), so the response
is an `N x n_species` presence matrix.

## Usage

``` r
simulate_jsdm(
  N = 100,
  n.species = 10,
  beta.comm.mean = c(0, 0.5),
  beta.comm.sd = c(0.5, 0.3),
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- n.species:

  Number of species (default 10).

- beta.comm.mean:

  Community-mean occupancy coefficients, length `1 + n_occ_covs`
  (intercept first).

- beta.comm.sd:

  Between-species SD of each occupancy coefficient (same length as
  `beta.comm.mean`).

- seed:

  Random seed.

## Value

A list with `y` (an `N x n_species` 0/1 presence matrix), `data`, and
`truth` (per-species coefficients and the community hyperparameters).

## Examples

``` r
sim <- simulate_jsdm(N = 60, n.species = 5, seed = 1)
dim(sim$y)   # 60 sites x 5 species presence matrix
```
