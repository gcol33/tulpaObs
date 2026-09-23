# Simulate single-species occupancy data

Simulate single-species occupancy data

## Usage

``` r
simulate_occu(
  N = 100,
  J = 4,
  n_occ_covs = 2,
  n_det_covs = 1,
  beta_occ = NULL,
  beta_det = NULL,
  n_visit_groups = 0L,
  sigma_visit = 1,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Number of visits (default 4).

- n_occ_covs:

  Number of occupancy covariates (default 2).

- n_det_covs:

  Number of detection covariates (default 1).

- beta_occ:

  Occupancy coefficients (auto-generated if NULL).

- beta_det:

  Detection coefficients (auto-generated if NULL).

- n_visit_groups:

  Number of groups a visit can belong to (an observer, say), each adding
  an effect to the detection logit of the visits it made. Visits are
  assigned to groups uniformly at random. `0` (default) simulates no
  such effect.

- sigma_visit:

  Standard deviation of the per-group detection effects, drawn from
  `N(0, sigma_visit^2)`.

- seed:

  Random seed.

## Value

A list with `y`, `data`, and `truth`. With `n_visit_groups > 0` it also
carries `visits`, a data frame with one row per site-visit in site-major
order and the factor column `visit_group`, and `truth` gains `b_visit`
(the group effects) and `sigma_visit`.

## Examples

``` r
sim <- simulate_occu(N = 50, J = 3, seed = 1)
dim(sim$y)

obs <- simulate_occu(N = 50, J = 3, n_visit_groups = 5, seed = 1)
table(obs$visits$visit_group)
```
