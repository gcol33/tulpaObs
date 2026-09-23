# Simulate integrated (multi-source) occupancy data

Simulate integrated (multi-source) occupancy data

## Usage

``` r
simulate_int_occu(
  N_total = 150,
  n_data = 2,
  J = c(4, 3),
  n_shared = 20,
  beta_occ = c(0.5, 0.3),
  beta_det = list(c(0.2, -0.4), c(-0.1, 0.3)),
  seed = NULL
)
```

## Arguments

- N_total:

  Total number of unique sites (default 150).

- n_data:

  Number of data sources (default 2).

- J:

  Vector of visits per source (default c(4, 3)).

- n_shared:

  Number of sites shared across sources (default 20).

- beta_occ:

  Occupancy coefficients (default c(0.5, 0.3)).

- beta_det:

  List of detection coefficient vectors per source.

- seed:

  Random seed.

## Value

A list with `y` (list of matrices, each row named by the site of `data`
it measures), `data`, `site_maps`, and `truth`.

## Examples

``` r
sim <- simulate_int_occu(N_total = 60, n_shared = 10, seed = 1)
lapply(sim$y, dim)
```
