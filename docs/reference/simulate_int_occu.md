# Simulate integrated (multi-source) occupancy data

Simulate integrated (multi-source) occupancy data

## Usage

``` r
simulate_int_occu(
  N.total = 150,
  n.data = 2,
  J = c(4, 3),
  n.shared = 20,
  beta.occ = c(0.5, 0.3),
  beta.det = list(c(0.2, -0.4), c(-0.1, 0.3)),
  seed = NULL
)
```

## Arguments

- N.total:

  Total number of unique sites (default 150).

- n.data:

  Number of data sources (default 2).

- J:

  Vector of visits per source (default c(4, 3)).

- n.shared:

  Number of sites shared across sources (default 20).

- beta.occ:

  Occupancy coefficients (default c(0.5, 0.3)).

- beta.det:

  List of detection coefficient vectors per source.

- seed:

  Random seed.

## Value

A list with `y` (list of matrices, each row named by the site of `data`
it measures), `data`, `site_maps`, and `truth`.

## Examples

``` r
sim <- simulate_int_occu(N.total = 60, n.shared = 10, seed = 1)
lapply(sim$y, dim)
```
