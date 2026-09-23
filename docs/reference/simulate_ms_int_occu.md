# Simulate integrated multi-species occupancy data

Simulate integrated multi-species occupancy data

## Usage

``` r
simulate_ms_int_occu(
  N = 100,
  J = c(3, 4),
  n_species = 5,
  n_data = 2,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Vector of visits per source (default c(3, 4)).

- n_species:

  Number of species (default 5).

- n_data:

  Number of data sources (default 2).

- seed:

  Random seed.

## Value

A list with `y` (list of 3D arrays), `data`, and `truth`.

## See also

[`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md),
the family this simulates for.

## Examples

``` r
sim <- simulate_ms_int_occu(N = 40, n_species = 3, seed = 1)
lapply(sim$y, dim)
```
