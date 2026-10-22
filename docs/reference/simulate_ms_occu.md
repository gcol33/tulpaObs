# Simulate multi-species occupancy data

Simulate multi-species occupancy data

## Usage

``` r
simulate_ms_occu(
  N = 100,
  J = 4,
  n.species = 10,
  beta.comm.mean = c(0, 0.5),
  beta.comm.sd = c(0.5, 0.3),
  alpha.comm.mean = c(0),
  alpha.comm.sd = c(0.5),
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Number of visits (default 4).

- n.species:

  Number of species (default 10).

- beta.comm.mean:

  Community mean for occupancy (default c(0, 0.5)).

- beta.comm.sd:

  Community SD for occupancy (default c(0.5, 0.3)).

- alpha.comm.mean:

  Community mean for detection (default c(0)).

- alpha.comm.sd:

  Community SD for detection (default c(0.5)).

- seed:

  Random seed.

## Value

A list with `y` (3D array), `data`, and `truth`.

## Examples

``` r
sim <- simulate_ms_occu(N = 40, J = 3, n.species = 4, seed = 1)
dim(sim$y)
```
