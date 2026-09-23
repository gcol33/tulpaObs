# Simulate multi-species occupancy data

Simulate multi-species occupancy data

## Usage

``` r
simulate_ms_occu(
  N = 100,
  J = 4,
  n_species = 10,
  beta_comm_mean = c(0, 0.5),
  beta_comm_sd = c(0.5, 0.3),
  alpha_comm_mean = c(0),
  alpha_comm_sd = c(0.5),
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 100).

- J:

  Number of visits (default 4).

- n_species:

  Number of species (default 10).

- beta_comm_mean:

  Community mean for occupancy (default c(0, 0.5)).

- beta_comm_sd:

  Community SD for occupancy (default c(0.5, 0.3)).

- alpha_comm_mean:

  Community mean for detection (default c(0)).

- alpha_comm_sd:

  Community SD for detection (default c(0.5)).

- seed:

  Random seed.

## Value

A list with `y` (3D array), `data`, and `truth`.

## Examples

``` r
sim <- simulate_ms_occu(N = 40, J = 3, n_species = 4, seed = 1)
dim(sim$y)
```
