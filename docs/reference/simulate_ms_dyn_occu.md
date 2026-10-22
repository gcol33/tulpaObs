# Simulate temporal multi-species occupancy data

Simulate temporal multi-species occupancy data

## Usage

``` r
simulate_ms_dyn_occu(
  N = 50,
  J = 3,
  n.species = 5,
  n.seasons = 4,
  beta.comm.mean = c(0),
  beta.comm.sd = c(0.5),
  gamma = 0.15,
  epsilon = 0.1,
  field = NULL,
  trend = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 50).

- J:

  Visits per season (default 3).

- n.species:

  Number of species (default 5).

- n.seasons:

  Number of seasons (default 4).

- beta.comm.mean:

  Community mean for occupancy (default c(0)).

- beta.comm.sd:

  Community SD for occupancy (default c(0.5)).

- gamma:

  Colonization probability (default 0.15).

- epsilon:

  Extinction probability (default 0.1).

- field:

  Optional per-site shared areal field (length `N`) added to the
  first-season occupancy logit of every species – the shared field of
  the community dynamic-spatial model (stMsPGOcc). Default `NULL` (no
  field).

- trend:

  Optional per-site shared varying-coefficient field (length `N`), added
  to the first-season occupancy logit as `trend_i * x_i` (weighted by
  the site covariate `x`) – the varying-coefficient field of the
  community dynamic-spatial model (svcTMsPGOcc). Default `NULL`.

- seed:

  Random seed.

## Value

A list with `y` (4D array), `data`, and `truth`.

## See also

[`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
the family this simulates for.

## Examples

``` r
sim <- simulate_ms_dyn_occu(N = 30, J = 2, n.species = 3, n.seasons = 3,
                            seed = 1)
dim(sim$y)
```
