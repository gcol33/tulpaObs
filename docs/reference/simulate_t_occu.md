# Simulate multi-season occupancy with an AR1 year effect (tPGOcc)

Simulate multi-season occupancy with an AR1 year effect (tPGOcc)

## Usage

``` r
simulate_t_occu(
  N = 150,
  T_seasons = 8,
  J = 3,
  beta_occ = c(0.2),
  p = 0.4,
  rho = 0.6,
  sigma = 0.7,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 150).

- T_seasons:

  Number of seasons / years (default 8).

- J:

  Visits per season (default 3).

- beta_occ:

  Occupancy coefficients (intercept first; default `c(0.2)`).

- p:

  Detection probability (default 0.4).

- rho:

  AR1 correlation of the year effect (default 0.6).

- sigma:

  AR1 innovation SD of the year effect (default 0.7).

- seed:

  Optional random seed.

## Value

A list with `y` (3D array `[N x T x J]`), `data`, and `truth`.

## See also

[`t_occu()`](https://gillescolling.com/tulpaObs/reference/t_occu.md),
the family this simulates for.

## Examples

``` r
sim <- simulate_t_occu(N = 30, T_seasons = 4, J = 2, seed = 1)
dim(sim$y)
```
