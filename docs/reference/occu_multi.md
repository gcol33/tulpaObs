# Multi-species co-occurrence occupancy family

Joint occupancy of `S` species where the occupancy state `z in {0,1}^S`
follows a log-linear model with first-order (per species) and
second-order (per species pair) natural parameters (Rota et al. 2016;
unmarked `occuMulti`). The second-order parameters are the species
interactions (positive = co-occur more than independent; negative =
avoid). The state `formula` is the shared occupancy covariate design
(each natural parameter carries its own coefficients); `detection` the
shared site-level per-species detection design. The latent state
integrates out by enumerating the `2^S` states, so the exact marginal is
maximised with an observed-information vcov.

## Usage

``` r
occu_multi()
```

## Value

A `tobs_family` object for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Details

`y` is a length-`S` list of `N x J` 0/1/NA detection matrices (or a 3D
`[sites x visits x species]` array); `species` names the arms.

## Examples

``` r
# \donttest{
sim <- simulate_occu_multi(S = 2, N = 300, seed = 1)
fit <- tobs(~ scov1, data = sim$data, family = occu_multi(),
            detection = ~ 1, y = sim$y, species = sim$species,
            control = list(verbose = FALSE))
coef(fit)
# }
```
