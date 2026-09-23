# Estimate species richness from community model

Per-site expected species richness `sum_s psi_{s,i}` with a posterior
credible interval propagated from the community-mean occupancy draws.
For the dynamic community family this is the first-season richness (from
`psi1`).

## Usage

``` r
tobs_richness(object)
```

## Arguments

- object:

  A `tobs_fit` object from a community occupancy model
  ([`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
  [`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
  or
  [`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md)).

## Value

A data.frame with site-level richness estimates.

## Examples

``` r
# \donttest{
sim <- simulate_ms_occu(N = 40, J = 3, n_species = 5, seed = 1)
fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
            y = sim$y, species = paste0("sp", 1:5),
            control = list(verbose = FALSE))
head(tobs_richness(fit))
# }
```
