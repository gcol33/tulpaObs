# Multispecies (community) single-season occupancy family

Per-species occupancy and detection with Gaussian community hyperpriors
on the per-species coefficients of each arm (the spOccupancy `msPGOcc`
model). The occupancy and detection per-species deviations are
independent, each with its own community covariance, fit by a community
Laplace-EM.

## Usage

``` r
ms_occu()
```

## Value

A `tobs_family` object.

## See also

[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
[`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
[`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md)

## Examples

``` r
# \donttest{
sim <- simulate_ms_occu(N = 80, J = 4, n_species = 8, seed = 1)
fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
            y = sim$y, species = paste0("sp", seq_len(8)), method = "laplace")
summary(fit)
# }
```
