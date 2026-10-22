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

## Engines

- `method = "laplace"` – the shared community Laplace-EM, the default;
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors alone (lfMsPGOcc) run here through the block-coordinate
  driver.

- `method = "nuts"` – the non-spatial community sampler over the
  closed-form per-(species, site) two-state marginal, sampling the
  community means, per-species deviations and both per-arm community
  covariances jointly, warm-started at the Laplace-EM mode; an areal
  term errors with a pointer to `"nested_laplace"`.

- `method = "pg_gibbs"` – the hierarchical Polya-Gamma Gibbs sampler
  (spOccupancy `msPGOcc`), giving a calibrated community-variance
  posterior; non-spatial.

- `method = "nested_laplace"` – a shared areal field
  ([`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md))
  on the occupancy arm: a plain intercept field runs the in-tree C++
  community-spatial nested Laplace-EM, a varying-coefficient bar
  (svcMsPGOcc) or
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors beside the field (sfMsPGOcc) run the block-coordinate fitter;
  the field is required on this route.

## See also

[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
[`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
[`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md)

## Examples

``` r
# \donttest{
sim <- simulate_ms_occu(N = 80, J = 4, n.species = 8, seed = 1)
fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
            y = sim$y, species = paste0("sp", seq_len(8)), method = "laplace")
summary(fit)
# }
```
