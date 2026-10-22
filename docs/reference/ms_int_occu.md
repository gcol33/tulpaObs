# Community (multispecies) integrated occupancy family

Multiple detection data sources share one latent occupancy state per
species, with Gaussian community hyperpriors on the per-species
coefficients of the occupancy arm and of every per-source detection arm.
Fit by a shared community Laplace-EM.

## Usage

``` r
ms_int_occu()
```

## Value

A `tobs_family` object.

## Engines

- `method = "laplace"` – the shared community Laplace-EM, the default:
  the shared occupancy mean and the per-source detection community
  components recover near nominal across seeds and more than one source
  (see `tests/testthat/test-ms-int-occu.R`). The binary community-mean
  intervals carry the mild Laplace under-dispersion typical of occupancy
  data (measured 95% CI coverage ~0.89 at small per-species n).

- `method = "pg_gibbs"` – the community integrated Polya-Gamma Gibbs
  sampler (one latent state per species-site observed by every source),
  giving a calibrated community-variance posterior; site-level
  per-source detection.

- `method = "nuts"` – the community integrated sampler over the exact
  multi-source two-state per-(species, site) marginal (in-tree C++),
  sampling the community means, per-species deviations and the `D + 1`
  per-arm community covariances jointly, warm-started at the Laplace-EM
  mode.

No route takes a structured term: a field or random effect on either
formula errors with a pointer, and there is no nested-Laplace path for
this family.

## See also

[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md),
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)

## Examples

``` r
f <- ms_int_occu()
f
```
