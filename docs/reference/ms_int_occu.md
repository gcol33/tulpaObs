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

## Scope

The Laplace engine is the supported route: the shared occupancy mean and
the per-source detection community components recover near nominal
across seeds and more than one source (see
`tests/testthat/test-ms-int-occu.R`). A NUTS sampler and an areal-field
path are a deliberate follow-up, not part of this family's working
surface; `method = "nuts"` / `"nested_laplace"` error from the
dispatcher with a pointer rather than silently downgrading. The binary
community-mean intervals carry the mild Laplace under-dispersion typical
of occupancy data (measured 95% CI coverage ~0.89 at small per-species
n).

## See also

[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md),
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)

## Examples

``` r
f <- ms_int_occu()
f
```
