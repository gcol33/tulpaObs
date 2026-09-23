# Community (multispecies) dynamic occupancy family

Per-species dynamic (HMM) occupancy with Gaussian community hyperpriors
on the per-species first-season occupancy and detection coefficients,
and shared community-wide colonisation / extinction transition
coefficients.

## Usage

``` r
ms_dyn_occu()
```

## Value

A `tobs_family` object.

## Scope

The Laplace engine is the supported route: the shared colonisation /
extinction dynamics and the per-species first-season occupancy /
detection components recover across seeds (see
`tests/testthat/test-ms-dyn-occu.R`, community-mean 95% CI coverage
measured ~0.98). A shared areal field on the first-season occupancy
formula (`~ 1 + icar(graph = adj)`) fits the `spOccupancy` `stMsPGOcc`
model under `method = "nested_laplace"`: the field is shared across
species, and because the first-season occupancy `psi1` only sets the
initial mixing weight of each species' HMM, the block-coordinate driver
alternates the community EM (field as a `psi1` offset) with an areal
field Newton, the field recovering cleanly (`cor` ~0.94). A
spatially-varying-coefficient bar
(`~ spatial(~ 1 + w || cell, graph = adj)`) adds a shared
covariate-weighted field alongside the intercept field, fitting the
`svcTMsPGOcc` model through the same K-field weighted-ICAR solve as the
community count SVC (both fields recover, `cor` ~0.90 / ~0.89).
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
only; a NUTS sampler and
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
fields are a deliberate follow-up; `method = "nuts"` errors from the
dispatcher with a pointer rather than silently downgrading.

## See also

[`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md),
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)

## Examples

``` r
f <- ms_dyn_occu()
f
```
