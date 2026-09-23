# Integrated occupancy family

Multiple data sources informing a shared latent occupancy state, each
with its own detection process.

## Usage

``` r
int_occu()
```

## Value

A `tobs_family` object.

## Details

A continuous-mesh
[`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term on the `detection` formula loads a spatially-structured detection
field, broadcast to every source, under `method = "laplace"`. The
per-source fields are reported as a named list in
`fit$spatial_field_det`. All sources must share one detection formula;
areal terms belong on the state arm.

## Examples

``` r
f <- int_occu()
f
```
