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

## Engines

- `method = "laplace"` – EM + Laplace over the multi-source two-state
  marginal, the default.

- `method = "laplace_sla"` – the Laplace fit with simplified-Laplace
  skew-corrected marginals.

- `method = "laplace_gibbs"` / `method = "laplace_mi"` – the Laplace fit
  with a post-EM Gibbs / multiple-imputation variance correction.

- `method = "nested_laplace"` – an areal field on the state arm through
  the nested-Laplace EM, its hyperparameters integrated on an outer
  grid.

- `method = "nuts"` – the unified C++ NUTS sampler over the same
  marginal.

## Examples

``` r
f <- int_occu()
f
```
