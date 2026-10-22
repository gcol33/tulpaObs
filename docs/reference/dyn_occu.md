# Dynamic (multi-season) occupancy family

Bernoulli occupancy state with colonisation + extinction transitions
across seasons (the MacKenzie et al. dynamic model). Colonisation /
extinction are constant across a site's seasons by default; a
season-varying rate is supplied as a `[n_sites x (T - 1)]` matrix column
of `data` (one column per transition interval) on `colonization` /
`extinction`. Detection is site-level by default; a `[n_sites x T]`
matrix column of `data` (one column per primary season) on the
`detection` formula gives per-season detection. All season-varying paths
run under `method = "laplace"` (the exact HMM-forward marginal refine
calibrates the coefficients); they are gated under `"nuts"` with a
pointer.

## Usage

``` r
dyn_occu()
```

## Value

A `tobs_family` object.

## Engines

- `method = "laplace"` – EM + Laplace over the exact HMM-forward
  marginal, the default and the route every season-varying path runs on.

- `method = "laplace_sla"` – the Laplace fit with simplified-Laplace
  skew-corrected marginals.

- `method = "laplace_gibbs"` / `method = "laplace_mi"` – the Laplace fit
  with a post-EM Gibbs / multiple-imputation variance correction.

- `method = "nested_laplace"` – an areal field on the occupancy arm
  through the nested-Laplace EM, its hyperparameters integrated on an
  outer grid.

- `method = "nuts"` – the unified C++ NUTS sampler over the HMM
  marginal; season-varying colonisation / extinction / detection error
  with a pointer.

## Examples

``` r
f <- dyn_occu()
f
```
