# Compute SLA skewness coefficients for a dyn_occu (HMM) fit

Drives the generic finite-difference d3 primitive against the R-side HMM
forward log-likelihood. The joint Sigma is assembled block-diagonally
with per-block precision matrices (see file-header note for the
justification).

## Usage

``` r
.sla_compute_dyn_occu(model, em_result, spatial = NULL, prior_spec = NULL)
```

## Arguments

- model:

  A `tobs_model` (model_type = "dynamic").

- em_result:

  The EM-Laplace return list with `$fits$occ/det/col/ext` and `$weights`
  (n_sites x n_seasons matrix).

- spatial:

  Optional `tobs_spatial`. A state-arm field reaches this evaluator
  (`.validate_spatial_laplace()` blocks only a detection-arm SPDE for
  the dynamic model), and the correction is intentionally NOT applied
  under one – see `.sla_spatial_reason()`, the same decision the single-
  season and integrated paths take. Returns `valid = FALSE`.

- prior_spec:

  Optional prior spec; passed to the init-block Louis info so the
  penalty is included.

## Value

List(gamma, valid, reason).
