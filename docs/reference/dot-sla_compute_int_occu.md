# Compute SLA skewness coefficients for an integrated occu fit

Builds joint block-diagonal Sigma (Louis I_obs on the occ block, raw
H_beta on each det_s block), then calls the generic FD primitive against
the original-likelihood evaluator `.loglik_int_occu`.

## Usage

``` r
.sla_compute_int_occu(model, em_result, spatial = NULL, prior_spec = NULL)
```

## Arguments

- model:

  A `tobs_model` (model_type = "integrated").

- em_result:

  EM-Laplace return list (with `$fits$occ`,
  `$fits$det1`,...,`$fits$det<S>`, `$weights`).

- spatial:

  Optional `tobs_spatial`. When set, the skewness correction is
  intentionally NOT applied (Gaussian marginals retained) for the same
  reason as the single-season path: the simplified-Laplace
  third-cumulant correction captures hyperparameter-free fixed-effect
  marginal skewness, but a spatial field's marginal skewness is
  dominated by hyperparameter-marginalisation, which the correction does
  not capture (validated against NUTS). The Gaussian fallback is the
  correct conservative behaviour, not a stub.

- prior_spec:

  Optional prior spec; passed to `.louis_info_psi_single()` so the
  penalty enters Louis I_obs.

## Value

`list(gamma, valid, reason)`.
