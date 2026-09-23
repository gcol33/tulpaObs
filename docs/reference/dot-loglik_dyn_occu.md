# R-side multi-season HMM log-likelihood for dyn_occu

Computes the marginal observation log-likelihood for the dynamic
occupancy model under the C++ forward-algorithm semantics in
`src/dyn_occ_likelihood.h`. Used by
[`.sla_gamma_fd()`](https://gillescolling.com/tulpaObs/reference/dot-sla_gamma_fd.md)
to evaluate the third-cumulant correction; not exported.

## Usage

``` r
.loglik_dyn_occu(beta, model)
```

## Arguments

- beta:

  Numeric joint coefficient vector, length p_psi+p_p+p_col+p_ext.

- model:

  A `tobs_model` of `model_type = "dynamic"`.

## Value

Scalar log-likelihood.

## Details

Joint coefficient vector layout (matches `model$process_info`): beta =
c(beta_psi1, beta_p, beta_gamma, beta_epsilon)

Site-level or season-varying detection: the
[`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md)
detection design is site-level, optionally season-indexed (a
`[n_sites x n_seasons]` covariate, \#124); there is no visit-level
(within-season) detection covariate.

No priors. No pseudo-binomial encoding. Returns the original-data
log-likelihood, matching the convention of
[`.loglik_occu_single()`](https://gillescolling.com/tulpaObs/reference/dot-loglik_occu_single.md).
