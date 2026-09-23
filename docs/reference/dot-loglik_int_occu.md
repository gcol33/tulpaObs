# Integrated multi-source occupancy log-likelihood (R-side)

Computes the marginal observation log-likelihood
`log P(y | beta) = Sum_i log( psi_i * Prod_s P(y_is | z=1, p_is) + (1 - psi_i) * I(no detections at any source) )`
where `psi_i = plogis(X_occ[i,] beta_psi)`,
`p_is = plogis(X_det_s[i,] beta_det_s)`.

## Usage

``` r
.loglik_int_occu(beta, model)
```

## Arguments

- beta:

  Joint parameter vector in process_info order:
  `c(beta_psi, beta_det_1, beta_det_2, ..., beta_det_S)`.

- model:

  A `tobs_model` with `model_type = "integrated"`.

## Value

Scalar log-likelihood.

## Details

Reads `model$y_sources` (list of integer matrices, source-local row
indexing, NA encoded as -1), `model$site_maps` (0-indexed global site
per source row), `model$X_processes` (occ design + one det design per
source, all n_sites x p), `model$n_sources`, `model$n_sites`. No priors,
no pseudo-binomial encoding.
