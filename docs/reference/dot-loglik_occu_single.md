# R-side single-season occupancy log-likelihood

Computes the marginal observation log-likelihood
`log P(y | beta) = Sum_i log[ psi_i * P(y_i | z=1) + (1 - psi_i) I(no detections) ]`
with psi_i = plogis(X_occ_i beta_occ), p_i = plogis(X_det_i beta_det).

## Usage

``` r
.loglik_occu_single(beta, model)
```

## Details

Used to cross-check
[`.sla_gamma_fd()`](https://gillescolling.com/tulpaObs/reference/dot-sla_gamma_fd.md)
against the analytical
[`.sla_gamma_diag()`](https://gillescolling.com/tulpaObs/reference/dot-sla_gamma_diag.md)
on single-season fits.
