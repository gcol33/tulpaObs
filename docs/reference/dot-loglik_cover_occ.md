# Occurrence-arm log-likelihood (Bernoulli)

Computes the per-arm log-likelihood for the cover-hurdle occurrence
indicator at a candidate `beta_occ`, with no prior or pseudo-binomial
encoding. Used by
[`.sla_gamma_fd()`](https://gillescolling.com/tulpaObs/reference/dot-sla_gamma_fd.md).

## Usage

``` r
.loglik_cover_occ(beta_occ, enc)
```

## Arguments

- beta_occ:

  Length-`p_occ` numeric coefficient vector.

- enc:

  Encoded data list from
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md)
  (uses `enc$occ_data$y` and `enc$occ_data$X`).

## Value

Scalar log-likelihood.
