# Beta-arm log-likelihood (positive subset only)

Beta log-density at the encoded positive responses with logit link,
summed over the positive subset. See `dev_notes/` for the closed form;
phi is treated as a fixed nuisance parameter (held at the fitted value)
and SLA gamma is computed for beta_pos only.

## Usage

``` r
.loglik_cover_pos_beta(beta_pos, phi, enc)
```

## Arguments

- beta_pos:

  Length-`p_pos` numeric coefficient vector.

- phi:

  Beta precision (positive scalar; held fixed at the fitted value).

- enc:

  Encoded data list from
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md)
  (uses `enc$pos_data$y` and `enc$pos_data$X`).

## Value

Scalar log-likelihood.
