# Lognormal-arm log-likelihood (positive subset only)

Note: the encoded `enc$pos_data$y` is already on the log scale (i.e.
`log(cover[occur == 1])`; see
[`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md)).
On the log scale the lognormal density reduces to a Gaussian on `log y`
with mean `X beta_pos` and SD `sigma_pos`. The Jacobian `- log y` is
constant in `beta_pos` and so does not affect the third derivative; we
drop it. (For a finite-difference d3 along the beta_pos direction,
beta-independent terms cancel.)

## Usage

``` r
.loglik_cover_pos_lognormal(beta_pos, sigma_pos, enc)
```

## Arguments

- beta_pos:

  Length-`p_pos` numeric coefficient vector.

- sigma_pos:

  Positive scalar (held fixed at the fitted value).

- enc:

  Encoded data list from
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md).

## Value

Scalar log-likelihood.
