# Coefficients for a tobs_fit

The fixed-effect estimates as one named numeric vector, named
`<arm>_<term>` as [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html)
and [`stats::confint()`](https://rdrr.io/r/stats/confint.html) name
them: `psi_(Intercept)`, `p_det_cov1`, `lambda_elev`, ... The
visit-level detection coefficients (`p_visit_<cov>`) belong to the
detection arm. A coordinate that enters no single arm, such as the
`log_r` overdispersion, keeps its own name.

## Usage

``` r
# S3 method for class 'tobs_fit'
coef(object, arm = NULL, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- arm:

  Optional name of one linear predictor (`"psi"`, `"p"`, ...). When
  given, only that arm's coefficients are returned, named by their
  design column without the arm prefix.

- ...:

  Ignored.

## Value

A named numeric vector.
