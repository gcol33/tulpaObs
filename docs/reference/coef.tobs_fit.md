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

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
coef(fit)
coef(fit, arm = "psi")
# }
```
