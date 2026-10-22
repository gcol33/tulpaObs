# Fixed-effect covariance and credible intervals for tobs_fit

On a nested-Laplace fit that keeps the per-cell pieces of its outer grid
(for example
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
or [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) with
an areal field), the covariance is the grid-marginalized one and the
bounds invert the Gaussian mixture the grid defines, rather than being
read off the Monte Carlo draws the fit also carries;
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) read the same
mixture. Every other fit reads as
[`tulpa::vcov.tulpa_fit()`](https://gillescolling.com/tulpa/reference/vcov.tulpa_fit.html)
/
[`tulpa::confint.tulpa_fit()`](https://gillescolling.com/tulpa/reference/confint.tulpa_fit.html)
describe.

## Usage

``` r
# S3 method for class 'tobs_fit'
vcov(object, ...)

# S3 method for class 'tobs_fit'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Ignored.

- parm:

  Parameter names or indices (default: all fixed effects). Indices
  address the fixed effects; a name that is not a fixed effect is read
  off the posterior draws when they carry a column of that name.

- level:

  Interval level (default 0.95).

## Value

[`vcov()`](https://rdrr.io/r/stats/vcov.html): the fixed-effect
covariance matrix. [`confint()`](https://rdrr.io/r/stats/confint.html):
a two-column matrix of lower and upper bounds, one row per requested
parameter, carrying the `interval_source`, `interval_declined`,
`retained_mass` and `skew_applied` attributes
[`tulpa::confint.tulpa_fit()`](https://gillescolling.com/tulpa/reference/confint.tulpa_fit.html)
documents.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
vcov(fit)
confint(fit, level = 0.9)
# }
```
