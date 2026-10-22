# Number of observations

The count of observed response entries the fit's likelihood reads, which
is the sample size [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html)
and [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) form an
information criterion at. What one entry is follows the family's own
response layout: a surveyed (site, visit) cell for the occupancy
families, a (site, band) count for distance sampling, a (site, visit,
source) cell for the integrated ones.

## Usage

``` r
# S3 method for class 'tobs_fit'
nobs(object, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Ignored.

## Value

Integer count of the fit's observed response entries. A model type with
no registered count is an error naming that type.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
nobs(fit)
# }
```
