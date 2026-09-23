# Log-likelihood of a tobs fit

The value and degrees of freedom are the engine's (the `tulpa_fit`
method); the `"nobs"` attribute is
[`nobs.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/nobs.tobs_fit.md)
on the same fit, so [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html)
and anything else reading the attribute counts the observations
[`nobs()`](https://rdrr.io/r/stats/nobs.html) reports.

## Usage

``` r
# S3 method for class 'tobs_fit'
logLik(object, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Passed on to the next method.

## Value

A `logLik` object.
