# Summary for tobs_fit, with skewness column when simplified Laplace is used

Extends `summary.tulpa_fit` with an extra `skew` column populated from
`object$skew` when the fit was produced with an SLA method
(`method = "laplace_sla"` / `"nested_laplace_sla"`). All quantile / mean
/ sd columns come from `object$draws`, which under simplified Laplace
are skew-normal samples – so 2.5%/97.5% quantiles are already
SLA-corrected.

## Usage

``` r
# S3 method for class 'tobs_fit'
summary(object, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Forwarded to `summary.tulpa_fit`.

## Value

Data frame as for `summary.tulpa_fit`, one row per coefficient named
`<arm>_<term>` as
[`coef.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/coef.tobs_fit.md)
names it, with columns `estimate`, `std.error` and the two interval
bounds, plus a `skew` column when `$skew` is present and `rhat`,
`ess_bulk`, `ess_tail` on a sampled fit. Every family, and every method
within a family, returns this layout.
