# Per-species coefficients from a batched occu_cover fit.

Per-species coefficients from a batched occu_cover fit.

## Usage

``` r
# S3 method for class 'tobs_batch'
coef(object, ...)
```

## Arguments

- object:

  A `tobs_batch` returned by
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) on
  multi-response `y`.

- ...:

  Passed to each per-species
  [`coef()`](https://rdrr.io/r/stats/coef.html).

## Value

A named list of per-species coefficient vectors.
