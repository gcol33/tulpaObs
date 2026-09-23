# Residuals from occupancy model

One contract for every family: a list with `occ`, the STATE-level
residual series (per site, or per site and season for a multi-season
family), and `det`, the OBSERVATION-level series (per visit / pass /
distance bin). A family carrying only one of the two levels fills that
one and leaves the other `NULL`; nothing returns a bare matrix, so a
caller reading `$occ` gets the state residual or `NULL`, never an error.

## Usage

``` r
# S3 method for class 'tobs_fit'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- type:

  One of `"deviance"` (default), `"pearson"`, or `"response"`.

- ...:

  Ignored.

## Value

A list with `occ` (state-level) and `det` (observation-level) residuals,
either of which may be `NULL` for a family that has no such level. A
model type with no registered handler is an error naming it.
