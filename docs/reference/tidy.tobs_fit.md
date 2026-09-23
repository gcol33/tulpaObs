# Coefficient table for a tobs_fit

Coefficient table for a tobs_fit

## Usage

``` r
# S3 method for class 'tobs_fit'
tidy(x, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `tobs_fit`.

- conf.level:

  Interval level (default 0.95).

- ...:

  Ignored.

## Value

A data frame with one row per coefficient: `arm` (the linear predictor
it enters, `NA` for a coordinate that enters none), `term` (its design
column, without the arm prefix), `estimate`, `std.error`, `conf.low` and
`conf.high`. Every family returns this layout.
