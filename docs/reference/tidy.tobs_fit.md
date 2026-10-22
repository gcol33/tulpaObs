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

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
tidy(fit)
# }
```
