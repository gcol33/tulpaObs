# Plot an occupancy aggregation scan

Heatmap of the identifiability score over the scanned (cell size x
temporal block) grid: information smallest-eigenvalue for
`score = "info"`, fraction of replicated units for `score = "count"`.
Identifiable cells are outlined.

## Usage

``` r
# S3 method for class 'tobs_aggregation_scan'
plot(x, ...)
```

## Arguments

- x:

  A `tobs_aggregation_scan` object.

- ...:

  Ignored.

## Value

Invisible `NULL`.
