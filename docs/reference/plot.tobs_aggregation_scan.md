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

## Examples

``` r
set.seed(1)
d <- data.frame(east = runif(200), north = runif(200),
                year = sample(2001:2006, 200, replace = TRUE),
                detected = rbinom(200, 1, 0.3))
scan <- occu_aggregation_scan(d, response = "detected",
                              coords = c("east", "north"), year = "year",
                              cell.sizes = c(0.2, 0.4),
                              block.lengths = c(2, 3))
plot(scan)
```
