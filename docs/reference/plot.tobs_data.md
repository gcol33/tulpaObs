# Plot detection history patterns

2x2 panel: detection frequency histogram, per-visit detection rates,
visit completeness, and spatial detection map (if coordinates
available).

## Usage

``` r
# S3 method for class 'tobs_data'
plot(x, ...)
```

## Arguments

- x:

  An `tobs_data` object.

- ...:

  Ignored.

## Value

Invisible `NULL`.

## Examples

``` r
sim <- simulate_occu(N = 30, J = 3, seed = 1)
long <- data.frame(site = rep(seq_len(30), times = 3),
                   visit = rep(1:3, each = 30),
                   det = as.vector(sim$y),
                   occ_cov1 = sim$data$occ_cov1,
                   effort = runif(90))
dat <- tobs_data(long, y = "det", site = "site", visit = "visit",
                 occ.covs = "occ_cov1", det.covs = "effort")
plot(dat)
```
