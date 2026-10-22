# Summarise occupancy data

Returns detection statistics including naive occupancy, naive detection,
per-visit rates, and detection frequency table.

## Usage

``` r
# S3 method for class 'tobs_data'
summary(object, ...)
```

## Arguments

- object:

  An `tobs_data` object.

- ...:

  Ignored.

## Value

An `tobs_data_summary` object (printed automatically).

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
summary(dat)
```
