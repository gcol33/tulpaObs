# Print method for tobs_fit

Print method for tobs_fit

## Usage

``` r
# S3 method for class 'tobs_fit'
print(x, ...)
```

## Arguments

- x:

  a `tobs_fit` object.

- ...:

  forwarded to underlying print methods.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
fit
# }
```
