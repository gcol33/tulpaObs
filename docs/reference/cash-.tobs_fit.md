# Access spOccupancy-compatible fields from tobs fits

Allows accessing spOccupancy-style fields (e.g., `$beta.samples`,
`$psi.samples`) on tobs_fit objects. Since tobs stores actual posterior
draws, this is a thin remapping layer.

## Usage

``` r
# S3 method for class 'tobs_fit'
x$name
```

## Arguments

- x:

  A `tobs_fit` object.

- name:

  Field name to access.

## Value

The requested field value.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
dim(fit$beta.samples)
dim(fit$alpha.samples)
# }
```
