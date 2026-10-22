# Fitted values (occupancy and detection probabilities)

Fitted values (occupancy and detection probabilities)

## Usage

``` r
# S3 method for class 'tobs_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Ignored.

## Value

A list with `psi` (occupancy probabilities), `p` (detection
probabilities), and `z` (posterior state probability at posterior-mean
parameters). For an
[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
fit `p` is a named list, one detection probability vector per source
(matching `y`'s names), since each source carries its own detection
design and, when present, its own detection-arm spatial field. For
single-season,
[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
and community models `z` is the per-site Bayes posterior `P(z=1 | y)`
(an
[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
site is occupied with certainty if ANY source detected it, and otherwise
pools every source's non-detection evidence); for dynamic models `z` is
the forward-backward (HMM smoothing) state posterior
`P(z_t=1 | y_{1:T})` as an `[n_sites x n_seasons]` matrix.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
str(fitted(fit))
# }
```
