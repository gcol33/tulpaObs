# PIT residuals for a tobs fit

The
[`tulpa::pit_residuals()`](https://gillescolling.com/tulpa/reference/pit_residuals.html)
method for `tobs_fit`: the family's own predictive CDF, with the
randomisation step a discrete response needs.

## Usage

``` r
# S3 method for class 'tobs_fit'
pit_residuals(
  object,
  n.samples = 250,
  nsim = NULL,
  seed = NULL,
  observed = NULL,
  ...
)
```

## Arguments

- object:

  A `tobs_fit` object.

- n.samples:

  Number of posterior samples (default 250).

- nsim:

  Alias for `n.samples`, and the name
  [`tulpa::test_uniformity()`](https://gillescolling.com/tulpa/reference/test_uniformity.html)
  forwards under. An explicit `nsim` wins, so a uniformity test run at
  `nsim = 2000` draws 2000 samples rather than the default 250.

- seed:

  Optional RNG seed for the posterior-draw selection, forwarded by
  [`tulpa::test_uniformity()`](https://gillescolling.com/tulpa/reference/test_uniformity.html).
  Set it and repeated calls return the same residuals; the caller's RNG
  stream is restored afterwards.

- observed:

  Must be `NULL`, which is what
  [`tulpa::test_uniformity()`](https://gillescolling.com/tulpa/reference/test_uniformity.html)
  forwards by default: a latent-state fit scores the response it was
  fitted to.

- ...:

  Must be empty; an unrecognised argument is an error.

## Value

Numeric vector of PIT residuals.

## See also

[`tulpa::test_uniformity()`](https://gillescolling.com/tulpa/reference/test_uniformity.html)
for the uniformity test on the result.
