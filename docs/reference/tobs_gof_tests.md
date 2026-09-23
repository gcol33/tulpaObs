# Goodness-of-fit tests for a tobs fit

The tulpa goodness-of-fit generics for `tobs_fit`: dispersion,
zero-inflation, and outlier counts compared against the fitted model's
simulated replicates. For the Kolmogorov-Smirnov uniformity test, hand
[`pit_residuals()`](https://gillescolling.com/tulpa/reference/pit_residuals.html)
to
[`tulpa::test_uniformity()`](https://gillescolling.com/tulpa/reference/test_uniformity.html),
which takes a PIT vector directly.

## Usage

``` r
# S3 method for class 'tobs_fit'
test_dispersion(
  object,
  n.samples = 250,
  nsim = NULL,
  seed = NULL,
  observed = NULL,
  alternative = c("greater", "two.sided", "less"),
  ...
)

# S3 method for class 'tobs_fit'
test_zero_inflation(
  object,
  n.samples = 250,
  nsim = NULL,
  seed = NULL,
  observed = NULL,
  ...
)

# S3 method for class 'tobs_fit'
test_outliers(
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

  A fitted `tobs_fit`.

- n.samples:

  Number of posterior-predictive replicates to simulate.

- nsim:

  The same budget under the name the tulpa default methods use;
  overrides `n.samples` when given.

- seed:

  Optional random seed for the replicates; the caller's RNG stream is
  restored afterwards.

- observed:

  Must be `NULL`: a tobs fit scores the response it was fitted to.
  Accepted so a call spelled for the tulpa default reaches the method.

- alternative:

  Tail of the dispersion p-value: `"greater"` (default,
  over-dispersion), `"two.sided"` or `"less"`.

- ...:

  Must be empty; an unrecognised argument is an error.

## Value

A list with the observed statistic (`observed`), its
posterior-predictive expectation (`expected`), their ratio (`ratio`) and
a tail p-value (`p.value`). The names are the same on every family, so a
caller can read the statistic without branching on the model type.
