# Posterior predictive check

The discrepancy-statistic check (the spOccupancy `ppcOcc` construction):
latent occupancy is drawn from its full conditional, a detection
replicate from the fitted model, and the two are compared through a
discrepancy statistic to give a Bayesian p-value. This is a test, not a
plot, so it is its own generic rather than a method on
[`tulpa::pp_check()`](https://gillescolling.com/tulpa/reference/pp_check.html),
which draws the graphical check.

## Usage

``` r
ppc(object, ...)

# S3 method for class 'tobs_fit'
ppc(object, fit.stat = c("freeman-tukey", "chi-squared"), n.samples = 500, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Passed to methods.

- fit.stat:

  `"freeman-tukey"` (default) or `"chi-squared"`.

- n.samples:

  Number of posterior samples (default 500).

## Value

A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 60, J = 3, seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE))
ppc(fit, n.samples = 100)$bayesian.p
# }
```
