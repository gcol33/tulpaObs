# Simulate replicate datasets from the posterior

Each replicate draws a parameter vector from the fit's posterior draws
and simulates the latent state and the response from it, so the
replicates carry parameter uncertainty.

## Usage

``` r
# S3 method for class 'tobs_fit'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- nsim:

  Number of simulated datasets (default 1).

- seed:

  Optional random seed. As for
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html), a
  supplied seed is set for the call and the caller's random number
  stream is restored afterwards.

- ...:

  Must be empty.

## Value

A list of length `nsim`, one simulated response per element (named
`sim_1`, `sim_2`, ...), each in the layout the family's `y` takes: a
detection matrix, a `[sites x visits x seasons]` array, a list of source
or species matrices, and so on. The list carries the `"seed"` attribute
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) specifies:
the `.Random.seed` in force when `seed` is `NULL`, otherwise `seed` with
its [`RNGkind()`](https://rdrr.io/r/base/Random.html) as the `"kind"`
attribute. Every family returns this layout, whatever `nsim` is.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
y_rep <- simulate(fit, nsim = 2, seed = 1)
lapply(y_rep, colMeans)
# }
```
