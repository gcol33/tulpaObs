# Compute marginal effect of a covariate

Returns the fitted response over the range of one covariate, holding the
others at their means, on the process's response scale (occupancy /
detection probability, or abundance intensity).

## Usage

``` r
tobs_marginal_effect(
  object,
  covariate,
  process = c("occupancy", "detection", "abundance"),
  n_points = 100L
)
```

## Arguments

- object:

  A `tobs_fit` object.

- covariate:

  Name of covariate.

- process:

  `"occupancy"` (default), `"detection"`, or `"abundance"` (the state
  process of a count family, on the intensity scale).

- n_points:

  Number of prediction points (default 100).

## Value

A data.frame with covariate values and the predicted response.

## Examples

``` r
# \donttest{
sim <- simulate_abun(N = 100, J = 4, n_abund_covs = 2, n_det_covs = 1,
                     seed = 1)
fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data,
            family = abun(K_max = 50), detection = ~ det_cov1, y = sim$y,
            control = list(verbose = FALSE))
me <- tobs_marginal_effect(fit, "abund_cov1", process = "abundance")
head(me)   # estimate is on the abundance (lambda) scale, not a probability
# }
```
