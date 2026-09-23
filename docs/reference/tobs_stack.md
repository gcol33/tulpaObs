# Stack fitted models into a LOO-weighted ensemble

Combines two or more fitted models into one predictive distribution with
leave-one-out predictive stacking. Each member contributes its
Pareto-smoothed importance-sampling LOO predictive density; the weights
maximize the stacked leave-one-out predictive density and so favour
members that predict held-out observations best. Works across every
family (single-season / community / dynamic / integrated occupancy,
JSDM, and the cover hurdle), each scored by the marginal likelihood in
`.tobs_pointwise_loglik()`.

## Usage

``` r
tobs_stack(..., method = c("stacking", "pseudobma"))

# S3 method for class 'tobs_stack'
print(x, ...)

# S3 method for class 'tobs_stack'
fitted(object, ...)

# S3 method for class 'tobs_stack'
predict(
  object,
  X.0 = NULL,
  quantiles = c(0.025, 0.5, 0.975),
  n.draws = 4000L,
  ...
)
```

## Arguments

- ...:

  ignored.

- method:

  weighting scheme passed to
  [`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html):
  `"stacking"` (default) or `"pseudobma"`.

- x:

  a `tobs_stack` object.

- object:

  a `tobs_stack` object.

- X.0:

  optional occupancy design matrix for out-of-sample prediction. When
  supplied, all members must share the same occupancy design.

- quantiles:

  credible-interval quantile levels.

- n.draws:

  size of the pooled stacked-predictive sample (design-matrix mode);
  each member contributes a share proportional to its weight.

## Value

An object of class `tobs_stack`: a list with `weights` (named numeric,
summing to 1), `fits` (the members), `loo` (per-member
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) objects),
`comparison` (a data frame of `elpd_loo`, `weight`, and worst Pareto-k
per member), and `method`.
[`stats::predict()`](https://rdrr.io/r/stats/predict.html) /
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html) on the
result return the weight-combined predictive.

## Details

The members must be fit to the same observation set. They may otherwise
differ freely – different covariates, priors, or inference `method`s.
(An ensemble of pure seed-variants, e.g. from
`tobs(..., control = list(n.seeds = K))`, is a valid but degenerate
case: the members are statistically identical, so the weights come out
roughly uniform.)

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n_occ_covs = 2, n_det_covs = 1,
                     seed = 1)
ctrl <- list(verbose = FALSE, progress = FALSE)
f1 <- tobs(~ occ_cov1, data = sim$data, family = occu(),
           detection = ~ det_cov1, y = sim$y, method = "laplace",
           control = ctrl)
f2 <- tobs(~ occ_cov1 + occ_cov2, data = sim$data, family = occu(),
           detection = ~ det_cov1, y = sim$y, method = "laplace",
           control = ctrl)
ens <- tobs_stack(simple = f1, full = f2)
ens$weights
# weight-combined in-sample psi / p / z
str(predict(ens))
# }
```
