# Predict cover from a cover_fit

Occurrence probability is always `p = plogis(X * beta_occ)`. The
conditional positive cover `mu` depends on the positive-part family:

## Usage

``` r
# S3 method for class 'cover_fit'
predict(
  object,
  newdata = NULL,
  type = NULL,
  include.RE = FALSE,
  times = NULL,
  time.col = NULL,
  level = 0.95,
  nsim = 1000L,
  draws = TRUE,
  weights = NULL,
  aggregate = FALSE,
  mc.tol = 0.05,
  nsim.max = 10000L,
  mc.floor = 0.001,
  ...
)
```

## Arguments

- object:

  A `cover_fit`.

- newdata:

  A data frame of covariates matching the original formula. For the
  nested-Laplace fit, one row per spatial unit (or a `cell` column).

- type:

  Separate-Laplace fit: one of `"expected"`, `"occupancy"`,
  `"conditional"`. Nested-Laplace fit: `"occurrence"`, `"cover_cond"`,
  `"cover_exp"`, or `"change"`; the separate-Laplace names `"expected"`,
  `"occupancy"` and `"conditional"` map to `"cover_exp"`, `"occurrence"`
  and `"cover_cond"`.

- include.RE:

  Ignored for the separate-Laplace fit (no spatial projection); the
  nested-Laplace fit always projects the shared field.

- times, time.col, level, nsim, draws:

  Nested-Laplace fit only: `times = c(t1, t2)` and `time.col` drive the
  `"change"` map; `level` is the credible level, `nsim` the draw count
  (or `"auto"` for a trajectory, as in
  [`predict.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md)),
  `draws` whether to attach the draw matrices.

- weights:

  Nested-Laplace fit only: pool several `newdata` rows into one
  prediction per `cell`, as in
  [`predict.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md).
  A numeric vector (one weight per row) or the name of a `newdata`
  column.

- aggregate:

  Nested-Laplace fit, `type = "trajectory"` only: `TRUE` averages each
  quantity over the prediction units per posterior draw and returns one
  row per time, as in
  [`predict.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md).

- mc.tol, nsim.max, mc.floor:

  Nested-Laplace fit, trajectory only: the Monte Carlo precision target,
  the draw cap and the SD floor the error is measured against, as in
  [`predict.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md).

- ...:

  Unused.

## Value

Separate-Laplace fit: a numeric vector. Nested-Laplace fit: a
`tobs_prediction`.

## Details

- `positive = "lognormal"`: `mu = exp(eta_pos + sigma_pos^2 / 2)`
  (lognormal back-transform on log-cover).

- `positive = "beta"`: `mu = plogis(eta_pos)` (mean of the beta on the
  natural cover scale with logit link).

Expected cover is `E[y] = p * mu` under both positive parts.

The separate-Laplace fit (`method = "laplace"`) returns a
fixed-effects-only numeric vector. The nested-Laplace shared-field fit
(`method = "nested_laplace"`) instead projects the shared
occupancy-cover field and returns a `tobs_prediction` of posterior draws
– the same tidy / `change` contract as
[`predict.tobs_fit()`](https://gillescolling.com/tulpaObs/reference/predict.tobs_fit.md)
for
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md):
pass `type = "change"` with `times = c(t1, t2)` and `time.col` for a
per-cell delta map. Each prediction unit is a row of `newdata` (or a
`cell` column indexing the field cells), and every quantity is
marginalized per draw over the grid-integrated joint posterior.

## Examples

``` r
# \donttest{
fit <- tobs(cover ~ moisture, data = meadow_cover,
            family = cover(response = "beta"))
nd <- data.frame(moisture = c(-1, 0, 1))
predict(fit, newdata = nd, type = "occupancy")
predict(fit, newdata = nd, type = "conditional")
predict(fit, newdata = nd, type = "expected")
# }
```
