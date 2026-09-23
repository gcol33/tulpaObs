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
  include_RE = FALSE,
  times = NULL,
  time_col = NULL,
  level = 0.95,
  nsim = 1000L,
  draws = TRUE,
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
  `"cover_exp"`, or `"change"` (the legacy aliases are accepted and
  mapped).

- include_RE:

  Ignored for the separate-Laplace fit (no spatial projection); the
  nested-Laplace fit always projects the shared field.

- times, time_col, level, nsim, draws:

  Nested-Laplace fit only: `times = c(t1, t2)` and `time_col` drive the
  `"change"` map; `level` is the credible level, `nsim` the draw count,
  `draws` whether to attach the draw matrices.

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
pass `type = "change"` with `times = c(t1, t2)` and `time_col` for a
per-cell delta map. Each prediction unit is a row of `newdata` (or a
`cell` column indexing the field cells), and every quantity is
marginalized per draw over the grid-integrated joint posterior.
