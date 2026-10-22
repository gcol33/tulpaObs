# Predict from occupancy model

Four modes:

- **In-sample**: `predict(fit)` returns fitted values.

- **State posterior / NA-response**: `predict(fit, type = "state")`
  returns the marginalised per-site occupancy posterior of a
  `"nested_laplace"` fit: the per-row eta posterior is a Gaussian
  mixture over the hyperparameter grid (per-cell fitted linear predictor
  and predictive variance), so `psi` is its Gauss-Hermite mean and
  `psi_lower` / `psi_upper` are the mixture-CDF quantiles – a calibrated
  95% credible interval, exact for every latent prior. Rows with
  `heldout = TRUE` are the INLA-style NA-response prediction targets
  (single-season sites whose detection history was all-missing), where
  occupancy is interpolated from the latent field rather than informed
  by detections.

- **Design-matrix**: `predict(fit, X.0 = ...)` predicts at new covariate
  values.

- **Terms-based**: `predict(fit, terms = "elev")` varies one covariate,
  others at mean.

- **Joint occu_cover**: for a `occu_cover` joint-coupled fit,
  `predict(fit, newdata, type = "occurrence" | "cover_cond" | "cover_exp" | "change")`
  samples the joint latent from the grid-integrated posterior (the
  outer-grid mixture via
  [`tulpa::tulpa_posterior_draws()`](https://gillescolling.com/tulpa/reference/tulpa_posterior_draws.html))
  and marginalises every derived quantity per draw. `type = "change"`
  with `times = c(t1, ..., tK)` (K \>= 2) returns a per-cell change
  table. Each quantity – occupancy `psi`, conditional cover
  `cover_cond`, expected cover `cover_exp` – has a level column per step
  (`psi_T1..psi_TK`, ...) carrying the posterior mean with `.sd` /
  `.lwr` / `.upr` at `level`, and a delta per step against `t1`
  (`delta_psi`, `delta_cover_cond`, `delta_cover_exp`) with `.lwr` /
  `.upr` and `.prob_pos`, the posterior probability `P(delta > 0)` per
  cell. The expected-cover change is split exactly into its occupancy-
  and abundance-driven parts (`delta_cover_from_occ`,
  `delta_cover_from_ab`), reported the same way. With two times the
  delta columns are unsuffixed; with K \> 2 each carries `_T<k>`. Every
  step is evaluated on ONE draw set, so levels and deltas share a
  posterior and are jointly valid. A standalone
  [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) joint
  fit returns the same layout with occupancy only. The result is a
  `tobs_prediction` table (one row per cell) carrying per-unit
  `[cell x nsim]` draw matrices in `attr(, "draws")`; map it yourself,
  e.g. `left_join(cents, pr, by = "cell")` then
  `geom_tile(aes(x, y, fill = delta_psi))` (or `geom_sf()` on polygon
  cells). `type = "trajectory"` with `times = c(t1, ..., tK)` predicts
  occupancy `psi`, conditional cover `cover_cond` and expected cover
  `cover_exp` at each time, for any times the user chooses, as a long
  table with the posterior mean, median and `level` interval of each
  (one row per cell and time); `aggregate = TRUE` averages over the
  cells per draw and gives the overall trajectory, one row per time. All
  times share one draw set.

- **Community families**:
  [`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
  [`ms_dyn_occu()`](https://gillescolling.com/tulpaObs/reference/ms_dyn_occu.md),
  [`ms_int_occu()`](https://gillescolling.com/tulpaObs/reference/ms_int_occu.md),
  [`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md),
  [`ms_distance()`](https://gillescolling.com/tulpaObs/reference/ms_distance.md)
  and
  [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)
  predict per-species matrices `[rows x species]`, in-sample from
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and at
  `newdata` from the per-species coefficients. Each family answers to
  its own response types: `"occupancy"` / `"detection"` (`ms_occu`,
  `ms_dyn_occu`, `ms_int_occu`), `"abundance"` / `"detection"`
  (`ms_abun`), `"lambda"` / `"sigma"` (`ms_distance`), and `"occupancy"`
  / `"detection"` / `"cover_cond"` / `"cover_exp"` (`ms_occu_cover`). A
  fit carrying an areal field or a latent-factor loading is in-sample
  only, the field having no value at an unseen site.

- **Spatial-factor community**: for a reduced-rank spatial-factor
  [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)
  fit, `predict(fit, type = "occupancy" | "cover_cond" | "cover_exp")`
  returns the per-species per-cell posterior as a long table (one row
  per cell x species): `"occupancy"` gives `psi`, `"cover_cond"` the
  conditional cover mean `E[cover | present]`, `"cover_exp"` the
  unconditional expected cover `psi * E[cover | present]`, each with a
  `_median` and a `_lower` / `_upper` interval. The maps are
  marginalised over the loading + field posterior, so a rare species
  borrows strength across the shared factors for a calibrated map. The
  latent fields are tied to the cell graph, so `X.0` / `newdata`
  prediction is not supported.

## Usage

``` r
# S3 method for class 'tobs_fit'
predict(
  object,
  X.0 = NULL,
  type = c("occupancy", "detection", "both", "state"),
  quantiles = c(0.025, 0.5, 0.975),
  terms = NULL,
  n.points = 50L,
  newdata = NULL,
  times = NULL,
  level = 0.95,
  nsim = 1000L,
  draws = TRUE,
  time.col = NULL,
  X.det.0 = NULL,
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

  A `tobs_fit` object.

- X.0:

  Optional design matrix for occupancy prediction.

- type:

  `"occupancy"` (default), `"detection"`, `"both"`, or `"state"`
  (nested-Laplace marginalised per-site psi, incl. held-out sites). For
  an `occu_cover` fit: `"occurrence"`, `"cover_cond"`, `"cover_exp"`,
  `"change"` or `"trajectory"`. `"detection"` / `"both"` on a
  `single`/[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  fit needs `X.det.0`.

- quantiles:

  Quantile levels for credible intervals: a strictly increasing length-3
  vector (lower, middle, upper), each in (0, 1). The levels name the
  columns they fill, so `c(0.05, 0.5, 0.95)` reports `q5` / `q50` /
  `q95`, and the levels used travel on the returned table in
  `attr(, "quantiles")`.

- terms:

  Name of a single term to vary: one column of the predicted process's
  design matrix, every other column held at its column mean. A vector of
  more than one term is an error, since a second term would be held at
  its mean rather than grouped over its levels. For a grid over more
  than one covariate, build the design matrix and pass it as `X.0`.

- n.points:

  Number of prediction points per continuous term.

- newdata:

  data.frame of prediction units. For `occu_cover`: one row per field
  cell (or carrying a `cell` column mapping rows to field cells),
  defaulting to the training data. For the families predicting from a
  design matrix
  ([`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
  [`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md),
  [`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md),
  [`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
  [`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md),
  [`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md),
  [`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md),
  [`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md)):
  the covariates of the arm `type` selects, expanded through that arm's
  fitted formula (factor levels and transforms as fitted) in place of
  `X.0` / `X.det.0`; an
  [`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  detection prediction expands it once per source. A random-effect term
  on the predicted arm whose grouping variable is a column of `newdata`
  adds each row's fitted group effect (its BLUP / posterior mean, times
  the slope covariate for a random slope); a level the fit did not see
  adds 0, the population mean. Without that column, and always through
  `X.0` / `X.det.0`, the prediction is population-level. Give `newdata`
  or `X.0` / `X.det.0` / `terms`, not both.

- times:

  `occu_cover` `type = "change"` / `"trajectory"`: numeric values of the
  time covariate. For `"change"`, at least two: `c(t1, t2)` differences
  one against the other, `c(t1, ..., tK)` differences every step against
  `t1`. For `"trajectory"`, one or more times to predict at, any values
  of the covariate (between or beyond the fitted times; beyond them the
  prediction follows the fitted time effect).

- level:

  `occu_cover` only: credible level for the interval columns (default
  0.95).

- nsim:

  `occu_cover` only: number of joint posterior draws (default 1000), or
  `"auto"` (`type = "trajectory"` on a joint nested-Laplace fit) to
  choose the count by Monte Carlo precision: draws are taken until the
  Monte Carlo standard error of every reported median and interval
  bound, over every unit, time and quantity, is at most `mc.tol` times
  that quantity's posterior SD, or `nsim.max` draws are reached. Each
  round redraws the whole table at the count the \\1/\sqrt{n}\\ scaling
  projects, so memory is that of one run at the final count and all
  quantities share one draw set.

- draws:

  `occu_cover` only: if `TRUE` (default), carry the per-unit
  `[cell x nsim]` draw matrices in `attr(, "draws")`.

- time.col:

  `occu_cover` only: name of the time covariate weighting the trend
  field / driving the change map; auto-resolved from the fit's stored
  trend weight when omitted.

- X.det.0:

  Optional detection design for out-of-sample `"detection"` / `"both"`
  prediction on a `single`-season or
  [`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  fit: a plain design matrix for a single-season fit, or a list of one
  design matrix per source (named by source, or in source order) for an
  [`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  fit – matching `fitted()$p`'s per-source shape. In-sample locations
  only: like the state-arm
  [`svc()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)/[`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  surfaces, this does not krige a spatial detection field to unseen
  points.

- weights:

  `occu_cover` (and the joint
  [`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) SVC
  route) only: pool several `newdata` rows into one prediction per cell.
  `newdata` then carries more than one row per `cell` (for example one
  row per habitat level), and `weights` gives each row's share of its
  cell, as a numeric vector with one entry per row or the name of a
  `newdata` column. Weights are normalised within each cell. Every
  quantity is pooled per posterior draw before it is summarised or
  differenced: occupancy, detection and expected cover are the weighted
  means of the rows, and conditional cover is the pooled expected cover
  divided by the pooled occupancy (the weighted mean of the rows'
  conditional cover when occupancy does not vary across them). The
  result has one row per cell, in order of first appearance, and a
  change map's occupancy and cover columns share one draw set. `NULL`
  (default) predicts each row on its own.

- aggregate:

  `type = "trajectory"` only: `FALSE` (default) returns one row per
  prediction unit and time; `TRUE` averages each quantity over the units
  per posterior draw first and returns one row per time (the overall
  trend, with its own interval).

- mc.tol:

  `nsim = "auto"` only: the Monte Carlo standard error allowed on each
  reported median and bound, as a fraction of its posterior SD (default
  0.05, about 2,900 draws for a Gaussian-shaped quantity at
  `level = 0.95`).

- nsim.max:

  `nsim = "auto"` only: the most draws to take (default 10000).

- mc.floor:

  Trajectory only: the smallest posterior SD a Monte Carlo error is
  measured against, on the quantity's own scale (default 0.001). Where
  occupancy saturates near 0 or 1 its posterior SD goes to about 0, and
  an error relative to it would be unstable; the error is measured
  against `max(SD, mc.floor)` instead. `psi`, `cover_cond` and
  `cover_exp` are probabilities or cover fractions, so one floor serves
  all three.

- ...:

  Ignored.

## Value

Depends on mode. In-sample:
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) result.
`"state"`: a data.frame with `row`, `psi` (marginalised posterior mean),
`psi_lower` / `psi_upper` (equal-tailed 95% credible interval; `NA` when
the engine did not return per-cell predictive variance), and `heldout`.
Design-matrix/ terms: data.frame with estimate and CIs;
`type = "detection"` on an
[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
fit returns a named list of such data.frames (one per source), and
`type = "both"` returns `list(occupancy = , detection = )`.
`occu_cover`: a `tobs_prediction` table (one row per cell) with per-unit
draw matrices in `attr(, "draws")`; `type = "trajectory"` a long table
(`cell`, `time`, then `<q>_mean` / `_median` / `_lwr` / `_upr` for
`psi`, `cover_cond` and `cover_exp`), one row per cell and time, or per
time with `aggregate = TRUE`. A trajectory states its Monte Carlo
precision in `attr(, "mc_se_max")` (the largest MC standard error of any
reported median or bound, in posterior SDs) and `attr(, "nsim_used")`;
under `nsim = "auto"` also `attr(, "mc_tol")` and
`attr(, "nsim_capped")` (`TRUE` when `nsim.max` stopped the draws before
`mc.tol` was met). `attr(, "mc_binding")` names the bound that sets
`mc_se_max`: a one-row data frame with its `cell` (`NA` with
`aggregate = TRUE`), `time`, `quantity`, `bound` (`"median"`, `"lwr"` or
`"upr"`) and `mc_se`.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
pr <- predict(fit, type = "occupancy")
summary(pr$psi)
# }
```
