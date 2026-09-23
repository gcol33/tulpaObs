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
  with `times = c(t1, t2)` returns a per-cell change table (`delta_p`,
  `delta_cover_cond`, `delta_cover_exp`, the occupancy / abundance
  decomposition, and `.lwr` / `.upr` at `level`), plus the start / end
  occupancy `p_T1` / `p_T2` with their own `.sd` / `.lwr` / `.upr`, and
  a `.prob_pos` column per headline delta giving the directional
  posterior probability `P(delta > 0)` per cell. Passing more than two
  times, `times = c(t1, ..., tK)`, widens the same table into a
  trajectory: a level column per step (`p_T1..p_TK`, `cover_cond_T*`,
  `cover_exp_T*`) and a `_T<k>`-suffixed delta per step, each
  differenced against `t1` and each carrying the same decomposition and
  interval columns. Every step is evaluated on ONE draw set, so the
  steps share a posterior and their deltas are jointly valid; with two
  times the delta columns keep their unsuffixed names, there being only
  one step to name. The result is a `tobs_prediction` table (one row per
  cell) carrying per-unit `[cell x nsim]` draw matrices in
  `attr(, "draws")`; map it yourself, e.g.
  `left_join(cents, pr, by = "cell")` then
  `geom_tile(aes(x, y, fill = delta_p))` (or `geom_sf()` on polygon
  cells).

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
  n_points = 50L,
  newdata = NULL,
  times = NULL,
  level = 0.95,
  nsim = 1000L,
  draws = TRUE,
  time_col = NULL,
  X_det.0 = NULL,
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
  an `occu_cover` fit: `"occurrence"`, `"cover_cond"`, `"cover_exp"`, or
  `"change"`. `"detection"` / `"both"` on a
  `single`/[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
  fit needs `X_det.0`.

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

- n_points:

  Number of prediction points per continuous term.

- newdata:

  `occu_cover` only: data.frame of prediction units, one row per field
  cell (or carrying a `cell` column mapping rows to field cells).
  Defaults to the training data.

- times:

  `occu_cover` `type = "change"` only: numeric values of the time
  covariate, at least two. `c(t1, t2)` differences one against the
  other; `c(t1, ..., tK)` returns a trajectory, every step differenced
  against `t1`.

- level:

  `occu_cover` only: credible level for the interval columns (default
  0.95).

- nsim:

  `occu_cover` only: number of joint posterior draws (default 1000).

- draws:

  `occu_cover` only: if `TRUE` (default), carry the per-unit
  `[cell x nsim]` draw matrices in `attr(, "draws")`.

- time_col:

  `occu_cover` only: name of the time covariate weighting the trend
  field / driving the change map; auto-resolved from the fit's stored
  trend weight when omitted.

- X_det.0:

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
draw matrices in `attr(, "draws")`.
