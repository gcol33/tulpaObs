# Simulate joint occupancy-detection + cover data

Per-cell mixture: latent z_i ~ Bernoulli(psi_i), per-visit detection
y_ij \| z_i = 1 ~ Bernoulli(p_ij), per-visit cover y_pos_ij \| y_ij = 1
drawn from `positive` (`"beta"` or `"lognormal"`) on the cover-arm
linear predictor. Used by the recovery test and as the generator for the
joint occupancy + cover hurdle family (see
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)).

## Usage

``` r
simulate_occu_cover(
  N = 200L,
  J = 4L,
  n_occ_covs = 1L,
  n_det_covs = 1L,
  n_pos_covs = 1L,
  beta_occ = NULL,
  beta_p = NULL,
  beta_pos = NULL,
  positive = c("lognormal", "beta", "gaussian"),
  phi = 30,
  sigma_pos = 0.4,
  adj = NULL,
  sigma = 0.6,
  alpha = 1,
  trend = FALSE,
  sigma_trend = 0.6,
  alpha_trend = 1,
  pos_field = FALSE,
  sigma_pos_int = 0.5,
  sigma_pos_trend = 0.6,
  det_field = FALSE,
  sigma_p_int = 0.5,
  sigma_p_trend = 0.6,
  re_det_groups = NULL,
  sigma_re_p = 0.7,
  re_pos_groups = NULL,
  sigma_re_pos = 0.7,
  re_det = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (cells).

- J:

  Number of visits per site.

- n_occ_covs, n_det_covs, n_pos_covs:

  Number of covariates on each arm (drawn IID standard normal).

- beta_occ, beta_p, beta_pos:

  Coefficient vectors c(intercept, slopes). Defaults pick
  weakly-informative values: psi intercept at logit(0.4), p intercept at
  logit(0.5), cover intercept on the appropriate link.

- positive:

  `"beta"` or `"lognormal"`.

- phi:

  Beta precision when `positive = "beta"` (default 30).

- sigma_pos:

  Lognormal residual SD when `positive = "lognormal"` (default 0.4).

- adj:

  Optional N x N adjacency matrix. When supplied, generates the shared
  ICAR field; when NULL, the simulator is non-spatial (matches v1).

- sigma:

  Spatial field amplitude (used only when `adj` is supplied), in the
  geo-mean-marginal-SD (Sorbye-Rue) convention: the field is drawn with
  geometric-mean marginal variance 1 before scaling by `sigma`, so
  `sigma` is directly comparable to a fitted `field_sd`
  (`fit$spatial$field_sd_mean` on the `nested_laplace` path,
  `fit$hyper_draws[, "field_sd"]` on `nuts`) – never to the raw `sigma`
  a `nested_laplace` fit reports, which is the amplitude against the
  unscaled intrinsic precision and differs by `sqrt(scale_q)`.

- alpha:

  Cover-arm scaling on the shared field (used only when `adj` is
  supplied). 1.0 = arms see the field identically; positive = same sign,
  negative = opposite.

- trend:

  Logical; when `TRUE` (and `adj` is supplied) a SECOND shared ICAR
  field `f2` (a spatially-varying temporal trend) is generated on the
  same graph, weighted by a per-cell covariate `time` drawn IID standard
  normal. The trend enters the occupancy and cover predictors as
  `sigma_trend * time_i * f2[i]` (occupancy) and
  `alpha_trend * sigma_trend * time_i * f2[i]` (cover); the detection
  predictor is unaffected. The `time` covariate is per-cell and
  broadcast to every visit of that cell.

- sigma_trend:

  Trend-field amplitude (used only when `trend = TRUE`).

- alpha_trend:

  Cover-arm scaling on the trend field (used only when `trend = TRUE`).

- pos_field:

  Logical; when `TRUE` (and `adj` is supplied) draw an INDEPENDENT areal
  field on the cover (positive) arm only – an intercept field plus a
  time-weighted trend field, each unrelated to the occupancy field and
  with no alpha copy. Adds a `time` column to the returned `data` and
  reports `g0` / `g1` (the two fields) and their SDs in `truth`. Fit by
  placing `spatial(~ 1 + time || cell, graph = adj)` in the `positive`
  formula.

- sigma_pos_int:

  Cover-arm intercept-field SD (used only when `pos_field = TRUE`).

- sigma_pos_trend:

  Cover-arm trend-field SD (used only when `pos_field = TRUE`).

- det_field:

  Logical; when `TRUE` (and `adj` is supplied) draw an INDEPENDENT areal
  field on the detection arm only – an intercept field plus a
  time-weighted trend field, unrelated to the occupancy and cover fields
  and with no alpha copy, so detection varies spatially on its own. Fit
  by placing
  `spatial(~ 0 + time || cell, graph = adj, to = "detection")` in the
  `detection` formula.

- sigma_p_int:

  Detection-arm intercept-field SD (used only when `det_field = TRUE`).

- sigma_p_trend:

  Detection-arm trend-field SD (used only when `det_field = TRUE`).

- re_det_groups:

  Optional integer `>= 2`: the number of levels of a per-visit detection
  random intercept (a `habitat` factor on `visit_data`, levels
  `hab1..K`), drawn `b_g ~ N(0, sigma_re_p^2)` and centred. `NULL`
  (default) adds no detection RE. Recover it with
  `detection = ~ ... + (1 | habitat)`.

- sigma_re_p:

  SD of the `re_det_groups` random intercept (default 0.7).

- re_pos_groups:

  Optional integer `>= 2`: the number of levels of a per-visit COVER-arm
  random intercept (a `habitat` factor on `visit_data`, levels
  `hab1..K`), drawn `b_g ~ N(0, sigma_re_pos^2)` and centred, added to
  the positive-cover linear predictor. `NULL` (default) adds no cover
  RE. Recover it with `positive = ~ ... + (1 | habitat)` under
  `cover_aggregate = "none"`. Truth in `truth$b_pos_re` /
  `truth$sigma_re_pos`.

- sigma_re_pos:

  SD of the `re_pos_groups` random intercept (default 0.7).

- re_det:

  Optional named list of FURTHER per-visit detection random effects, for
  crossed / nested / slope designs. Each element
  `list(K =, sigma =, prefix =, nested_in =, slope_cov =, sigma_slope =, rho =)`
  adds a factor column (levels `<prefix>1..K`). Without `slope_cov` it
  is a centred `N(0, sigma^2)` random intercept; `nested_in = "<name>"`
  nests its codes within a previously listed grouping (matching
  `(1 | parent/child)`), otherwise crossed. With
  `slope_cov = "<column>"` (a per-visit covariate, generated
  `N(0, slope_sd)` if absent, `slope_sd` default 1) it is a random
  slope: a slope-only uncorrelated block when `rho` is unset, or a
  correlated intercept + slope block (covariance from `sigma` /
  `sigma_slope` / `rho`) when `rho` is given. Truth is returned in
  `truth$re_det[[name]]` (named by the level label a fit reconstructs):
  `b` / `b_slope` BLUP vectors, or the `B` BLUP matrix plus `s0` / `s1`
  / `rho` for a correlated slope.

- seed:

  Optional integer seed.

## Value

A list with `y` (N x J detection matrix), `y_pos` (N x J cover matrix,
NA where not detected), `data` (per-site covariate frame, gaining a
`time` column when `trend = TRUE`), `visit_data` (per-visit covariate
frame, N\*J rows in site-major order), and `truth` (the coefficients,
dispersion, and field(s) if generated; `f2`, `sigma_trend`,
`alpha_trend`, and `time` when `trend = TRUE`).

## Details

When `adj` is supplied (a square adjacency matrix), an ICAR field
`f[1..N]` is drawn from `MVN(0, Q^-)` (with sum-to-zero constraint), and
the linear predictors become

    eta_psi_i = X_psi[i, ] %*% beta_occ + sigma * f[i]
    eta_pos_ij = X_pos[i, ] %*% beta_pos + alpha * sigma * f[i]

matching the v2 nested-Laplace fit's parameterisation.

## Examples

``` r
sim <- simulate_occu_cover(N = 50, J = 3, seed = 1)
dim(sim$y)
```
