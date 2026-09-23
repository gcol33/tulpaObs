# Fit cover_hurdle as a joint binomial+(gaussian\|beta) model with shared spatial field via [`tulpa::tulpa_nested_laplace_joint()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace_joint.html).

For both positive parts the dispersion scalar is integrated on the outer
joint hyperparameter grid (per-arm `phi_pos` axis):

## Usage

``` r
fit_cover_hurdle_joint_nested(
  enc,
  data,
  positive = enc$positive,
  control = list(),
  temporal = NULL,
  re = NULL,
  priors = NULL
)
```

## Arguments

- enc:

  Output of
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md).

- data:

  The original (un-subsetted) data frame – required to resolve the
  spatial spec (group_var lookup, n_spatial_units check).

- positive:

  `"lognormal"` or `"beta"`.

- control:

  List with optional `max_iter`, `tol`, `n_threads`, `sigma_grid`,
  `rho_grid`, `rho_car_grid`, `alpha_grid`, `phi_init`, `phi_bounds`
  (the last two are forwarded to the beta pre-fit when
  `positive = "beta"`). For ICAR / CAR_proper backends `tau_grid` is
  also accepted and translated to `sigma_grid` as
  `sigma = 1 / sqrt(tau)`. The cover arm sees the shared field at
  amplitude `alpha * sigma`, so `alpha_grid` is the cover-arm amplitude
  axis and `sigma_grid` the donor's. Regularizing hyperpriors on the
  joint (sigma, alpha) axes can be set via `prior_sigma` (donor
  amplitude) and `prior_alpha` (copy coefficient) – each a length-2 list
  `list(family, params)` matching tulpa's `prior_sigma` / `prior_alpha`
  args. The prior on alpha directly regularizes the copy scalar at small
  `n_pos`, replacing the per-arm `prior_sigma_pos` of the pre-reparam
  API. `prior.phi` puts the same kind of regularizing hyperprior on the
  cover-arm dispersion grid (the beta precision under
  `positive = "beta"`, the log-scale SD under `lognormal`), re-weighting
  the `phi.grid` axis by the chosen density in place of the engine's
  default prior on it; same `list(family, params)` form, forwarded to
  tulpa's `prior_phi`.

- temporal, re:

  Structured
  [`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  / [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  blocks from the formula, stacked onto the shared spatial block via the
  multi-block joint engine.

- priors:

  Optional
  [`cover_priors()`](https://gillescolling.com/tulpaObs/reference/cover_priors.md)
  object (or coercible list). When supplied, the per-arm fixed-effect
  prior reaches the joint engine as a `beta_prior_mean` /
  `beta_prior_prec` on the occurrence and positive responses, mirroring
  the separate-Laplace path. `NULL` / `FALSE` / `"none"` leave both arms
  unpenalised.

## Value

List shaped like the single-Laplace fit output but with extra `joint`
field carrying the raw `tulpa_nested_laplace_joint` result.

## Details

- `positive = "lognormal"` (gaussian arm): the residual SD is the
  per-grid phi. The default 7-point log-spaced grid is centred on the
  non-spatial prefit from `.prefit_lognormal_sigma()` and spans
  `[sigma_hat / 3, sigma_hat * 3]`. The posterior mean and SD across
  that axis are surfaced as `sigma_pos` / `sigma_pos_sd` on the returned
  `cover_fit`.

- `positive = "beta"`: the beta precision is the per-grid phi. The
  default 7-point log-spaced grid spans `[2, 300]`; posterior mean and
  SD are surfaced as `phi_pos` / `phi_pos_sd`.

Override the per-arm phi grid via `control$phi.grid`.
