# Cover hurdle family (vegetation cover, MOTIVATE pattern)

Latent presence (Bernoulli) plus conditional positive cover (beta or
lognormal). Does not share the replicate-detection assumption of the
other families – see `vignette("families")` for the conceptual caveat.

## Usage

``` r
cover(
  response = c("beta", "beta_oi", "lognormal", "lognormal_trunc", "ordinal", "gaussian"),
  breaks = NULL
)
```

## Arguments

- response:

  likelihood for the positive cover part:

  - `"beta"` – cover in (0, 1) with a logit link;

  - `"beta_oi"` – a one-inflated Beta: a point mass at cover = 1 (plots
    recorded at full cover, a genuine boundary mass rather than a near-1
    continuous value) plus a Beta on the interior (0, 1). The inflation
    probability is a constant, estimated as the share of positive plots
    at the ceiling; the interior Beta is fit on the (0, 1) plots.
    Reported as `pi_one`. (With the zero hurdle this is the
    zero-one-inflated Beta.);

  - `"lognormal"` – a Gaussian on log-cover (unbounded above);

  - `"lognormal_trunc"` – a Gaussian on log-cover upper-truncated at
    `log(1) = 0`, the bounded-support form for cover in (0, 1\] (it
    cannot place mass above cover = 1 the way `"lognormal"` can).
    Requires `method = "nested_laplace"`;

  - `"ordinal"` – an interval-censored Gaussian on log-cover with KNOWN
    class thresholds (an ordered probit for Braun-Blanquet-style class
    data), set with `breaks`. Requires `method = "nested_laplace"`;

  - `"gaussian"` – an identity-Gaussian magnitude (the delta-normal
    hurdle): Bernoulli presence times a Gaussian on the raw response.
    For a pre-transformed or otherwise unbounded response on a real
    scale (it permits negative fitted values), NOT for cover fractions
    in (0, 1); those stay on `"beta"` / `"lognormal"` / `"ordinal"`.
    Presence is the nonzero sentinel (`y != 0`) rather than `y > 0`.

- breaks:

  for `positive = "ordinal"` only, the interior cover-class boundaries
  on the (0, 1) cover-fraction scale (strictly ascending, all in (0,
  1)); the open outer classes are added automatically. `NULL` for the
  other families.

## Value

A `tobs_family` object.

## Response on the formula left-hand side

The cover response is a single length-N vector, so it may sit on the top
formula left-hand side and `y =` is dropped. The two calls are
equivalent:

    # response on the LHS (y = omitted)
    tobs(cover.flat ~ time.sc + habitat +
           spatial(~ 1 + time.sc || cell_idx, graph = adj),
         data = dat, family = cover(response = "beta"),
         method = "nested_laplace")

    # the same fit with a one-sided formula and an explicit y =
    tobs(~ time.sc + habitat +
           spatial(~ 1 + time.sc || cell_idx, graph = adj),
         y = dat$cover.flat, data = dat, family = cover(response = "beta"),
         method = "nested_laplace")

Naming the response makes the per-arm spatial labels read naturally:
`cover()` splits `cover.flat` into a `presence` arm and a `positive`
arm, the arm names that per-arm formulas and share() address. The LHS is
evaluated against `data` (then the calling environment), so it may be a
bare column or an expression.

## Joint nested-Laplace engine – spatial-prior parameterisation

When fitted with `method = "nested_laplace"` and an areal spatial term
in the latent-presence formula (`bym2(graph = adj)`, or `car()` /
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)),
the engine identifies a single latent field `z` per region and
parameterises the two arms as `eta_occ = X beta_occ + sigma * z` and
`eta_pos = X beta_pos + alpha * sigma * z`. Both BYM2 sub-blocks (`phi`,
`theta`) are subject to hard sum-to-zero constraints (see
`.joint_inner_var()` for the math): the cover-arm intercept
`beta_pos[1]` is identified as the population mean of `eta_pos` under
`mean(z) = 0`. Any simulator generating data for this engine must demean
both `phi_f` and `theta_f` before scaling, otherwise the estimator
targets `beta_pos_0_truth + alpha * mean(w_s_sim)` and coverage of the
*population* truth collapses with alpha. See
[`simulate_cover_joint()`](https://gillescolling.com/tulpaObs/reference/simulate_cover_joint.md)
for a ready-made demeaned simulator.

## Spatially varying trend (a weighted areal term)

A spatially varying trend is model structure, so it is declared in the
formula as a SECOND, weighted areal term on the same graph as the
intercept field, not through `control`:

    ~ time.sc +
      icar(graph = adj, group_var = "cell_idx") +
      icar(graph = adj, weight = time.sc, group_var = "cell_idx")

The unweighted term is the shared intercept field; the weighted term
`icar(..., weight = col)` is the trend field, whose contribution to each
arm's predictor is `weight_i * z[cell_i]`, coupled onto the cover arm
with its own scale (`alpha_trend`, reported in `fit$alpha_trend` /
`fit$sigma_trend`) integrated over the outer grid. The umbrella spelling
`spatial(graph = adj, model = "icar", weight = col)` resolves
identically. The trend block inherits the intercept block's amplitude;
give it its own with
`share(spatial(), terms = list(intercept = ..., trend = ...))`. Requires
`method = "nested_laplace"`. A coupled trend cannot currently combine
with
[`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/ [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
blocks in the same fit.

The axis the amplitude rides carries prior structure – a point mass at
`alpha = 0` ("no coupling") and a log-spaced slab above it – so it is
set in one of two ways. `share(alpha = grid(...))` STATES its nodes, and
with them that structure. `share(alpha = grid(n = 9))` states a
RESOLUTION: the engine re-reads its own axis with that many slab nodes,
point mass and slab bounds unchanged, so sharpening the axis never
restates its structure. `terms =` gives either, per block. The
resolution is the only way to raise this axis, because it does not
densify when the donor `control$sigma.grid` does: on an informative data
set the outer grid's quadrature effective sample size saturates on the
copy amplitude while every other axis tracks the request (measured
engine-side, `NOTES_measurements.md`). One block takes one of the two;
giving both is an error. `control$alpha.grid[.trend]` and
`control$alpha.n[.trend]` are the wire format these compile into and are
refused as user input.

Set it when the fit's hyperparameter intervals are reported. On the
coupled SBC fixture the declared resolution leaves the field SD
miscalibrated once the data are informative – at 10 visits per site its
uniformity p-value is 9.1e-05, against 0.17 at 3 visits – and
`control$alpha.n = 21` returns every scored quantity to nominal at both.
Thirteen slab nodes do not (the field SD still reads 4.5e-03), and the
price is a 2.3-3x longer fit, since the outer grid is a tensor. Point
estimates are not affected the way the intervals are; the measurement is
in `NOTES_measurements.md`.

## Varying-coefficient spatial bar (the compact single-term form)

The intercept field plus its weighted trend field can also be written as
one
[`spatial()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term carrying an lme4-style coefficient formula:

    ~ time.sc + habitat +
      spatial(~ 1 + time.sc || cell_idx, graph = adj)

The bar left-hand side spells the coefficient fields: the intercept
column (`1`) is the unweighted field; each covariate column (`time.sc`)
is a weight-scaled coefficient field (`weight_i * z[cell_i]`). The bar
right-hand side (`cell_idx`) is the graph node index (the areal
`group_var`); `||` requests independent intercept and slope fields, a
single `|` makes them correlated. This desugars to exactly the two-term
weighted-areal form above, so the two spellings give the same fit.

## Choosing a field's arm

placement and share(): The cover hurdle's two arms are `presence` (the
`y > 0` Bernoulli arm) and `positive` (the `y | y > 0` arm);
[`summary()`](https://rdrr.io/r/base/summary.html) and the coefficient
output print these same labels. A field is placed on an arm by writing
it in that arm's per-arm formula (`presence = ~ ...`,
`positive = ~ ...`), and shared across arms with share(). A field in the
single shared `formula` reaches both arms.

*Both arms (shared / copied).* A field in the shared `formula` is one
presence-anchored latent copied onto the positive arm with an estimated
coupling per coefficient field (`alpha` for the intercept field,
`alpha_trend` for the trend field), marginalized on the outer grid:

    # one latent field, shared by both arms
    tobs(cover.flat ~ time.sc + habitat +
           spatial(~ 1 + time.sc || cell_idx, graph = adj),
         data = dat, family = cover(response = "beta"), method = "nested_laplace")

    presence: eta_presence = ... + u_cell + time.sc * s_cell
    positive: eta_positive = ... + alpha * u_cell + alpha_trend * time.sc * s_cell

Per-arm formulas make the shared field explicit: place it on `presence`
and copy it onto `positive`. `share(spatial())` estimates the coupling
on the default axis; `share(spatial(), alpha = grid(c(...)))` integrates
it over supplied nodes, `alpha = grid(n = 9)` over the engine's own axis
read at `n` slab nodes, and `alpha = 0.5` fixes it.
`share(spatial(), prior = list(...))` regularizes the coefficient
itself; a fit copying both the intercept and a weighted trend block is
refused it, since one prior reaches the engine per fit.

    tobs(presence = ~ time.sc + habitat +
                      spatial(~ 1 + time.sc || cell_idx, graph = adj),
         positive = ~ time.sc + habitat + share(spatial()),
         data = dat, family = cover(response = "beta"), method = "nested_laplace")

*One arm only (free / separate).* A field written in one arm's formula
is a separate latent on that arm alone, with its own precision and no
cross-arm coupling. A field in each arm's formula gives two independent
latents:

    tobs(presence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
         positive = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
         data = dat, family = cover(response = "beta"), method = "nested_laplace")

    presence: eta_presence = ... + u_presence_cell + time.sc * s_presence_cell
    positive: eta_positive = ... + u_positive_cell + time.sc * s_positive_cell

The free fit reports `sigma_armspecific`. Arm-specific fields are their
own spatial structure: they do not combine with a shared field, a
weighted trend, or
[`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/ [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
in the same formula, and at most one targets each arm.

The `||` and `|` axis is separate from shared / free: `||` makes the
intercept and slope fields independent, while a single `|` makes them
correlated (a free cross-covariance, MCAR). A correlated `|` bar shared
across both arms (in the shared `formula`, or on `presence` with
`share()` on `positive`) is copied with one estimated amplitude; placed
on one arm's formula it is a free-Sigma correlated field on that arm
alone, no cross-arm copy.

    field in the shared formula, no share()   one presence-arm latent, pinned at alpha = 0
    field + share() (either spelling)         one shared / copied latent (presence anchor, coupling estimated)
    a field in each arm's formula             separate / free latents, no coupling
    ||                                        independent intercept and slope coefficient fields
    |                                         correlated (MCAR) coefficient fields, copy-only

Without a `share()` (and without `control$alpha.grid`) the presence
field is NOT carried onto the cover arm: the amplitude is pinned at
`alpha = 0` and the field rides the presence arm alone. Structure nobody
wrote is not in the model, the reading
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
and
[`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md)
also have. Write `+ share(spatial())` in the shared formula, or place
the field on `presence` and `share(spatial())` on `positive`, to couple
the two arms. A correlated `|` bar is the exception: that spelling is
copy-only by definition, so the bar itself states the coupling.

## A formula bar is a random effect, not a spatial field

A bare lme4 bar in the formula – `(1 | cell)`, `(1 + x | cell)`,
`(x || cell)` – is a grouped random effect, not a spatial field, even
when the grouping factor names the areal graph nodes. The engine's
inline-MCAR call `tulpa::spatial(graph, ~ 1 + x | cell)` reads `| cell`
as a separable spatial field, but the same spelling inside a `cover()`
formula is parsed as an IID random effect on `cell`. For a spatial field
on the cells write either the compact
[`spatial()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
bar of the section above, `spatial(~ 1 + x || cell, graph = adj)`, or
the two-term weighted-areal form,
`icar(graph = adj, group_var = "cell") +`
`icar(graph = adj, weight = x, group_var = "cell")`. When a formula
carries a bar whose grouping factor is also an areal term's `group_var`,
`cover()` emits a one-time message noting the bar is fitted as a random
effect; suppress it with
[`base::suppressMessages()`](https://rdrr.io/r/base/message.html).

## Checkpoint / resume

A spatial cover-hurdle fit integrates over a large outer hyperparameter
grid and can run for hours.
`control$checkpoint = list(path = "fit.ckpt", resume = TRUE)` makes the
joint engine append each completed grid cell to `path`; a
`resume = TRUE` run loads the finished cells and solves only the rest,
reproducing the from-scratch fit, so a killed or rebooted fit resumes
instead of restarting. `resume = FALSE` starts a fresh file. Forwarded
to
[`tulpa::tulpa_nested_laplace_joint()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace_joint.html).

## Examples

``` r
# \donttest{
sim <- simulate_cover(N = 200, seed = 1)
fit <- tobs(~ x, data = sim$data, family = cover("beta"), y = sim$y)
summary(fit)
# }
```
