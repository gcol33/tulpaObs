# Joint occupancy-detection + cover hurdle family

Cell-level latent presence `psi`, per-visit binomial detection `p`, and
per-visit positive cover `f_pos` (beta or lognormal) on a third linear
predictor. The N-mixture analogue for vegetation cover: where
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) couples
occupancy to counts via Royle's binomial-N marginal, `occu_cover()`
couples it to a positive-cover observation at each detected visit. The
latent presence z marginalises out in closed form (two states), so the
marginal log-likelihood is exact.

## Usage

``` r
occu_cover(
  response = c("beta", "lognormal", "gaussian"),
  cover_aggregate = NULL
)
```

## Arguments

- response:

  likelihood for the positive cover arm. `"beta"` (cover in (0, 1)),
  `"lognormal"` (log-cover Gaussian), or `"gaussian"` (an
  identity-Gaussian magnitude on a real, unbounded scale – the
  delta-normal hurdle; not for cover fractions, which stay on
  `"beta"`/`"lognormal"`).

- cover_aggregate:

  how the cover arm collapses per occupancy unit on the shared-field
  spatial path: `"mean"`, `"median"`, `"latent"` (a per-unit cover
  random effect integrated out), or `"none"` (per-visit). `NULL`
  (default) is `"mean"` on the spatial path and `"none"` on the
  non-spatial path. See the *Cell-aggregated cover* section.

## Value

A `tobs_family` object.

## Details

Per-cell likelihood:

    any_det_i : L_i = psi_i * prod_j h_ij
    no_det_i  : L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)

    h_ij = (1 - p_ij) * 1{y_ij = 0}
         + p_ij       * f_pos(y_pos_ij; eta_pos_ij, ...) * 1{y_ij = 1}

A detected visit (`y_ij = 1`) with a missing cover (`y_pos_ij = NA`)
keeps its detection term but drops the `f_pos` factor: cover is taken
missing-at-random, so the cover likelihood runs over the detected visits
with an observed cover.

Reduces to
[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) when
the cover arm is degenerate, and to the plot-level cover hurdle
([`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md))
when J = 1 and detection is perfect.

## Engines

The non-spatial fit is a direct Laplace approximation
(`method = "laplace"`) or a NUTS sampler over the same exact two-state
marginal (`method = "nuts"`, beta or lognormal cover), giving calibrated
intervals and a per-draw pointwise likelihood for WAIC / LOO. A shared
areal field across the occupancy and cover arms is the
`method = "nested_laplace"` path (a structured
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term on the psi formula); a structured term on a non-spatial route
errors from the dispatcher with a pointer to it.

`method = "nuts"` also samples the coupled areal field(s)
([`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/ [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)),
written either as areal terms or as the bar
`spatial(~ 1 || cell, graph = adj)`, TOGETHER WITH their
hyperparameters: the field SD, the mixing (`bym2`) or
spatial-correlation (`car_proper`) parameter, and - where the formula
asks for one - the cover-arm copy amplitude. Each is bounded to the span
the `nested_laplace` path's outer quadrature integrates that axis over,
in that grid's own coordinate, and the same `control$sigma.grid` /
`alpha.grid` / `rho.car.grid` knobs set it on both routes - so the
sampler is an independent reference for the hyperparameter layer rather
than a fit conditioned on that layer's point estimate. The measure over
that span is the one the outer grid declares for the axis: flat in the
grid's coordinate for the field SD and the mixing parameter, and for the
copy amplitude the penalized-complexity Exponential slab, with
`control$copy.slab = "flat"` asking for the flat alternative the grid
also accepts. The grid weighs that slab against a point mass at
`alpha = 0` ("no coupling"), which a gradient sampler cannot visit, so
the sampled copy amplitude is the posterior conditional on a coupled
field; the point mass's share of the prior is read on the
`nested_laplace` route. `fit$nuts$sampled_hyper` and
`fit$nuts$fixed_hyper` name which is which per fit (`icar` pins `rho` at
1: the intrinsic precision has no mixing parameter; a field with no
`share()` pins the copy amplitude at 0), `fit$hyper_draws` carries their
posterior alongside `field_sd`, the geometric-mean marginal SD the field
implies. `control$fixed.hyper = TRUE` conditions on the warm
nested-Laplace estimate instead.

The grid-integrated `nested_laplace` path reports `field_sd` too
(`fit$spatial$field_sd_mean` / `field_sd_sd`, alongside `sigma_mean`),
in the SAME geo-mean-marginal-SD convention as the NUTS path's
`field_sd`. It is a fixed multiple of `sigma` (not an independent grid
axis), so it is NOT in `fit$means` / `fit$vcov` – folding it in there
would make the joint parameter vcov exactly singular. Do not compare
fits, or a fit against
[`simulate_occu_cover()`](https://gillescolling.com/tulpaObs/reference/simulate_occu_cover.md)'s
`sigma`, by reading the raw `sigma` on this path: it is the field's
amplitude against the unscaled intrinsic ICAR precision `Q = D - W`,
which differs from `field_sd` by `sqrt(scale_q)`, a graph-size-dependent
factor (about 2.1 for a 30-node chain graph). `field_sd` is the number
comparable across fits and to a simulation truth; `sigma` is the raw
amplitude the engine's grid axis is spelled in.

A spatially-varying coefficient - a *weighted* areal term beside the
unweighted intercept field, or a bar column
`spatial(~ 1 + year || cell, graph = adj)` - is sampled as a SECOND
field block, with its own whitened surface and its own (SD, mixing, copy
amplitude) coordinates: two fields share no hyperparameter. The surfaces
are `fit$spatial_field` and `fit$trend_field` / `fit$trend_fields`
(named by the weight column), and each block's hypers carry that block's
suffix in `fit$hyper_draws` and in `fit$nuts$sampled_hyper`
(`sigma_trend`, `alpha_trend`, indexed when there are several). Each
field's amplitude is addressed on its own through
`share(spatial(), terms = list(...))` – with stated nodes, or
`grid(n = )` for the engine's axis read more finely – or through
`control$alpha.grid[.trend]` / `control$alpha.n[.trend]`, the same knobs
the grid-integrated route reads. A second field multiplies the warm
fit's outer grid, so a defaulted axis is thinned to three nodes over the
same span - the sampled prior is unchanged, the warm fit stays
affordable. A correlated bar (`|`, one free-Sigma MCAR block across the
fields) stays on the grid-integrated `nested_laplace` path; the
community spatial-factor route is
[`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md).

On the shared-field `nested_laplace` path the field-coupled occupancy
slope Wald interval is mildly anti-conservative at small N (the
grid-integrated Laplace under-disperses that coefficient; the pooled
coverage across the three arms stays near nominal). The non-spatial
`nuts` path gives fully calibrated occupancy intervals.

## What a bare areal term means for the cover arm

A term's process is the formula it sits in. An areal term written in the
occurrence formula therefore loads on the OCCURRENCE arm alone: the
cover arm sees no field, and the copy amplitude `alpha` is 0. This holds
under `method = "nested_laplace"` and `method = "nuts"` alike, so the
same input fits the same model on either engine.

To carry the occurrence field onto the cover arm as well, write the copy
in the `positive` formula:

    positive = ~ pos_cov + share(spatial())                      # amplitude estimated
    positive = ~ pos_cov + share(spatial(), alpha = grid(c(...))) # over given nodes
    positive = ~ pos_cov + share(spatial(), alpha = grid(n = 9))  # engine axis, finer
    positive = ~ pos_cov + share(spatial(), alpha = 0.5)          # fixed
    positive = ~ pos_cov + share(spatial(), prior = list("pc.prec", c(4, 0.01)))

The amplitude is integrated over that axis on the `nested_laplace` path
and sampled over its span on the `nuts` path; a scalar `alpha =` pins it
on both. `prior =` regularizes the coefficient itself, in the joint
driver's `list(<family>, <params>)` shape (see `prior_alpha` in
[`tulpa::tulpa_nested_laplace_joint()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace_joint.html));
one reaches the engine per fit, so a fit copying both an intercept and a
trend block is refused it rather than given it on the first.

`control$alpha.grid[.trend]` and `control$alpha.n[.trend]` are the WIRE
FORMAT these compile into, not user surface: a fit that sets one is
refused and told which formula form to write. `control$prior.alpha` is
still the lower-level spelling of `prior =`; set it in one place, not
both.

Without a `share()` the occurrence field is NOT carried onto the cover
arm – there is no implicit coupling – so the amplitude is not a
parameter of the model, and neither engine reports one. Under
`nested_laplace` `alpha` is absent from
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and the `n_params` count;
under `nuts` it is absent from `fit$nuts$sampled_hyper`,
`fit$nuts$fixed_hyper` and the `fit$hyper_draws` columns.
[`predict()`](https://rdrr.io/r/stats/predict.html) returns a cover arm
the occurrence field does not enter. Add `share(spatial())` to couple
the two arms.

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

## Coupled fields and spatially-varying trends

The spatial engines (`method = "nested_laplace"`, default
`control$engine = "joint"`, and `method = "nuts"`) share one areal field
(the cell intercept, an unweighted
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term) across the occupancy and cover arms, when the `positive` formula
copies it. ADDITIONAL coupled fields - spatially-varying coefficients,
e.g. a temporal trend - are added as *weighted* areal terms in the
formula:

    ~ elev + icar(graph = adj) + icar(graph = adj, weight = year)

Each weighted term `icar(graph, weight = col)` is a second field whose
contribution to a predictor row is `weight_i * z[cell_i]`. All such
fields share the same graph; a `share()` carries one onto the cover arm
with its own scale (`alpha` for the intercept field, `alpha_trend` for a
trend field), integrated over the outer grid, and a field the `positive`
formula does not copy stays on occurrence.
`share(spatial(), terms = list(...))` gives a per-block amplitude and
must address every block. The intercept field is reported in
`fit$spatial_field`, trend fields in `fit$trend_field` /
`fit$trend_fields`. The trend coupling grid defaults to the intercept
block's amplitude; give the trend its own with
`share(spatial(), terms = list(intercept = ..., trend = ...))`. Under
`method = "nuts"` each field is its own sampled block, with the scales
sampled over those same spans rather than integrated over the grid.

A single trend field may also be requested out-of-band with
`control = list(trend = list(weight = "<col>"))`, naming a numeric
per-cell covariate in the cell `data`; this is the equivalent of one
weighted formula term. Specify the trend field one way or the other, not
both.

## Sites larger than cells (`group_var`)

By default each site (one row of `y` / `data`, one latent occupancy
state) is its own field node, so the graph must have one node per site.
Passing `group_var = "<col>"` to the
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term maps each site to a field node named by that integer column, so
several sites can share one node. The field then stays length `n_cells`
(the graph) while occupancy, detection, and cover run over `n_sites`.
The motivating layout is a site = cell x time-period: plots in a
cell-period are the detection replicates, occupancy is per cell-period,
and a per-site time weight on a coupled trend field
(`icar(graph, weight = time, group_var = "cell")`) gives a
detection-corrected occupancy trend on a shared cell field. Supported on
the default `joint` engine; the `v2_joint` / `v3_nested` escape hatches
bind the field 1:1 to sites and reject `group_var`.

## Per-group random intercept on the shared-field path ([`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) / `(1 | g)`)

On the spatial `nested_laplace` path a single random INTERCEPT term on
the psi formula – `re(g)` or the `(1 | g)` bar – layers a per-group
occupancy offset on top of the shared field. It joins the joint fit as
one `iid` prior block whose latent rides the occupancy arm only; its
variance integrates on the outer grid alongside the field sigma / alpha
and is reported as the `sigma_re` hyperparameter, with the per-group
BLUPs in `fit$re` and via
[`ranef()`](https://gillescolling.com/tulpa/reference/ranef.html). A
default `re.sigma.grid` (log-spaced) sets the grid; pass
`control$re.sigma.grid` to override. Scope: one random-intercept term (a
scalar per group maps onto the one iid block); a random slope, a
correlated multi-coefficient block, or an RE without a shared field is
rejected (the joint engine integrates every variance component on its
grid, so multiple variance components do not scale – richer RE is the
non-spatial cover-hurdle EM's route). This is the single-species
analogue of the community spatial occu_cover
([`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md))
and the tulpaObs consumer of tulpa's field + per-group RE engine
composition.

A bar is therefore always a random effect, never a spatial field, even
when the grouping factor names the areal graph nodes. The engine's
inline-MCAR call `tulpa::spatial(graph, ~ 1 + x | cell)` reads `| cell`
as a separable spatial field, but the same spelling in an `occu_cover()`
formula is an IID random effect on `cell`. For a spatial field on the
cells use an areal term – `icar(graph = adj, group_var = "cell")` (plus
`icar(graph = adj, weight = x, group_var = "cell")` for a
spatially-varying trend) – not a bar. When a formula carries a bar whose
grouping factor is also an areal term's `group_var`, `occu_cover()`
emits a one-time message noting the bar is fitted as a random effect;
suppress it with
[`base::suppressMessages()`](https://rdrr.io/r/base/message.html).

## Independent field on the cover arm (placement)

A `share(spatial())` in the `positive` formula shares the occupancy
field: the cover arm sees it as `alpha * (occupancy field)`, the
coregionalization copy on the outer `(sigma, alpha)` grid. When the
cover trend is spatially structured but is not a scalar multiple of the
occupancy field, `alpha` collapses toward 0 and the cover arm inherits
no field, so per-cell conditional cover (and its change over time,
`delta_cover_cond`) comes out flat. A
[`spatial()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
bar placed in the `positive` formula adds an INDEPENDENT, non-copied
areal field on the cover arm alone:

    tobs(occurrence = ~ x + icar(graph = adj, group_var = "cell"),
         detection  = ~ 1,
         positive   = ~ time + spatial(~ 1 + time || cell, graph = adj),
         data = cell_dat, y = y, y_pos = y_pos,
         family = occu_cover("lognormal"), method = "nested_laplace")

The occupancy field still drives psi and, via the alpha copy,
`delta_cover_exp`; the independent cover field carries a cover-specific
structure the alpha copy cannot express. Its intercept + coefficient
columns become separate ICAR blocks on the cover arm, each with its own
precision integrated on the outer grid and reported as `sigma_pos_field`
(intercept) / `sigma_pos_field_<col>` (a covariate column). The field SD
grid defaults to `control$sigma.grid`; set
`control$sigma.grid.pos.field` to integrate it over its own (usually
coarser) grid, which keeps the added axis from multiplying the
outer-grid cost. The per-cell field posterior is in `fit$pos_field` /
`fit$pos_field_table` (and `fit$pos_field_tables` per column). Scope: it
composes with the shared occupancy field but uses per-visit cover
(`cover_aggregate = "none"`) and does not combine with the correlated
`|` MCAR field, the latent cover RE, or the batched fused path; like the
shared field it is fitted as ICAR (bym2/car read as ICAR). A static
intercept cover field is only weakly identified against the occupancy
field's alpha copy (both are per-cell cover offsets); a time-weighted
trend cover field is identified separately and is what makes
`delta_cover_cond` spatially varying.

The same placement works on the detection arm: a spatial-field term in
the `detection` formula
(`detection = ~ 1 + spatial(~ 0 + time || cell, graph = adj)`) fits an
independent field on the detection predictor, for a spatially-structured
detection probability. Its SD is reported as `sigma_p_field` /
`sigma_p_field_<col>`.

## Checkpoint / resume

A full-field `occu_cover()` fit integrates over a large outer
hyperparameter grid and can run for hours, so a reboot or OOM kill
otherwise loses the whole run.
`control$checkpoint = list(path = "fit.ckpt", resume = TRUE)` makes the
joint engine append each completed grid cell to `path`; a
`resume = TRUE` run loads the finished cells and solves only the
remaining ones, reproducing the from-scratch fit. `resume = FALSE`
starts a fresh file (any stale checkpoint at `path` is removed first). A
torn tail from a killed write is discarded and re-solved, and a
checkpoint written for different data or settings is rejected rather
than resumed onto. Forwarded to
[`tulpa::tulpa_nested_laplace_joint()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace_joint.html).

## Cell-aggregated cover (`cover_aggregate`)

On the shared-field spatial path the cover arm has one observation per
detected visit, so a cell with many detected plots informs the shared
field far more than the single occupancy observation for that cell; the
field is then driven almost entirely by the cover arm and the
detection-corrected occupancy surface flattens. `cover_aggregate`
collapses the cover arm to a single response per occupancy unit
(cell-period) so the two arms inform the field with comparable weight:

- `"mean"` – model the mean cover over the unit's detected visits;

- `"median"` – model the median cover;

- `"latent"` – model a per-unit cover random effect \\u_i \sim N(0,
  \sigma_u^2)\\ shared across the unit's detected visits and integrated
  out per unit, so the cover arm contributes one marginal observation
  per unit while keeping every detected visit (the principled
  counterpart of mean / median aggregation). The lognormal arm
  integrates in closed form (compound-symmetry); the beta arm uses
  adaptive Gauss-Hermite quadrature (`control$n.quad`, default 15). The
  within-unit dispersion is pre-fit from the within-unit spread and held
  fixed; \\\sigma_u\\ is integrated on the outer grid
  (`control$sigma.u.grid`).

- `"none"` – the per-visit cover arm (one cover observation per detected
  visit).

`NULL` (the default) selects `"mean"` on the shared-field spatial
`nested_laplace` path and `"none"` on the non-spatial `laplace` path
(which has no shared field to over-weight). Aggregation needs a
cell-level positive design: a visit-level covariate in the `positive`
formula cannot be collapsed to one value per cell and errors.
Aggregation is currently wired on the spatial path only; requesting
`"mean"` / `"median"` on the non-spatial `laplace` fit errors rather
than silently using per-visit cover.

## Random effects on the observation arms

Under `method = "nested_laplace"` (shared-field spatial) and
`method = "nuts"` (non-spatial, which samples the group SD), random
intercepts may be written on the **detection** or **positive-cover**
formula with the usual `lme4` bar or
[`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
spelling, e.g. `detection = ~ effort + (1 | habitat)`. The grouping is
per visit – one code per `(site, visit)` entry, distinct from the
per-site occupancy-arm grouping – so a many-level categorical visit
covariate (an EUNIS habitat class, an observer) enters as a partially
pooled random intercept. **Crossed** (`(1 | habitat) + (1 | observer)`)
and **nested** (`(1 | region/site)`, which expands to
`re(region) + re(region:site)`) groupings are supported: each term joins
the joint fit as its own `iid` latent block. Each block's variance
integrates on the outer grid (`sigma_re_p` / `sigma_re_pos` for a lone
term on the detection / cover arm, suffixed by the grouping variable –
`sigma_re_p_habitat` – when several terms share an arm; tune with
`control$re.sigma.grid.p` / `re.sigma.grid.pos`). The per-group BLUPs
are reported in `fit$re` (one entry per term, keyed by arm or
`"<arm>:<var>"`) and via
[`ranef()`](https://gillescolling.com/tulpa/reference/ranef.html).
[`predict()`](https://rdrr.io/r/stats/predict.html) sums every term's
BLUP offset on the predicted arm when `newdata` carries the grouping
columns and `type = "detection"` (or `"occurrence"` / `"cover_cond"`);
an unseen / held-out level shrinks that term to the population mean.
Several crossed terms multiply the outer grid, so set
`control$integration = "ccd"` (and/or coarsen the RE grids) at scale.

Random **slopes** are also supported (tulpa engine \>= 0.0.39): an
uncorrelated slope (`(x || g)`, `(0 + x | g)`) is one per-row weighted
`iid` block per coefficient, and a correlated slope (`(1 + x | g)`) one
multivariate free-Sigma `miid` block. A slope term's `fit$re` entry
carries an `[n_groups x n_coefs]` BLUP matrix, a per-coefficient
`sigma`, and (correlated) a `cor` matrix; the hyperparameters are
`sigma_re_p_<coef>` and `cor_re_p_<ci>_<cj>`, marginalized over the
grid. [`predict()`](https://rdrr.io/r/stats/predict.html) weights each
coefficient by its covariate column in `newdata` (intercept = 1). Each
slope covariate is **standardized** to unit SD before fitting (mirroring
the fixed-effect design), so the variance grid is scale-invariant – the
reported slope BLUPs and `sigma` are back-transformed to the covariate's
natural units (correlation is scale-free). A correlated slope adds
`p(p+1)/2` Sigma axes to the outer grid, so its free-Sigma grid uses a
compact principled default (symmetric correlation nodes including 0);
widen it with `control$re.logchol.grid.p` / `re.logchol.grid.pos`.

The positive-cover RE needs per-visit cover
(`cover_aggregate = "none"`). As with the occupancy-arm RE, each
grid-integrated variance carries the binary / small-cluster
inner-Laplace attenuation and is a lower bound on the truth; the BLUPs
recover the per-group structure.

`method = "nuts"` samples an observation-arm random **intercept**
instead of integrating it on a grid: each grouping factor is a
non-centered block (`b_g = sigma_re z_g`, `z_g ~ N(0, 1)`) whose group
SD is a coordinate of the sampled vector under a `N(0, 1.5^2)` prior on
`log sigma_re` – the same weakly informative width the other
observation-family samplers use, and the reason the sampled SD carries
no small-cluster attenuation to correct for. Crossed and nested
groupings are simply several blocks, on either or both observation arms.
`fit$re` keys and
[`ranef()`](https://gillescolling.com/tulpa/reference/ranef.html) match
the `nested_laplace` path; each term additionally reports `sigma_median`
and the per-draw `sigma_draws`, the summaries to quote for a
right-skewed variance component, with per-SD convergence in
`fit$nuts$re_sigma_rhat` / `re_sigma_ess`. The sampled coefficient
surface ([`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), `fit$draws`) is
unchanged. Still on `nested_laplace` only: a random slope, and an
observation-arm RE composed with the coupled areal field.
[`waic()`](https://mc-stan.org/loo/reference/waic.html) /
[`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
on a sampled fit read the coefficient draws, so they score the RE at
zero – the same limitation the sampled coupled field has.

## See also

[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) (no
cover),
[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
(plot-level hurdle, no detection),
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) (counts
not cover).

## Examples

``` r
# \donttest{
N <- 120; J <- 4
sim <- simulate_occu_cover(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                           n_pos_covs = 1, positive = "beta", seed = 1)
long <- data.frame(site_id = rep(seq_len(N), each = J),
                   visit    = rep(seq_len(J), times = N),
                   y        = as.vector(t(sim$y)),
                   det_cov1 = sim$visit_data$det_cov1,
                   pos_cov1 = sim$visit_data$pos_cov1)
od  <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                 det.covs = c("det_cov1", "pos_cov1"))
cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
fit <- tobs(~ occ_cov1, data = cell_dat, family = occu_cover("beta"),
            detection = ~ det_cov1, positive = ~ pos_cov1,
            y = od$y, y_pos = y_pos, visits = od$det.covs, method = "laplace")
summary(fit)
# }
```
