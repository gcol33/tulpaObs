# Simulation-based calibration for a fitted tobs model

Runs the POSTERIOR simulation-based calibration experiment (Talts et al.
2018; Sailynoja, Schmitt, Buerkner and Vehtari 2026, Algorithm 2) on a
fitted `tobs_fit`, through the engine's
[`tulpa::sbc()`](https://gillescolling.com/tulpa/reference/sbc.html)
front door. The truth is drawn from the fit's own posterior at the
observed data, a replicate is simulated at that truth on FRESH cells,
the model is refitted on the two data sets together, and the truth's
rank under that augmented posterior is recorded. Under exact inference
those ranks are Uniform(0, 1); the departure from uniformity is a
verdict on the approximation.

## Usage

``` r
# S3 method for class 'tobs_fit'
sbc(
  object,
  n.sim = 100L,
  n.draws = 1000L,
  n.ref = 200L,
  quantities = NULL,
  controls = character(),
  bad.factor = 1.25,
  level = 0.95,
  seed = 0L,
  model.only = FALSE,
  fit.control = list(),
  control = list(),
  ...
)
```

## Arguments

- object:

  A fitted `tobs_fit`. Its family must be registered; see Details for
  the roster.

- n.sim:

  Number of simulations. Each costs one refit on the pooled data.

- n.draws:

  Posterior draws per fit backing the reported predictives.

- n.ref:

  Reference draws behind the `log_lik` rank. `0` drops that arm.

- quantities:

  Optional character vector restricting what is scored.

- controls:

  Mis-scaled control arms to add: any of `"wide"`, `"narrow"`.

- bad.factor:

  The control arms' SD multiplier.

- level:

  Simultaneous band level.

- seed:

  Seed offset. Simulation `i` draws its truth at `seed + i`.

- model.only:

  Return the callback list instead of running the experiment, for
  inspection or for passing to
  [`tulpa::sbc()`](https://gillescolling.com/tulpa/reference/sbc.html)
  directly.

- fit.control:

  Merged into the `control` list of every refit.

- control:

  Passed to
  [`tulpa::sbc()`](https://gillescolling.com/tulpa/reference/sbc.html)
  (`progress`, `rand_seed`).

- ...:

  Unused.

## Value

An object of class `sbc` from tulpa, carrying the per-quantity rank
ECDFs, the exact simultaneous band, and the uniformity tests. When
`model.only = TRUE`, the callback list instead.

## Details

Registered families: `occu_cover` (the coupled occupancy + cover hurdle
on the joint nested-Laplace engine), and `occu`, `abun`, `count`,
`removal`, `distance`, `fp_occu`, `royle_nichols`, `occu_ttd` and
`double_observer`, whose site marginals multiply. For those the
replicate is the family's own
[`simulate()`](https://rdrr.io/r/stats/simulate.html) kernel at the
drawn theta, the rank arm is the family's exact marginal log-likelihood,
and two kinds of fit are refused rather than approximated: a structured
term (its field is a latent quantity shared across sites that theta does
not hold, which is what the coupled `occu_cover` route handles by
drawing on fresh cells) and a visit-level observation design (the
replicate kernel assembles the observation arm from the site-level
design, so the two sides would carry different models).

`dyn_occu` (constant-rate only) pools on the SITE axis and leaves the
season axis alone, since its response is 3D
`[n_sites x max_visits x n_seasons]`: a bespoke replicate generator
forward-simulates the colonization/extinction HMM, but the posterior
draws and the rank statistic are the same shared machinery every other
registered family uses.

`int_occu` (full source-site overlap only) pools on the SITE axis and
leaves the per-source axis alone, since its response is a list of one
detection matrix per source; the rank statistic is the shared machinery
unchanged. Registering it surfaced a real bug (since fixed): the
detection intercept was silently left in standardized-covariate units on
every fit with an autoscaled detection covariate, independent of SBC.

`gdistremoval` pools two response matrices (`yDist` band counts, `yRem`
period counts) together on the SITE axis, reusing the same
list-of-matrices `pool` `int_occu` uses; its replicate generator is its
own [`simulate()`](https://rdrr.io/r/stats/simulate.html) handler.

`dyn_abun` (constant-rate only) shares `dyn_occu`'s 3D
`[n_sites x max_visits x n_seasons]` response and site-axis pooling, but
unlike `dyn_occu` it already has a working
[`simulate()`](https://rdrr.io/r/stats/simulate.html) handler (the
Dail-Madsen open N-mixture forward is not a two-state HMM that needs a
bespoke generator), so its replicate is the same
family-[`simulate()`](https://rdrr.io/r/stats/simulate.html) route as
`occu`/`abun`/etc.

`occu_categorical` is a `tobs_multiarm_fit` with two independent
Laplace-Gaussian blocks (presence, class) rather than a joint MVN draw
matrix; `draws()` samples the two blocks independently and
`loglik_many()` scores the two-arm likelihood directly, since this
family has no `.tobs_pointwise_loglik` dispatch to reuse.

`distsamp_open` (constant-dynamics, Poisson only) shares `dyn_abun`'s 3D
response and site-axis pooling, and its `fit$means`/`fit$draws` are the
standard `.tobs_bfgs_marginal_fit()` shape, so `draws`/`simulate`/
`loglik_many` are all the shared generic ones.

`occu_multi`'s response is a list of S per-species matrices, the same
shape
[`int_occu()`](https://gillescolling.com/tulpaObs/reference/int_occu.md)
pools; its [`simulate()`](https://rdrr.io/r/stats/simulate.html) is
custom, since species share one joint multi-species state rather than
being independently observed.

`dyn_int_occu` is the product of the multi-season and multi-source
shapes: a named list of S per-source 3D arrays, pooled by a new
`.tobs_sbc_pool_named_3d` that composes both existing pooling rules.

`t_occu` is a `pg_gibbs` family whose `fit$draws` is already the real
pooled posterior sample; it has no
[`simulate()`](https://rdrr.io/r/stats/simulate.html) handler and no
`.tobs_pointwise_loglik` dispatch, so both are custom – `simulate` draws
a fresh AR1 year effect at the theta's own `(sigma, rho)`, and
`loglik_many` is a Laplace approximation to that year effect's marginal
(FD-validated, brute-force cross-checked at T = 2).

`cover` is a `tobs_multiarm_fit` with the same two-independent-block
shape as `occu_categorical`; `positive = "lognormal"` only for v1, with
dispersion held fixed (no SE is reported for it anywhere in the
package).

`occu_multiscale_cover` is a standard single-block fit (unlike `cover`'s
two-block one); the exchangeable unit is the CELL rather than the plot,
so pooling and site labeling track cell indices.

A family is registered by adding one entry to the internal registry – a
replicate generator, a refit call and optionally a joint statistic; the
pooling, the grouping labels, the arm construction and the controls are
shared.

## What is scored

The arm fixed effects, the areal field SD on each arm, the copy scale
`alpha`, and the dispersion when the fit estimates one. The per-cell
field itself is not scored – it is integrated out by the fit and carries
no truth theta could hold. The cover-arm field SD and `alpha` are read
per draw off the outer hyperparameter grid, whose cells are sampled by
their own normalized weight, so what is ranked is the grid-marginalized
posterior of each rather than a function of component modes.

A quantity whose posterior has no spread – a dispersion the joint engine
holds fixed in its cell-coupling spec, for instance – is dropped rather
than scored, and named in `attr(x, "fixed")`.

An extra `log_lik` quantity ranks a joint statistic of the whole
parameter vector, which catches a posterior getting each marginal right
while getting the dependence between them wrong.

## Fresh cells

The replicate is drawn on new cells whose areal graph is block-diagonal
against the observed one. `occu_cover` couples its occurrence and cover
arms through one shared field, and theta carries no per-cell field
value, so a replicate re-observing the observed cells would depend on
the observed data given theta and break the factorization posterior SBC
rests on. `group_ids` is supplied to
[`tulpa::sbc()`](https://gillescolling.com/tulpa/reference/sbc.html),
which verifies the observable half of that (disjoint group labels)
rather than assuming it.

## Controls

`controls = c("wide", "narrow")` adds arms reporting the same draws
rescaled about their mean by `bad.factor` and its reciprocal. They are
deliberately mis-scaled posteriors and are expected to leave the band: a
calibration read that nothing can fail is not evidence.

## References

Talts, S., Betancourt, M., Simpson, D., Vehtari, A. and Gelman, A.
(2018). Validating Bayesian inference algorithms with simulation-based
calibration. arXiv:1804.06788.

Sailynoja, T., Schmitt, M., Buerkner, P.-C. and Vehtari, A. (2026).
Posterior SBC: simulation-based calibration checking conditional on
data. Statistics and Computing 36:78.

## See also

[`tulpa::sbc()`](https://gillescolling.com/tulpa/reference/sbc.html),
[`waic()`](https://mc-stan.org/loo/reference/waic.html)

## Examples

``` r
# \donttest{
N <- 20L; J <- 3L
adj <- matrix(0L, N, N)
for (s in seq_len(N)) {
  if (s > 1L) adj[s, s - 1L] <- 1L
  if (s < N)  adj[s, s + 1L] <- 1L
}
sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                           adj = adj, sigma = 0.7, alpha = 1, seed = 1L)
long <- data.frame(site_id = rep(seq_len(N), each = J),
                   visit = rep(seq_len(J), times = N),
                   y = as.vector(t(sim$y)),
                   det_cov1 = sim$visit_data$det_cov1,
                   pos_cov1 = sim$visit_data$pos_cov1)
od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                det.covs = c("det_cov1", "pos_cov1"))
y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
# The dispersion goes on the outer grid so it is estimated and scored, and
# the field-SD grid is pinned so both stages integrate the same support.
ctl <- list(engine = "joint", verbose = FALSE,
            sigma.grid   = exp(seq(log(0.15), log(2.0), length.out = 9)),
            phi.grid.pos = exp(seq(log(0.20), log(0.90), length.out = 7)))
fit <- tobs(~ occ_cov1 + icar(graph = adj),
            data = cbind(data.frame(site_id = seq_len(N)), sim$data),
            family = occu_cover("lognormal"),
            detection = ~ det_cov1,
            positive = ~ pos_cov1 + share(spatial()),
            y = od$y, y_pos = y_pos, visits = od$det.covs,
            method = "nested_laplace", control = ctl)
sbc(fit, n.sim = 20L, controls = "narrow", fit.control = ctl)
# }
```
