# Structured terms inside a tobs() formula

Spatial fields, random effects, temporal structure, varying coefficients
and latent factors are written inside the model formula, in the process
they belong to, rather than passed as
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
arguments. A term entered in the state formula acts on occupancy or
abundance; the same term written in `detection = ~ ...` acts on
detection.

## Details

These are formula specials, not exported functions: they are resolved
against an internal registry while the formula is parsed, and calling
them outside a
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) formula
is an error. This page documents the call shapes and their arguments.

## Areal spatial fields

`icar(graph)`, `bym2(graph, scale_factor)`, `car(graph)` and
`car_proper(graph)` put a field over the nodes of a symmetric adjacency
matrix `graph`. `bym2()` computes its Riebler `scale_factor` from the
graph when it is not supplied. All four accept:

- `group_var`:

  Column naming each observation's graph node, when the graph is over
  regions rather than rows.

- `weight`:

  Per-node numeric column turning the field into a spatially varying
  coefficient (`weight_i * z[node_i]`) instead of an intercept field.
  `car()` does not take it.

- `id`:

  Name for the realization, so `share("id")` can share it with another
  process.

The umbrella form
`spatial(~ 1 + w || cell, graph = adj, model = "icar")` forwards to the
same terms and additionally reads a bar, which is how areal varying
coefficients are usually written (see below).

An intrinsic areal field (`icar`, `bym2`, `car`) has one constant null
direction per CONNECTED COMPONENT of its graph, and each of them is
identified separately: the engine augments the precision with one
sum-to-zero per component (INLA's `adjust.for.con.comp = TRUE`), and the
components share one field variance. A graph with islands or with
separate survey regions is therefore fitted as intended. Because the
component count is also what an adjacency looks like when a spatial join
has quietly dropped edges, and because it is not visible in the fitted
object, a graph with more than one component reports its component count
and sizes at fit time.

A component of a single node has no neighbours to borrow from. The areal
cover paths reject one outright, with an error naming the node; where
such a node is accepted it keeps a proper `N(0, 1/tau)` effect, because
the identification augments the precision rather than imposing a hard
constraint, which on a lone node would pin it to exactly zero.

## Continuous spatial fields

`gp(lon, lat, ...)` is an NNGP-approximated Gaussian process over
coordinates, `multiscale_gp(lon, lat, ...)` a two-scale (local plus
regional) version, and `spde(lon, lat, ...)` a Matern field over a
triangular mesh. Coordinates come either from bare column names or from
`coords = `, a two-column matrix.

- `cov`:

  Covariance function: `"exponential"` (default), `"matern"`,
  `"gaussian"`, `"spherical"`.

- `nu`:

  Matern smoothness.

- `nn`:

  Number of nearest neighbours in the Vecchia approximation, capped at
  `n - 1`. `multiscale_gp()` splits this into `nn_local` and
  `nn_regional`.

- `prior_range`:

  Length-2 `c(r0, alpha)` PC prior on the range, read as
  `P(range < r0) = alpha`. Required by `gp()` and `svc()`, which ship no
  default; `spde()` defaults to `c(0.5, 0.5)` alongside `prior_sigma`.

- `mesh`, `max_edge`, `cutoff`:

  `spde()` mesh controls, passed to
  [`tulpa::spatial_spde()`](https://gillescolling.com/tulpa/reference/spatial_spde.html).

## Spatially varying coefficients

A coefficient that varies over space is a spatial field like any other,
so it is written under `spatial()` and the `model = ` choice says
whether space is a graph or a set of coordinates. In both forms you name
the coefficients that vary; the difference is only what they vary over.

**Areal**, over graph nodes, written as a bar:

    ~ 1 + w + spatial(~ 1 + w || cell, graph = adj)

with `method = "nested_laplace"`. The bar's left-hand side names the
varying coefficients (one field per column, the intercept included), its
right-hand side the graph node. Broad family support (single-season,
dynamic, integrated and the community routes) and recovery-tested
throughout. This is what most spatially-varying-coefficient questions
want. `model = ` selects the areal structure (`"icar"` by default, or
`"bym2"` / `"car"` / `"car_proper"`).

**Continuous**, over coordinates, as NNGP surfaces:

    ~ elevation + spatial(lon, lat, model = "svc",
                          coefficients = "elevation", prior_range = c(50, 0.05))

There is no node index to group on, so the coefficients are named
directly instead of through a bar.

- `coefficients`:

  Names of the design columns whose coefficients vary, matched against
  the design of the arm the surfaces load on. One surface per entry;
  `"(Intercept)"` gives a plain intercept surface.

- `indices`:

  The same selection by column position, for a design column with no
  usable name. Give `coefficients` or `indices`, not both.

- `coords`:

  Two-column coordinate matrix, as an alternative to bare `lon, lat`
  column names. One coordinate pair per site.

- `cov`:

  `"exponential"` (default), `"matern"` or `"gaussian"`.

- `nn`:

  Nearest neighbours in the Vecchia approximation, default 15, capped at
  `n - 1`.

- `prior_range`:

  Required, no default: `c(r0, alpha)` with `P(range < r0) = alpha`.

- `sigma2_prior_scale`:

  Scale of the prior on the surface's marginal variance.

`svc(lon, lat, coefficients = )` is the direct constructor and stays
available, the way `icar()` does alongside `spatial(model = "icar")`.

The continuous form fits on single-season
[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md) under
`"laplace"`, `"nested_laplace"` and `"nuts"`, and on
[`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md),
[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md),
[`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md)
and
[`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md)
under `"laplace"` and `"nested_laplace"`. The surfaces load on the state
arm (the occupancy logit, or log lambda for the count families) and come
back as `fit$svc_field`, with the integrated hyperparameters in
`fit$svc_hyper` and the resolved design columns in `fit$svc_indices`. On
any other family, method or arm the term errors and points at the areal
bar rather than being silently dropped.

## Random effects

`re(group, ...)` is a grouped random effect, and `lme4` bars are sugar
over it: `(1 | g)`, `(x | g)`, `(x || g)`, `(0 + x | g)`, `(1 | g:h)`
and `(1 | g/h)` are rewritten to `re()` calls before the formula is
evaluated.

- `type`:

  `"intercept"` (default), `"slope"` or `"iid"`.

- `covariate`:

  Slope column(s): a column name, several names, a
  [`cbind()`](https://rdrr.io/r/base/cbind.html) matrix or a one-sided
  formula. Required for `type = "slope"`.

- `model`:

  Structure over the group levels: `"iid"` (default), `"ar1"`, `"rw1"`,
  `"rw2"`.

- `correlated`:

  Whether slopes share a correlated covariance with the group intercept.

- `intercept`:

  `FALSE` drops the implicit group intercept, giving a slope-only block.

## Temporal structure

`temporal(time, ...)` puts a field over the levels of `time`.

- `type`:

  `"ar1"` (default), `"rw1"`, `"rw2"` or `"iid"`.

- `group`:

  Optional grouping, for one series per group.

- `cyclic`:

  Wraps the last level onto the first.

- `tau_shape`, `tau_rate`:

  Gamma prior on the field precision.

## Latent factors and sharing

`latent(n_factors)` adds per-site latent factors with per-species
loadings to a community model, giving residual co-occurrence.
`constraint` and `sigma_prior_rate` control the identification anchor
and the loading prior.

`share("id")` shares one realization of the term named `id` across
processes, optionally scaled: `alpha = 0.5` pins the coupling amplitude,
`alpha = grid(c(0.25, 0.5, 1))` marginalizes over it on the outer
integration.

`residual =` carries an arm-specific deviation beside the copy, so the
sharing arm's surface is `alpha * w + delta` rather than `alpha * w`
alone. Use it where the two arms respond to related but not identical
surfaces: a plain copy has no way to bend, so a pattern the shared field
does not have is not fitted at all. `residual = "full"` gives the
deviation one latent per field node, spanning the same surfaces as a
field placed on the arm itself. `residual = r` gives it `r` basis
functions on a basis orthogonalized against the shared field, so the
latent problem is `K + r` rather than `2K`, and the deviation cannot
reproduce the copied component. Because the two halves are orthogonal,
the reported amplitude is identified: it is the projection of the fitted
arm surface onto the occurrence surface, taken on the linear-predictor
scale. A rank-`r` deviation is pinned at the shared field's own SD,
which costs one warm fit of the model without it. Available on
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
under `method = "nested_laplace"`, one per fit, and it needs a shared
field to deviate from.

[`grid()`](https://rdrr.io/r/graphics/grid.html) states the axis's
NODES, and the span the fit works over is wider than the nodes
themselves, for two separate reasons.

First, each node represents the cell around it, so the axis runs half a
node step past the outermost node at each end, in the axis's own log
coordinate. That is what makes the quadrature weights sum over a whole
axis rather than stopping on the last node, and under `method = "nuts"`
the sampled amplitude takes the same span as the support of its flat
prior, so one stated grid moves both backends over one region. For `k`
equally log-spaced nodes this is exact: `k - 1` steps between the nodes
plus a half step at each end, so the span is `k / (k - 1)` times the
stated range in the log coordinate. It is therefore largest for the
shortest grid – 2x at two nodes, 1.125x at nine. With refinement off,
`c(0.2, 0.5)` gives `[0.127, 0.791]`.

Second, refinement. `control$adaptive.grid` is `TRUE` by default, and
what it does depends on who placed the axis. An axis YOU stated is
DENSIFIED: further nodes are placed inside the range you wrote, never
past its ends, so the nodes you give are a bound and the half step above
is the whole of the widening. An axis the fit placed for you – the
default copy axis, a field SD axis – is EXTENDED when a boundary cell
still holds weight (`adaptive.grid.edge.thresh`, default 0.02), which is
what lets a fit whose default grid misses the mode go and find it.

Refining an axis changes how densely its nodes sit, not the region they
cover: the reported span stays the declared span unless a refinement
point lands past it. A stated `c(0.2, 0.5)` is only densified, so it
spans `[0.127, 0.791]` with refinement on as with it off. Holding an
axis at exactly the nodes stated takes BOTH
`control$adaptive.grid = FALSE` and
`control$var.of.means.consistency = FALSE`: refinement reaches an axis
from two independent passes, and the consistency pass runs whatever the
adaptive grid is set to, so a fit that turns off only the first still
integrates nodes placed inside the range it wrote. Under
`method = "nuts"` a fit reports the span each sampled hyperparameter was
worked over in `fit$nuts$hyper_support`, which is the region that
hyperparameter's flat prior is supported on.

## See also

[`tobs`](https://gillescolling.com/tulpaObs/reference/tobs.md) for the
model front door.

## Examples

``` r
# `occu()` takes the detection history as `y` and a ONE-sided state formula:
# a two-sided formula would be resolved against a response the family does
# not carry. `data` is site-level, `visits` visit-level.
sim <- simulate_occu(N = 30, J = 3, seed = 1)
d   <- transform(sim$data, cell = seq_len(30),
                 lon = runif(30), lat = runif(30))

# A chain graph over the 30 cells: neighbours are consecutive indices.
adj <- matrix(0L, 30, 30)
adj[cbind(1:29, 2:30)] <- 1L
adj[cbind(2:30, 1:29)] <- 1L

# \donttest{
ctrl <- list(verbose = FALSE, progress = FALSE)

# Areal field on occupancy
fit <- tobs(~ occ_cov1 + icar(graph = adj, group_var = "cell"),
            data = d, family = occu(), detection = ~ 1, y = sim$y,
            method = "nested_laplace", control = ctrl)
fit$spatial_field[1:5]

# Detection random effect by observer, who changes between visits
obs <- simulate_occu(N = 60, J = 3, n_visit_groups = 4, seed = 1)
fit_re <- tobs(~ occ_cov1, data = obs$data, family = occu(),
               detection = ~ (1 | visit_group), y = obs$y,
               visits = obs$visits, method = "laplace", control = ctrl)
coef(fit_re)

# Areal spatially varying coefficient on `occ_cov2`
fit_svc <- tobs(~ occ_cov2 + spatial(~ 1 + occ_cov2 || cell, graph = adj),
                data = d, family = occu(), detection = ~ 1, y = sim$y,
                method = "nested_laplace", control = ctrl)
coef(fit_svc)

# Continuous NNGP varying coefficient on the `occ_cov2` slope
fit_nngp <- tobs(~ occ_cov2 + spatial(lon, lat, model = "svc",
                                      coefficients = "occ_cov2",
                                      prior_range = c(0.5, 0.05)),
                 data = d, family = occu(), detection = ~ 1, y = sim$y,
                 method = "laplace", control = ctrl)
coef(fit_nngp)
# }
```
