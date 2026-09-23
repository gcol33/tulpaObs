# Model criteria for occupancy / cover models

[`waic()`](https://mc-stan.org/loo/reference/waic.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html),
[`dic()`](https://gillescolling.com/tulpa/reference/criteria_doors.html),
and
[`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
on a `tobs_fit` build the family pointwise log-likelihood matrix once
and hand it to the engine's single criteria layer
[`tulpa::tulpa_criteria()`](https://gillescolling.com/tulpa/reference/tulpa_criteria.html),
which derives WAIC / DIC / CPO / LPML / PSIS-LOO. DIC additionally
evaluates the deviance at the posterior mean of the parameters, supplied
by the family-specific `.tobs_loglik_at_mean()`.

## Usage

``` r
# S3 method for class 'tobs_fit'
waic(x, n.draws = 1000L, loo.unit = c("obs", "cell"), n.threads = NULL, ...)

# S3 method for class 'tobs_fit'
loo(x, n.draws = 1000L, loo.unit = c("obs", "cell"), n.threads = NULL, ...)

# S3 method for class 'tobs_fit'
dic(object, n.draws = 1000L, n.threads = NULL, ...)

# S3 method for class 'tobs_fit'
cpo(
  object,
  n.draws = 1000L,
  loo.unit = c("obs", "cell"),
  n.threads = NULL,
  ...
)

# S3 method for class 'tobs_fit'
pointwise_loglik(object, ndraws = NULL, ...)
```

## Arguments

- x, object:

  A `tobs_fit` object.

- n.draws:

  Posterior draws used to build the pointwise log-likelihood (the cover
  / occu_cover paths sample this many; the draw-matrix families use the
  first `n.draws` stored draws). Default 1000.

- loo.unit:

  The cross-validation unit for
  [`waic()`](https://mc-stan.org/loo/reference/waic.html) /
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) /
  [`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html).
  `"obs"` (default) is the family's pointwise unit – one column of the
  log-likelihood per plot (cover) or site (occu_cover) – and is
  byte-identical to the call without the argument. `"cell"` switches to
  leave-one-group-out cross-validation (LOGO-CV): the fit's own
  per-observation cell map folds the log-likelihood columns of a spatial
  cell into one, so each cell is a fold instead of each plot / site.
  [`waic()`](https://mc-stan.org/loo/reference/waic.html) and
  [`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
  hand the map to
  [`tulpa::tulpa_criteria()`](https://gillescolling.com/tulpa/reference/tulpa_criteria.html)
  as `group`; [`loo()`](https://mc-stan.org/loo/reference/loo.html)
  applies the same fold itself and runs PSIS on the resulting
  `[n_draws x n_cells]` matrix, whose column is the cell's joint
  conditional log-likelihood per draw, so the importance ratio inverts
  the whole cell's likelihood. A whole cell is a larger perturbation of
  the posterior than a single row, so the Pareto k values are
  correspondingly higher and are the diagnostic to read before trusting
  the cell-level number. Implemented for
  [`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
  (the areal field node, when sites are grouped via `group_var`) and
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  (the `site_cell` map); a non-spatial fit has no cells, so `"cell"`
  errors there. Equivalent to passing `group =` the cell map directly,
  without hand-building it.
  [`dic()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
  has no cross-validation unit – it is a plug-in deviance over all
  observations – and rejects `loo.unit`.

- n.threads:

  Threads for the parallel
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  pointwise log-likelihood (the compact / ragged path). The draw loop is
  embarrassingly parallel, so this is the WAIC / LOO analogue of the
  fit's `n.threads.outer`. `NULL` (default) uses all but four logical
  cores, at most two under `R CMD check`; other families ignore it.

- ...:

  For
  [`dic()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
  and
  [`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html),
  forwarded to
  [`tulpa::tulpa_criteria()`](https://gillescolling.com/tulpa/reference/tulpa_criteria.html)
  (e.g. `chunk_size`, or an explicit `group =` for a custom
  leave-one-group-out unit).
  [`waic()`](https://mc-stan.org/loo/reference/waic.html) and
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) accept only
  `group =`, which folds the pointwise matrix before the loo call, and
  reject any other argument.

- ndraws:

  For
  [`pointwise_loglik()`](https://gillescolling.com/tulpa/reference/criteria_doors.html),
  the posterior draws used to build the matrix, as `n.draws` above.
  `NULL` (default) uses 1000.

## Value

[`dic()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
and
[`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
return a `tulpa_criteria` object.
[`waic()`](https://mc-stan.org/loo/reference/waic.html) returns a `waic`
object and [`loo()`](https://mc-stan.org/loo/reference/loo.html) a
`psis_loo` object from the loo package; read the estimates from
`$estimates` (e.g. `waic(fit)$estimates["waic", "Estimate"]`). Both
carry the cross-validation unit they scored in their `"loo_unit"`
attribute.

## Details

[`waic()`](https://mc-stan.org/loo/reference/waic.html) and
[`loo()`](https://mc-stan.org/loo/reference/loo.html) are the loo
package's generics and build their objects through
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) and
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) on the same
pointwise matrix ([`loo()`](https://mc-stan.org/loo/reference/loo.html)
via PSIS with relative effective sample sizes), so
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
and the rest of that ecosystem read them directly.
[`dic()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
and
[`cpo()`](https://gillescolling.com/tulpa/reference/criteria_doors.html)
are tulpa's, and return a `tulpa_criteria` object.

## See also

[`tulpa::tulpa_criteria()`](https://gillescolling.com/tulpa/reference/tulpa_criteria.html),
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
