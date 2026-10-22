# Convergence record for a fitted model

The public accessor for whether a fit converged, with one return shape
across every
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) family.
Every family's optimiser / EM / sampler verdict is recorded under
`fit$convergence`, and a Laplace-family fit that stops without meeting
its criterion warns once, naming `control$max.iter`. `convergence()`
returns the full record (`converged`, `n_iter`, `sla_status` when the
simplified-Laplace marginals were used, and the per-group AGHQ solve
status on the families fitted by that engine), and `converged()` returns
the single logical.

## Usage

``` r
convergence(object, ...)

# S3 method for class 'tobs_fit'
convergence(object, ...)

converged(object, ...)

# S3 method for class 'tobs_fit'
converged(object, ...)
```

## Arguments

- object:

  A fitted `tobs_fit` (occupancy / abundance / cover / ...).

- ...:

  Ignored.

## Value

`convergence()`: a list with `converged` (logical), `n_iter` (integer,
`NA` for grid / closed-form fits with no iteration count), and
`sla_status` (character, when present). A fit whose random effects were
integrated by adaptive Gauss-Hermite quadrature also carries `group_ok`
(logical, one entry per group: `FALSE` where that group's posterior
solve failed and its estimates are `NA`), `groups_failed` (their
indices) and `groups_failed_names` (their labels where the family names
its groups – species, for the community models). `converged` is `FALSE`
whenever any group failed, since the quantities that group contributes
to are not estimates of anything; filter on these rather than on the
text of the warning the engine raises. A single-species fit whose
grouped random effects were integrated that way carries `re_boundary`:
one record per random-effect SD (`sigma_<g>_<coef>`) and, on a
correlated block, per correlation (`cor_<g>_<ci>_<cj>`), with the
boundary Wald `statistic`, its `critical` value and `distinguishable`;
`FALSE` means the SD cannot be told from zero or the correlation from
+-1, and the fit says so in one message ("boundary (singular) fit").
`converged()`: a single `TRUE` / `FALSE`.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
converged(fit)
convergence(fit)$n_iter
# }
```
