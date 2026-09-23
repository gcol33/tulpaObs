# One-row model summary for a tobs_fit

Every family, and every method within a family, glances to one column
set: `nobs`, `df`, `logLik`, `n_fixed`, `n_samples`, `n_divergent`,
`mean_accept` and `converged`, the sampler columns `NA` on a fit with no
draws. The joint nested-Laplace fits add the outer grid placement and,
when requested, the outer Pareto-k diagnostic. The joint-coupled
families
([`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md),
[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md),
[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md)
spatial,
[`occu_multiscale_cover()`](https://gillescolling.com/tulpaObs/reference/occu_multiscale_cover.md))
carry both at the fit top level.

## Usage

``` r
# S3 method for class 'tobs_fit'
glance(x, ...)
```

## Arguments

- x:

  A fitted `tobs_fit`.

- ...:

  Ignored.

## Value

The one-row `glance` data frame. A joint-coupled fit adds two outer-grid
placement columns, reported whether or not the grid moved so an inert
auto-recenter is visible in a batch summary:

- `outer_grid_placement`:

  `"fixed"` (the grid the fit was given) or `"auto_recentered"` (the
  engine moved it onto the hyperparameter mode).

- `outer_grid_recenter_declined`:

  On a `"fixed"` placement, why the recenter did not apply – for
  instance a user-pinned axis, or `"auto_recenter_disabled"`. `NA` when
  the grid was recentered.

A fit that also requested the Pareto-k diagnostic
(`control$diagnose.k = TRUE`, off by default) adds three more:

- `pareto_k`:

  The outer importance-sampling \\\hat{k}\\ for the hyperparameter
  Gaussian summary; `< 0.7` indicates a reliable summary, `NA` when the
  diagnostic did not run or the proposal was degenerate.

- `pareto_k_is_ess`:

  The importance-sampling effective sample size on the PSIS-smoothed
  weights (numeric); `pareto_k_is_ess / control$k.samples` is the
  relative IS efficiency.

- `pareto_k_proposal_source`:

  How the importance proposal was built: `"mode_hessian"` from the
  Laplace curvature at the hyperparameter mode – curvature-backed, so
  the \\\hat{k}\\ stays trustworthy even when a sharp posterior
  collapses the integration grid to ~1 cell; `"grid_moment"` from the
  grid-weighted covariance of the integration nodes (with
  `"moment_matched"` its refinement); or `"grid_mixture"`, the
  local-bump-per-cell mixture matching what the engine samples on a
  spread grid.
