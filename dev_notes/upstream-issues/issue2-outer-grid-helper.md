## Summary

Expose tulpa's outer hyperparameter-grid integration (currently baked into the cover hurdle's `theta_grid` machinery) as a reusable callback-driven helper, so downstream tulpaObs families can build outer-grid integration over hyperparameters without rebuilding the grid+weights+law-of-total-covariance plumbing.

## Motivation

tulpaObs has at least three places that want the same pattern:

1. **`cover()` joint engine** -- already has it, baked into `tulpa_nested_laplace_joint`. Returns `fit$joint$theta_grid` + `theta_weights` so the user can post-process.
2. **`abun()` + areal spatial** -- has its own grid integration via `nmix_laplace_{icar,bym2,car_proper}()` returning per-grid `cov_blocks`; assembled in tulpaObs's `.nmix_grid_vcov()` via the law of total covariance.
3. **`occu_cover()` v3 nested-Laplace** (just shipped 2026-05-28, `gcol33/tulpaObs/R/occu_cover_nested.R`) -- currently uses Laplace at the joint MAP for the outer hyperparameters (point + observed-Fisher SE). A proper outer grid integration would give marginal posteriors for `(alpha, sigma)` instead of just a point estimate.

Each of these reinvents the same wheel: per-grid inner fit + log-marg evaluation + weight normalisation + posterior moment assembly.

## Proposed API sketch

```r
tulpa_hyper_grid(
  hyper_specs = list(
    list(name = "alpha", prior = function(x) dnorm(x, 0, 5, log = TRUE),
         grid = c(-2, -1, 0, 1, 2)),
    list(name = "sigma", prior = function(x) dexp(x, 1, log = TRUE),
         grid = c(0.25, 0.5, 1, 2, 4))
  ),
  inner_fit = function(hypers) {
    # returns list(log_marg = , beta_mean = , beta_cov = )
  },
  combine = c("law_of_total_cov", "weighted_mean_only"),
  ...
)
# -> returns: grid points, weights, posterior beta_mean / beta_cov,
#    marginal posterior summaries per hyperparam (mean, sd, quantiles)
```

The `inner_fit` callback is the new bit -- it's whatever per-grid mode-finder the family uses (e.g. tulpa's `tulpa_nested_laplace()` for the cover hurdle, or a custom Newton for occu_cover). The grid + weights + posterior assembly is shared.

## Downstream impact

- `cover()` joint engine: refactor to use the helper (no behaviour change, deduplication).
- `abun()` spatial path: same.
- `occu_cover()` v4: gets proper outer-grid integration without rebuilding the plumbing.
- Any future family that wants outer-grid integration over hypers gets it as a one-callback drop-in.

## Repro / current implementations

- Cover hurdle grid: `gcol33/tulpa/.../tulpa_nested_laplace_joint.cpp` (and the R wrapper around `fit$joint$theta_grid`).
- N-mixture grid: `gcol33/tulpaObs/R/abun.R` (`.nmix_grid_vcov()`).
- occu_cover v3 (no grid yet, would benefit): `gcol33/tulpaObs/R/occu_cover_nested.R`.

## Priority

Medium. Not a blocker for any current family, but every new family that wants nested-Laplace pays the duplication tax until this exists. The cleanest time to do it is alongside the 3-arm engine in the companion issue (could share API surface).

## Related

- Companion request: 3-arm nested-Laplace joint engine (would naturally use this helper).
