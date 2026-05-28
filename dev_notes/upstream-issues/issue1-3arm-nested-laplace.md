## Summary

Extend `tulpa_nested_laplace_joint()` (currently a 2-arm engine for tulpaObs's cover hurdle) to a generic N-arm joint nested-Laplace engine, so likelihoods with more than two linear predictors can use the same C++ machinery instead of pure-R outer/inner loops.

## Motivation

tulpaObs landed a new family `occu_cover()` (2026-05-28, see `gcol33/tulpaObs/R/occu_cover_nested.R`) that has **three** linear predictors:

- `psi` (cell-level latent state)
- `p` (per-visit detection)
- `pos` (per-visit positive cover, on a shared cell-level field with scaling `alpha`)

The v3 spatial path runs a pure-R inner Newton on `z` profiled out per outer `(alpha, sigma)` candidate, with an outer BFGS over ~10 fixed-effect + hyperparameter params. It works (recovery is clean: alpha 1.11 vs truth 1.0, sigma 0.82 vs truth 1.0, field cor 0.98 at N=100, J=6, single seed) but each fit is ~2 min on the smoke and would be hours on real EVA data with ~1000+ cells.

The cover hurdle's `tulpa_nested_laplace_joint()` already does the right thing for 2 arms. Generalising it to N arms would let tulpaObs delete `R/occu_cover_nested.R`'s pure-R fitter and pick up tulpa's C++ speed.

## Proposed API sketch

```r
tulpa_nested_laplace_joint(
  arms = list(
    list(name = "psi", link = "logit", X = X_psi, y = ..., field_coef = 1.0),
    list(name = "p",   link = "logit", X_visit = X_p_visit, y_visit = ...,
         field_coef = 0.0),     # detection doesn't see the field
    list(name = "pos", link = "identity", X_visit = X_pos_visit, y_visit = y_pos,
         family = "lognormal", field_coef_param = "alpha")  # cover-arm scaling free
  ),
  field = list(type = "icar", graph = adj, scaled = TRUE),
  hyperparams = list(alpha = list(prior = "normal", mean = 0, sd = 5),
                     sigma = list(prior = "halfnormal", scale = 2)),
  outer = list(method = "grid", n.outer = NULL),   # or "laplace_at_map"
  inner = list(method = "newton", max_iter = 50, tol = 1e-6)
)
```

The current 2-arm engine collapses to N = 2 with no `field_coef_param` (the cover hurdle's alpha is the second arm's field_coef). The cleanest extension is probably "arm list + per-arm field coefficient (constant OR free hyperparam)".

## Downstream impact

- tulpaObs's `occu_cover()` v3 can drop pure-R for tulpa C++ -- expected 10-100x speedup on real EVA scale.
- Any future tulpaObs family that wants a shared latent field across >2 arms (e.g., joint occu + cover + community, or occu + cover + counts) inherits the engine for free.
- The existing 2-arm cover hurdle keeps working with the same call surface.

## Repro / current implementation pointer

- v3 pure-R fitter: `gcol33/tulpaObs/R/occu_cover_nested.R` (~350 lines).
- Recovery test: `gcol33/tulpaObs/tests/testthat/test-occu-cover-spatial.R`.
- v2 joint-Laplace alternative (for comparison; has a (z, alpha, sigma) ridge that v3 breaks): `gcol33/tulpaObs/R/occu_cover_spatial.R`.

## Priority

Nice-to-have. v3 works as pure R; the upstream lift is a speed/reuse play, not a correctness blocker.
