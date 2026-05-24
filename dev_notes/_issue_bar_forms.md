# lme4 bar syntax: slope-only / multi-slope / nested (gcol33/tulpaObs#10)

**Status: RESOLVED.** All three previously-rejected bar forms are implemented
and covered by recovery + parse tests.

## What `.tobs_desugar_bars()` now maps (R/formula_parse.R)

```
(1 | g)          -> re(g, type = "intercept")
(x | g)          -> re(g, type = "slope", covariate = cbind(x))
(1 + x + z | g)  -> re(g, type = "slope", covariate = cbind(x, z))
(x || g)         -> re(g, type = "slope", covariate = cbind(x), correlated = FALSE)
(0 + x | g)      -> re(g, type = "slope", covariate = cbind(x), intercept = FALSE)
(1 | g:h)        -> re(interaction(g, h, drop = TRUE), type = "intercept")
(1 | g/h)        -> re(g) + re(interaction(g, h, drop = TRUE))   [LHS distributed]
```

## Where the work landed

* `re()` takes a covariate matrix / multiple names / one-sided formula, plus an
  `intercept =` flag (R/formula_terms.R).
* `build_re_spec()` derives `n_coefs` from the slope-matrix width and passes a
  per-term `re_has_intercept` vector (R/occu_fit.R).
* `.tobs_bar_group_terms()` / `.tobs_flatten_op()` / `.tobs_interaction_call()`
  expand crossed/nested grouping (R/formula_parse.R).
* Cholesky size corrected `k*(k+1)/2` -> `k*(k-1)/2` (src/populate_helpers.h)
  to match tulpa's tanh-Cholesky prior.

## Engine note (correction to the original issue's scope)

The issue assumed no tulpa change was needed. Two things required it:
1. The correlated Cholesky count was wrong on the tulpaObs side (above).
2. Slope-only needed a per-term `re_has_intercept` flag in the engine
   (tulpa ABI 21 -> 22): `slope_at()` / `obs_re_contrib()` (laplace_spec.cpp)
   and the autodiff RE contribution (log_post_generic_impl.h) hardcoded the
   implicit intercept. See tulpa NEWS.

Random slopes fit under NUTS only (the EM-Laplace path drops formula random
effects; nested-Laplace rejects slopes). Recovery tests use NUTS:
`tests/testthat/test-re-bar-recovery.R`.
