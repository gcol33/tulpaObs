## Summary

Investigate whether `tulpa::tulpa_re_aghq()` (currently the variance-component debias engine for community / multispecies N-mixture, plus tulpaObs's RE-AGHQ debias) can be repurposed as the **outer hyperparameter integration engine** for nested-Laplace family fits, replacing per-family bespoke grid loops.

## Motivation

`tulpa_re_aghq()` already does the right kind of thing:

- Quadrature grid over a low-dimensional parameter
- Per-grid evaluation via a `make_site` / oracle callback
- Log-Cholesky covariance parameterisation
- LKJ penalty
- Joint optimization + marginal Hessian

The community N-mixture case treats the per-species coefficient deviations as the "RE" being debiased. But the SAME pattern fits nested-Laplace hyperparameter integration:

- The "RE" is the small hyperparameter block `(sigma, alpha, ...)` instead of per-species coefficients
- The "oracle" returns the inner-Laplace log-marginal evaluated at the candidate hyperparameter
- The "quadrature grid" is the outer hyperparameter grid
- The "marginal Hessian" is exactly the hyperparameter observed-Fisher information we want

If the abstraction holds, an occu_cover() v4 outer-grid integration can be written as:

```r
tulpa::tulpa_re_aghq(
  oracle = function(theta) {
    # theta = c(alpha, log_sigma); inner Newton on z, return log_marg
    list(log_lik = log_marg, grad = ..., hess = ...)
  },
  init = c(alpha = 0, log_sigma = 0),
  n.quad = 5,
  ...
)
```

vs. the bespoke outer BFGS that the v3 fitter currently runs (`gcol33/tulpaObs/R/occu_cover_nested.R`, `.tobs_fit_occu_cover_nested()`).

## Open questions

This is more of a design question than a feature request. Worth answering before pursuing it:

1. **Does the AGHQ engine actually handle the "outer hyperparameters" use case cleanly?** The current callers (`ms_abun()` via `cpp_nmix_community_oracle()` and tulpaObs's `R/re_aghq.R`) all treat the AGHQ block as "per-group RE deviations". Hyperparameter integration is structurally different (one block, no per-group structure) -- maybe AGHQ degenerates to a trivial case (n_groups = 1) and it just works, maybe it doesn't.
2. **What does `make_site` / oracle look like for the hyperparameter case?** The current callback signature is per-site marginal; the analogue for hypers is "given this hyperparameter value, return the log marginal of the data after profiling z out".
3. **Does the LKJ prior / log-Cholesky covariance still make sense?** For hyperparameter integration, the natural priors are univariate (half-normal on sigma, normal on alpha) -- not a covariance matrix. The AGHQ engine might need a "univariate priors per coord" mode.

## If the answer is "yes, reuse AGHQ"

Then companion issue (generic outer-grid helper) becomes: "thin wrapper around `tulpa_re_aghq()` for the hyperparameter use case", and the implementation is much smaller. That's the win.

## If the answer is "no, AGHQ is genuinely a different engine"

Then the companion outer-grid helper is the right path -- build a separate outer-grid helper rather than overload AGHQ.

## Downstream impact

If reusable: meaningfully less new code for occu_cover v4 and any future nested-Laplace family. If not: at least we know the answer.

## Priority

Low. Best done as a design pass before starting either companion request -- if AGHQ already covers the outer integration case, both of those become smaller. If not, those two issues are the path.

## Related

- Companion: 3-arm nested-Laplace joint engine
- Companion: generic outer-grid integration helper
