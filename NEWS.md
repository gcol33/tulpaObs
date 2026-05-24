# tulpaObs NEWS

## Unreleased

* feat(re): formula random effects are now fit by the default `engine =
  "laplace"` instead of being silently dropped (gcol33/tulpaObs#11). A
  variance-component EM (`R/em_laplace_re.R`) wraps tulpa's fixed-sigma Laplace
  in the occupancy missing-data EM (feeding the random-effect mode back into
  psi) and an EM/REML update for the variance components, fitting iid intercept
  RE (`(1 | g)`) and uncorrelated random slopes (`(x || g)`, `(0 + x | g)`,
  `(1 + x || g)`) on the occupancy predictor of a single-season model.
  Deterministic Laplace variance estimates for binary occupancy carry the usual
  small-cluster (PQL) bias; `engine = "nuts"` remains the calibrated route.
  Forms the deterministic path cannot fit -- correlated slopes (a
  Cholesky-factored covariance, `(1 + x | g)`), random effects on the detection
  predictor, RE combined with a spatial term, RE with visit-level detection
  covariates, or RE on a non-single family -- now error with a pointer to
  `engine = "nuts"` rather than being silently dropped.

* feat(re): random-effect parameters are named and summarised. Under NUTS the
  `log_sigma` / `chol` / `z` columns are labelled (replacing `param[i]`), and
  `ranef()` reconstructs the per-group BLUPs on the natural scale
  (`b_{g,c} = sigma_c (L z_g)_c`, marginalising over the draws). The
  deterministic path reports the variance-component sigma and the Schur
  posterior standard errors. `ranef()` is re-exported from tulpa, and
  `coef()` now includes the visit-level detection coefficients (`p_visit_*`)
  that `summary()` already reported.

* feat(formula): lme4 bar syntax now supports the slope-only, multi-slope, and
  nested/crossed grouping forms (gcol33/tulpaObs#10):
  - `(1 + x + z | g)` stacks several random slopes into one correlated block
    (`re()` accepts a covariate matrix; `build_re_spec()` derives `n_coefs`
    from the slope-matrix width instead of hardcoding two).
  - `(0 + x | g)` is a slope-only block with no group intercept (threaded
    through tulpa's new `re_has_intercept` flag; requires tulpa >= ABI 22).
  - `(1 | g:h)` groups over the interaction factor; `(1 | g/h)` expands to one
    `re()` per implied grouping factor (`g`, `g:h`), distributing the LHS
    slopes across each. `||` makes the block's covariance diagonal.
  Recovery tests cover a 3x3 correlated intercept+2-slope block and a
  slope-only block under NUTS (`tests/testthat/test-re-bar-recovery.R`).

* fix(re): the correlated-slope Cholesky factor is now sized `k*(k-1)/2`
  (strictly-lower triangle) to match tulpa's tanh-Cholesky prior, replacing an
  oversized `k*(k+1)/2` that left `k` unconstrained parameters per term.

* fix(data): `tobs_data()` output composes with `tobs(visit_data = ...)`
  (gcol33/tulpaObs#8). A `det.covs` named list of `[N, J]` matrices is now
  reshaped internally to the visit-level detection design instead of erroring
  with `object '<covariate>' not found`, giving the visit-level-detection path
  a clean public route. Recovery test in
  `tests/testthat/test-issue8-visit-detection.R`.
