# tulpaObs NEWS

## Unreleased

* feat(re): correlated random slopes `(1 + x | g)` now fit under the default
  `method = "laplace"` (gcol33/tulpaObs#11), not only NUTS. The variance-
  component EM (`R/em_laplace_re.R`) carries a full per-term RE covariance
  `Sigma` (diagonal for `(x || g)`, full for `(1 + x | g)`) and updates it with
  `Sigma_k <- mean_g [b_g b_g' + Cov(b_g | y)]`, consuming the per-group
  posterior covariance from `tulpa::tulpa_laplace(return_re_cov = TRUE)`. Because
  the M-step fits the package's M-inflated pseudo-binomial (`n = M`, prior
  `Sigma/M`), the natural-scale covariance is `M` times the block tulpa returns
  (the inflated Hessian is `M` times the natural one). The estimated off-diagonal
  is reported as a `cor_<g>_<ci>_<cj>` correlation alongside the `sigma_` marginal
  SDs. This removes the previous `.validate_re_laplace()` rejection that routed
  `(1 + x | g)` to NUTS; the duplicate R-side Schur for the RE posterior variance
  is dropped in favour of the engine's `cov_blocks`. Deterministic Laplace still
  carries the small-cluster (PQL) bias -- NUTS is the calibrated route. Recovery
  test in `tests/testthat/test-re-laplace-recovery.R` (sigmas, correlation, and
  BLUPs vs simulated truth).

* feat(priors): `cover_priors()` adds opt-in Gaussian fixed-effect priors to the
  cover hurdle, penalising *both* arms -- the occurrence (Bernoulli) intercept /
  slopes and the positive-part (beta or lognormal) intercept / slopes. The
  penalty threads through `tulpa::tulpa_laplace()` for the occurrence and
  lognormal arms and through `tulpa::tulpa_laplace_beta(beta_prior = )` for the
  beta arm. Priors are opt-in (`priors = NULL`/`FALSE` fits unpenalised); an
  `sd = Inf` component is a no-op. `occu_priors()` is rejected for `cover()` with
  a pointer to `cover_priors()`, and the prior errors (no silent drop) when
  combined with a spatial term or with `method = "nested_laplace"`. Recovery
  tests in `tests/testthat/test-cover-priors.R`.

* break(api): `abun()`, `ms_abun()`, and `dyn_abun()` rename the latent-mixture
  argument `family =` to `mixture =` (`"poisson"` / `"negbin"`, after
  `unmarked::pcount()`), removing the collision with the family-object concept
  that `tobs(family = )` already owns.

* break(api): the exported `tobs_priors()` constructor (and its print method)
  are removed -- it was wired to no fitting path. Use the family-group prior
  builders `occu_priors()` (occupancy group) and `cover_priors()` (cover).

* break(control): cover-hurdle `control = list(...)` keys are renamed from
  underscores to the package's dotted convention (`max_iter` -> `max.iter`,
  `prior_sigma` -> `prior.sigma`, `prior_alpha` -> `prior.alpha`,
  `sigma_pos_grid` -> `sigma.pos.grid`, etc.), matching every other
  user-facing `control` key.

* fix(print): `print.tobs_fit` and `print.tobs_family` label the default route
  "default method" instead of "engine", matching the `method = ` argument users
  actually type.

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
