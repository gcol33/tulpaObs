# tulpaObs NEWS

## Unreleased

* refactor(re): the AGHQ quadrature engine moved to tulpa (`tulpa::tulpa_re_aghq`,
  requires tulpa >= 0.0.2). `R/re_aghq.R` previously reimplemented the
  Gauss-Hermite nodes, the log-Cholesky covariance parametrization, and the LKJ
  penalty -- all of which tulpa already owns -- which duplicated generic
  inference machinery across the package boundary. It is now a thin wrapper that
  supplies only the family-specific occupancy / detection per-site marginal (the
  `make_site` callback) and delegates the quadrature, mode-finding, joint
  optimization, and marginal Hessian to the engine. Behaviour is unchanged
  (the recovery numbers are identical); the engine is reusable by other families
  and is recovery-tested in tulpa. A failed refine now warns (it kept the EM
  result silently before).

* feat(re): random effects on the **detection** predictor now fit under
  `method = "laplace"` (gcol33/tulpaObs#11 follow-up). Previously a `(1 | g)` /
  `(x | g)` term on the detection formula errored toward NUTS; the
  variance-component EM (`R/em_laplace_re.R`) now partitions the RE terms into an
  occupancy and a detection arm (`.tobs_re_split_arms()`), fits each arm's RE
  block in its own M-step (occupancy as the M-inflated pseudo-binomial, detection
  as a genuine weighted binomial at natural prior scale), and folds both RE
  posterior modes back into psi and p at the E-step. The AGHQ debias
  (`R/re_aghq.R`) is generalized to either arm: when the random effect enters the
  detection predictor `b` moves `p` (not `psi`), so the per-group marginal uses
  the binomial-in-`p` site derivatives (finite-difference verified,
  `dev_notes/probe_re_det_aghq_deriv.R`). Detection-arm RE parameters are named
  for the detection process (`sigma_p<t>_*`, `re_p<t>_*`). Measured recovery
  (`dev_notes/probe_re_det_*.R`, 24 seeds, per-group n ~ 10): the raw `nAGQ = 1`
  EM attenuates the detection `sigma` by ~70% (only occupied sites inform `p`);
  the AGHQ refine restores it to within ~1% of truth (0.805 at truth 0.80) with
  88-96% fixed-effect CI coverage. Supports the same iid-intercept /
  uncorrelated-slope / correlated-slope forms as the occupancy arm. A single RE
  shared across BOTH predictors, RE + visit-level detection, RE + spatial, and
  non-single families still point to NUTS.

* feat(re): adaptive Gauss-Hermite (AGHQ) debias of the random-effect variance
  components under `method = "laplace"`, on by default (`R/re_aghq.R`,
  `.tobs_re_aghq()`). The variance-component EM integrates the RE block `b` by
  Laplace (the glmer `nAGQ = 1` regime), which attenuates `sigma` / the RE
  correlation toward zero for binary occupancy at small per-group sample size.
  After the EM converges, the per-group marginal `int prod_i L_i(eta_i + Z_i b)
  N(b; 0, Sigma) db` -- reusing the exact closed-form occupancy site marginal
  (`z` integrated out) -- is refined by `n.quad`-point adaptive Gauss-Hermite
  quadrature centred at the EM mode, and `(beta, chol Sigma)` are re-optimized on
  it; the fixed-effect SEs are read from the exact-marginal Hessian. Measured
  recovery (`dev_notes/probe_re_aghq*.R`): at per-group `n = 8` the EM attenuates
  `sigma` by ~18% (bias -0.16 at truth 0.9), AGHQ cuts that to ~4% (matching
  NUTS); correlated-slope `sigma`s recover to ~1% on the seed average. Controls:
  `re.aghq` (default `TRUE`; `FALSE` keeps the raw `nAGQ = 1` EM), `n.quad`
  (default 9; `n.quad = 1` is the plain Laplace marginal), and `re.lkj` (default
  1.5). Applies to a single grouping factor with RE dimension <= 3 (the
  recovery-tested forms); crossed / nested groupings fall back to the EM. A
  weakly-identified *correlated* random slope's correlation is regularized off
  the `+-1` boundary by a default LKJ(`re.lkj = 1.5`) penalty on the block's
  correlation matrix -- `(eta - 1) log det R`, maximized at independence,
  leaving the marginal SDs untouched and `O(1)` against the `O(n_groups)`
  likelihood; `re.lkj = 1` disables it. On the recovery sim (truth rho = 0.61)
  this removes every `+-1` boundary hit while keeping rho near-unbiased (bias
  -0.00 at per-group n = 12, +0.02 at n = 25), where unregularized ML
  over-estimates rho and hits the boundary. The fit is then a MAP and the
  reported SEs come from the penalized (posterior-precision) Hessian. NUTS
  remains available for a full posterior treatment of the correlation. Recovery
  tests in `tests/testthat/test-re-laplace-recovery.R`.

* feat(nested): calibrated credible intervals for NA-response prediction.
  `predict(fit, type = "state")` now returns `psi_lower` / `psi_upper` (95%)
  alongside `psi` for single-season nested-Laplace fits. The intervals are
  calibrated by refining the EM field with one exact-marginal pass
  (`.tobs_occu_state_marginal_fit()`): integrating out the latent occupancy
  state makes each site a Bernoulli on `D = 1{>=1 detection}` with mean
  `q * plogis(eta)`, where `q = 1 - (1 - p)^J` is the per-site detection
  probability (a held-out site has no visits, so `q = 0` and it is interpolated
  by the field). This is fit through tulpa's generic `family = "bernoulli"` with
  a per-observation probability scale `det_prob = q`, so the converged Hessian is
  the marginal curvature and the per-cell predictive variance
  (`tulpa_nested_laplace()$fitted_eta_var`) is calibrated directly. The per-row
  eta posterior is a Gaussian mixture over the hyperparameter grid; `psi` is its
  Gauss-Hermite mean and the interval is the mixture-CDF quantiles, both per the
  marginalise-derived-quantities rule. The EM's M-inflated pseudo-binomial
  (which weights the data ~M times the prior and whose unit-trial Hessian is the
  complete-data information) is kept only for the mode and detection estimate;
  reading variance or grid weights off it under-covers and collapses the
  hyperparameter grid. Held-out coverage measures ~1.0 (calibrated, slightly
  conservative) with cor ~0.88 / MAE ~0.11 on a 10x10 icar/bym2 grid; recovery
  test in `tests/testthat/test-nested-laplace-prediction.R`. Older tulpa without
  `fitted_eta_var` reports `psi` with `NA` interval columns. Other model types
  keep the EM occ fit (their NA-response mapping is not yet wired).

* feat(nested): nested Laplace generalised beyond single-season occupancy.
  `method = "nested_laplace"` now fits `int_occu()`, `ms_occu()`, and
  `dyn_occu()` as well as `occu()` -- the multi-block latent prior (spatial /
  temporal / iid) is attached to the state ("occ") M-step block of the same
  per-model-type callbacks the single-Laplace path uses, so there is one set of
  callbacks rather than a `build_*_callbacks_nested` duplicate. The block's
  per-row `spatial_idx` maps each state row to its site, so a community model
  shares one site-level field across the species at a site, and integrated /
  dynamic models carry a spatial / temporal field on the shared psi / psi1
  predictor. The driver `.tobs_em_nested_laplace()` is now a thin wrapper over
  `.tobs_laplace(latent_prior = )`. The registry (`.tobs_family_methods`) lists
  `nested_laplace` for these families; calibrated recovery for the multi-block
  engine remains the same follow-up tracked for single-season occupancy. Smoke
  tests in `tests/testthat/test-nested-laplace-families.R`.

* feat(nested): INLA-style NA-response prediction. A single-season occupancy
  site whose detection history is all-missing (all `NA`) is held out of the
  likelihood (`n_trials = 0`) but kept in the latent field, so its occupancy is
  interpolated from neighbours (`.tobs_heldout_sites()`). `predict(fit, type =
  "state")` returns the marginalised per-site psi posterior -- the weighted mean
  over the hyperparameter grid of `plogis(eta)`, integrated rather than plugged
  in at the mode -- with the held-out rows flagged. The E-step is now
  field-aware (`psi_i = plogis(X_i beta + field[idx(i)])`) so the field tracks
  the data instead of converging to the fixed-effect-only fixed point. The
  marginalised psi is read from the engine's per-cell fitted linear predictor
  (`tulpa_nested_laplace()$fitted_eta`), so it is exact for every latent prior
  -- including `bym2`, whose predictor mixes structured and unstructured
  components with hyperparameter-dependent scales that the engine reconstructs
  with the right `d_fac` (older tulpa without `fitted_eta` falls back to
  mode reconstruction, exact for the d_fac = 1 priors only). Calibrated
  predictive intervals need the latent field's per-cell posterior variance and
  are deferred to engine support -- only the marginalised mean is reported.
  Recovery tests (icar + bym2) in
  `tests/testthat/test-nested-laplace-prediction.R`.

* fix(api): backend coverage is now enforced from a single source of truth.
  `.tobs_family_methods` (in `R/tobs.R`) declares the `method`s each working
  family supports, and `tobs()` validates the resolved method against it,
  erroring with a pointer to the supported set. This removes the silent
  `nested_laplace` -> single-Laplace downgrade that `dyn_occu()` / `ms_occu()` /
  `int_occu()` / `jsdm()` previously hit (`.map_engine()` emitted only a
  `message()` and then stamped `fit$method <- "nested_laplace"` on a fit that was
  actually single-Laplace -- a provenance bug). The nested-Laplace engine is
  wired only for single-season occupancy and the cover hurdle joint path; the
  cover hurdle has no NUTS likelihood or EM-correction engine, so its
  `nuts` / `laplace_gibbs` / `laplace_mi` rejections (previously scattered across
  `.dispatch_cover()`) now flow through the same central check. Tests in
  `tests/testthat/test-family-method-coverage.R`.

* feat(formula): single-verb `spatial(..., model = ...)` umbrella over the
  areal (`icar`, `bym2`, `car`, `car_proper`) and continuous (`gp`,
  `multiscale_gp`, `spde`) spatial terms, mirroring `temporal(time, type = ...)`
  and `INLA`'s `f(i, model = ...)`. `spatial(graph = adj, model = "bym2")` is
  identical to `bym2(graph = adj)` and `spatial(lon, lat, model = "spde")` to
  `spde(lon, lat)`; the specific constructors still work. Dispatches through the
  term registry (single source of truth), forwarding coords / `graph =` /
  per-model arguments and `id` unchanged. Named arguments are validated against
  the target constructor's formals, so a typo'd or wrong-model argument
  (`spatial(lon, lat, model = "gp", graph = adj)`) errors at the call site
  rather than being silently absorbed as a coordinate by the continuous terms'
  `...`.

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
  is dropped in favour of the engine's `cov_blocks`. The variance components
  still carry the Laplace approximation's small-cluster bias for binary data
  (the glmer nAGQ=1 regime, not Breslow-Clayton PQL) -- NUTS is the calibrated
  route for the covariance. Recovery
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
  `(1 + x || g)`) on the occupancy predictor of a single-season model. The
  variance-component estimates carry the Laplace approximation's small-cluster
  bias for binary data (the glmer nAGQ=1 regime, not Breslow-Clayton PQL);
  `engine = "nuts"` remains the calibrated route for the covariance.
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
