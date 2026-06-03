# tulpaObs NEWS

## 0.0.12 (2026-06-03)

* docs: clean up two roxygen warnings surfaced on `R CMD Rd2pdf` / `document()`.
  `simulate_cover()` had a bare `%*%` in its generative-model `\describe` block
  (the `%` opened an Rd comment and mangled the `\item` entries); the criteria
  page (`tobs_waic()` / `tobs_dic()` / `tobs_cpo()`) linked to the unexported
  `.tobs_loglik_at_mean()`, which has no Rd topic. Both now render as inline
  code.

## 0.0.11 (2026-06-03)

* Bump the tulpa dependency to (>= 0.0.9) / gcol33/tulpa@v0.0.9, which carries
  the finite-guarded outer-grid weight normalisation (gcol33/tulpa#65). The
  defensive softmax fallback in the joint-coupled fitter is retained, but the
  upstream NaN-weight path it guards against is now fixed in the engine itself.

## 0.0.10 (2026-06-03)

* fix(occu_cover): the beta latent cover spec honours the engine's Expected
  (Fisher) curvature request, so the outer-grid corner cells converge instead of
  returning a non-finite `log_marginal` (tulpaObs#35). The latent marginal's
  observed information `E_pi[sneg] - Var_pi(s)` can go indefinite at extreme
  hyperparameter cells (large `sigma_u` driving the beta mean toward 0/1), where
  it stalled the inner Newton; under `hessian = "fisher"` (the beta default) the
  inner step now uses the always-positive marginal Fisher information
  `E_pi[sum_j fisher_beta(eta + u)]`. The reported `log_marginal` / SEs are
  unchanged (the final mode-pass re-evaluates with observed curvature); the fix
  recovers the ~20% of grid mass previously discarded at the corners and removes
  the fragile NaN path for beta latent. The lognormal latent path is exactly
  quadratic (observed == expected) and is unchanged. Ground-truth check (Expected
  curvature == the brute-force Fisher marginal integral, PSD) in
  `test-occu-cover-latent.R`.

* fix(occu_cover): GOF (`tobs_waic` / `tobs_dic` / `tobs_cpo`) now scores the
  cover term at the granularity the fit used (tulpaObs#34). The pointwise
  log-likelihood evaluated the positive-arm density at every detected visit even
  for a `cover_aggregate = "mean"` / `"median"` / `"latent"` fit, so it scored a
  likelihood the model was never fit to: `p_waic` grew super-linearly in the
  visits-per-site count and the LOO Pareto-k fraction went pathological. The
  cover term is now evaluated once per occupancy unit at the aggregated cover
  (mean / median) or via the per-unit cover-RE marginal (latent), matching the
  fitter. The pointwise log-likelihood also reads the dispersion the spec held
  fixed (`sigma_pos` / beta precision) instead of a bare unit default, so the
  WAIC / DIC / LOO of every spatial `occu_cover()` fit is on the fitted
  dispersion scale. Regression test in `test-occu-cover-aggregate.R`.

* fix(occu_cover): `predict()` / WAIC grid sampling no longer fails when the
  outer grid has non-converged cells. `tulpa_posterior_draws()` samples the grid
  mixture by `fit$weights`; when a corner of the grid carries a non-finite
  `log_marginal` (e.g. the beta latent spec's Gauss-Hermite arm not converging)
  the engine's normalized weights could collapse to all-zero, leaving the sampler
  with no positive-weight cell ("nothing to sample"). The joint-coupled fitter now
  falls back to the same finite-cell softmax weights it uses for the reported
  posterior moments, so `predict()` and WAIC stay consistent with the point
  estimates. Covered by the beta latent predict test in
  `test-occu-cover-latent.R`.

* feat(occu_cover): latent cover-per-unit (`cover_aggregate = "latent"`). The
  principled counterpart of mean / median aggregation: instead of collapsing a
  unit's detected covers to one number, the cover arm carries a per-unit cover
  random effect `u_i ~ N(0, sigma_u^2)` shared across the unit's detected visits
  and integrates it out per unit, keeping every detected visit. The lognormal
  arm integrates in closed form (compound-symmetry sufficient statistics); the
  beta arm uses adaptive Gauss-Hermite quadrature (`control$n.quad`, default 15).
  Because the cover predictor is unit-level the per-unit marginal is a scalar
  function of one eta, so it slots into the one-row-per-unit layout with no
  within-arm Hessian coupling. The within-unit dispersion is pre-fit from the
  within-unit spread and held fixed; `sigma_u` is integrated on the outer grid
  (reported as the `phi_pos` hyperparameter; `control$sigma.u.grid`). Same gates
  as aggregation (cell-level positive design, shared-field spatial path,
  `joint_coupled` engine). New stateful `_latent` cell-coupling specs
  (FD-checked vs brute-force numerical integration in `test-occu-cover-latent.R`);
  `sigma_u` + field + coefficient recovery in the same file.

* feat(occu_cover): cell-aggregated cover (`cover_aggregate`, tulpaObs#33). On
  the shared-field spatial path the cover arm carried one observation per
  detected visit, so a cell with many detected plots drove the shared ICAR field
  far more than the single occupancy observation for that cell and the
  detection-corrected occupancy surface flattened. `cover_aggregate = "mean"`
  (the new default on the spatial path) / `"median"` collapses the cover arm to
  one response per occupancy unit (the mean / median cover over its detected
  visits), so occupancy and cover inform the field with comparable weight;
  `"none"` keeps the per-visit arm. Aggregation reads a cell-level positive
  design from `data`; a visit-level `positive` covariate keeps the per-visit arm
  (the bare default falls back, an explicit request errors). New `_agg`
  cell-coupling specs evaluate the cover density once per cell (FD-checked in
  `test-occu-cover-coupling.R`); recovery + gates in
  `test-occu-cover-aggregate.R`.

* fix(occu_cover): regularise the cover (pos) arm intercept by default on the
  joint spatial path (tulpaObs#32). The cover intercept was left at the engine's
  flat 1e-4 ridge while occupancy / detection carried the `occu_priors()`
  defaults, so on a shared field it could float along the field-level confound to
  a huge posterior SD (occupancy stayed tight) and blow up `predict()`'s
  conditional cover via Jensen. It now carries the weakly-informative
  `cover_priors()` intercept prior by default, like the load-bearing detection
  prior; `priors = FALSE` / `"none"` still disables it.

## 0.0.9 (2026-06-02)

* feat(occu_aggregation_scan): suggest a spatial cell size and yearly clustering
  that make a single-season occupancy model identifiable. Single-visit plot data
  carries no within-unit replication, so psi and p are confounded until records
  are pooled into (cell, year-block) buckets; the scan scores candidate (cell
  size x year block) pairs by structural replication (`"count"`) or the curvature
  of the constant-model occupancy likelihood (`"info"`, smallest eigenvalue /
  posterior SEs of the 2x2 (logit psi, logit p) information), with a `plot()`
  method.

* feat(occu_cover, cover): forward the `integration` control key to the tulpa
  joint backend (gcol33/tulpaObs#31). `control$integration = "ccd"` / `"grid"`
  now reaches `tulpa_nested_laplace_joint()`, so the coupled cover-hurdle fit
  can select CCD outer integration over the latent + phi hyperparameter axes.

* feat(ms_dyn_occu, ms_int_occu): community (multispecies) dynamic and integrated
  occupancy families. `ms_dyn_occu()` is the community version of `dyn_occu()`
  (per-species first-season occupancy + detection coefficient random effects,
  shared community-wide colonisation / extinction); `ms_int_occu()` is the
  community version of `int_occu()` (multiple detection sources share one latent
  occupancy state per species, per-species occupancy + per-source detection RE).
  Both fit by a shared community Laplace-EM (`R/community_em.R`): the latent state
  marginalizes in closed form (HMM forward for dynamic, two-state mixture for
  integrated), the per-species coefficient deviations are integrated by a
  joint-Newton mode-find with the RE blocks Schur-folded, and a closed-form
  M-step updates the per-arm community covariance. Parameter-recovery + 95% CI
  coverage tests (`test-ms-dyn-occu.R`, `test-ms-int-occu.R`).

* fix(ms_occu): community single-season occupancy is now a correct community
  model (gcol33/tulpaObs#30). The previous `ms_occu()` did not fit per-species
  random effects in either engine: the Laplace route collapsed to a pooled GLM
  over the stacked species rows (no species RE), and the NUTS route forced one
  shared species intercept onto both the occupancy and detection arms. `ms_occu()`
  now uses the shared community Laplace-EM with independent per-arm Gaussian
  community covariances (the spOccupancy `msPGOcc` model), recovering the
  per-species occupancy and detection coefficients its own `simulate_ms_occu()`
  generates. `ranef()`, per-species `fitted()`, and `tobs_richness()` read the
  per-species structure; recovery + coverage tests in `test-ms-occu.R`. The
  legacy generic-engine community path (`build_community_callbacks`,
  `.tobs_build_community`, the `community` model_type in `src/occu_fit.cpp`, and
  the community entries in the Laplace / nested-Laplace switches) is removed.
  `ms_occu()` is Laplace-only; a correct community NUTS / areal-spatial path
  needs independent per-arm RE blocks in the sampler and is deferred.

## 0.0.8 (2026-06-02)

* feat(occu_multiscale_cover): three-level occupancy + cover hurdle family
  (gcol33/tulpaObs#29). A cell-level occupancy gate (`psi`), a plot-level
  availability gate (`theta`), per-visit detection (`p`) and the cover hurdle
  (`pos`), for vegetation data where a site's "visits" are spatially distinct
  plots aggregated into a `(cell, period)` rather than temporal revisits
  (Nichols et al. 2008; Mordecai et al. 2011). Where `occu_cover()` conflates
  within-cell prevalence into the detection arm (Kendall & White 2009), the
  explicit middle level separates them. Both `z` (over cells) and `a` (over
  plots) marginalize in closed form, so the joint marginal log-likelihood is
  exact -- a new four-arm cell-coupling spec
  (`src/cell_coupling_occu_multiscale_cover.h`, the nested two-state mixture)
  drives `tulpa::tulpa_nested_laplace_joint()` over the shared `(sigma, alpha)`
  field grid. Spatial-only (`method = "nested_laplace"`); a single shared areal
  field. The no-detection occupancy-mixture math is now shared with
  `occu_cover()` via `nodet_mixture_block` (`src/occu_coupling_shared.h`). Adds
  `occu_multiscale_cover()` and `simulate_occu_multiscale_cover()`. Tests:
  `test-occu-multiscale-cover-coupling.R` (FD-checks every closed-form
  derivative, both families, branches A/B, Expected curvature, 2-level
  reduction), `test-occu-multiscale-cover-recovery.R` (parameter recovery + CI
  coverage + field shape). `fitted()` / `predict()` for the family are pending.

* feat(ms_occu_cover): community (multispecies) joint occupancy-detection +
  cover family, the community version of `occu_cover()`. Per-species coefficient
  random effects with Gaussian community hyperpriors on all three arms
  (occupancy `psi`, detection `p`, positive cover), so rare species borrow
  strength from common ones through the shared community means and covariances.
  The latent presence `z` integrates out per species-cell in closed form (the
  same two-state mixture as `occu_cover()`); the per-species deviations are
  integrated by a Laplace-EM (arrowhead joint Newton with the per-species RE
  blocks Schur-folded, closed-form community-covariance M-step, Louis 1982
  Schur-complement community-mean SEs). Beta + lognormal positive arms,
  non-spatial Laplace only (`method = "nested_laplace"` errors: the per-group RE
  on a shared coupled field needs upstream engine support). Adds
  `ms_occu_cover()` and `simulate_ms_occu_cover()`. Tests: `test-ms-occu-cover.R`
  (community-mean recovery + 15-seed CI coverage + per-species coefficient
  recovery). Status `"experimental"`; NUTS / negbin / per-species dispersion RE /
  AGHQ variance-component debias pending.

## 0.0.7 (2026-06-01)

Requires tulpa (>= 0.0.7) and tulpaMesh (>= 0.1.2).

* feat(occu_cover/cover): forward grid-cell checkpoint/resume into the joint
  nested-Laplace engine (gcol33/tulpa#50). `control$checkpoint = list(path =,
  resume =)` is passed verbatim to `tulpa::tulpa_nested_laplace_joint()` from
  both the `occu_cover()` joint-coupled path and the `cover()` hurdle path, so a
  full-field fit killed by a reboot or OOM resumes from the last completed outer
  grid cell instead of restarting. `"checkpoint"` is on the `occu_cover` + `cover`
  control allowlist and documented as a Checkpoint/resume section on both
  families. Tests: `test-occu-cover-checkpoint.R`.

## 0.0.6 (2026-06-01)

Requires tulpa (>= 0.0.6) and tulpaMesh (>= 0.1.2).

* perf(occu_cover): forward `control$diagnose.k` / `control$k.samples` to the
  joint engine. The outer Pareto-k diagnostic re-solves the inner Laplace at
  `k.samples` sampled hyperparameters; at field scale the draws stall at extreme
  sigma and the diagnostic costs ~50x the grid integration while returning NA
  for the multi-block ICAR config (gcol33/tulpa#51). `control$diagnose.k = FALSE`
  skips it for a production fit; small fits keep the engine default.
* feat(occu_cover): `group_var` on the `icar()` / `bym2()` term decouples the
  occupancy units (sites) from the field nodes (cells), so many sites share one
  areal field node. A site = cell x time-period then carries a per-site trend
  weight, giving a detection-corrected occupancy trend on a shared cell field
  (the field stays length n_cells while psi / p / cover run over n_sites).
  Joint_coupled engine only; the v2/v3 escape hatches reject it. Recovery in
  `test-occu-cover-group-var.R`.
* feat(cover): the joint cover-hurdle predict substrate now handles the coupled
  multi-block case (an ICAR intercept field plus one or more SVC trend fields)
  under the per-block `(sigma, alpha)` copy convention -- the occupancy arm
  scales block `k` by `sigma`, the positive arm by `alpha * sigma`
  (gcol33/tulpaObs#15).
* feat(cover): `cover()` accepts the `trend`, `alpha.grid`, and
  `alpha.grid.trend` control knobs for the trend-field integration.

## 0.0.5 (2026-06-01)

Requires tulpa (>= 0.0.5): the joint cover-hurdle path links against the
engine's `cell_coupling.h` / `model_data.h` (ABI 32).

* perf(occu_cover): speed up the beta positive arm in the joint cover-hurdle
  cell-coupling spec (gcol33/tulpa#46, lever 3). The per-observation
  `digamma`/`trigamma` now use tulpa's portable, inlinable, OpenMP-safe
  implementations instead of the `R::` math-library calls; the score and
  curvature share their `digamma` terms in a single pass; and the spec honours
  the engine's `CellDerivs::grad_only` request, skipping the `trigamma` entirely
  on a factor-reuse inner-Newton step. The `joint_coupled` engine now also
  exposes `control$n.threads.outer` (the engine's parallel sparse outer grid)
  and `control$force.sparse`. End-to-end invariance of the beta cover fit to
  `inner.refresh` and `n.threads.outer` is covered in
  `tests/testthat/test-occu-cover-joint-reuse.R`.
* feat(occu_cover): joint occupancy-detection + cover-hurdle family
  (`occu_cover("lognormal")` / `occu_cover("beta")`, gcol33/tulpa#32). A site's
  occupancy/detection arm and a positive-cover arm (lognormal or beta) are fit
  jointly, with the two linear predictors sharing a spatial field through a
  cell-coupling spec. `method = "nested_laplace"` routes through the
  `joint_coupled` engine by default; the engine integrates the shared field's
  hyperparameters and the coupling coefficient on the outer grid. `coef()`,
  `vcov()`, `ranef()`, `fitted()` and `simulate()` carry both arms.

* feat(occu_cover): `predict()` for the joint fit (gcol33/tulpaObs#22). Samples
  the grid-integrated joint latent via `tulpa::tulpa_posterior_draws()` and
  marginalizes each derived quantity per draw, returning a `tobs_prediction`
  with per-unit draw matrices and the change-column contract (`delta_p`,
  `delta_cover_cond`, `delta_cover_exp`) with `.lwr` / `.upr` at the requested
  level.

* feat(ms_abun): per-species negative-binomial dispersion (gcol33/tulpaObs#14).
  `ms_abun(mixture = "negbin")` gives each species its own overdispersion
  `log r_s ~ N(mu_log_r, sigma_log_r)`, partially pooled across the assemblage;
  `fit$ms_dispersion` reports the per-species `r_s` with the community
  `mu_log_r` / `sigma_log_r`, and `ranef()` carries a `logr` arm.

* feat(spde): native SPDE (Matern-via-mesh) continuous spatial fields on the
  occupancy state and detection arms and on the single-species / community
  N-mixture abundance arm, with the (range, sigma) hyperparameters integrated on
  the outer grid under a PC prior. Requires `tulpaMesh` for mesh construction.

* feat(ms_abun): opt-in exact-Newton inner solver for the areal shared-field
  community N-mixture (`control$inner_solver = "newton"`, default `"em"`,
  gcol33/tulpaObs#12) -- an accuracy/validation alternative to the EM M-step on
  the same outer field-hyperparameter grid; both return the same `tobs_fit`
  shape and `ms_community$optimizer` records which ran.

* feat(abun): grouped random effects on the single-species N-mixture
  (gcol33/tulpaObs#13). `abun()` now accepts `(1 | g)` (and the slope /
  uncorrelated / correlated variants the formula RE machinery already
  understands) on either the abundance or the detection predictor of a
  single-species fit. `.tobs_fit_nmix_re()` warm-starts the betas with the
  no-RE Laplace fit and refines through `.tobs_nmix_re_aghq()`, a thin
  `make_site` callback over `nmix_site_marginal()` driven by
  `tulpa::tulpa_re_aghq()`. `control$n.quad = 1` is the joint Laplace
  (production); `n.quad > 1` is the AGHQ debias of the small-cluster sigma
  attenuation. Poisson and negbin (the global `log_r` is carried as the
  trailing theta coordinate, jointly estimated with the betas). Gated: RE
  combined with an areal spatial term, RE with visit-level detection
  covariates, and RE shared across both arms (each errors with a pointer
  rather than silently dropping the requested structure). `coef()`,
  `vcov()`, and `ranef()` carry the RE component. `tests/testthat/test-abun-re.R`
  covers the structural surface, the capability gates, and (under
  `NOT_CRAN=true`) sigma recovery, fixed-effect CI coverage, and the NB+RE
  path.

## 0.0.2 (2026-05-28)

* fix(build): clean-slate compile against tulpa restored. The
  in-tree N-mixture move (commit c8b6912) left two casing mismatches
  (`using tulpaObs::NMix_*` for functions defined as `tulpaObs::nmix_*`) in
  `src/nmix_spatial.cpp` / `src/nmix_spatial_bym2.cpp`, and a broken include
  guard in `src/nmix_spatial_assemble.h` (`#ifndef TULPA_NMIX_SPATIAL_ASSEMBLE_H`
  vs `#define TULPAOBS_NMIX_SPATIAL_ASSEMBLE_H`) that re-included the file and
  redefined its templates on the second pass. All three fixed; cold parallel
  build is ~13 s (`R CMD INSTALL -j8`, rtools45). Requires tulpa >= 0.0.3.

* refactor: the N-mixture observation model (single-species, areal-spatial,
  community) has moved from the tulpa engine into tulpaObs as a
  consumer-side `LikelihoodSpec`, restoring the principled
  engine/model-package boundary. No user-facing API change: `abun()` and
  `ms_abun()` continue to be the public surface. The native
  `NMixCommunityOracle` is now a tulpaObs `XPtr<tulpa::REGroupOracle>` that
  reaches the engine through `<tulpa/aghq_oracle.h>`; the community fit
  drives it through `tulpa::tulpa_re_aghq()`. Bumps the tulpa pin to
  v0.0.3.

* feat(ms_abun): community / multispecies N-mixture (`ms_abun()`, the
  spAbundance `msNMix` model) now fits under `method = "laplace"`. Per-species
  abundance and detection coefficients are random effects with Gaussian
  community hyperpriors (`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`,
  `beta_p_s ~ N(mu_p, Sigma_p)`); the latent abundances integrate out in closed
  form per species-site, and the fit is a C++ Laplace-EM (`nmix_laplace_re()`,
  driving a native `NMixCommunityOracle`) -- per-species coefficient
  modes, a closed-form covariance M-step, and fixed-effect SEs from the marginal
  observed-information Schur complement (with the `Var[N|y]` rank-1 correction).
  `y` is a 3D array `[n_sites x max_visits x n_species]` or a named list of
  count matrices; pass `species =`. `coef()` returns the community means;
  `ranef()` the per-species coefficient deviations; `vcov()` / `confint()` the
  community-mean covariance; `fitted()` / `simulate()` the per-species
  `lambda` / `p` / counts. `simulate_ms_abun()` +
  `tests/testthat/test-ms-abun.R` cover community-mean recovery, 95% CI
  coverage over 20 seeds, per-species coefficient recovery, and the S3 surface.
  Poisson only for now (a global negative-binomial size and an areal-spatial
  community field are upstream-pending). Requires tulpa >= 0.0.2.

* feat(abun): `abun(mixture = "negbin")` now fits. The negative-binomial
  abundance mixture (`Var(N) = lambda + lambda^2 / r`) is wired through to
  tulpa's N-mixture kernel (`mixture = "NB"`) on both the non-spatial Laplace
  path and the areal-spatial nested-Laplace path (`icar()` / `bym2()` /
  `car_proper()`). Non-spatially the log size `log_r` is estimated jointly with
  the coefficients and reported with a standard error (the trailing `vcov`
  coordinate); spatially the size `r` is integrated over the outer
  hyperparameter grid and reported as a posterior mean / sd. `coef()` /
  `confint()` / `vcov()` report the abundance and detection arms; the dispersion
  is surfaced via `fit$nmix_dispersion`. `simulate()` and `simulate_abun()` draw
  `N ~ NegBin(mu, r)` under the NB mixture. Matches `unmarked::pcount(mixture =
  "NB")`; `test-abun.R` adds NB point recovery, dispersion recovery, 95% CI
  coverage, and a spatial NB fit. Requires tulpa >= 0.0.2.

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
