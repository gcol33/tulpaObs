# tulpaObs NEWS

## 0.0.180 (2026-08-03)

* **Pinned to `tulpa (>= 0.0.115)` (gcol33/tulpa#270).** The 0.0.114 pin
  below (needed for #169's `set_spatial_adjacency()` fix) exposed a
  separate engine bug: `occu_cover()`'s coupled ICAR joint Hessian dropped
  10592-124160 nonzero contributions whenever the field met either a
  correlated detection-arm random slope or its own detection-arm beta
  under the rank-1 s2z fold path (fields above 256 units) --
  `HessianPatternGuard` (new since 0.0.101) caught the pattern builder
  missing a cross-arm entry for a latent block reached by only one side of
  a coupled arm pair. Fixed in tulpa 0.0.115; no tulpaObs-side change
  needed.

* **ICAR/BYM2 occupancy NUTS fits work again (bugfix, #169).**
  `populate_spatial()` (the shared `ModelData` builder `cpp_occu_fit`'s
  unified NUTS entry point uses) assigned `n_spatial_units` /
  `adj_row_ptr` / `adj_col_idx` / `n_neighbors` directly for ICAR/BYM2,
  leaving `spatial_partition` at its default (0 nodes) and
  `n_spatial_components` at its default (1) -- the exact case tulpa's
  ABI 39->40 `set_spatial_adjacency()` was added to prevent
  (`tulpa_compute_param_layout()` now rejects a partition that doesn't
  describe its adjacency instead of silently pinning nothing). Every
  ICAR/BYM2 occupancy NUTS fit through this path errored at first use.
  Now calls `set_spatial_adjacency()`. Requires `tulpa (>= 0.0.114)`.

* **Distance-sampling marginal sums truncate per site instead of to a
  shared ceiling (bugfix, #168).** Each site's binned-distance
  marginal-count sum now caps at `K_hi = min(K_max, K_lo + headroom)`
  instead of always summing to the shared `K_max`, mirroring
  `nmix_kernel.h`'s per-site truncation. A verified-widening safety
  layer ports nmix's approach: after fitting under a capped headroom,
  the score at the fitted coefficients is compared under both the
  capped and uncapped truncation, and disagreement widens the headroom
  4x and re-fits, escalating to fully uncapped if needed. Wired into
  `distance_laplace()`, `.tobs_fit_distance_spatial()`,
  `.tobs_fit_ms_distance()`, and `.tobs_fit_distance_nuts()`.

* `distance_laplace()`'s `quad_xptr` parameter is now documented
  (codoc fix; the parameter itself shipped in 0.0.179 for #165).

## 0.0.179 (2026-07-27)

* **macOS `R CMD check` builds again** (`tulpa (>= 0.0.101)`,
  `Remotes: gcol33/tulpa@v0.0.101`). `tulpa/portable_math.h` returns
  `std::pair` while including only `<cmath>` and `<limits>`; libstdc++ supplies
  `<utility>` transitively, Apple clang's libc++ does not, so
  `count_grouped_oracle.cpp` failed with `no template named 'pair' in namespace
  'std'` and the whole macOS job died before producing a check directory. The
  header now includes `<utility>` itself -- a header-only fix with no behaviour
  on any platform, so nothing in this package moves with it.

  Worth recording how long it hid: `R-CMD-check.yaml` runs ubuntu alone on push
  and adds windows and macOS only on the weekly cron, so a macOS-only compile
  error is invisible to every push and surfaces once a week at most. The two
  platforms that pass are the two that share libstdc++.

## 0.0.178 (2026-07-27)

* **The cover-hurdle occurrence reduction no longer shifts the reported
  log-marginal (bugfix).** `control$aggregate.occ` (default ON) collapses the
  occurrence arm to one Binomial row per group of exchangeable Bernoulli plots.
  That leaves the gradient, the Hessian, the mode, the SEs and every posterior
  weight pointwise unchanged, but the binomial log-density carries
  `lchoose(n, y)`, which is zero on a Bernoulli row and non-zero on the group
  replacing it -- so `fit$log_marginal` came back on a different additive scale
  depending on an internal performance switch (measured 74.7 and 72.4 nats on an
  18-cell x 10-plot fixture). The constant is now subtracted, putting the
  reported value back on the scale of the observed per-plot data. Estimates are
  unaffected in either direction; only the reported marginal moves, and only
  when the reduction runs. A model comparison between two fits made with the same
  setting was already valid and is unchanged.

* **`TULPAOBS_PROGRESS=0` silences the console progress bar process-wide.** The
  bar defaults on and is set per fit from `control$progress`, which `tobs()`
  writes to the scoped `tulpa.nl_progress` option on every call -- so a caller
  that cannot reach the individual fits could not turn it off at all. A test
  suite is exactly that caller: the weekly recovery job's log is mostly
  per-iteration bars from hundreds of fits, which buries the failure output the
  log exists to show. Setting the environment variable (also `false` / `no` /
  `off`) flips the default for the whole process. An explicit
  `control$progress` still wins in both directions, and the `progress.file`
  heartbeat is untouched -- it is the only liveness signal on a detached run,
  so silencing the console must not take it along.

* **Two recovery blocks that had never run are now green, and one of them was
  hiding a real bug.** Both sit behind `skip_if_fast()` and `skip_on_cran()`, so
  neither smoke nor `R CMD check` reaches them and the weekly full-recovery job
  has never completed a run -- the cover-hurdle aggregation block was byte
  identical to the commit that introduced it and had simply never executed. It
  found the log-marginal shift above. The second, `test-dyn-int-occu-recovery.R`,
  failed on the test runner rather than the model: the runner passes
  `package = "tulpaObs"` so that oracle files can reach internals through `:::`,
  which also puts the namespace on the symbol-lookup chain ahead of the search
  path and shadows `unmarked`'s S4 `logLik()` with the S3 generic. Both
  `unmarked` calls are now namespace-qualified.

* **The cover trend-field test asserts the amplitude ratio, not the amplitude.**
  `sigma` and `sigma_trend` are ICAR conditional-scale hyperparameters, and the
  map onto a field's marginal SD is an engine parameterisation: both moved by the
  same factor (-33% and -32% over 12 seeds, sd 0.031) while the occurrence slope,
  its coverage and `alpha_trend` stood still. Two independently simulated fields
  with different truths shifting together is a convention change, not a recovery
  regression. The test now asserts the convention-free estimand -- the ratio of
  the two amplitudes, 1.129 against a simulated 1.143 -- and keeps an absolute
  bracket only as a gross-regression guard, labelled as such.

* **Engine floor moved to `tulpa (>= 0.0.100)`** (`Remotes: gcol33/tulpa@v0.0.100`),
  which carries exact dispersion and closed-form outer curvature for
  zero-inflated and hurdle mixtures. `gcol33/tulpaMesh` stays at `@v0.1.3`.

## 0.0.177 (2026-07-19)

* **`n.iter` is the post-warmup sample count on every NUTS path (bugfix).** The
  observation-family and community samplers already ran `n.iter + n.warmup`
  iterations and kept `n.iter` post-warmup draws, but the single-season
  occupancy path (`occu()` / `dyn_occu()` / `int_occu()`) and the cover paths
  (`cover()` / `occu_cover()`) instead treated `n.iter` as the *total* including
  warmup, so they kept only `n.iter - n.warmup` draws -- and returned zero draws
  (all-`NaN` summaries) when `n.iter == n.warmup`. All paths now agree:
  `n.iter` = kept samples, total run per chain = `n.iter + n.warmup`. The `?tobs`
  control docs are corrected to match. Existing calls keep more draws (and run
  `n.warmup` more iterations) than before on the occu/cover paths.

## 0.0.176 (2026-07-19)

NUTS calibration coverage and areal-field exposure from the read-only
test-coverage audit (closes #139, #142; #137/#140 test gaps also closed).

* **The single-season occupancy NUTS path now exposes its areal field
  (#142).** `occu()` with an `icar()` / `bym2()` term under `method = "nuts"`
  sampled the field but left the field parameters unnamed (`param[k]`) and
  `fit$spatial_field` `NULL`, so it could only be smoke-tested. The engine
  already exports the block offsets on its parameter layout; `occu_fit.cpp` now
  emits them as `spatial_layout` and names the columns
  (`spatial_field[i]` / `spatial_theta[i]` / `log_tau_spatial` /
  `log_sigma_spatial` / `logit_rho_spatial`), and `fit$spatial_field` carries
  the centred per-cell surface for `icar` / `car_proper` (the intrinsic level is
  confounded with the intercept, so it is centred as the nested-Laplace summary
  is). The `bym2` field is the Riebler rho-mix of two blocks scaled by the graph
  scale factor, so its columns are named but its reconstruction is left to the
  named draws. Recovery-tested: `icar` field correlation ~0.81 across seeds
  (`test-occu-areal-nuts-recovery.R`).

* **Observation-family, abundance, community, and single-season NUTS paths get
  20-seed CI-coverage tests (#139).** The family paths (`distance` / `removal` /
  `fp_occu` / `dyn_abun` / `abun`) and community paths (`ms_count` /
  `ms_int_occu` / `ms_dyn_occu` / `ms_abun`) had only single-seed recovery;
  `dyn_occu` / `int_occu` NUTS had no recovery test at all. All now carry a
  >=20-seed 95% CI-coverage study on the calibrated interval the NUTS path
  exists to get right (`test-family-nuts-coverage.R`). The single-season
  `cpp_occu_fit` path reports no per-parameter SD, so those blocks score the
  credible interval off the draws.

## 0.0.174 (2026-07-18)

Consistency, single-source, and WAIC-boundary fixes from the read-only audit
(closes #141, #135, #133).

* **The Polya-Gamma Gibbs conjugate machinery is a single shared engine
  (#135).** The six PG-Gibbs fitters (`occu()`, `t_occu()`, `ms_occu()`,
  `ms_int_occu()`, `ms_dyn_occu()`, `ms_count()` / `jsdm()`) re-inlined the same
  conjugate coefficient draw, the community mean + Inverse-Gamma variance loop,
  the hyperprior constants, and the post-sampling summary / fit assembly. These
  now live once in `R/pg_gibbs_shared.R` (`.tobs_pg_draw_beta`,
  `.tobs_pg_community_update`, `.tobs_pg_hyper`, `.tobs_pg_summarize`,
  `.tobs_pg_finalize_fit`); each fitter keeps only its latent-state step and
  design layout. Behaviour is byte-identical (verified: every family's fit is
  bit-for-bit unchanged at a fixed seed).

* **The `#116` observation families give the friendly method error (#141).**
  `royle_nichols()`, `occu_ttd()`, `occu_multi()`, `double_observer()`,
  `gdistremoval()`, and `distsamp_open()` already sat in `.tobs_family_methods`,
  so an unsupported `method =` was already rejected with a "Supported: laplace"
  pointer; a regression test now pins that contract.

* **`occu_cover()` WAIC / LOO pointwise density matches the fit kernel at the
  cover boundary (#133).** The WAIC / LOO pointwise log-density
  (`src/occu_cover_ploglik.cpp`) is now the fit kernel
  `src/occu_coupling_shared.h::pos_log_density` itself, and the R positive-arm
  density (`.occu_cover_pos_logdens`, via the new `.tobs_log_safe`) drops the eta
  clamp and guards `log(0)`. A detected visit with cover exactly 0 or 1 now
  yields a finite pointwise density (the number the model was fit with) instead
  of `-Inf` / `NaN`, and the density no longer clamps the predictor. The
  three-level `src/occu_ms_cover_ploglik.cpp` gains the missing
  `isfinite(cover)` guard at a detected plot, so a missing-at-random cover drops
  out instead of poisoning the pointwise sum. Non-spatial fits are numerically
  unchanged.

## 0.0.173 (2026-07-18)

Three inline-variant scope notes graduated to working features.

* **`distsamp_open()` gains alternative population dynamics
  (`dynamics = "notrend" / "trend" / "autoreg" / "ricker" / "gompertz"`).**
  Following `unmarked::distsampOpen()`, the transition `N_t | N_{t-1}` can now be
  the density-independent `"constant"` (default, `Binomial(N_{t-1}, omega) +
  Poisson(gamma)`), the stationary `"notrend"` (recruitment tied to
  `(1 - omega) * lambda`, no free `gamma` arm), the exponential `"trend"`
  (`Poisson(N_{t-1} * gamma)`, no survival), the `"autoreg"` density-dependent
  recruitment (`Binomial(N_{t-1}, omega) + Poisson(gamma * N_{t-1})`), or the
  carrying-capacity forms `"ricker"` / `"gompertz"` (an estimated `K` reported as
  `K`, a growth rate reported as `r`). `"constant"` / `"notrend"` fit with the
  exact analytic gradient; the density-dependent forms use the same exact
  forward-HMM marginal through a value-only kernel
  (`compute_dyn_abun_site_dyn()` in `src/dyn_abun_kernel.h`, validated against a
  brute-force reference to `< 4e-15`) with a numeric gradient. The abundance and
  distance-scale arms recover on every dynamics; survival / recruitment / K / r
  sit on the usual short-series ridge, and the density-regulated forms grow
  toward the carrying capacity, so bound `K_max` (the forward is cubic in the
  truncation) when abundance is high. `simulate_distsamp_open()` gains
  `dynamics`; NUTS / areal / season-varying detection stay follow-ups.
* **`double_observer()` gains the dependent (removal-style) protocol
  (`type = "dependent"`).** A primary observer records what it detects and a
  secondary observer records only the primary's misses (two cell counts:
  primary-detected and secondary-only). A single fixed primary observer gives
  two cell means for three parameters (`lambda`, `p1`, `p2`) -- a ridge -- so the
  two detections are not separately identified; observer role-swapping breaks it,
  with observer 1 primary at some sites and observer 2 primary at others (a
  per-site `primary` indicator in `{1, 2}`), giving four cell means that identify
  all three. The marginal is again a product of independent thinned Poissons, fit
  by `optim` BFGS with an observed-information vcov. A single-observer `primary`
  warns (no role-swapping). `simulate_double_observer(type = "dependent")` draws
  the two cells with an alternating primary; the independent three-cell protocol
  is unchanged.
* **`count()` gains a continuous NNGP Gaussian-process field via `gp()`.**
  `count(...) + gp(lon, lat, prior_range = c(r0, alpha))` on the abundance formula
  routes to tulpa's single-block `nngp` nested-Laplace kernel: the GP marginal
  variance and range are integrated on the kernel's own outer grid and the field
  is Schur-folded out, so the fit reports grid-integrated fixed effects plus the
  GP hyperparameter posterior (`fit$gp_hyper`, `sigma2` / `phi_gp`). The per-cell
  field is integrated out on this path (`fit$spatial_field` is `NULL`); use
  `spde()` for a reconstructed continuous field map. Poisson / binomial only (a
  negbin size / gaussian residual variance is not jointly identified against a
  per-node field, gated with a pointer), as for the areal `icar()` / `car_proper()`
  path; `multiscale_gp()` (two-scale) is not hosted by the nested-Laplace engine
  and errors with a pointer to `gp()` / `spde()`.

## 0.0.172 (2026-07-18)

* **`distsamp_open()` gains negative-binomial and zero-inflated initial
  abundance (`mixture = "negbin" / "zip" / "zinb"`).** The open-population
  distance marginal already runs through the `dyn_abun()` forward-HMM kernel,
  which carries `use_nb` / `eta_logr`, so the negative-binomial initial abundance
  is threaded straight through with an analytic `log_r` score (FD-validated,
  `cor = 1`, `max|an - fd| = 6e-8`). Zero inflation is a pure-R additive layer
  over the composed per-site marginal (the `abun()` / `dyn_abun()` ZIP recipe): a
  structural-zero site is never occupied across any primary period, so all its
  band counts are zero and the marginal is
  `omega * 1{all zero} + (1 - omega) * L_open`, with `omega` an intercept-only
  structural-zero probability reported as `zi_logit` / `zi_omega` (distinct from
  `omega`, the survival arm). The ZI structural-zero share recovers cleanly (a
  single-seed probe recovered 0.30 against a truth of 0.30); the abundance and
  distance-scale arms recover on every mixture. `simulate_distsamp_open()` gains
  `mixture` / `size` / `zi`. Poisson fits are byte-identical. As with the base
  family the survival / recruitment ridge (and the NB size) is weakly identified
  at short series; the forward is cubic in the truncation `K_max`, so a high
  abundance is slow. NUTS / areal / season-varying detection stay a follow-up.
* **Fix: `distsamp_open()` WAIC / pointwise log-likelihood.** The per-site
  `.dso_site_loglik()` read `ev$log_lik` (the summed dataset scalar) instead of
  `ev$log_lik_site` (the per-site vector), adding the whole-dataset HMM total to
  each site's band term. WAIC / LOO / `tobs_ploglik()` are now per-site correct.

## 0.0.171 (2026-07-18)

* **New family `distsamp_open()` -- open-population distance sampling (#116, the
  final item; closes #116).** The unmarked `distsampOpen` model: a Dail-Madsen
  open N-mixture (as [dyn_abun()]) observed by distance sampling at each primary
  period -- the open-population counterpart of the single-season
  `gdistremoval()`. `N_1 ~ Poisson(lambda)`, `N_t = Binomial(N_{t-1}, omega) +
  Poisson(gamma)` (apparent survival `omega`, recruitment `gamma`); at each period
  the detected birds are distance-sampled into bins. The band allocation is
  conditional on the period total detected, so it factors out of the abundance
  HMM (the `gdistremoval()` trick): the marginal reuses the validated `dyn_abun`
  forward kernel with the detection probability set to the overall distance
  detection `pdist` (`eta_p = logit(pdist)`), plus the per-period band
  multinomials -- no new HMM kernel. Four site-level arms: log abundance
  (`formula`), log distance scale (`detection`), logit survival (`omega = ~ ...`),
  log recruitment (`gamma = ~ ...`). Maximised by `optim` BFGS with an ANALYTIC
  gradient (one kernel call per evaluation -- the `dyn_abun` kernel returns
  `grad_eta_*`, and the sigma block chains through `pdist` via a distance-integral
  finite difference; FD-validated to 4e-7) and an observed-information vcov. The
  `K_max` truncation defaults to the detection-corrected abundance scale
  (`max(period total) / pdist + headroom`), since the forward recursion is cubic
  in `K_max`. `simulate_distsamp_open()`, full S3 (`fitted` / `predict` /
  `residuals` / WAIC), `y` = a `[n_sites x n_bins x n_seasons]` array of
  distance-band counts per primary period. Recovery-tested: 20/20 seeds converge,
  every parameter's 95% coverage >= 0.90 (pooled 0.933), max coefficient bias
  ~0.068 at `n_seasons = 4` (`test-distsamp-open.R`). Scope for the first ship:
  half-normal key, line / point transect, Poisson initial abundance, constant
  Dail-Madsen dynamics, site-level arms; NB / ZIP, other dynamics, and
  season-varying detection are follow-ups.

  This closes #116 -- every enumerated model family (Royle-Nichols, ZIP/ZINB,
  continuous SPDE fields on the count families, generalized multinomial /
  double-observer, time-to-detection occupancy, community distance sampling,
  multi-species co-occurrence occupancy, joint distance+removal, and now
  open-population distance sampling) is shipped and recovery-tested.

## 0.0.170 (2026-07-18)

* **New family `gdistremoval()` -- joint distance + removal sampling (#116, the
  last open item).** The unmarked `gdistremoval` model (Amundson et al. 2014): a
  single-season point count where the detected birds are cross-classified by a
  distance band (distance sampling) AND the removal period of first detection
  (removal sampling). Site abundance `N ~ Poisson(lambda)`; the total detected is
  a binomial thinning of `N`, and Poisson is closed under binomial thinning, so
  the marginal is closed-form -- `ysum ~ Poisson(lambda * pdist * prem)` with the
  band and period counts as two conditional multinomials (the `double_observer()`
  Poisson-multinomial pattern, here with a half-normal distance multinomial and a
  depleting-removal multinomial `pi_k = r(1-r)^{k-1}`). Three site-level arms: log
  abundance (`formula`), log distance scale (`detection`), logit per-period
  removal capture (`removal = ~ ...`). Maximised by `optim` BFGS with an
  observed-information vcov. The half-normal band integrals are closed-form for
  both line and point transects (matching numeric integration to machine
  precision). `simulate_gdistremoval()`, full S3 (`fitted` / `predict` /
  `residuals` / WAIC), `y` = distance-band counts + `y_rem` = removal-period
  counts (per-site totals must match). Recovery-tested: 20/20 seeds converge,
  coefficients unbiased (max bias 0.011), pooled 95% coverage ~0.975
  (`test-gdistremoval.R`). Scope for the first ship: half-normal key, line / point
  transect, Poisson abundance, constant removal capture, availability fixed at 1
  (a single primary period does not identify availability separately). Hazard-rate
  key, NB / ZIP abundance, and a multi-period availability arm are follow-ups
  (each closed-form under the same thinning). This closes #116.

## 0.0.169 (2026-07-18)

* **Community-covariance CI coverage validation for `ms_occu_cover()` NUTS
  (#115, B7 DoD).** The community-model NUTS Definition of Done asks for community
  *covariance* recovery at nominal 95% coverage over >= 20 seeds, not only the
  mean and the de-attenuation direction. Added a 20-seed test
  (`test-ms-occu-cover-nuts.R`) that reconstructs, per fit, the posterior for the
  six community SDs (`sqrt(diag(Sigma_arm))` for the occ / p / pos arms, 2 coords
  each, from the sampled log-Cholesky coordinates via `.ms_ocs_chol_unpack`) and
  checks each central 95% interval against the simulated truth. Measured: pooled
  coverage 0.908 over the six SD components x 20 seeds (every component >= 0.85),
  0 divergences on every seed. This is the joint-cover analogue of the `ms_occu`
  covariance-coverage sweep (0.0.167); the covariance sampling machinery
  (`.ms_ocs_*` log-Cholesky pack / unpack + prior) is byte-shared across all four
  community NUTS families, and the two direct measurements (a two-state
  single-season marginal and the heavier three-arm joint-cover marginal) bracket
  the family likelihoods. Closes the last #115 item -- both B2 (community
  occupancy NUTS for `ms_occu` / `ms_dyn_occu` / `ms_int_occu`) and B7
  (`ms_occu_cover` NUTS + dispersion-RE) are shipped and validated to the DoD.

## 0.0.168 (2026-07-18)

* **Continuous SPDE field on the count families (`count()` + `spde()`; #116).**
  A continuous-mesh Matern field via `spde(lon, lat, ...)` on the abundance
  formula now fits under `method = "nested_laplace"`, alongside the areal
  icar / car_proper / bym2 fields. The latent lives on the mesh nodes
  (`fit$spatial_field`, length `n_mesh`) and the barycentric projector
  `fit$spatial$tulpa_spec$A` (`n_sites x n_mesh`) maps it onto the sites, so
  `fitted()` / `predict()` add the projected per-site contribution
  `A %*% mesh_field`. Three touch points: the count dispatch gate accepts
  `spde`; `.tobs_nested_attach_field_summary` reconstructs the mesh field from
  the grid-averaged block modes (the same icar-style summary, one field node per
  mesh node); `.count_spatial_field_offset` A-projects the mesh field to sites.
  The areal (icar/car_proper/bym2) path is byte-identical -- those fits carry no
  projector (`tulpa_spec$A` is `NULL`), so the A-projection branch is skipped.
  Poisson only (as the areal count field: a per-node field is not jointly
  identified with a negbin size / gaussian residual variance under the
  fixed-dispersion nested loop). Recovery-tested: the mesh field projects to the
  sites with `cor > 0.9` against the simulated truth and the abundance slope
  recovers (`test-count-spatial.R`, gated `skip_if_no_tulpamesh()`).

## 0.0.167 (2026-07-18)

* **Community-covariance CI coverage validation for `ms_occu()` NUTS (#115).**
  The community-model NUTS Definition of Done asks for community-*covariance*
  recovery, not only the mean. Added a 20-seed test that reconstructs, per fit,
  the posterior for the community SDs (`sqrt(diag(Sigma_arm))` from the sampled
  log-Cholesky coordinates) and checks their 95% CIs against the simulated truth.
  Measured coverage: `sd_psi[1]` 0.85, `sd_psi[2]` 0.95, `sd_p[1]` 0.90 (pooled
  ~0.90), at/above the 0.85 rubric floor -- the calibration the Laplace-EM
  variance lower bound cannot provide (`test-ms-occu-nuts.R`).

## 0.0.166 (2026-07-18)

* **bym2 areal field on the count families (#116).** `count(...) + bym2(graph)`
  under `method = "nested_laplace"` now fits and reconstructs its field, joining
  the icar / car_proper areal support. The generic nested-Laplace field summary
  (`.tobs_nested_attach_field_summary`) gained a bym2 branch: it reads the block's
  `sigma` / `rho` grid axes and the (phi | theta) mode slices and reconstructs the
  rho-mixed unit field `z = sqrt(rho / scale) * phi + sqrt(1 - rho) * theta`
  (Riebler 2016), grid-averaged. Measured field correlation ~0.93 with the
  simulated (structured) truth over three seeds; the abundance slope recovers
  (`test-count-spatial.R`). The same summary path drives `occu()` etc., so a bym2
  areal field on those families is now reconstructed (`fit$spatial_field`) rather
  than dropped. icar / car_proper reconstruction is unchanged. spde() / gp() on
  the count families remain follow-ups.

## 0.0.165 (2026-07-18)

* **Per-species dispersion random effect for `ms_occu_cover()` NUTS (#115 B7).**
  `ms_occu_cover(..., method = "nuts", control = list(dispersion.re = TRUE))`
  promotes the shared community cover log-dispersion to a per-species random
  effect `log_disp_s = mu_ld + sigma_ld * z_ld_s` -- a fourth 1-D community arm
  (the `ms_abun` `log_r_s` analogue), sampled jointly with the coefficient arms.
  It reuses the same non-centered machinery (the 1x1 log-Cholesky factor is
  `sigma_ld`) and the extracted per-cell sweep, so the shared-dispersion path is
  byte-identical. The C++ `FullGradFn` is byte-exact against the FD-validated R
  oracle (lp diff 0); 0 divergences. `fit$ms_dispersion` carries
  `sigma_log_disp` (the community dispersion SD) and `log_disp_species` (the
  per-species posterior-mean log-dispersions). `test-ms-occu-cover-nuts.R`. (A
  simulated per-species dispersion spread for a full `sigma_ld` recovery sweep is
  a simulator follow-up; on shared-dispersion data `sigma_ld` is small, as
  expected.)

## 0.0.164 (2026-07-18)

* **NUTS for the community joint occupancy + cover family (`ms_occu_cover()`,
  #115 part B7).** `ms_occu_cover(..., method = "nuts")` (non-spatial) now samples
  the exact joint posterior -- the community means, the per-species
  occupancy / detection / cover deviations, the three independent per-arm
  community covariances (occ / p / pos), AND the shared community log-dispersion
  -- over the per-(species, cell) two-state occu_cover marginal via an in-tree C++
  `FullGradFn` (`src/ms_occu_cover_nuts.cpp`), warm-started at the Laplace-EM mode.
  The joint-cover analogue of the `ms_occu` / `ms_int_occu` community NUTS targets
  (per-arm non-centered blocks `b = C z` + a shared log-dispersion scalar); it
  wraps the existing cover-density kernels (`occu_coupling_shared.h`,
  `pos_log_density` / `pos_grad_eta` / `pos_grad_logdisp`) in the community loop
  rather than reimplementing them. Non-centered, so 0 divergences, and it
  **removes the documented Laplace community-variance attenuation** (the
  `fit$ms_community$var_attenuation` caveat) -- the sampled per-arm community SD
  de-attenuates toward truth where the Laplace estimate collapses at small
  per-species n. The R oracle is FD-validated and the C++ port is byte-exact
  against it (`test-ms-occu-cover-nuts.R`: R-oracle vs FD, C++ vs oracle,
  `b_from_z` round-trip, community-mean recovery + de-attenuation). With this the
  community-model NUTS family (#115 B2) is joined by the joint-cover route (B7).
  Non-spatial lognormal / beta / gaussian; negbin / dispersion-RE remain
  follow-ups. The shared spatial-factor NUTS route (a field on the occupancy arm)
  is unchanged.

## 0.0.163 (2026-07-18)

* **NUTS for the community integrated occupancy family (`ms_int_occu()`, #115).**
  `ms_int_occu(..., method = "nuts")` now samples the exact joint posterior --
  the community means, the per-species occupancy / per-source detection
  deviations, and the D + 1 independent per-arm community covariances -- over the
  multi-source two-state per-(species, site) marginal via an in-tree C++
  `FullGradFn` (`src/ms_int_occu_nuts.cpp`), warm-started at the community
  Laplace-EM mode. The multi-source generalisation of the `ms_occu` / `ms_dyn_occu`
  NUTS targets (per-arm non-centered blocks `b = C z`, no shared globals), reusing
  the shared `.ms_ocs_` epilogue helpers (#128). Non-centered, so 0 divergences,
  and it removes the Laplace-EM small-cluster attenuation of the community
  variance (the sampled community SD de-attenuates vs the Laplace lower bound).
  The R oracle `.tobs_ms_int_occu_nuts_logpost` is FD-validated and the C++ port
  is byte-exact against it (`test-ms-int-occu-nuts.R`: R-oracle vs FD, C++ vs
  oracle, `b_from_z` round-trip, community-mean recovery + de-attenuation).
  Non-spatial only.

## 0.0.162 (2026-07-18)

* **Deduplicated the community-NUTS R epilogue into shared `.ms_ocs_*` helpers
  (#128).** The tail half of every in-tree community-NUTS fitter (`ms_occu`,
  `ms_dyn_occu`, `ms_abun` non-spatial + spatial, `ms_count`) carried
  near-byte-identical copies of the warm-start metric, the multi-chain harness,
  and the posterior-summary code. These now live once in the `.ms_ocs_`
  namespace: `.ms_ocs_b_idx()` (per-species z-block indices), `.ms_ocs_pd()`
  (symmetrise + jitter to PD), `.ms_ocs_fd_metric(grad_fn, theta)` (FD-Hessian
  inverse mass, taking the per-family gradient closure), `.ms_ocs_sig_mean()`
  (community-covariance posterior mean), and `.ms_ocs_run_chains(run_chain,
  n_chains)` (the multi-chain runner). The five families call them instead of
  private twins; net ~180 fewer duplicated lines. Behaviour-preserving: the
  refactor only relocates identical code and factors the gradient closure out of
  the metric, so the NUTS draws are byte-identical on a fixed seed (verified per
  family) and every recovery test is unchanged. The per-family `_b_from_z`
  reconstruction stays local (the scalar dispersion / Gaussian-variance arms make
  it genuinely family-specific), and `ms_count`'s log-phi prior likewise.

## 0.0.161 (2026-07-18)

* **bym2 / proper-CAR areal recovery tests for the observation families (#131).**
  The `nested_laplace` areal path on `distance`, `removal`, `fp_occu`, and
  `dyn_abun` was recovery-tested for `icar` but only smoke-tested (finite vcov,
  non-null field) for `bym2` / `car_proper`. `bym2` is a distinct code path -- it
  reconstructs the rho-mixed unit field `z = sqrt(rho) * phi + sqrt(1 - rho) *
  theta`, which the `icar` recovery never exercises. Each family now has a field
  correlation + FE slope recovery block for `bym2` (all four) and `car_proper`
  (`removal` upgraded from smoke, `fp_occu` / `dyn_abun` newly covered; `distance`
  proper-CAR is already recovered via the hazard-key variant). Measured field
  correlations: count families recover strongly (distance ~0.92, removal ~0.89,
  dyn_abun ~0.82), the weakly-identified occupancy `fp_occu` field lower (~0.37 --
  one binary site per node) but well above chance. Reuses the existing `icar`
  fixtures; only the term and the field reconstruction differ.

## 0.0.160 (2026-07-18)

* **Parameter-recovery tests for three grouped-RE paths (#130).** The detection
  (`p`) arm RE on `abun()`, the cover (`pos`) arm RE on `occu_cover()`, and the
  uncorrelated detection random slope `(0 + x | g)` on `occu_cover()` were
  shipped with label/shape smoke tests only. Each now has a multi-seed recovery
  block asserting the RE SD recovers (band, measured from the sibling paths) and
  the per-group BLUPs correlate with the simulated offsets. To support the cover
  arm, `simulate_occu_cover()` gained `re_pos_groups` / `sigma_re_pos`, which
  inject a per-visit random intercept on the positive-cover linear predictor
  (mirroring the existing `re_det_groups` on the detection arm) and return
  `truth$b_pos_re` / `truth$sigma_re_pos`.

## 0.0.159 (2026-07-18)

* **Maintenance: dead-helper removal and test-suite hygiene (#129, #132).**
  Removed four internal helpers with no call sites (`.cover_arm_to_slot`,
  `.occu_cover_joint_coef_names`, `.occu_cover_ppc_cover`, `.se_from_hessian`);
  the `.tobs_ppc_occu_cover()` PPC computes its cover discrepancy inline, so the
  extracted twin was orphaned. Kept the third-derivative helper `.l3_poisson_log`
  for family symmetry with its exercised `.l3_binomial_logit` /
  `.l3_gaussian_identity` siblings and added its unit test (#129).
  `test-sla-int-occu.R` gained the `skip_if_fast()` / `skip_on_cran()` guards its
  four integrated-occupancy SLA fits were missing, so the fast pass no longer runs
  them (#132). The SPDE recovery suite's `skip_if_not_installed("tulpaMesh")`
  gates became `skip_if_no_tulpamesh()`, which fails loudly under
  `TULPAOBS_REQUIRE_SPDE=1` instead of silently skipping the whole SPDE surface
  when `tulpaMesh` is absent on a CI runner (#132). Also refreshed the stale
  `ms_dyn_occu` method-coverage assertions to reflect its new NUTS route (0.0.158).

## 0.0.158 (2026-07-18)

* **NUTS for the community dynamic occupancy family (`ms_dyn_occu()`, #115).**
  `ms_dyn_occu(..., method = "nuts")` now samples the exact joint posterior --
  the community means, the per-species first-season / detection deviations, the
  two independent per-arm community covariances, and the shared colonisation /
  extinction transition globals -- over the per-(species, site) HMM-forward
  marginal via an in-tree C++ `FullGradFn` (`src/ms_dyn_occu_nuts.cpp`,
  `R/ms_dyn_occu_nuts.R`), warm-started at the community Laplace-EM mode with a
  diagonal Laplace metric. Non-centered per-species blocks
  `b_{s,arm} = C_arm z_{s,arm}` break the covariance/deviation funnel; the C++
  gradient is byte-exact against the R oracle. Removing the Gaussian
  approximation de-attenuates the community variance components the Laplace-EM
  under-reports on the binary arms. Non-spatial only (a structured term still
  routes to `nested_laplace`). `test-ms-dyn-occu-nuts.R` (R-oracle vs
  finite-difference gradient, C++ byte-exact, community-mean recovery + 0
  divergences, S3 methods, variance de-attenuation, spatial gate).

## 0.0.157 (2026-07-17)

* **Regularize the community zero-inflated N-mixture structural-zero variance
  (#116).** `ms_abun(mixture = "zip" / "zinb")` now puts a weak Penalized-
  Complexity prior on the per-species structural-zero random-effect SD
  `sigma_omega` (new `control$omega.sigma.prior = c(U, alpha)`, default
  `c(1, 0.05)`; set `NULL` to restore pure ML). `sigma_omega` is the softest
  AGHQ direction and at few species can collapse to the boundary, flattening the
  marginal Hessian and attenuating the recovered SD; the prior adds curvature
  there without biasing an identified fit (measured: a collapse-prone seed's
  `sigma_omega` recovers 0.07 -> 0.18 toward a 0.3 truth, identified seeds
  unmoved). Consumes tulpa 0.0.85's `tulpa_re_aghq(sigma_prior = ...)`.
* **Actionable error when `K_max` is below the largest observed count.** The
  community N-mixture marginal sums the latent N only to `K_max`, so a
  user-supplied `K_max < max(y)` made the per-(species,site) marginal
  structurally `-Inf` and the joint optimum singular -- previously surfaced as an
  opaque "singular marginal Hessian". It now fails early with a message pointing
  at `K_max`. The default (`max(y) + 100`) never trips it.
* **Faster community N-mixture AGHQ fits**, recovery preserved. Three changes to
  the joint AGHQ integration of the community count families: the per-node
  log-likelihood takes a fast path that skips the count-moment (digamma /
  trigamma) pass (byte-identical `log_lik`); the default zero-inflated `n_quad`
  drops to 3 (from 5), making a default `zinb` fit tractable; and the scalar
  nuisance blocks (the NB dispersion `log_r`, the structural-zero `omega`)
  integrate at fewer quadrature nodes than the correlated abundance / detection
  blocks via tulpa 0.0.85's per-block `n_quad` (`control$n.quad.scalar`, default
  2). Measured on the test fixtures: NB ~24% faster, ZINB ~33% faster, with
  `lambda` / `r` / `omega` unchanged.
* Require `tulpa (>= 0.0.85)` and pin `Remotes: gcol33/tulpa@v0.0.85`.

## 0.0.156 (2026-07-17)

* Require `tulpa (>= 0.0.84)` and pin `Remotes: gcol33/tulpa@v0.0.84` so an
  install pulls the engine build in which a joint `occu_cover` fit with a
  checkpoint and `diagnose_k = TRUE` no longer aborts with a fingerprint
  mismatch from the Pareto-k diagnostic re-solve (tulpa#161).

## 0.0.155 (2026-07-17)

* Zero-inflated counts for the community N-mixture `ms_abun(mixture = "zip" /
  "zinb")` (spAbundance-style structural absence; #116). Structural zeros are a
  per-species random effect `logit_omega_s ~ N(mu_omega, sigma_omega)`: a share
  `omega_s` of a species' sites is never occupied. It mirrors the per-species NB
  dispersion `log_r_s` design, integrated jointly with the abundance / detection
  (and, under `zinb`, dispersion) coefficients by the community AGHQ oracle. The
  per-site marginal wraps the Royle marginal in a structural-zero mixture in
  `NMixCommunityOracle` -- an all-zero site is `log(omega_s + (1 - omega_s)
  exp(L_i))`, a detection site rules out `N = 0` -- with the score / observed-
  information / Fisher composed in closed form from the plain-marginal kernel
  outputs (no marginal-kernel change). `mu_omega` joins the coefficient surface
  as `logit_omega`; `ms_zi` carries the community-mean `omega` and the
  per-species `omega_s`; `ranef()` gains the structural-zero deviation. A shared
  field / `latent()` factor on a `zip` / `zinb` fit errors rather than silently
  dropping the structural-zero share. `simulate_ms_abun(mixture = "zip" /
  "zinb")` generates the truth; 20-seed community-mean recovery holds at the
  coverage floor. This completes the ZIP / ZINB axis across `abun()` /
  `ms_abun()` / `dyn_abun()`.

* Fixed `simulate()` on a zero-inflated `abun()` fit dropping the structural
  zeros: the posterior-predictive count draws ignored `zi_omega`, so a `zip` /
  `zinb` fit's replicates had no excess zeros. The structural-zero Bernoulli is
  now drawn per site (byte-identical stream when no zero-inflation).

## 0.0.154 (2026-07-17)

* Partial season overlap for `dyn_int_occu()` (spOccupancy `tIntPGOcc`; #122).
  A detection source that does not observe a `(site, season)` marks it `NA`, so a
  staggered survey where sources cover different seasons is expressed by
  NA-padding each source to the common `[n_sites x max_visits x T]` grid: the
  forward drops an absent source from that season's emission (`nvalid = 0`) and
  marginalises a `(site, season)` unobserved by every source (`e0 = e1 = 1`).
  `simulate_dyn_int_occu(source_seasons = list(1:4, 3:6))` generates staggered
  coverage; 20-seed recovery holds. Two boundary anchors pin the reduction to its
  neighbours: one source (the others all `NA`) reproduces `dyn_occu()`, and one
  season of data (`T = 2`, season 2 all `NA`) reproduces `int_occu()`, both to
  ~1e-3 -- the #122 definition-of-done anchors.

* Fixed the non-spatial `dyn_int_occu()` variance-covariance sign: the observed
  information was built from the Jacobian of the *log-likelihood* gradient (`-I`),
  so `solve()` returned the negative vcov and every standard error was clamped to
  zero, silently breaking `confint()` and interval coverage. It is now the
  Jacobian of the negative-log-likelihood gradient (`+I`); SEs are proper and
  95% Wald intervals cover at ~0.90 over 20 seeds. The spatial fitter was
  unaffected.

## 0.0.153 (2026-07-17)

* Visit-varying detection for `royle_nichols()` (unmarked `occuRN` with
  observation covariates; #116). Detection was site-level only; passing `visits`
  now lets `logit r_ij` carry visit-level covariates. The latent-`N` marginal
  generalises from the `(k_i, n_i)` sufficient-statistic form to the full
  per-visit product inside the same Poisson sum,
  `prod_j [1 - (1 - r_ij)^N]^{y_ij} [(1 - r_ij)^N]^{1 - y_ij}`, which reduces
  exactly to the site-level form when `r_ij` is constant (verified to machine
  precision). `fitted()` / `predict(type = "detection")` / `simulate()` return the
  per-visit `[n_sites x max_visits]` detection grid; the site-level path is
  byte-identical. `simulate_royle_nichols(beta_r_visit =)` draws a visit-level
  detection covariate and returns `visits` for `tobs(..., visits =)`. Recovery
  over 20 seeds: max rel-bias ~0.04, min 95% CI coverage ~0.95
  (`test-royle-nichols.R`).

## 0.0.152 (2026-07-17)

* Varying-coefficient areal field on multi-season integrated occupancy:
  `dyn_int_occu()` + a `spatial(~ 1 + w || cell, graph)` bar (spOccupancy
  `svcTIntPGOcc`; #122). Alongside the shared intercept field (`stIntPGOcc`), a
  covariate-weighted (spatially-varying-coefficient) field loads on the
  first-season occupancy logit. The areal-BFGS nested-Laplace driver
  (`R/areal_bfgs.R`) already accepts a list of field blocks and scatters the
  per-site psi1 score `w1 - psi1` to each; the weighted block is the ICAR field
  with a `w`-weighted loading (`eta += w * z[map]`, the field `z` reported
  unweighted). `.areal_field_car()` gained an optional `weight` (byte-identical
  when absent), and `.tobs_areal_field_blocks()` expands a bar into an intercept
  block plus one weighted block per covariate. `fit$spatial_field` is the
  intercept surface, `fit$trend_field(s)` the varying-coefficient surface(s);
  both recover by correlation (intercept ~0.85, trend ~0.72;
  `test-dyn-int-occu-areal-recovery.R`). `simulate_dyn_int_occu(trend =)`. icar
  only in v1.

## 0.0.151 (2026-07-17)

* Multi-season occupancy with an AR1 year random effect: `t_occu()`
  (spOccupancy `tPGOcc`; #124). Per-`(site, season)` Bernoulli occupancy with a
  shared AR1 year effect on the occupancy logit and NO colonization / extinction
  transition (that is `dyn_occu()`); use it for an occupancy trend over years
  with a temporal random effect. Given the year effects the seasons factorise,
  so the fit is an exact Polya-Gamma Gibbs sampler (`method = "pg_gibbs"`, the
  engine spOccupancy uses): draw the latent state, then the joint
  `(beta_occ, eta)` as one Gaussian Markov random field update with the AR1
  precision as the year-effect prior, then `beta_p`, then the AR1
  hyperparameters (`sigma^2` conjugate Inverse-Gamma, `rho` on a grid). The
  year-effect surface (`fit$temporal_field`), the occupancy / detection
  coefficients, and the innovation SD `sigma` recover cleanly (20-seed
  `test-t-occu.R`, year-effect correlation ~0.97, `Rhat` ~1.00); the AR1
  correlation `rho` is a weakly-identified parameter of a short time series that
  recovers only as the number of seasons grows (an identifiability property of
  the AR1, not the sampler -- given even the true year effects `rho_hat` climbs
  from ~0.08 at 8 seasons to ~0.59 at 200, `dev_notes/_probe_tpg_rho.R`), so it
  is reported but not asserted tightly. `y` is a 3D array
  `[n_sites x n_seasons x max_visits]` or a list of per-season matrices;
  `simulate_t_occu()`. v1 = site-level occupancy / detection covariates,
  `pg_gibbs` only.

## 0.0.150 (2026-07-17)

* Polya-Gamma Gibbs for community multi-source integrated occupancy:
  `ms_int_occu(method = "pg_gibbs")` (#115, #126). The community integrated
  extension of msPGOcc: one latent occupancy state per (species, site) observed by
  D detection sources, with per-species occupancy and per-source detection
  coefficients under Gaussian community hyperpriors. Conditional on the
  Polya-Gamma auxiliaries each sweep samples, per species, the single latent z
  (occupied if any source detects, else Bernoulli on the pooled
  occupied-undetected mass across sources), the PG-conjugate `beta_psi_s`, and the
  D per-source `beta_pd_s` (each at that species' occupied, covered sites), then
  the conjugate community mean and Inverse-Gamma community variance per coordinate
  per arm. This gives a calibrated community-variance posterior (`sd_psi` ~0.53 vs
  a 0.5 truth) where the community Laplace-EM attenuates it; the community means
  recover and split-Rhat ~ 1.0. `ms_community` matches the Laplace fit's layout
  (`Sigma_`/`sd_`/`coef_`/`blup_<arm>`), so `coef()` / `ranef()` / `fitted()`
  work unchanged. `test-ms-int-occu-pg-gibbs.R`. Completes the community-occupancy
  Polya-Gamma trio (`ms_occu` / `ms_dyn_occu` / `ms_int_occu`).

## 0.0.149 (2026-07-17)

* Shared areal field on multi-season integrated occupancy: `dyn_int_occu()` +
  `~ 1 + icar(graph = adj)` under `method = "nested_laplace"` (the `spOccupancy`
  `stIntPGOcc` model, #122). The field sits on the first-season occupancy (`psi1`)
  arm; because `psi1` only sets the initial mixing weight of each site's HMM, the
  exact per-site field gradient is the Fisher-identity score `w1 - psi1` (the
  smoothed season-1 occupancy). Fit through the shared areal-BFGS nested-Laplace
  driver (the same recipe as `fp_occu` / `dyn_abun`), one field unit per site; the
  colonization / extinction / per-source detection arms carry fixed effects only.
  The interior field recovers (`cor` ~0.8), transition rates recover.
  `simulate_dyn_int_occu(field =)`; `test-dyn-int-occu-areal-recovery.R`. `icar()`
  only in v1 (bym2 / car_proper, a varying-coefficient bar, and NUTS are
  follow-ups).

* `dyn_int_occu()` now fits with analytic forward-backward gradients (the
  Fisher-identity score of the multi-source colext HMM: `w1 - psi1` on occupancy,
  smoothed pairwise transition joints on colonization / extinction, occupancy-
  weighted binomial score per detection source) instead of a finite-difference
  BFGS with a numeric Hessian. Faster (a fit that took seconds of numeric
  differencing now completes in well under a second) and the observed-information
  vcov comes from the FD-Jacobian of the analytic gradient; the non-spatial fit is
  otherwise unchanged (recovery regression tests pass).

## 0.0.148 (2026-07-17)

* Season-varying detection on `dyn_occu()` (#124): a detection covariate supplied
  as a `[n_sites x T]` matrix column of `data` (one column per primary season)
  drives per-season detection, `logit p_it = X_it beta_p`, fit with
  `detection = ~ det_cov`. The E-step emission reads the season's detection
  probability, the M-step encodes one detection row per `(site, season)` at that
  season's covariate, and the exact season-varying HMM-forward marginal refine
  calibrates the coefficients (recovered with ~95% coverage over 20 seeds). A
  plain per-site detection covariate keeps the byte-identical constant-detection
  path. Detection unrolling shares the period-agnostic
  `.tobs_period_arm_design()` with the season-varying colonization / extinction
  path (the `#80` recipe); gated under `method = "nuts"` with a pointer.
  `simulate_dyn_occu(beta_det_season =)`; `test-dyn-occ.R`. `fitted()$z` is now
  period-varying-aware (also corrects the smoothed state for season-varying
  colonization / extinction, which the previous per-site indexing got wrong).

## 0.0.147 (2026-07-17)

* Zero-inflated open N-mixture: `dyn_abun(mixture = "zip" / "zinb")` (#116). A
  structural-zero site is never occupied across any season, so all its counts are
  zero; the observed per-site marginal is the two-component mixture
  `omega * 1{all y = 0} + (1 - omega) * L_DailMadsen`, the same additive
  structural-zero layer the static `abun()` ZIP uses (`nmix_zip.R`) but over the
  exact forward-HMM open-population marginal. Pure-R over the C++ per-site
  marginal; the Poisson / negbin paths are byte-identical. Analytic-gradient BFGS
  (the Dail-Madsen per-site eta gradients weighted by the structural-zero
  posterior, plus the mixture score for the ZI logit), observed-information vcov
  from the FD-Jacobian of the analytic gradient. The structural-zero share
  recovers with nominal coverage over 20 seeds; the ZI logit is named `zi_logit`
  (dyn_abun's `omega_*` is the SURVIVAL arm). `simulate_dyn_abun(zi =)` injects
  structural zeros. Non-spatial laplace, intercept-only omega; a field / RE / NUTS
  stay Poisson / negbin with a pointer. `test-dyn-abun-zip.R`. This completes the
  ZIP / ZINB mixture axis across `abun()` and `dyn_abun()`.

* Spatially-varying coefficient on community dynamic occupancy:
  `ms_dyn_occu()` + `spatial(~ 1 + w || cell, graph)` (the `svcTMsPGOcc` model,
  #123). An intercept field plus a shared covariate-weighted field on the
  first-season occupancy arm, both recovering (`cor` ~0.90 / ~0.89). No new code
  was needed: the psi1 oracle already returns per-site/per-species score+curv, and
  the K-field weighted-ICAR block-coordinate solve is the same machinery the
  community count SVC (`svcMsAbund`) uses -- the weighted bar flows through the
  existing dynamic-spatial dispatch unchanged. `simulate_ms_dyn_occu(trend =)`
  adds the varying-coefficient surface; `test-ms-dyn-occu.R`. Completes the
  spatial multi-season community occupancy line (`stMsPGOcc` + `svcTMsPGOcc`);
  `bym2()` / `car_proper()` fields and a NUTS sampler remain follow-ups.

## 0.0.146 (2026-07-17)

* Spatial Polya-Gamma Gibbs for single-season occupancy: `occu() + icar()` under
  `method = "pg_gibbs"` (#126; the spOccupancy `spPGOcc` engine). An intrinsic
  areal (ICAR) field on the occupancy logit, `logit(psi_i) = X_i beta + f_i`,
  `f ~ ICAR(tau)`. Conditional on the Polya-Gamma auxiliaries the joint
  `(beta, f)` update is a single Gaussian Markov random field draw (dense
  `(p + n) x (p + n)` precision `t(W) Omega W + blkdiag(B0inv, tau Q)`,
  `W = [X | I]`); `tau` has a conjugate Gamma full conditional; the field is
  centred to sum-to-zero each sweep with its level moved into the intercept (so
  `eta` is preserved -- omitting that leaks a systematic shift). The field
  recovers (`cor ~ 0.74`), the intercept and detection recover, split-Rhat
  ~ 1.00. As with any spatial occupancy model, a site-level occupancy covariate
  spatially-confounds with the saturated one-node-per-site field, so the
  recovery target is the field + intercept + detection. icar only in v1
  (bym2 / car_proper are follow-ups). `test-occu-pg-gibbs-spatial.R`.

  This completes the Polya-Gamma engine (#126): PGOcc + spPGOcc (`occu`),
  msPGOcc (`ms_occu`), tMsPGOcc (`ms_dyn_occu`), and the community Bernoulli /
  binomial GLMM (`jsdm()` / `ms_count("binomial")`).

## 0.0.145 (2026-07-17)

* Polya-Gamma Gibbs for the community Bernoulli / binomial GLMM: `jsdm()` and
  `ms_count(response = "binomial")` with `method = "pg_gibbs"` (#126). The
  community logistic GLMM has no latent state and no detection sub-model -- `y`
  is the observed k-of-n (n = 1 for jsdm's presence/absence) -- so this is the
  simplest of the PG engines: per-species coefficients with Gaussian community
  hyperpriors, each an exactly conjugate Gaussian update conditional on the
  Polya-Gamma auxiliaries, plus the conjugate community mean and near-Jeffreys
  Inverse-Gamma community variance. One fitter serves both front doors (jsdm =
  bernoulli, `ms_count("binomial")` = k-of-n). Samples the exact community
  posterior, so the community variance recovers where the Laplace-EM attenuates
  (jsdm `sd_mu` 0.674 / 0.439 vs truth 0.7 / 0.5, vs Laplace 0.635 / 0.401);
  community means + per-species coefficients recover, split-Rhat ~ 1.00. Only the
  logistic responses are routed here -- Poisson / negbin / gaussian reject
  `pg_gibbs` with a pointer. Non-spatial in v1 (the sfMsPGBinom / lfJSDM
  spatial-factor PG variants are follow-ups). `R/ms_count_pg_gibbs.R`;
  `test-ms-count-pg-gibbs.R`.

## 0.0.144 (2026-07-17)

* Community dynamic occupancy Polya-Gamma Gibbs: `ms_dyn_occu()` with
  `method = "pg_gibbs"` (#115, #126; the spOccupancy `tMsPGOcc` engine). Combines
  the community PG machinery (msPGOcc) with a two-state HMM forward-filter
  backward-sample (FFBS) latent-occupancy step: each sweep FFBS-samples the
  per-species occupancy path z, then does the PG-augmented conjugate coefficient
  updates -- per-species season-1 occupancy `psi1` and detection `p` (with
  Gaussian community hyperpriors), and the SHARED community-level colonization
  `gamma` and extinction `eps` from the aggregated 0-origin / 1-origin
  transitions across all species -- plus the conjugate community mean and
  Inverse-Gamma community variance for the psi1 / p arms. Samples the exact
  community posterior, so the community variance recovers where the Laplace-EM
  attenuates. The shared gamma / eps recover tightly (informed by every species),
  the community season-1 occupancy mean and SD recover, and split-Rhat ~ 1.00. v1:
  constant transitions, site-level detection, no structured terms.
  `R/ms_dyn_occu_pg_gibbs.R`; `test-ms-dyn-occu-pg-gibbs.R`.

## 0.0.143 (2026-07-17)

* Community occupancy Polya-Gamma Gibbs: `ms_occu()` with `method = "pg_gibbs"`
  (#115, #126; the spOccupancy `msPGOcc` engine). The hierarchical extension of
  the single-species PGOcc engine: per-species occupancy / detection coefficients
  with Gaussian community hyperpriors, each PG-augmented conjugate-Gaussian
  update wrapped in conjugate community-mean and near-Jeffreys Inverse-Gamma
  community-variance draws (a diagonal per-arm community covariance, as
  spOccupancy uses). Conditional on the Polya-Gamma auxiliaries every coefficient
  update is exactly conjugate, so this samples the exact community posterior --
  and, unlike the community Laplace-EM (whose variance components carry a
  documented small-cluster attenuation), it recovers the community VARIANCE:
  across seeds the PG `sd_psi` / `sd_p` land on the true community SDs (e.g.
  0.667 / 0.417 vs truth 0.6 / 0.4) where the Laplace-EM attenuates below them
  (0.595 / 0.342). The SD is reported as the posterior median (robust to the
  variance-component skew at moderate species counts). Community means and
  per-species coefficients recover; split-Rhat ~ 1.00. v1: single-season
  community occupancy, site-level detection, no structured terms (the shared-field
  sfMsPGOcc is a follow-up). `R/ms_occu_pg_gibbs.R`; `test-ms-occu-pg-gibbs.R`.

## 0.0.142 (2026-07-17)

* Polya-Gamma Gibbs engine for single-season occupancy: `method = "pg_gibbs"`
  (#126, the spOccupancy `PGOcc` engine). A REAL Gibbs chain over the exact
  occupancy posterior via Polya-Gamma data augmentation (Polson, Scott & Windle
  2013) -- conditional on the PG auxiliaries, both the occupancy and detection
  logistic coefficient updates are exactly conjugate Gaussian. This is distinct
  from `method = "laplace_gibbs"`, which is a stochastic-EM variance correction
  (mode-finds + Rubin pooling, no PG augmentation, stationary distribution not
  the posterior). Each sweep samples the latent occupancy `z`, then draws
  `omega ~ PG` and the conjugate Gaussian coefficients on each arm, using tulpa's
  tested Polson-Scott-Windle sampler (`tulpa:::cpp_rpg`). Reports split-Rhat and
  bulk-ESS (a real MCMC). Validated: on a well-identified model the PG-Gibbs
  posterior matches the Laplace observed-Fisher fit (means within 1 SE, SDs
  within 20%), Rhat ~ 1.00 / ESS ~ 700-900, and 95% credible intervals cover at
  the nominal rate over 15 seeds. v1: single-season `occu()`, site-level
  detection, no structured terms (the PG-spatial extensions --
  `pg_binomial_{icar,bym2,rsr,gp}` in tulpa -- are the documented follow-up).
  `R/occu_pg_gibbs.R`; `test-occu-pg-gibbs.R`.

## 0.0.141 (2026-07-17)

* Multi-season integrated occupancy: `dyn_int_occu()` (#122, spOccupancy
  `tIntPGOcc`). The product of the two shipped families -- a dynamic
  (multi-season HMM) occupancy state (`psi1`, colonization `gamma`, extinction
  `eps`) whose per-season emission pools SEVERAL detection sources (integrated
  occupancy). Pooling sources across seasons is how colonization / extinction
  estimates come out of individually-sparse opportunistic data. The latent
  occupancy sequence integrates out by the two-state HMM forward recursion (the
  `dyn_occu` exact-marginal forward generalised with a per-season emission pooled
  over sources), maximised by `optim` with an observed-information vcov -- pure R,
  no new C++. `y` is a list of `S` `[sites x visits x seasons]` detection arrays;
  the state `formula` models `logit psi1`, `colonization = ~` / `extinction = ~`
  the transitions (required, as in `dyn_occu()`), `detection` the shared
  per-source design (each source carries its own coefficients). Full `fitted` /
  `predict` / `residuals` / `simulate` / WAIC surface; `simulate_dyn_int_occu()`.
  Recovers `psi1` / `gamma` / `eps` and per-source detection over 12 seeds. v1:
  full site/season overlap, constant transitions, site-level detection,
  non-spatial laplace; partial overlap, season-varying rates, an areal `psi1`
  field (`stIntPGOcc`), a weighted bar (`svcTIntPGOcc`), and NUTS are documented
  follow-ups. `test-dyn-int-occu.R`.

## 0.0.140 (2026-07-17)

* Double-observer abundance: `double_observer()` (#116, \pkg{unmarked}
  `multinomPois` with the independent double-observer pi-function). Site
  abundance `N ~ Poisson(lambda)` surveyed by two independent observers with
  detection `p1` / `p2`; each individual is recorded by observer 1 only, 2 only,
  or both. By Poisson-multinomial thinning the three observable cell counts are
  independent Poissons (`n_c ~ Poisson(lambda * pi_c)`), so the marginal is
  closed form with NO latent-abundance summation -- the simplest of the
  self-contained pure-R families. `y` is an `N x 3` cell-count matrix
  (observer-1-only, observer-2-only, both); the state `formula` models
  `log lambda`, `detection` the shared per-observer design (observers carry
  separate coefficients). `optim` BFGS mode-find with observed-information vcov
  (`R/double_observer.R`); full `fitted` / `predict` / `residuals` / `simulate` /
  WAIC surface; `simulate_double_observer()`. Recovers abundance + both detection
  probabilities with calibrated slope coverage. Site-level detection, non-spatial
  laplace in v1. `test-double-observer.R`.

## 0.0.139 (2026-07-17)

* Multi-species co-occurrence occupancy: `occu_multi()` (#116, \pkg{unmarked}
  `occuMulti`; Rota et al. 2016). The joint occupancy state of `S` species at a
  site, `z in {0,1}^S`, follows a log-linear model whose natural parameters are
  first-order (per species) and second-order (per species pair -- the
  interaction: positive = co-occur more than independent, negative = avoid); each
  species is then detected conditional on presence. The latent state integrates
  out by enumerating the `2^S` states, so -- like `royle_nichols()` / `occu_ttd()`
  -- the exact marginal is maximised by `optim` with an observed-information vcov
  (`R/occu_multi.R`, no new C++). `y` is a list of `S` `N x J` detection matrices
  (or a 3D array) with `species =`; the state `formula` is the shared occupancy
  covariate design (each natural parameter carries its own coefficients),
  `detection` the shared site-level per-species design. Full `fitted`
  (`$psi` = marginal per-species occupancy) / `predict` / `residuals` / `simulate`
  / WAIC surface; `simulate_occu_multi()`. The log-linear natural-parameter
  covariate slopes are only weakly identified individually (they trade off), but
  the marginal occupancy and the interaction recover cleanly -- the recovery test
  targets those (interaction sign + magnitude over 15 seeds, positive and
  negative). Shared design, first + second order, site-level detection, laplace
  in v1 (per-parameter formulas, higher-order terms, visit-level detection, areal
  fields are follow-ups). `test-occu-multi.R`.

## 0.0.138 (2026-07-17)

* Time-to-detection occupancy: `occu_ttd()` (#116, \pkg{unmarked} `occuTTD`). A
  survey records the TIME to first detection rather than a 0/1 outcome; at an
  occupied site the time-to-detection is exponential with rate `lambda` (constant
  hazard) over a survey of length `surveyLength`, censored if it reaches the end.
  The state `formula` models `logit psi`, `detection` the `log lambda` rate. The
  latent occupancy state integrates out in closed form (two states), so -- like
  `royle_nichols()` -- the exact marginal is maximised by `optim` with an
  observed-information vcov (`R/occu_ttd.R`, `.ttd_site_loglik` over per-site
  sufficient statistics). `y` is an `N x J` matrix of detection times (a value
  `>= surveyLength` is a non-detection, `NA` a survey not conducted).
  `simulate_occu_ttd()`, full `fitted` / `predict` / `residuals` / `simulate` /
  WAIC surface. Site-level rate, exponential TTD, non-spatial laplace in v1
  (Weibull shape, visit-varying rate, areal fields, NUTS are follow-ups).
  `test-occu-ttd.R` (constructor + gates, S3, 15-seed recovery with calibrated
  rate-slope coverage).

## 0.0.137 (2026-07-17)

* Zero-inflated N-mixture abundance: `abun(mixture = "zip")` / `"zinb"` (#116).
  A structural-zero share `omega` of sites have `N = 0` regardless of `lambda`
  (`L_i = omega * 1{all y_i = 0} + (1 - omega) * L_royle_i`). Implemented as a
  pure-R additive layer over the shared per-site Royle marginal
  (`nmix_site_marginal()`), so no marginal-kernel change and the plain
  Poisson / negbin paths are byte-identical; `.tobs_fit_nmix_zip`
  (`R/nmix_zip.R`) mode-finds by L-BFGS-B over
  `[beta_lambda | beta_p | logit_omega | (log_r)]` with box bounds only on the
  pathological `logit_omega` / `log_r` corners (betas unbounded, no interior
  bias). `logit_omega` is a model coordinate in `coef`/`vcov`/`summary` (like the
  negbin `log_r`); `fit$zi_omega` reports the structural-zero probability. ZIP
  recovers betas + omega cleanly with calibrated slope coverage; ZINB is weakly
  identified (structural zeros vs NB overdispersion both absorb zeros) and
  recovers in a higher-count regime. Non-spatial laplace only in v1 -- a
  nuts / nested_laplace / structured-term request errors at the dispatcher with a
  pointer rather than silently dropping the zero-inflation. `simulate_abun(mixture
  = "zip"/"zinb", omega =)`. `test-abun-zip.R`.

## 0.0.136 (2026-07-16)

* Identity-Gaussian positive arm on the community and multiscale cover-hurdle
  families (#127): `ms_occu_cover(response = "gaussian")` and
  `occu_multiscale_cover(response = "gaussian")`, the delta-normal magnitude
  (`mu = eta`, shared residual SD `sigma_pos`, no log Jacobian) for a
  pre-transformed / unbounded positive response. The density
  (`GaussianPositive`) already shipped for the single-species arms (#112); this
  threads the `response` enum, the `OccuMultiscaleCoverGaussianCoupling` typedef,
  and a shared positive-code (`0` lognormal / `3` beta / `4` gaussian) through
  the community Laplace-EM, the multiscale nested-Laplace coupling, the
  non-spatial Laplace / NUTS marginal, the WAIC ploglik kernel, and the simulate
  draws (`.occu_cover_pos_code`, replacing the binary `is_beta` branches). The
  gaussian residual SD recovers to lognormal parity (0.392 vs 0.392 at a matched
  config); the shared-field spatial-factor `ms_occu_cover()` variant stays
  lognormal-only and now rejects the new arm explicitly.

## 0.0.135 (2026-07-16)

* Detection-arm areal fields on the observation families (#114). A
  `detection = ~ icar()` (or `car_proper()` / `bym2()`) term now loads a
  spatially-varying detection surface on the detection arm of `removal()`,
  `distance()`, `fp_occu()`, and `dyn_abun()`, alongside the existing
  abundance/occupancy-arm areal fields. Each routes through the shared areal-BFGS
  driver with a `det_arm` flag that sends the field offset to the detection scale
  and reads the detection-arm gradient: the capture logit for `removal()` (per-
  pass design, so the per-pass score is `rowsum`-summed to a per-site field
  gradient), the true-positive logit `eta_p11` for `fp_occu()`, the per-site
  detection logit applied across every season for `dyn_abun()`, and the detection
  scale `eta_sigma` for `distance()` (working under both the half-normal and the
  hazard key, with the log-shape threaded as a global in the same evaluation).
  `fit$spatial_field_arm` labels which arm carries the field. Poisson and negbin;
  the detection-arm field is `nested_laplace` only (NUTS carries a field on the
  abundance/occupancy arm alone). Recovery-tested in `test-obs-det-field.R`.
* Distance-sampling grouped random effects under the hazard key (#114).
  `distance()` with a grouped random effect on the abundance arm now fits under
  the hazard-rate key as well as the half-normal. The scalar log-shape is not a
  per-site design column, so instead of a second global parameter it is profiled
  over the AGHQ log-marginal -- an outer optimisation over the shape, each
  candidate a full AGHQ fit at a fixed shape the `DistanceGroupedOracle` carries
  -- and the profile log-shape row/column is inserted into the AGHQ variance-
  covariance matrix (profile standard error, zero off-diagonal). Poisson and
  negbin.
* Temporal-only fixed-hyper fields under NUTS for the observation families
  (#114). A `temporal()` term (AR1 / RW1 / RW2 / iid) on the
  abundance/occupancy arm of `removal()`, `distance()`, `fp_occu()`, and
  `dyn_abun()` now samples under `method = "nuts"` via the same non-centered
  fixed-hyper field block the areal NUTS path uses, with the field map keyed by
  period and a temporal whitened loading. Recovery-tested (AR1 field + slope, 0
  divergences) in `test-obs-nuts-temporal.R`; a simultaneous areal-and-temporal
  field under NUTS remains gated to `nested_laplace` with a pointer.

## 0.0.134 (2026-07-16)

* Shared areal field on community dynamic occupancy (#123), the `spOccupancy`
  `stMsPGOcc` model. A shared `icar()` field on the first-season occupancy
  formula of `ms_dyn_occu()` -- `~ 1 + icar(graph = adj)` under
  `method = "nested_laplace"` -- fits a spatial multi-season community occupancy
  model. Because the first-season occupancy only sets the initial mixing weight
  of each species' HMM, the per-(species, site) marginal is linear in it, the
  same two-component mixture the single-season community field uses, so the field
  routes through the shared block-coordinate driver with the HMM-forward
  conditional likelihoods in place of the single-season emission. The forward-
  backward kernels are vectorised over sites and the community EM is driven by an
  analytic Fisher-identity gradient (validated against finite differences), so
  the fit is fast. The shared field recovers cleanly (`cor` ~0.94 on an interior
  field) and the community transition-dynamics coverage clears the 0.85 floor.
  `simulate_ms_dyn_occu()` gained a `field =` argument. `icar()` only; a NUTS
  sampler and `bym2()` / `car_proper()` fields remain follow-ups (gated with a
  pointer).

## 0.0.133 (2026-07-16)

* Season-varying colonization / extinction on `dyn_occu()` (#124). A covariate
  supplied as an `[n_sites x (T-1)]` matrix column of `data` now drives an
  interval-indexed colonization or extinction rate --
  `colonization = ~ season_cov` / `extinction = ~ season_cov` -- the recipe
  `dyn_abun()` already uses for season-varying survival / recruitment (#80),
  ported to the colonization-extinction HMM forward. The forward-backward E-step
  uses a per-interval transition matrix (constant rates broadcast, so an existing
  site-level or constant fit is byte-identical to before), the transition M-step
  becomes a weighted logistic on the smoothed per-interval transitions, and the
  exact-marginal refine gains a season-varying HMM-forward marginal that
  calibrates the standard errors. `simulate_dyn_occu()` gained `beta_gamma` /
  `beta_epsilon` to draw season-varying truth. Recovery and ~95% Wald coverage
  hold across 20 seeds; the season-varying rate is gated under `method = "nuts"`
  (the compiled forward reads one rate per site) with a pointer to `laplace`.

## 0.0.132 (2026-07-16)

* Binomial `k`-of-`n` count GLMM without replicates (#125), the spOccupancy
  `svcPGBinom` family. `count(response = "binomial")` with a per-site `trials =`
  argument (default 1, i.e. Bernoulli) fits detection-free binomial data through
  the same `count()` front door as the Poisson / negbin / Gaussian responses --
  logit link, no dispersion (the variance is pinned by `n`). A plain areal field
  (`icar()` / `car_proper()`) is supported (`svcPGBinom`): unlike the negbin /
  Gaussian responses, the binomial is identified against a per-node field, at
  `trials = 1` too. `ms_count(response = "binomial")` extends the community
  Laplace-EM the same way (community `svcPGBinom`; NUTS and a shared field /
  latent factor are gated to laplace non-spatial for now). `simulate_count()` /
  `simulate_ms_count()` gained a binomial draw with `trials =`; fitted values are
  expected successes `n * p`, `predict(newdata=)` the per-trial probability, and
  WAIC / residuals score the binomial pmf. The non-spatial fit matches
  `glm(..., family = binomial)` to ~1e-3 and runs ~100x faster than
  `spOccupancy::svcPGBinom()`. `jsdm()` stays the Bernoulli (`trials = 1`) alias.

## 0.0.131 (2026-07-16)

* Identity-Gaussian positive arm for the cover hurdle (#112). `cover(response =
  "gaussian")` and `occu_cover(response = "gaussian")` fit the delta-normal
  hurdle -- a Bernoulli presence process times a Gaussian magnitude on the
  positive scale -- the lognormal arm without the log transform or its Jacobian
  (residual `y - eta`, mean `mu = eta`, draw `y = eta + sigma z`). This is for a
  signed / unbounded positive magnitude, NOT cover fractions on `(0, 1)` (use
  `beta` / `beta_oi` for those). Recovery and 95% interval coverage hold at the
  lognormal bar (20 seeds, pooled coverage 0.94-0.95). WAIC/LOO carry through; the
  Freeman-Tukey PPC (non-negativity assumed) and NUTS on the joint path are gated
  with a pointer. `simulate_cover()` / `simulate_cover_joint()` /
  `simulate_occu_cover()` gained the Gaussian draw.

* Intrinsic ICAR / BYM2 areal fields now sample under NUTS for the observation
  families and `ms_abun` (#113), via the #71 sum-to-zero reparameterisation
  (whitened loading drops the precision null-space, so `z` is auto-centred and
  `n_raw = n_units - 1` for ICAR / `2n - 1` for BYM2). Previously
  car_proper-only. The shared field-block C++ (`nuts_field_block.h`) generalised
  to a non-square loading byte-identically to the square path; a single R helper
  `.tobs_nuts_field_loading()` is the one source of the whitened loading. Per
  family the field posterior mean is sum-to-zero centred and correlates with
  truth at 0 divergences (recovery in `test-count-spatial-nuts.R`).

* Observation-family spatial/temporal breadth (#114). Four structured-term
  corners on the observation families:
  - Temporal-only fields under `nested_laplace` on `removal()`, `distance()`,
    `fp_occu()`, `dyn_abun()`: a `temporal()` term on its own (no companion areal
    field) runs the shared areal-BFGS driver with a single AR1/RW1/RW2/iid block
    on the family's structured arm. Recovers a known AR1 truth (cor 0.90-0.99),
    slope on truth; per-family multi-seed recovery tests.
  - `dyn_abun()` NUTS + `temporal()`: a fixed-hyper non-centered AR1/RW1/RW2/iid
    field on the initial-abundance arm rides the SAME NUTS field block as the
    areal field (`nuts_field_block.h`), with `field_map` = period index and a
    temporal whitened loading fixed at the nested-Laplace temporal-only estimate.
    0 divergences, temporal field cor ~1.0. Areal + temporal under NUTS stays
    gated to `nested_laplace`.
  - `distance()` detection-arm areal field: a field in the `detection=` formula
    (`detection = ~ icar(graph)`) loads on the per-site detection scale (log
    sigma) instead of the abundance arm -- a spatially-varying detection scale,
    via the per-site sigma gradient the kernel already exposes. Half-normal key;
    field recovers (cor ~0.97); recovery + gate tests.
  - `distance()` hazard-key areal field on the abundance arm was already wired and
    recovery-tested (`test-distance.R`), confirming that corner.
  A temporal term under `nuts` on the other observation families, and a temporal
  term on the `count()` family, remain gated with a pointer.

## 0.0.130 (2026-07-16)

* `temporal()` now errors on the Laplace engine instead of being dropped.
  `.tobs_laplace()` consumes the `spatial` and `re` terms but has no temporal
  channel, and `.tobs_fit_model()` never passed one, so
  `tobs(~ x + temporal(year), family = occu(), method = "laplace")` -- the
  DEFAULT method for `occu()` / `dyn_occu()` / `int_occu()` -- silently fitted a
  model without the field the formula asked for, with no warning. The term is
  carried by `nested_laplace` (grid-integrated, composing with an areal field and
  `re()` blocks) and by `nuts` (`populate_temporal`); the gate points at both.
  This is the failure the `svc()` gate directly above it already guarded against.

* Repaired 25 malformed `sprintf()` gate messages across 15 files. A multi-line
  message written as comma-separated literals inside `sprintf()` without a
  `paste0()` wrapper leaves only the FIRST literal as the format string; the rest
  silently become `...` arguments. 17 of the 25 errored at call time rather than
  reporting anything: where a string literal landed on a `%d` (every
  "spatial term has %d units but the model has %d sites" gate -- `abun`,
  `ms_abun`, `removal`, `distance`, `fp_occu`, `dyn_abun`, `ms_occu`, and the
  shared `areal_bfgs` driver -- plus the `distance()` / `ms_distance()`
  `cutpoints` length check and the community `icar` node-count check), and one
  invalid `%q` directive in the `.tobs_family_methods` dispatcher fallback. The
  remaining 8 truncated mid-sentence, dropping the value the user needed: the
  `abun()` / `ms_abun()` continuous-field gates stopped before naming `spde()`,
  the package's only continuous-field escape hatch for that family, and the
  unsupported-areal-type gate never interpolated the type it had rejected.
  Every gate now emits its full message with the offending value interpolated.

* `tests/testthat/test-gate-messages.R` covers both: the message text on the
  repaired gates, and a package-wide AST scan asserting no `sprintf()` call in
  `R/` has a format/argument arity mismatch or an invalid directive. Gates fire
  from error paths, so no fitting test exercised any of them -- which is how 25
  broken messages accumulated unnoticed.

## 0.0.129 (2026-07-16)

Requires tulpa >= 0.0.82 (ABI 34) and rebuilds against it.

* BREAKING: the continuous NNGP fields `gp()` and `svc()` now require
  `prior_range = c(r0, alpha)`, a PC prior on the spatial range encoding
  `P(range < r0) = alpha` -- the contract `spde()` has always had. tulpa 0.0.82
  ships the range anchors unset and refuses a NNGP block without them
  (gcol33/tulpa#144), and neither package invents a default: the range is in
  coordinate units, so no value suits every dataset, and the range is weakly
  identified by the likelihood, so a default would be a prior doing real work
  on the posterior rather than a convenience. `phi_prior_lower` /
  `phi_prior_upper` are gone (they parameterized the Uniform-behind-a-wall
  prior that #144 removed). On unit-square coordinates
  `prior_range = c(0.1, 0.05)` reads as "a 5% chance the range is under 0.1".

* `svc()` on single-season `occu()` under NUTS is now recovery-validated rather
  than shape-validated (gcol33/tulpaObs#118, #119). This is the measurement the
  engine fixes could not make upstream: `svc()` is a tulpaObs term, so no fit
  reachable from tulpa's own R side exercises the NNGP SVC path. Over seeds
  1/2/3/11 at N = 150, J = 6, p = 0.6, truth `phi = 0.25`, `sigma = 1.3`:
  divergences fell from 72-83% of post-warmup draws to **0 on every seed**, and
  `phi` from ~4 (the old Uniform prior's mean) to 0.14-0.23. Two upstream causes
  contributed -- the range prior (tulpa#144) and the SVC marginal-SD prior,
  which was improper on the coordinate it is sampled on and left `sigma`
  unbounded above; fixing the latter is also what pulled `phi` onto its truth,
  the two being the two ends of the GP ridge.

  Surface correlation with the known truth did **not** improve (0.73 mean over
  seeds, against 0.66 before). The calibration test used to predict it would,
  and no longer asserts on it: sampler health and surface accuracy are separate
  axes, and the latter is bounded by the information at these settings. The
  `skip()`ped calibration assertions are enabled.

## 0.0.128 (2026-07-15)

* Occupancy with an areal field under `method = "nested_laplace"` estimates the
  field precision. The state block was encoded at `M = 1000` pseudo-trials per
  site whenever the latent block was areal, but a state row carries one binary
  occupancy observation, so this overstated its information roughly 1000-fold.
  The field prior was swamped, between-cell binomial noise was read as a real
  field, and the state slope inflated with the field through the logistic
  conditional-to-marginal factor. Any nested latent block now encodes the state
  arm at `M = 1`, the site's actual information content. On `dyn_occu()` +
  `icar()` this moves the fitted field from sd 1.35 to 0.164 where the truth is
  0, and from 2.46 to 0.917 where the truth is 1.0; the season-1 slope moves from
  0.64 to 0.53 against a truth of 0.5, and interval coverage from 0.83 to 0.95.
  The encoding is monotone in `M` and no arm regresses, single-season included.

* `dyn_occu()` and `int_occu()` under `method = "nested_laplace"` report state
  standard errors. Both fell through to a path that finds no `H_beta` on a
  nested-Laplace fit and returns `NA`, so every such fit reported `NA` state SEs
  (and therefore no intervals) while the detection, colonization and extinction
  arms were finite. The Louis identity that the single-season arm already used is
  not family-specific -- the state arm's complete-data score is
  `x_i (z_i - psi_i)` in all three -- so it now covers them, reading the season-1
  smoothed weight column for dynamic and the shared per-site weight for
  integrated.

* `int_occu()` + an areal field additionally reached neither of the two fixes
  above: its M-step had no latent-block branch (so a nested areal block, which
  arrives with an empty spatial slot, took the `M = 1000` path), and the driver
  never handed it the latent prior, leaving its field-aware E-step unreachable.
  Its field now recovers at 0.976 against a truth of 1.0 and shrinks to 0.178
  against a truth of 0, with finite SEs on every seed.

* These combinations are now recovery-tested rather than smoke-tested
  (`test-occu-areal-recovery.R`, `test-int-occu-areal-recovery.R`,
  `test-dyn-occu-areal-recovery.R`, `test-dyn-occu-svc-recovery.R`). The previous
  tests asserted a class, a type string and a field length on fixtures containing
  no field, and passed throughout. Note that a null-field fixture cannot test this
  path at all: the outer grid runs over `tau` in [0.3, 30], so a true field sd of
  0 lies outside it and the marginal can only report the smallest field the grid
  expresses -- "correctly shrunk" and "pinned at the grid edge" give the same
  number. The new tests use a field the grid represents in its interior and assert
  the grid mode is off both boundaries.

* Community N-mixture latent factors: `ms_abun()` now takes `latent(n)` on the
  abundance formula (the spAbundance `lfMsNMix` analogue), giving residual
  species co-occurrence through per-site factors with per-species loadings, and
  composes with a shared field (`icar()` / `car_proper()` / `bym2()` / `spde()`)
  for the spatial-factor case. The latent abundance still marginalises in closed
  form per species-site, so the whole latent structure sits on `log lambda` and
  the family reduces to one working oracle over the Royle marginal, driving the
  shared block-coordinate engine added in 0.0.127. That oracle is the Louis
  (1982) block `nmix_site_marginal()` already exposes -- no new kernel. The
  residual species correlation recovers at ~0.99, and with a shared field the
  field and factors separate cleanly (field ~0.98, residual ~0.94). Poisson only:
  a negative-binomial size is a second per-site dispersion and is not identified
  against a per-site latent structure. A plain shared field with no factors keeps
  its existing C++ path.

* Community distance sampling: new `ms_distance()` family (the spAbundance `msDS`
  analogue) -- per-species binned distance sampling with Gaussian community
  hyperpriors on the abundance and detection-scale coefficients, with `latent(n)`
  factors (`lfMsDS`) and a shared field (`sfMsDS`). `y` is a
  `[n_sites x n_bins x n_species]` array or a named list of per-bin count
  matrices. The latent abundance marginalises per species-site exactly as for
  `distance()`, so the family adds no C++: the existing distance kernel already
  returns the per-site marginal pieces the community EM and the latent driver
  need, and the driver's oracle is the same Louis formula the community N-mixture
  uses. Under the hazard-rate key the log-shape is shared across species. The
  community means are unbiased over seeds with ~0.9 Wald coverage, and the
  residual species correlation recovers at ~0.99. Poisson only; no NUTS path yet.
  `simulate_ms_distance()` draws through the same quadrature the likelihood
  integrates against, so simulated data come from the `pi` the model is fit
  against.

* `abun()` and `ms_abun()` now reject a non-numeric `K_max`. `K_max` is the first
  argument, so `abun("negbin")` -- reaching for the mixing distribution -- bound
  the string to `K_max`, coerced it to `NA`, and surfaced much later as an
  unrelated comparison error inside a kernel. Write `abun(mixture = "negbin")`.

* The shared community factor update now rejects a non-finite Newton step,
  holding the previous iterate instead. A marginal that sums over a latent count
  can return a non-finite curvature once an iterate wanders far enough out, and
  the resulting `NaN` propagated into the factors and surfaced only later, as a
  non-finite standard deviation in the rescale. The guard is a no-op when every
  step is finite, so existing factor fits are unchanged.

* Caught up with the current `tulpa`: the nested-Laplace fitters no longer pass
  the removed `verbose` control, the joint paths pass `k_samples` (renamed from
  `diagnose_draws`), and the AGHQ calls pass `max_iter` (renamed from `maxit`).
  Against the current `tulpa` these had broken the areal `count()` path, the
  varying-coefficient abundance fit, and every `occu_cover()` nested-Laplace
  joint fit. The user-facing `control` names are unchanged.

* The shared community Laplace-EM gained two optional arguments, both defaulting
  to the previous behaviour byte-identically: `sp_info` (an analytic per-species
  observed information, for families whose kernel already exposes it -- the
  finite-difference fallback spends `2(P + G)` full marginal sweeps per species
  per Newton step rediscovering it) and `init_b` / `init_Sigma` (a warm start, so
  a block-coordinate caller re-entering the EM once per outer pass resumes
  instead of restarting cold).

## 0.0.127 (2026-07-15)

* One shared latent-structure engine for the community families
  (`R/community_latent.R`). The block coordinate ascent, the areal Newton, the
  latent-factor update and the field hyperparameter grids used to be duplicated
  across `ms_count_spatial.R`, `ms_count_factor.R` and `ms_occu_field.R`; they
  differed only in the per-(site, species) working score and curvature with
  respect to an additive offset on the structured arm. A family now supplies one
  callback -- `working(eta) -> list(score, curv)` (Poisson gives `(y - mu, mu)`,
  the occupancy two-state marginal its own pair, Bernoulli `(y - psi,
  psi (1 - psi))`) -- and gets the field, the factors and the SVC generalisation
  from the shared driver. The refactor is behaviour-preserving: every existing
  recovery test passes unchanged (the inner factor loops are exactly, not
  approximately, equivalent, since site `i`'s working weights depend only on
  `zeta[i, ]` and species `s`'s only on `lambda[s, ]`).

* Community latent-factor and spatial-factor occupancy (`ms_occu()` +
  `latent(n)`, with or without a shared field; the spOccupancy `lfMsPGOcc` /
  `sfMsPGOcc` analogues, gcol33/tulpaObs#119). Residual species co-occurrence on
  the occupancy arm via Q per-site factors with per-species loadings. Composed
  with a shared field the loadings are centred across species, so the field owns
  the shared spatial mean and the factors the between-species residual. The
  loadings / factors are identified only up to rotation, so the reported target
  is the residual species covariance `lambda lambda'`.

* `jsdm()` is now the community GLMM (gcol33/tulpaObs#121): per-species
  coefficients under a Gaussian community covariance -- the spOccupancy `lfJSDM`
  / `sfJSDM` model class -- rather than shared fixed effects with a scalar
  per-species intercept. A JSDM is exactly the `ms_count()` community model with
  a logit link, so it now shares that binder, community Laplace-EM, latent
  driver, NUTS target and S3 surface. `latent(n)` gives lfJSDM; a shared areal
  field alongside it gives sfJSDM. `method = "nuts"` samples the exact joint
  community posterior over the Bernoulli response through the family-aware
  in-tree C++ target (byte-exact against the R oracle). The former single-block
  correction routes (`laplace_sla` / `laplace_gibbs` / `laplace_mi`) do not apply
  to the community EM and are no longer offered.

* Continuous Matern (`spde()`) fields on the community latent driver. The mesh
  nodes carry the field and the barycentric projector `A` maps them onto sites,
  which is the same linear-map slot the areal `group_var` incidence uses, so the
  driver's `t(A) diag(w) A + tau Q` solve applies unchanged. The precision is
  `Q(kappa) = kappa^4 C0 + 2 kappa^2 G1 + G1 C0^-1 G1` scaled by the driver's
  `tau`, with `kappa` (the Matern range) chosen on a grid by the field marginal
  exactly as proper-CAR's `rho` is; it is proper, so the field carries no
  sum-to-zero constraint. `ms_occu()` and `ms_count()` gain continuous fields
  alongside `jsdm()`.

* Single-species spatially-varying-coefficient abundance (`count()` +
  `spatial(~ 1 + w || cell, graph)`; the spAbundance `svcAbund` analogue,
  gcol33/tulpaObs#120). No new engine: the multi-block prior already emitted one
  weighted latent block per resolved field for the count model, and the nested
  field summary already looped them into `spatial_field` + `trend_fields`. The
  fit records the per-site contribution `sum_k W[i,k] f_k[i]` so `fitted()` adds
  the full weighted field rather than the intercept field alone.

* `ms_count()` gains `predict()` and `residuals()`, which it never had (only
  `coef` / `ranef` / `fitted` / `simulate`).

* Bug fix: `tobs_waic()` / `tobs_loo()` / `tobs_cpo()` scored the community
  occupancy families **without a shared areal field**. `community_ploglik.R`
  built the occupancy predictor from the coefficients alone, so the shipped
  `ms_occu()` + `icar()` path was compared as a fixed-effect-only model. The
  field (and now the factor) offsets enter the pointwise log-likelihood.

* Bug fix: `control = list(max.outer = ...)` was read by the `ms_count()` /
  `ms_occu()` dispatchers but rejected by control validation, so the block
  coordinate outer cap was unreachable.

* Bug fix: the areal count fit recorded the intercept field as its per-site
  offset, which is the whole contribution only when there is no SVC field.

## 0.0.126 (2026-07-15)

* Community count NUTS now covers every response (`ms_count()` +
  `method = "nuts"`; gcol33/tulpaObs#117): the Poisson sampler generalises to
  negative-binomial and Gaussian. The negative binomial adds a per-species
  dispersion random effect `log_r_s ~ N(mu_log_r, sigma_log_r^2)` as a second
  community arm (the same structure the negbin Laplace-EM fits); the Gaussian
  adds S free per-species residual variances `log_phi_s` (matching the Laplace
  outer loop, which estimates each `phi_s` with no community prior), each with a
  weakly-informative prior. One family-aware in-tree C++ FullGradFn
  (`src/ms_count_nuts.cpp`) drives all three, byte-exact against the R oracle
  `.tobs_ms_count_nuts_logpost`; warm-started at the community Laplace-EM mode,
  0 divergences, NUTS == Laplace on the means, and the dispersion recovered
  (`fit$ms_dispersion` carries `r_s` / `variance`). This removes the Laplace-EM's
  mild negbin slope attenuation and returns calibrated non-Gaussian community
  intervals. `test-ms-count-nuts.R` adds the negbin / Gaussian oracle
  cross-check and recovery.
* Fixed a latent bug in the community count NUTS fitter that had gone unnoticed
  because its recovery test is `skip_if_fast`: the tulpa NUTS engine takes
  `n_iter` as the total (warmup + sampling) count, so passing `n.iter` alone
  returned zero post-warmup draws; the fitter now passes `n.iter + n.warmup` and
  surfaces the real per-draw diagnostics (`divergent` / `treedepth` / `epsilon`)
  instead of NA placeholders, and calls `.tobs_nuts_rhat_ess()` with its correct
  single argument.
* Community count NUTS now accepts missing (`NA`) site x species entries
  (`ms_count()` + `method = "nuts"`), matching the Laplace-EM path. An `NA` in the
  response matrix marks that (species, site) as unobserved and drops it from the
  per-(species, site) data sum, so a species keeps only its observed sites (ragged
  coverage). The C++ FullGradFn and the R oracle mask the missing entries
  identically, so the joint log-posterior + gradient stay byte-exact with `NA`
  present; the reported `N` counts only observed entries. `test-ms-count-nuts.R`
  adds the byte-exact NA cross-check (all three responses) and a NUTS-with-NA run
  that agrees with the Laplace-EM fit on the same ragged data.

## 0.0.125 (2026-07-14)

* Community count NUTS (`msAbund` NUTS; gcol33/tulpaObs#117): `ms_count()` +
  `method = "nuts"` samples the exact joint posterior of the non-spatial
  community Poisson GLMM -- the community means, the per-species coefficient
  deviations, and the community covariance -- via a new in-tree C++ FullGradFn
  (`src/ms_count_nuts.cpp`) over tulpa's NUTS engine, warm-started at the
  community Laplace-EM mode. The reduced counterpart of the `ms_abun` NUTS: no
  detection arm, no latent-N marginalisation, so the per-(species, site)
  contribution is a plain Poisson log-likelihood. NON-CENTERED (`b_s = C z_s`,
  `z_s ~ N(0, I)`) so the community covariance enters only the data term. The
  joint log-posterior + gradient are byte-exact vs the R oracle
  (`.tobs_ms_count_nuts_logpost`, checked to 1e-7); NUTS recovers the community
  means and agrees with the Laplace-EM mode with 0 divergences. Poisson;
  negbin / gaussian community NUTS are follow-ups.

## 0.0.124 (2026-07-14)

* Occupancy community SVC (`svcMsPGOcc` analogue; gcol33/tulpaObs#118): a
  varying-coefficient bar `spatial(~ 1 + w || cell, graph)` on the `ms_occu()`
  occupancy formula now fits an intercept field plus one spatially-varying-
  coefficient field per covariate,
  `logit psi_{s,i} = X_i (mu + b_s) + f0_{u(i)} + w_i f1_{u(i)}`. Fit by block
  coordinate ascent (`R/ms_occu_field.R`): the community occupancy Laplace-EM
  with the shared field as a psi offset, alternated with a two-state-marginal
  occupancy field update (`.ms_occu_field_solve`, a Newton over the per-(species,
  site) occupancy marginal + the ICAR prior). Both fields recover cleanly (field
  correlations ~0.96). `icar()` only; the plain single intercept field keeps the
  in-tree C++ community-spatial path (`sfMsPGOcc`), so that route is unchanged.
  Recovery-tested in `test-ms-occu-field.R`.

## 0.0.123 (2026-07-14)

* Community count field: `bym2()` (the Riebler scaled BYM2) now joins `icar()` and
  `car_proper()` as a shared-field kind for `ms_count()` (gcol33/tulpaObs#117). The
  field is the combined `phi = a v + b u` of a structured ICAR part `v` and an
  unstructured iid part `u`, with `a = sigma sqrt(rho/scale)`,
  `b = sigma sqrt(1 - rho)` and the Riebler scale factor; `(sigma, rho)` are chosen
  over a small grid by the field marginal (reported as `fit$spatial_hyper$sigma`
  / `rho`). Fit by a two-component joint Newton in the block-coordinate field step
  (`.ms_count_bym2_solve`). The single shared intercept field only (no
  varying-coefficient bar / group_var; use `icar()`/`car_proper()` there).
  Recovery-tested in `test-ms-count-spatial.R`. With this, the community count
  field supports the intrinsic (`icar`), proper (`car_proper`), and scaled BYM2
  (`bym2`) areal priors.

## 0.0.122 (2026-07-14)

* Community count field: `group_var` (sites > field cells) now works for
  `ms_count()` shared fields (gcol33/tulpaObs#117). When several sites share a
  spatial cell -- `icar(graph = cell_adj, group_var = "cell")` -- the field has
  one node per graph cell and a site->cell incidence aggregates the per-site
  working residuals; `fit$spatial_field` is the per-cell field. The shared field
  solve carries an optional incidence map `M` (identity when one node per site,
  so the one-node-per-site routes are byte-compatible). Recovery-tested in
  `test-ms-count-spatial.R`.

## 0.0.121 (2026-07-14)

* Community count field: `car_proper()` (a proper CAR) now joins `icar()` as a
  shared-field kind for `ms_count()` (gcol33/tulpaObs#117). The field is
  sum-to-zero deviations (the intercept owns the level); `icar()` fixes the
  dependence at the intrinsic limit, `car_proper()` estimates the dependence
  strength `rho` over a small grid by the field marginal likelihood (reported as
  `fit$spatial_hyper$rho`). Recovery-tested in `test-ms-count-spatial.R`. The
  field solve is shared with the shared-field / SVC / spatial-factor routes.

## 0.0.120 (2026-07-14)

* Spatial-factor `ms_count()`: a shared areal field `icar()` and latent factors
  `latent(n)` now compose on one formula (the spatial-factor `sfMsAbund`;
  gcol33/tulpaObs#117), `log mu_{s,i} = X_i (mu + b_s) + f_{u(i)} + sum_q
  lambda_{s,q} eta_{q,i}`. The block coordinate ascent runs the field update and
  the factor update in the same loop; when both are present the factor loadings
  are centred across species (`sum_s lambda_{s,q} = 0`) so the shared field owns
  the spatial mean and the factors own the between-species residual co-occurrence
  -- both recover cleanly (field correlation ~0.99, residual correlation ~0.99).
  The shared-field-only, factor-only, and combined routes are now one fitter
  (`.tobs_fit_ms_count_latent`, `R/ms_count_spatial.R`): single source of truth,
  with the previous fitters as thin special cases. Recovery-tested in
  `test-ms-count-factor.R`.

## 0.0.119 (2026-07-14)

* Community latent-factor `ms_count()` (the spAbundance `lfMsAbund` analogue,
  Poisson; gcol33/tulpaObs#117): a `latent(n_factors)` term on the abundance
  formula models residual species co-occurrence with Q per-site latent factors +
  per-species loadings, `log mu_{s,i} = X_i (mu + b_s) + sum_q lambda_{s,q}
  eta_{q,i}`. Fit by block coordinate ascent (`R/ms_count_factor.R`): the
  community Laplace-EM with the factor term as a per-observation offset,
  alternated with a Poisson factor update (alternating Newton on the factors and
  the loadings, with a unit-variance anchor). The loadings / factors are
  identified only up to rotation, so the recoverable target is the residual
  species covariance `Sigma_res = lambda lambda'` (reported as
  `fit$ms_factor$residual_cov` / `residual_cor`, with `loadings` and `factors`);
  the residual correlation recovers cleanly (correlation with truth ~0.95).
  `fitted()` / WAIC are factor-aware. Poisson; not composed with a shared areal
  field yet. Recovery-tested in `test-ms-count-factor.R`.

  With this, every tier of the spAbundance relative-abundance branch is covered:
  `abund` (`count()`), `spAbund` (`count()` + areal), `msAbund` (`ms_count()`),
  `sfMsAbund` (`ms_count()` + `icar()`), `svcMsAbund` (`ms_count()` + SVC bar),
  and `lfMsAbund` (`ms_count()` + `latent()`).

## 0.0.118 (2026-07-14)

* Community SVC for `ms_count()` (the spAbundance `svcMsAbund` analogue, Poisson;
  gcol33/tulpaObs#117, #118): a varying-coefficient areal bar
  `spatial(~ 1 + w || cell, graph)` on the abundance formula now fits an
  intercept field plus one shared spatially-varying-coefficient field per
  covariate, `log mu_{s,i} = X_i (mu + b_s) + f0_{u(i)} + w_i f1_{u(i)}`. The
  block-coordinate field step (`R/ms_count_spatial.R`) generalises to K
  covariate-weighted ICAR fields solved jointly (a K x K sparse block system with
  a per-field tau M-step); the intercept field is `fit$spatial_field`, the
  varying-coefficient field(s) `fit$trend_field(s)`. Both recover cleanly (field
  correlations ~0.98). Recovery-tested in `test-ms-count-spatial.R`. Together with
  the shared-field (`sfMsAbund`) case this closes the actionable community-SVC
  item of #118 for the abundance branch.

## 0.0.117 (2026-07-14)

* Community-spatial `ms_count()`: a shared areal field `icar()` on the abundance
  formula now fits under `method = "nested_laplace"` (the spAbundance `sfMsAbund`
  analogue, Poisson; gcol33/tulpaObs#117). One shared field across all species,
  `log mu_{s,i} = X_i (mu + b_s) + f_{u(i)}`. Fit by block coordinate ascent
  (`R/ms_count_spatial.R`): the community Laplace-EM with the field as a per-site
  offset, alternated with a self-contained Poisson-ICAR field update (an analytic
  sparse Newton + a closed-form tau M-step) -- pure R, no new C++. The field is
  informed by every species at each site, so it recovers cleanly (field
  correlation ~0.98) alongside the community means; `fitted()` and WAIC are
  field-aware. Poisson + `icar()` only (an overdispersed community count is not
  identified against a per-site field; `bym2()`/`car_proper()`, group_var, and
  the varying-coefficient bar are follow-ups). Recovery-tested in
  `test-ms-count-spatial.R` (community-mean recovery + pooled coverage + field
  recovery over 20 seeds).

## 0.0.116 (2026-07-14)

* Areal `count()` (negbin / Gaussian): the gate now states, and its error points
  out, that an overdispersed areal count with one field node per site is
  *fundamentally not identified* -- the latent field absorbs all extra-Poisson
  variation, so the field-integrated marginal likelihood is monotone in the
  dispersion toward the Poisson limit (verified over a dispersion grid). No
  estimator recovers the dispersion in this design; the message points to
  `abun()` (an N-mixture on replicated counts) for overdispersion with a spatial
  signal, or a Poisson areal `count()`.

## 0.0.115 (2026-07-14)

* New `ms_count()` family: the community / multispecies relative-abundance GLMM
  (the spAbundance `msAbund` model, gcol33/tulpaObs#117). Per-species GLMM on an
  observed count / continuous response with Gaussian community hyperpriors on the
  coefficients -- no detection, no latent state -- the abundance analogue of
  `ms_occu()` without the occupancy layer, and the community version of
  `count()`. Poisson / negative-binomial (per-species dispersion RE) / Gaussian.
  `y` is an `n_sites x n_species` matrix (or a named list of count vectors);
  `species` names the columns. Reuses the shared community Laplace-EM
  (`R/community_em.R`) -- pure R, no new C++ -- with a count `sp_ll` / `sp_grad`.
  Full S3 (`coef` / `vcov` / `confint` / `ranef` / `fitted` / `simulate` / WAIC)
  plus `simulate_ms_count()`. Recovery-tested in `test-ms-count.R`: community-mean
  recovery + pooled 95% coverage over 20 seeds (Gaussian and Poisson unbiased;
  the negbin slope carries a mild first-order-Laplace attenuation) + dispersion
  recovery. Non-spatial Laplace; the community-spatial (`sfMsAbund`) and NUTS
  variants are documented follow-ups.

## 0.0.114 (2026-07-14)

* Areal `count()`: a plain areal field -- `icar()` or `car_proper()` -- on the
  abundance formula now fits under `method = "nested_laplace"` (the `spAbundance`
  `spAbund` analogue, gcol33/tulpaObs#117). The count response is observed
  directly, so the fit is a single `tulpa::tulpa_nested_laplace()` call over the
  count block with the areal field as its latent GMRF prior: the field
  hyperparameters integrate on the nested outer grid, and the fixed effects plus
  their covariance are grid-integrated (law of total covariance). `fitted()` is
  field-aware in-sample. Poisson only -- with one field node per site a negbin
  size / gaussian residual variance is not jointly identified with the field
  under the fixed-dispersion loop, so those error with a pointer; `bym2()` (mixed
  structured/unstructured field) and the improper `car()` are likewise gated to
  `icar()` / `car_proper()`. Recovery-tested (fixed-effect coverage + field
  recovery over 20 seeds) in `test-count-spatial.R`.

## 0.0.113 (2026-07-14)

* New `count()` family: a GLMM on an observed count / continuous response
  directly, with no detection process and no latent state -- the abundance
  analogue of `jsdm()`, and the relative-abundance model of `spAbundance`
  (`abund`). Poisson / negative-binomial (log link) or Gaussian (identity), one
  value per site, supplied via `y =` or a two-sided `formula` left-hand side. The
  negative-binomial size and the Gaussian residual variance are estimated by an
  outer dispersion loop (reported in `fit$count_dispersion`). Full S3 surface
  (`coef` / `vcov` / `confint` / `fitted` / `predict` / `residuals` / `summary` /
  WAIC) plus `simulate_count()`. Non-spatial Laplace for this release; areal
  (`spAbund`), community (`msAbund`), and NUTS are the documented follow-ups
  (gcol33/tulpaObs#117). Recovery-tested (coefficient coverage + dispersion
  recovery over 20 seeds) in `test-count.R`.

* `svc()` on a family that does not consume it now errors with a pointer instead
  of silently dropping the term (gcol33/tulpaObs#118). The continuous NNGP `svc()`
  spatially-varying coefficient is wired only for single-season `occu()` under
  `method = "nuts"`; on the count families (and `occu()` under laplace /
  nested_laplace) it was extracted and discarded, fitting a model missing the term
  the user asked for. The recovery-tested route for a spatially-varying
  coefficient -- a weighted areal bar, `spatial(~ 1 + w || cell, graph)` with
  `method = "nested_laplace"` -- is unchanged.

## 0.0.112 (2026-07-14)

* Canonicalize the cover-hurdle direct-grid engine on the name `joint`
  everywhere, completing the rename the 0.0.111 message removal left partway.
  `control$engine = "joint"` is the documented value for the `occu_cover()` /
  `occu_multiscale_cover()` and standalone-occupancy spatial fits (the old
  `"joint_coupled"` string still routes to the same fitter), and the internal
  fitters, helpers, source files (`occu_cover_joint.R`, `occu_joint.R`,
  `occu_multiscale_cover_joint.R`), and the stored engine label all read `joint`
  now. No fitting behaviour changes.

## 0.0.111 (2026-07-14)

* Remove the last two deprecation nudges, completing the 0.0.110 clean-rename
  pass. `tobs()` no longer messages when the cover-hurdle state formula is given
  as `formula =` rather than `occurrence =` (both are accepted; `occurrence`
  reads symmetrically with `detection` / `positive`), and `occu_cover()` no
  longer messages on `control$engine = "joint_coupled"` (the documented engine
  name; the internal `"joint"` default routes to the same fitter). Nothing
  changes beyond dropping the two messages.

## 0.0.110 (2026-07-14)

* Drop the backward-compatibility shims added in 0.0.109. The extra-arm formula
  arguments of `tobs()` are the bare names only (`colonization` / `extinction`,
  `omega` / `gamma`, `p10` / `certainty`); the `<x>_formula` spellings and the
  `simulate_dyn_ms_occu()` / `simulate_int_ms_occu()` simulator aliases are gone,
  along with the deprecation-warning helper and the `tulpaObs-deprecated` page.
  There is no public release depending on the old names, so this is a clean hard
  rename rather than a deprecation cycle.

## 0.0.109 (2026-07-14)

* Unify the extra-arm formula arguments of `tobs()` on bare process / symbol
  names, matching `detection`, `positive`, and `availability`: `col_formula` ->
  `colonization`, `ext_formula` -> `extinction` (dynamic occupancy),
  `omega_formula` -> `omega`, `gamma_formula` -> `gamma` (open N-mixture), and
  `fp_formula` -> `p10`, `b_formula` -> `certainty` (false-positive occupancy;
  `certainty` is the `b` arm, renamed off `b` to avoid clashing with the `by`
  argument). The old spellings still work but emit a one-time deprecation warning.
* Rename the community simulators to match their family constructors:
  `simulate_dyn_ms_occu()` -> `simulate_ms_dyn_occu()` and
  `simulate_int_ms_occu()` -> `simulate_ms_int_occu()` (matching `ms_dyn_occu()`
  / `ms_int_occu()`). The old names remain as deprecated pass-throughs.
* Refresh the README: the model menu now lists all 19 family constructors, the
  install snippet no longer pins a stale tag, the documentation section links the
  full vignette set, and the lifecycle badge reflects the stable public API.

## 0.0.108 (2026-07-13)

* Fix the inherited S3 method surface on the multi-arm fit classes `cover_fit`
  (all `cover()` Laplace / nested-Laplace joint paths) and
  `occu_categorical_fit`. These fits carry independent per-arm coefficient blocks
  rather than a flat `$model` / `$draws` layout, so the inherited `tobs_fit`
  methods that branch on `object$model$model_type` errored (`nobs`, `residuals`,
  `simulate` crashed with "argument is of length zero"; `vcov`, `confint`,
  `logLik`, `tidy` errored in the tulpa summarizer; `fitted` returned
  silently-wrong values). A new intermediate class `tobs_multiarm_fit` supplies
  arm-aware `coef`, `vcov` (exact block-diagonal, since arms are fit
  independently), `confint`, `logLik`, `nobs`, `glance`, `tidy`, and `summary`,
  and turns `fitted` / `residuals` / `simulate` into informative refusals
  pointing at `predict()` / `simulate_<family>()`. The NUTS `cover()` path (which
  does carry `$draws`) defers to the tulpa draw-based methods unchanged.

## 0.0.107 (2026-07-12)

* Pin tulpa `>= 0.0.79`, which makes `integration = "grid_adaptive"` decline to
  the dense tensor before any inner solve on a small outer grid (default
  `control$adaptive_grid_min_cells = 48`), so the adaptive integrator is never
  slower than dense on the small hyperparameter grids typical of shared-trend
  occupancy fits.

## 0.0.106 (2026-07-12)

* Pin tulpa `>= 0.0.78` (Remotes `gcol33/tulpa@v0.0.78`), which adds the
  `integration = "grid_adaptive"` outer-grid integrator. Occupancy / cover joint
  fits (`occu_cover()`, `cover()`) forward `control$integration` unchanged, so
  `"grid_adaptive"` reaches the joint driver: on a sharply-peaked hyperparameter
  posterior (a strongly-identified field SD / Beta precision, the fine-grid
  regime) it evaluates only the mass-concentrated outer cells for a posterior
  that matches the dense tensor, and declines back to the tensor otherwise.

## 0.0.105 (2026-07-10)

* Pin the `tulpa` dependency to the memory-safe `v0.0.77` (the Imports floor and
  the Remotes ref), which sizes the joint outer-grid thread budget from the
  measured CHOLMOD factor rather than a fixed per-thread guess. The earlier
  Remotes pin (`@v0.0.70`) conflicted with installing `gcol33/tulpa` at HEAD, so
  a fresh `pak::pak(c("gcol33/tulpa", "gcol33/tulpaObs"))` could not solve the
  dependency graph; the two refs are consistent again.

## 0.0.104 (2026-07-08)

Missing-at-random cover on `occu_cover()`:

  A detected visit may now carry a missing cover value (`y_pos = NA`). The
  detection and occupancy arms keep the visit; only the Beta / lognormal cover
  factor drops out, so a species recorded as present but with the cover not
  measured contributes its detection without fabricating a cover. Missingness is
  taken as at-random given the model covariates: the cover likelihood is the
  product of the cover density over the observed detected visits only.

  The cover density gates on `detected & is.finite(cover)` uniformly across the
  builders (dense and compact), the joint nested-Laplace cell coupling, the
  non-spatial Laplace / NUTS fits, the aggregated and latent cover arms, and the
  WAIC / PSIS-LOO pointwise log-likelihood. With no missing cover the build and
  fit are byte-identical to before (the NA sentinel is never introduced).

## 0.0.103 (2026-07-08)

Performance (occu_cover joint fit):

  Detection-pattern compression on the coupled nested-Laplace fit. Within a
  site, all-undetected visits that share a detection design row enter the
  occupancy mixture only through prod_v (1 - p_v) = prod_u (1 - p_u)^{w_u}, so
  they collapse to one row of multiplicity w_u carried in the p arm's existing
  n_trials slot; detected visits stay individual (each keeps its own per-visit
  cover, so the cover arm and its alignment are untouched). This is exact
  sufficient statistics, not an approximation: the compressed fit matches the
  uncompressed fit to floating point (equivalence test in
  test-occu-cover-compact.R; on the MOT-scale detection design the 923k plots
  reduce to ~298k unique rows, and the large all-undetected cells collapse
  30-65x). On by default for the single-species path; off for the batched fused
  solve (per-species detection differs, so nodet rows are not exchangeable
  across species). getOption("tulpaObs.compress_nodet", TRUE) forces the
  uncompressed build for an equivalence check. No tulpa engine change (the
  weight rides the coupled arm's n_trials, already plumbed to CellResponse).

## 0.0.102 (2026-07-07)

New family:

* **`royle_nichols()`** -- Royle-Nichols occupancy (Royle & Nichols 2003;
  `unmarked::occuRN`): occupancy with abundance-induced detection heterogeneity,
  `N_i ~ Poisson(lambda_i)`, `detect ~ Bernoulli(1 - (1 - r_i)^{N_i})`. The latent
  `N` marginalises in closed form (a Poisson sum), so the exact marginal is
  maximised with an observed-information vcov. Site-level detection, `laplace`
  only for this first ship (`R/royle_nichols.R`). Parameter recovery validated
  over seeds (< 5% bias, ~0.92-1.00 CI coverage); full `fitted()` / `predict()` /
  `residuals()` / `tobs_waic()` / `tobs_dic()` / `tobs_cpo()` and a
  `simulate_royle_nichols()` generator.

Post-fit coverage: the S3 and diagnostic surface now spans the family set, closing
the gaps where whole family clusters had no prediction or goodness-of-fit.

* `tobs_waic()` / `tobs_dic()` / `tobs_cpo()` now work on the community-occupancy
  families (`ms_occu()`, `ms_dyn_occu()`, `ms_int_occu()`), which previously
  errored. The per-(species, site) marginal is scored over the community-mean
  pseudo-draws with the per-species BLUP deviations plugged in; the two-state
  marginal reproduces the C++ single-season kernel column-for-column
  (`R/community_ploglik.R`).
* `predict()` and `residuals()` now work on the community-occupancy families
  (per-species occupancy / detection, and per-species occupancy residuals)
  instead of raising a "not implemented" error.
* `jsdm()` gains `fitted()` / `predict()` / `residuals()` (per-species occupancy
  probability and presence residuals); previously `fitted()` errored.
* The count families (`abun()`, `removal()`, `distance()`, `dyn_abun()`) gain
  `tobs_test_dispersion()` / `tobs_test_zero_inflation()` / `tobs_test_outliers()`,
  scored on the per-site total count -- the natural overdispersion / excess-zero
  unit for an N-mixture-type model. They were previously gated to single-season
  occupancy.
* `tobs_predict_spatial()`, `predict(terms = )`, and `tobs_marginal_effect()` now
  apply the state process's inverse link, so a count fit returns the abundance
  intensity instead of a silently mis-linked logit-of-log-lambda.
  `tobs_marginal_effect(process = "abundance")` is now available.
* `tobs_richness()` now accepts all three community-occupancy families (not only
  `ms_occu()`); for the dynamic family it reports first-season richness.
* New simulator `simulate_jsdm()`, completing the per-family simulator set.

## 0.0.101 (2026-07-07)

* Adds three packaged example datasets (`foray_counts`, `meadow_cover`,
  `peatland_occu`) with a `data-raw/make_datasets.R` builder and `LazyData`,
  so the family help pages and vignettes can run on shipped data.
* Splits the large family and dispatch sources (`family_cover_hurdle.R`,
  `laplace.R`, `occu_cover.R`, `occu_cover_joint_coupled.R`, `tobs.R`,
  `ms_occu_cover_spatial.R`) into focused modules (dispatch, arms,
  postprocess, decode, diagnostics, helpers, callbacks, simulate). No
  user-facing behaviour change; the fast smoke suite is unchanged.
* Pins `LinkingTo: gcol33/tulpa@v0.0.70`.

## 0.0.100 (2026-07-06)

* `cover()` gains per-component `copy(spatial(), terms = list(intercept = ,
  trend = ))`, matching `occu_cover()`'s grammar: the presence field's intercept
  and weighted-trend blocks can now couple onto the positive arm with
  independent amplitude grids (`control$alpha.grid` / `control$alpha.grid.trend`)
  instead of one shared amplitude. The whole-field forms
  (`copy(spatial())`, `copy(spatial(), alpha = )`) are unchanged and stay
  byte-identical. Per-component `terms =` must address every field block, and a
  `trend` component named against an intercept-only field errors.
* Tests: the `cover()` copy-coupling amplitude handoff (the alpha ->
  `control$alpha.grid[.trend]` path) is now covered end-to-end -- whole-field
  `copy(alpha = )` byte-matches the shared-formula + control spelling,
  per-component `terms =` with equal grids byte-matches the whole-field form,
  and a decoupled `terms =` byte-matches the explicit control grids while
  differing from the coupled fit. The occu_cover cover-arm INTERCEPT field SD
  (the base per-cell cover-arm map, the direct deliverable) now has a
  parameter-recovery test, alongside the existing trend-field recovery.
* Docs: `?tobs` disambiguates the two `response` arguments -- the long-frame
  pivot `response =` (a column of `data`) versus the family constructor's
  `response =` (the positive-part distribution), which can coexist in one call.

## 0.0.99 (2026-07-06)

* `tobs()` documents the `positive` argument (the cover-hurdle positive-arm
  formula) under its own name. The roxygen entry carried the stale `response`
  label left over from the 0.0.95 family-constructor rename, so the generated
  help documented a `response` argument absent from the signature and left the
  real `positive` formal undocumented. Documentation only; no change to fitting.

## 0.0.98 (2026-07-03)

* `glance()` now reports `pareto_k_is_ess` as the numeric importance-sampling
  effective sample size the engine computes (matching `fit$pareto_k_is_ess` at
  the top level), instead of coercing it to a logical. The field was documented
  and consumed as a boolean "the k-hat column is a quad-ESS fallback" flag, but
  the engine has no such fallback: `pareto_k` is always the outer k-hat (or `NA`)
  and `pareto_k_is_ess` is always the IS-ESS on the PSIS-smoothed weights, so
  `as.logical()` discarded the number (a healthy IS-ESS of, e.g., 110 became
  `TRUE`). `pareto_k_is_ess / control$k.samples` is the relative IS efficiency.

## 0.0.97 (2026-07-03)

* `to =` is fully decoupled from the spatial-field call surface in both cover
  families. Placement now tags each field's arm on the evaluated spec directly,
  instead of round-tripping the arm through a deparsed formula, so `to =` is no
  longer an argument the bar form reads: the `spatial(~ ... || node)` bar takes
  only `graph` and `by`, and an explicit `to =` is an unknown argument. The arm
  is chosen by placement (write the field in that arm's formula) and shared
  across arms with `copy()`, exactly as in 0.0.96; fits are byte-identical (the
  cover placement-equals-shared and `occu_cover()` recovery suites are
  unchanged). The bespoke `to =` deprecation guard is removed.

## 0.0.96 (2026-07-02)

* `copy()` is now wired into the `cover()` hurdle engine, so a shared spatial
  field across the two cover arms is written the same way as in `occu_cover()`:
  place the field in the `presence` formula and `copy(spatial())` in the
  `positive` formula. `copy(spatial())` estimates the coupling amplitude on the
  default grid (byte-identical to the previous both-arm spelling);
  `copy(spatial(), alpha = grid(c(...)))` integrates it over a supplied grid and
  `copy(spatial(), alpha = 0.5)` fixes it.

* `to =` on a spatial field is retired from both cover families. An arm is now
  chosen by placement (write the field in that arm's formula) and a field is
  shared across arms with `copy()`; a field in the single shared `formula`
  reaches both arms. The `||` / `|` (independent / correlated MCAR) axis is
  unchanged. Every routing `to =` expressed is byte-identical under the new
  spelling: shared intercept/trend via a shared-formula field or `copy()`,
  arm-specific and correlated single-arm fields via placement, MCAR both-arm via
  a shared-formula `|` bar or `copy()`. A user-supplied `to =` now errors with a
  pointer to placement / `copy()`.

## 0.0.95 (2026-07-02)

* The positive-arm distribution argument of the cover families is renamed
  `positive` -> `response`: `cover(response = "beta")`,
  `occu_cover(response = "beta")`, `ms_occu_cover()`, `occu_multiscale_cover()`.
  This removes the name clash with the `positive = ~ x` arm formula on a
  `tobs()` call, where `positive` otherwise meant both the positive-cover arm
  formula and the positive-arm response distribution. Positional calls
  (`cover("beta")`) are unchanged; a named `positive =` now errors. The stored
  family field and the whole fitting engine are unchanged, so fits are
  byte-identical.

## 0.0.94 (2026-07-02)

* `cover()` per-arm formulas now carry spatial fields by placement, completing
  "arm = formula" for both hurdle families: a `spatial()` term in `presence` or
  `positive` becomes an arm-specific field on that arm, byte-identical to the
  shared-formula `to =` spelling (the field is routed through the same machinery
  and indexed onto the arm's rows by the fitter). Fixed effects and fields both
  follow placement; `copy()` reuses a named field across arms. `temporal()` /
  `re()` in a per-arm formula are still declared on the shared `formula`. Closes
  gcol33/tulpaObs#111.

## 0.0.93 (2026-07-02)

* `cover()` per-arm formulas (arm = formula): `cover(presence = ~ x, positive =
  ~ t)` gives the presence and positive hurdle arms their own fixed effects (two
  independent designs), matching `occu_cover()`'s per-arm formulas. The single
  shared `formula` stays the back-compat spelling (byte-identical; the full cover
  test suite is unchanged). First cut: per-arm formulas carry fixed effects only;
  declare fields on the shared `formula` with `to =` (per-arm field placement is
  gcol33/tulpaObs#111).

* Bugfix: an lme4 random-effect bar on the `occu_cover()` detection or
  positive-cover formula (`detection = ~ det_cov1 + (1 | habitat)`, and the
  `(x | g)` / `(0 + x | g)` slope spellings) is fitted again. The arm-field
  lifting introduced in 0.0.90 round-tripped each arm formula through
  `terms()` / `reformulate()`, which dropped the parentheses off a bar
  (`(1 | g)` -> `1 | g`); R then re-parsed `~ det_cov1 + 1 | g` as
  `(det_cov1 + 1) | g`, so the downstream RE parse saw no bar and silently
  dropped the random effect. The lift now re-parenthesizes bar terms it keeps.
  The `re(g)` spelling was unaffected. Restores gcol33/tulpaObs#102, #103.

## 0.0.92 (2026-07-02)

* `occu_cover()` detection-arm spatial field: a spatial-field term on the
  `detection` formula (`detection = ~ 1 + spatial(~ 0 + time || cell, graph =
  adj)`, or the explicit `to = "detection"` spelling) now fits a
  spatially-structured detection probability. The non-copied field block scatters
  onto the detection predictor by riding the detection arm with `field_coef = 1`
  while the shared occupancy field is kept off detection by the `spatial_idx = 0`
  sentinel -- the same mechanism the detection random effect uses
  (gcol33/tulpaObs#102), so no engine change was needed (closes gcol33/tulpa#140).
  Recovery of the detection field SD is tested across seeds.

## 0.0.91 (2026-07-02)

* `occu_cover()` field placement: a spatial-field term written in the `positive`
  formula (`positive = ~ t + spatial(~ 0 + time || cell, graph = adj)`) declares
  an independent cover-arm field by placement, byte-identical to the explicit
  `to = "positive"` spelling (which is retained). The intake lifts the term onto
  the occurrence formula with the arm tag; `copy()` stays on its arm's formula.
* The arm-specific field machinery (parse, block build, per-arm sigma naming) is
  now arm-generic, and `simulate_occu_cover(det_field = TRUE)` injects a known
  detection field. A detection-arm field (`to = "detection"`) is reserved but
  rejected at parse: the joint C++ substrate scatters fields onto the occupancy
  and cover arms only, so a detection field is unidentified until the substrate
  is wired (gcol33/tulpa#140).

## 0.0.90 (2026-07-02)

* `occu_cover()` now accepts `control$sigma.grid.pos.field`, the SD grid for the
  independent cover-arm field (`to = "positive"`, #110). The joint engine already
  read it and defaulted it to `control$sigma.grid`; it was missing from the
  family's control whitelist, so passing it errored. Setting a coarse field grid
  now keeps the added axis from multiplying the outer-grid cost.

## 0.0.82 (2026-07-01)

* The aggregated / latent-mode occu_cover() posterior predictive check moved its
  per-draw loop into C++ (`cpp_occu_cover_ppc_agg`): the detection replicate plus
  the aggregated (mean / median) or shared cover-RE (latent) cover replicate,
  drawing from R's RNG stream in the former order (byte-identical). With this,
  every posterior-predictive / PIT diagnostic that draws from posterior
  components -- occu_cover (all cover modes), cover, single-season -- is a C++
  kernel. The remaining R generators (`simulate()`) are left in R by design: they
  select a posterior draw with `sample.int`, whose index sampler is not
  byte-reproducible in C++, so porting them would change results for no runtime
  gain -- the same boundary the package already draws for Monte Carlo generation.

## 0.0.81 (2026-07-01)

* The leave-one-out PIT weighting (`.tobs_loo_pit_from_limits`) moved its
  per-observation PSIS loop into tulpa's C++ `cpp_psis_loo_pit` (PSIS columns
  parallel, the jitter in index order), byte-identical to the former R loop.
  Requires tulpa >= 0.0.64.

## 0.0.80 (2026-07-01)

* The cover() hurdle and single-season occupancy posterior diagnostics moved
  their per-draw / per-site loops into C++: cover PIT (`cpp_cover_pit_cdf`,
  deterministic) and PPC (`cpp_cover_ppc`); single-season PPC (`cpp_single_ppc`)
  and randomized PIT residuals (`cpp_single_pit`). The replicate draws come from
  R's RNG stream via the R:: samplers in the SAME order as the former R loops
  (the posterior-draw selection `sample.int` stays in R), so under a fixed seed
  every result is byte-identical. All four positive families for the cover PIT.

## 0.0.79 (2026-07-01)

* The `occu_cover()` posterior diagnostics moved their per-draw loops into C++:
  the deterministic detection-summary CDF limits (the randomized-PIT / LOO-PIT
  building block, `cpp_occu_cover_cdf_limits`, parallel over draws) and the
  posterior predictive check (`cpp_occu_cover_ppc`, cover_aggregate = "none").
  The PPC draws the latent state, detection replicate, and cover replicate from
  R's RNG stream via the R:: samplers in the SAME order as the former R loop, so
  under a fixed seed the discrepancy is byte-identical; the aggregated cover
  modes keep their R path.

## 0.0.78 (2026-07-01)

* The spatial-factor community occupancy + cover family
  (`ms_occu_cover_spatial`) WAIC / DIC / LOO pointwise log-likelihood moved its
  per-draw loop into C++ (`cpp_ms_ocs_ploglik`). The former R loop unpacked the
  NUTS draw (community mean, per-species deviation, shared fields `W`, occupancy
  loadings `L`, optional cover loadings `Lpos`, log-dispersion), assembled each
  species' predictors with the shared-factor offset `W L[s,]` on occupancy (and
  `W Lpos[s,]` on cover), and evaluated the dense occu_cover per-cell marginal;
  the kernel does all of that, parallel over draws. Both constrained and
  unconstrained loading parameterisations canonicalise to one packed layout.
  Byte-close (~1e-15) to the R oracle, thread-count invariant. With this, **every
  family's pointwise log-likelihood in the WAIC / DIC / LOO path is a C++
  kernel.**

## 0.0.77 (2026-07-01)

* The community N-mixture (`ms_abun`) WAIC / DIC / LOO pointwise log-likelihood
  moved its per-draw loop into C++ (`cpp_ms_nmix_ploglik_batch`). The former R
  loop reconstructed each species' deviation `b = C z` from the non-centered
  NUTS draw (a log-Cholesky factor per arm) and called the per-species Royle
  marginal in R; the kernel now does the log-Cholesky reconstruction and the
  per-(species, site) marginal (`compute_nmix_site`) internally, parallel over
  draws. Byte-identical to the former loop (same kernel), thread-count invariant.

## 0.0.76 (2026-07-01)

* The three-level multiscale occupancy + cover family
  (`occu_multiscale_cover`) builds its per-cell pointwise log-likelihood in a
  C++ OpenMP kernel (`cpp_occu_ms_cover_ploglik`), parallel over draws. This was
  the last pure-R marginal in the criteria path (the cell -> plot -> visit
  three-level mixture plus the per-detected-visit cover density, formerly an
  `apply()` over draws). The kernel mirrors the R reference
  `.occu_ms_cover_nonspatial_ll` draw for draw (~1e-14) and is thread-count
  invariant; both cover families and both visit-block layouts. WAIC / DIC / LOO
  for this family now honour `n.threads`.

## 0.0.75 (2026-07-01)

* The count / multistate families' WAIC / DIC / LOO pointwise log-likelihood
  moved its per-draw loop from R into C++: N-mixture (`cpp_nmix_ploglik_batch`),
  removal (`cpp_removal_ploglik_batch`), false-positive occupancy
  (`cpp_fp_occu_ploglik_batch`), binned distance sampling
  (`cpp_distance_ploglik_batch`), and open N-mixture / dyn_abun
  (`cpp_dyn_abun_ploglik_batch`). The per-site marginal was already C++
  (`compute_*_site`); each kernel now loops draws over it directly, parallel over
  draws, so the result is byte-identical to the former R loop (0 difference,
  same kernel) and thread-count invariant. `n.threads` reaches these through
  `.tobs_ploglik_from_draws`.

## 0.0.74 (2026-07-01)

* The remaining pure-R pointwise log-likelihood loops now run in C++ OpenMP
  kernels, parallel over the observation index, completing the WAIC / DIC / LOO
  compute port:
  - The dense (padded `[n_sites x max_visits]`) `occu_cover()` no-aggregation
    path is flattened to the ragged form and shares the 0.0.72 kernel, so every
    `occu_cover()` fit (compact or dense) builds its pointwise log-likelihood in
    parallel.
  - The draw-matrix occupancy families -- single-season (`cpp_occu_single_ploglik`),
    multi-season dynamic HMM (`cpp_occu_dynamic_ploglik`), and multi-source
    integrated (`cpp_occu_integrated_ploglik`) -- moved their per-observation
    marginal out of R. The R functions keep the (BLAS) linear-predictor
    assembly and the draw-invariant count gathering; the marginal runs in the
    kernel. `n.threads` is plumbed through `.tobs_ploglik_from_draws`.
  Each kernel mirrors its former R loop (reproduced as the test oracle),
  agreeing to libm rounding (~1e-15) and thread-count invariant.

## 0.0.73 (2026-07-01)

* The `cover()` hurdle pointwise log-likelihood (the WAIC / DIC / LOO input) now
  runs in a C++ OpenMP kernel parallel over posterior draws
  (`cpp_cover_hurdle_ploglik`), covering all four positive families (lognormal,
  lognormal_trunc, ordinal, beta). Same treatment as the `occu_cover()` ragged
  path in 0.0.72: `tobs_waic()` / `tobs_dic()` / `tobs_cpo()` on a `cover_fit`
  honour `n.threads`. The kernel mirrors the R reference `.tobs_cover_hurdle_ll`
  (retained for the posterior-mean plug-in and the tests), agreeing to libm
  rounding (~1e-14) and thread-count invariant.

## 0.0.72 (2026-07-01)

* WAIC / DIC / LOO for the compact (ragged) `occu_cover()` fit now build the
  pointwise log-likelihood in a C++ OpenMP kernel that parallelises over
  posterior draws (`cpp_occu_cover_ploglik_ragged`). This was the dominant
  serial cost on large single-species fits: the draw loop scales with
  `n.draws x total plots` and previously ran single-threaded in R. `tobs_waic()`,
  `tobs_dic()`, and `tobs_cpo()` gain an `n.threads` argument (default: all but
  four logical cores, matching the occu_cover fit's own outer-grid default). The
  kernel mirrors `.occu_cover_site_ll_ragged` draw for draw and accumulates each
  site's visit sums in visit order, so it agrees with the R path to libm rounding
  (~1e-13) and is thread-count invariant. On ~200k visits x 200 draws the pointwise
  build drops from 17.9 s (R) to 0.76 s (16 threads); the dense and aggregated
  (mean / median / latent) cover paths are unchanged.

## 0.0.71 (2026-07-01)

* `occu_cover()` can now give the cover (positive) arm its OWN independent
  spatial field, decoupled from the occupancy field's alpha copy
  (gcol33/tulpaObs#110). An arm-specific `spatial()` bar with a single
  `to = "positive"` on the occurrence formula --
  `spatial(~ 1 + time || cell, graph = adj, to = "positive")` -- adds a
  non-copied ICAR block on the cover arm (intercept + per-covariate trend
  fields, each with its own precision integrated on the outer grid). It composes
  with the shared occupancy field: occupancy still drives psi and, via the alpha
  copy, `delta_cover_exp`, while the independent cover field carries a
  cover-specific structure the alpha copy cannot express -- so `delta_cover_cond`
  is spatially varying instead of collapsing to a global slope when `alpha -> 0`.
  Reported as `sigma_pos_field` / `sigma_pos_field_<col>` with the per-cell
  posterior in `fit$pos_field` / `fit$pos_field_table`. `simulate_occu_cover()`
  gains `pos_field` / `sigma_pos_int` / `sigma_pos_trend` to simulate it. Scope:
  per-visit cover (`cover_aggregate = "none"`); ICAR only (bym2/car read as
  ICAR); not composed with the correlated `|` MCAR field, the latent cover RE,
  or the batched fused path.

## 0.0.70 (2026-06-30)

* A correlated (`|`, free-Sigma MCAR) intercept + slope spatial field can now sit
  on a single cover arm: `spatial(~ 1 + time.sc | cell, graph = adj, to =
  "presence")` puts a free 2x2 cross-coefficient Sigma on the occurrence arm
  alone (the occupancy intercept and time-slope fields covary), with no
  cross-arm copy. Previously the `|` form was copy-only and required both arms.
  The single-arm field uses the 0-sentinel `spatial_idx` on the other arm and
  carries no `alpha` (reported as NA); `sigma_mcar` / `rho_mcar` recover the
  field SDs and their correlation (gcol33/tulpaObs#109).

## 0.0.69 (2026-06-30)

* `cover(positive = "beta_oi")` -- a one-inflated Beta cover family. Plots
  recorded at exactly full cover (`y = 1`) are modelled as a genuine point mass
  instead of being clamped to `1 - 1e-6` (which biases the interior precision).
  With a constant inflation probability the likelihood factorizes: `pi` is the
  share of positive plots at the ceiling (a binomial proportion, reported as
  `pi_one` with its SE), and the interior Beta is fit on the `(0, 1)` plots.
  `predict()` returns the one-inflated conditional cover `pi + (1 - pi) * mu`.
  Works on the non-spatial and nested-Laplace (areal) paths
  (gcol33/tulpaObs#108).

## 0.0.68 (2026-06-30)

* Arm-specific cover-arm spatial fields accept `model = "bym2"`:
  `spatial(~ 1 || cell, graph = adj, to = "positive", model = "bym2")` now puts a
  BYM2 (structured ICAR + iid) field on the cover arm, previously restricted to
  icar / car / car_proper. The block fits as a non-copied length-2 latent over
  the paired (sigma, rho) grid, and the joint-draw projection reconstructs the
  rho-mixed unit field so `predict()` / WAIC see the full mix. Recovers the
  field and the structured fraction (gcol33/tulpaObs#107).

## 0.0.67 (2026-06-30)

* `cover()` / `occu_cover()` spatial fits accept `control$prior.phi`, a
  regularizing hyperprior on the cover-arm dispersion grid (the Beta precision
  under `positive = "beta"`, the log-scale SD under `lognormal`). Forwarded to
  tulpa's new `prior_phi`, it re-weights the `phi.grid` axis by the chosen
  density (`list("pc.prec", c(U, alpha))` / `list("half_normal", scale)`)
  instead of an implicit flat prior (gcol33/tulpa#139). Requires tulpa >= 0.0.62.

## 0.0.66 (2026-06-23)

* `tulpa` import floor raised to `>= 0.0.61`, the version the engine is built and
  tested against (the CCD outer-integration path and the `adjacency()` graph
  front door). The previous `>= 0.0.57` floor could pair this release with a
  tulpa too old for those, so an install that did not also upgrade tulpa would
  resolve to a broken combination. Metadata only; no code change.

## 0.0.65 (2026-06-23)

* `predict(type = "change")` for the joint cover-family (`occu_cover()`) and the
  rerouted standalone `occu()` SVC fit now reports the per-cell change certainty,
  not just the change. The change table gains, per cell: the start / end
  occupancy (`p_T1` / `p_T2`, or `psi_T1` / `psi_T2` for `occu()`) with their own
  `.sd` / `.lwr` / `.upr` interval, and a `.prob_pos` column per headline delta
  (`delta_p`, `delta_cover_cond`, `delta_cover_exp`; `delta_psi` for `occu()`)
  giving the directional posterior probability `P(delta > 0)`. All are taken over
  the grid-integrated draws, so they carry the joint posterior rather than a
  plug-in of the means, and they are pure additions (existing columns unchanged).
  This makes the two-point change prediction a complete drop-in for a per-cell
  occupancy-trend summary (start, end, change, direction certainty) without a
  hand-rolled trend pass downstream.

## 0.0.64 (2026-06-22)

* `tobs_data()` gains `type = "positive"`: a positive-real `(0, Inf)` response
  for the lognormal / gamma cover arm, alongside `type = "cover"` (a `[0, 1]`
  proportion, the beta arm). Both share the floor-to-`NA` absence policy; only
  `"cover"` enforces the upper bound (the shared validation lives in one
  `.tobs_floor_continuous()` helper). The long-frame builders pick the storage
  type from the family's positive distribution, so a lognormal / gamma cover
  that exceeds 1 now round-trips through `occu_cover_inputs()`, the single-fit
  long-frame path, and `by=` instead of being rejected by the `[0, 1]` check.
  `occu_cover_inputs()` exposes this as `positive = ` (default `"beta"`).
* A single `occu_cover()` fit now accepts a long / plot-level frame directly,
  the same contract the `by=` batch path already supported for many species:
  pass `site = `, `visit = `, `response = ` (the 0/1 detection column) and
  `y_pos = ` (the cover column), plus any visit-level `det.covs = `, with a long
  `data`, and `tobs()` builds the paired occurrence / cover arms and the
  site-level design internally. This removes the hand-built `tobs_data()` x2 plus
  the manual occurrence/cover alignment check from user scripts. The arms default
  to the compact (ragged) layout on the nested-Laplace route (no per-site visit
  cap); set `control$compact = FALSE` for the dense grid. The fit is byte-
  identical to the hand-built route (asserted in `test-occu-cover-long.R`).
* New `occu_cover_inputs()` exposes that builder: it returns the `y` / `y_pos` /
  `visits` / `site_data` bundle from a long frame so the arms can be inspected
  before fitting. The long -> arms construction (`.occu_cover_response_pair()`)
  is now single-sourced and shared by the `by=` batch loop and the single-fit
  path, rather than re-implemented per call site.
* `tobs(by = )` for `occu_cover()` now builds its per-species response arms and
  the shared visit grid compactly on the nested-Laplace route (its default
  looped backend), instead of allocating B dense `[n_sites x max_visits]`
  response matrices plus a dense visit grid. The padded-grid memory was the
  reason the batch could not run on the uncapped EVA data the single-fit path
  handles; the looped batch now scales the same way (memory is O(observations),
  not the padded grid). Each species' fit is unchanged to the compact-vs-dense
  tolerance (asserted in `test-occu-cover-by.R`). The non-joint `laplace` route
  stays dense (its engine reads the padded grid), as does the opt-in fused
  backend (`control$batch.backend = "fused"`), which stacks dense per-species
  columns by design; `control$compact` overrides the default.

## 0.0.63 (2026-06-22)

* `tobs_waic()` / `tobs_dic()` for `occu_cover()` now build the pointwise
  log-likelihood in memory-adaptive draw-chunks. The heavy transient was the two
  `[n_plots x n_draws]` per-visit predictor matrices (about 15 GB at 1000 draws on
  the full no-cap EVA data); `.occu_cover_ploglik_core()` now processes the draws
  in blocks sized to a fraction of free RAM (`/proc/meminfo` on Linux, the
  optional `ps` package elsewhere, a 4 GB default otherwise), so the peak stays
  bounded on a memory-tight node. WAIC is a sum over draws, so the result is
  byte-identical to the unchunked computation regardless of chunk size (asserted
  in `test-occu-cover-compact.R`). The full draw count is kept, so WAIC precision
  is unchanged -- callers no longer need to trade draws for memory.

## 0.0.62 (2026-06-22)

* The compact (ragged) `occu_cover()` path now carries an observation-arm random
  effect (`(1 | g)` on the detection or positive-cover formula), so a
  random-habitat detection fit runs uncapped like the fixed-effects spec. The RE
  group codes (and any random-slope design) are resolved over the V valid visits
  directly -- site-level groupings broadcast via the visit's site -- so the
  compact fit is byte-identical to the dense fit, BLUPs included (asserted in
  `test-occu-cover-compact.R`). Removes the earlier guard that errored on
  compact + observation RE.

## 0.0.61 (2026-06-22)

* DESCRIPTION `Remotes:` no longer pins `gcol33/tulpa@v0.0.50`. That tag predates
  the `tulpa (>= 0.0.57)` import floor (added in 0.0.60), so a fresh
  `pak::pak("gcol33/tulpaObs")` resolved tulpa to v0.0.50 and failed the version
  constraint. The `gcol33/tulpa` remote now tracks the default branch (where
  tulpa and tulpaObs are released together) so it cannot go stale against the
  import floor again; `gcol33/tulpaMesh` stays pinned at `@v0.1.3` to match
  tulpa's own remote (a differing tulpaMesh ref across the two would itself be a
  resolver conflict).

## 0.0.60 (2026-06-22)

* `tobs_data(compact = TRUE)` builds a compact (ragged) `tobs_data`: the response
  is stored as one row per valid site-visit (a `tobs_ragged` carrier) and each
  detection covariate as a length-V vector, instead of a padded
  `[n_sites x max_visits]` matrix. The joint nested-Laplace `occu_cover()` engine
  consumes the valid visits directly (it compacts the dense grid to exactly these
  rows anyway), so a fit on compact input is byte-identical to the dense fit on
  the same data -- now asserted in `tests/testthat/test-occu-cover-compact.R`
  (means, sds, `predict()`, and WAIC all match). Because the compact layout never
  materialises the padded grid, its memory is the number of observations rather
  than `n_sites x max_visits`, so a site with tens of thousands of visits no
  longer needs a per-site visit cap before the data can be built. Scoped to the
  joint nested-Laplace path; other engines, an observation-arm random effect, and
  cell-aggregated cover read the dense grid and error clearly on compact input.
* The `occu_cover` cell-coupling spec declares its dense cross-Hessian pairs
  (`dense_cross_pairs()`), so the joint engine (tulpa >= 0.0.57) no longer
  reserves a `J x J` slab for the all-undetected `(p, p)` block (it is the rank-1
  self-cross) or the always-zero `(p, pos)` / `(pos, pos)` blocks. A cell with `J`
  visits is now `O(J)`, not `O(J^2)`, which is what lets the uncapped compact path
  fit a grid with a cell holding tens of thousands of plots without running out of
  memory. The DESCRIPTION `tulpa` floor moves to `>= 0.0.57`.

## 0.0.59 (2026-06-19)

* `tobs()` exposes the single-batch + bootstrap outer Pareto-k controls for the
  joint-coupled families (`occu_cover()`, `cover()`, `occu()`, multiscale),
  forwarding them to the engine (tulpa >= 0.0.50, gcol33/tulpa#127):
  `control$diagnose.draws` (the precision knob; legacy `k.samples` accepted as an
  alias), `control$k.bootstrap`, `control$k.tail.points`, `control$k.conf.bands`.
  Replaces the removed `k.batches` / `k.adapt` / `k.batches.max` (#123/#124). The
  fit's `$joint_fit` carries `pareto_k_se_boot`, `pareto_k_ci_low` /
  `pareto_k_ci_high`, `pareto_k_se_formula`, `pareto_k_tail_points` /
  `pareto_k_tail_points_requested`, `pareto_k_band_confident`, and the top-level
  `diagnose_draws` / `diagnose_cost_ratio`. For a tighter k raise `diagnose.draws`,
  not `k.bootstrap`.

## 0.0.58 (2026-06-19)

* `tobs()` exposes `control$k.adapt` + `control$k.batches.max` for the
  joint-coupled families (`occu_cover()`, `cover()`, `occu()`, multiscale),
  forwarding them to the engine's adaptive batched outer Pareto-k
  (gcol33/tulpa#124). With `k.adapt = TRUE` the batch count grows from
  `k.batches` until the reliability band resolves (the k-hat band lands in one
  band) or the `k.batches.max` cap is reached; default off (`k.adapt = FALSE`).
  Requires tulpa >= 0.0.49.

## 0.0.57 (2026-06-19)

* `tobs()` exposes `control$k.batches` for the joint-coupled families
  (`occu_cover()`, `cover()`, `occu()`, multiscale), forwarding it to the engine's
  batched outer Pareto-k (gcol33/tulpa#123). With `k.batches > 1` the fit reports
  the median outer k-hat over that many independent importance batches plus the
  observed `pareto_k_lo` / `pareto_k_hi` range (the diagnostic's Monte Carlo
  spread, not a posterior CI); the band is classified off the median. Default
  `1L` (single batch), so existing fits are unchanged. Requires tulpa >= 0.0.48.

## 0.0.56 (2026-06-19)

* New `occu_categorical()` family (gcol33/tulpaObs#106): a presence + nominal
  (unordered) K-class hurdle. Each unit is absent (`y = 0`) or present in one of
  K classes (`y` in `1..K`); presence is a Bernoulli arm and the class given
  present is a baseline-category multinomial logit (the last class the baseline),
  the two arms factorising the likelihood exactly. This is the categorical
  counterpart of `cover()` (presence + magnitude) and the K-class generalisation
  of `fp_occu()` (its K = 2 confusion case): for an *ordered* class response use
  `cover(positive = "ordinal")`; this family is for classes with no ordering
  (colour morph, microhabitat use, classifier label). The multinomial math is
  the FD-validated tulpa kernel (`multinomial_logit.h`); the non-spatial Laplace
  fit is the vectorised R Newton over the same closed forms. Ships with
  `simulate_occu_categorical()`, a `predict()` method (presence, conditional, and
  joint class probabilities), and parameter-recovery tests. Non-spatial Laplace
  for this first ship; spatial fields / NUTS (the native multi-process
  likelihood) and the latent-class misclassification variant are documented
  follow-ups. Requires tulpa (>= 0.0.47).

* New `cover(positive = "lognormal_trunc")` positive arm (consumer of
  gcol33/tulpa#122): an upper-truncated lognormal for bounded cover -- a Gaussian
  on log-cover upper-truncated at `log(1) = 0`, so it cannot place mass above
  cover = 1 the way `"lognormal"` can. It rides the joint nested-Laplace path on
  the engine's new `truncated_gaussian` family via a per-plot truncation ceiling,
  threaded through encode / decode, the pointwise log-likelihood, the PIT CDF,
  and posterior-predictive replication (inverse-CDF truncated-normal draws), with
  the truncated-lognormal conditional cover mean in prediction. Requires
  `method = "nested_laplace"`. Recovery-tested against bounded cover data, and
  verified close to the plain lognormal fit under negligible truncation.
  Requires tulpa (>= 0.0.47).

## 0.0.55 (2026-06-18)

* `tobs_waic()` / `tobs_cpo()` gain `loo.unit = c("obs", "cell")`: the
  cross-validation unit. The default `"obs"` is the family's pointwise unit (one
  plot for `cover()`, one site for `occu_cover()`) and is byte-identical to the
  previous call. `"cell"` switches to leave-one-group-out cross-validation
  (LOGO-CV): the fit's own per-observation cell map is supplied to
  `tulpa::tulpa_criteria(group = )`, so each spatial cell is one fold instead of
  each plot / site, without the caller hand-building the map. Implemented for
  `cover()` (the areal field node via `spi_full`, when plots are grouped with
  `group_var`) and `occu_cover()` (the `site_cell` map); a non-spatial fit has no
  cells and errors with a pointer to `loo.unit = "obs"`. The per-family
  column -> cell map is the single `.tobs_loo_cell_map()` (each family's
  pointwise builder fixes its own column order). Equivalent to passing
  `group = ` the cell map directly. Requires tulpa (>= 0.0.45). (tulpaObs#105)

## 0.0.54 (2026-06-18)

* New `cover(positive = "ordinal", breaks = ...)`: an interval-censored Gaussian
  positive arm for Braun-Blanquet (ordinal cover-class) data, the
  measure-invariant counterpart of the `lognormal` / `beta` arms. The cover is
  recorded only as the ordinal class it falls in, so the latent log-cover is
  `Normal(eta, sigma^2)` censored to the observed class band and the likelihood
  is the class probability MASS -- an ordered probit with KNOWN thresholds, no
  free cutpoints and no change-of-variable Jacobian. `breaks` are the interior
  class boundaries on the (0, 1) cover-fraction scale (open outer classes added
  automatically); class bands are lower-closed so a representative sitting on its
  own boundary (the myscale class `a` rep at 0.03) maps to its own class. Wired
  on the joint nested-Laplace engine only (the per-observation `(lower, upper)`
  bounds are consumed by tulpa 0.0.44's built-in `interval_gaussian` family); the
  single-Laplace and NUTS paths error with a pointer, and the simplified-Laplace
  skew correction no-ops (`sla_status = "ordinal_unsupported"`). `sigma_pos`
  carries the latent log-cover SD, integrated on the same outer phi grid as the
  lognormal sigma. WAIC / PIT / PPC are the genuine discrete-class diagnostics.
  Recovery-tested from censored class data (`test-cover-ordinal.R`); shared
  positive-arm family / phi-grid logic extracted to `.cover_pos_family_grid()`.
  Requires tulpa (>= 0.0.44).

* Response / site / visit input handling now has a single source of truth in
  `R/inputs.R`, replacing the near-identical site-count cross-check each family
  binder hand-rolled (`occu`, `abun`, `removal`, `distance`, `fp_occu`,
  `dyn_abun`, `dyn_occu`, `jsdm`, the `ms_*` community families, and the cover /
  `occu_cover` / multiscale families):
    * `.tobs_check_site_count()` is the one check that the site dimension of `y`
      matches `nrow(data)`. The historical messages are preserved -- a 2D
      response counts "rows", a 3D / per-source response "sites", and the cover
      hurdle's response vector "values". (The integrated families' deliberate
      partial-coverage site maps are unchanged: they are not a mismatch.)
    * `tobs()` accepts a `tobs_data` frame in place of the `(data, y, visits)`
      triple. The frame's `occ.covs` / `y` / `det.covs` are unpacked through the
      same pipeline raw arguments take (`det.covs` is already the named-matrix
      shape `visits` consumes), so a frame fit is byte-identical to the
      equivalent raw fit. Passing `y =` / `visits =` alongside a frame errors.
    * `tobs()` records canonical response totals on the fit (`fit$dims`:
      `n_sites`, `max_visits`, `n_sources`) and, under `control$verbose`, reports
      them once before dispatch.
    * `test-inputs-frame.R`.

## 0.0.53 (2026-06-18)

* The `occu_cover()` outer Pareto-k diagnostic re-solves are validated against
  the faster tulpa 0.0.43 path (gcol33/tulpa#118, now the dependency floor):
    * `test-occu-cover-pareto-k.R` pins that the fast default (Shamanskii reuse,
      loosened inner tol, near-neighbour batch order, per-cell nearest-grid-mode
      warm start) reports the SAME k-hat as the byte-for-byte exact diagnostic
      (`tulpa.kdiag.refresh = 1`, `tol = 0`, no reorder, no per-cell warm), to a
      few 1e-4, and that the k-hat agrees with `loo::psis` on the diagnostic's
      actual importance ratios.
    * Measured speedup of the diagnostic re-solves is 2.6-2.8x at 144-256 cells
      with the k-hat byte-stable; the diagnostic still stays OFF by default
      (`control$diagnose.k`), since it reports k-hat only and does not move the
      betas, SDs, or field.
* Profiling note corrected in `occu_cover_joint_coupled.R` / `CLAUDE.md`: the
  binding per-solve cost is the per-Newton-iteration Hessian/gradient SCATTER
  (the beta arm's per-observation digamma/trigamma fill, 73-83%), NOT the sparse
  Cholesky factorize (a flat ~0.5 ms, 8-12%, not super-linear up to ~1100 cells).

## 0.0.52 (2026-06-18)

* The joint nested-Laplace outer Pareto-k diagnostic is now reachable on the
  `tobs_fit` itself (gcol33/tulpaObs#104). It was only at
  `fit$joint_fit$pareto_k`, so a diagnostic script reading `fit$pareto_k`
  directly got `NULL`.
    * `occu_cover()`, `occu()` spatial, and `occu_multiscale_cover()` joint-coupled
      fits promote `pareto_k`, `pareto_k_is_ess`, `pareto_k_scope`, and
      `pareto_k_proposal_source` to the fit top level (shared
      `.tobs_promote_pareto_k()`), and `glance.tobs_fit()` surfaces `pareto_k`,
      `pareto_k_is_ess`, and `pareto_k_proposal_source`.
    * `pareto_k_proposal_source` (the tulpa 0.0.41 mode-Hessian proposal,
      gcol33/tulpa#116, now the dependency floor) reads `"mode_hessian"` -- the
      importance proposal is curvature-backed, so the k-hat stays trustworthy
      even when a sharp posterior collapses the integration grid to ~1 cell -- or
      `"grid_moment"`, the regime to watch, when it comes from the grid-weighted
      node covariance as the grid concentrates.
    * Inert when the diagnostic was not requested: `control$diagnose.k` defaults
      OFF (gcol33/tulpaObs#101), so with it off no field or `glance()` column is
      added.

## 0.0.51 (2026-06-17)

* `occu_cover()` (and every other) nested-Laplace joint fit gets a usable
  outer-grid progress signal on a detached / redirected run, via the tulpa
  0.0.40 reporter fix (gcol33/tulpa#115, now the dependency floor):
    * the console line advances cell by cell instead of freezing at the serial
      pilot (`1/N`) -- the master thread emits from inside the parallel region;
    * the ETA rests on the realised per-cell throughput of completed cells
      rather than the cheap warm pilot (which projected ~10x optimistic), and is
      shown as a lower bound (`ETA >=`) until a parallel cell has finished.
  The detached-run heartbeat file (`control$progress.file`, gcol33/tulpaObs#43)
  is unchanged and remains the robust signal where a console flush is buffered
  away.
* Crossed / nested / correlated and uncorrelated random slopes on the
  detection / positive-cover arms of `occu_cover()` are verified end to end
  against simulated truth (crossed-intercept, nested, and correlated-slope
  parameter recovery; predict() shrink-to-mean for unseen levels), closing
  gcol33/tulpaObs#103.

## 0.0.50 (2026-06-17)

* `occu_cover()` random **slopes** on the detection / positive-cover arms now
  fit on a meaningful Sigma grid regardless of the covariate's scale. Each slope
  covariate is standardized to unit SD before fitting (mirroring the
  fixed-effect design autoscaling), and the reported slope BLUPs and SDs are
  back-transformed to the covariate's natural units (correlation is scale-free).
  The correlated-slope free-Sigma (`miid`) grid now uses a principled compact
  default -- symmetric correlation nodes including 0 and reaching strong +/-, and
  log-spaced SD nodes -- so the marginal correlation is no longer forced into a
  lop-sided range; widen it with `control$re.logchol.grid.p` /
  `re.logchol.grid.pos`. `simulate_occu_cover()` gains `re_det` slope-covariate
  scale control (`slope_sd`) for non-unit-scale recovery tests.

## 0.0.49 (2026-06-17)

* `occu_cover()` supports random effects -- intercepts AND slopes -- on the
  **detection** and **positive-cover** arms on the shared-field nested-Laplace
  path, with the usual `lme4` bar or `re()` spelling
  (gcol33/tulpaObs#102, #103). The grouping is per visit (one code per
  site-visit), distinct from the per-site occupancy-arm random intercept (#56),
  so a many-level categorical visit covariate (EUNIS habitat, observer) is
  partially pooled.
    - **Intercepts**, including **crossed** (`(1 | habitat) + (1 | observer)`)
      and **nested** (`(1 | region/site)`): each term is one `iid` latent block,
      named `sigma_re_p` / `sigma_re_pos` (suffixed by the grouping variable when
      several share an arm).
    - **Slopes** (needs the tulpa engine >= 0.0.39, gcol33/tulpa#114): an
      uncorrelated slope (`(x || g)`, `(0 + x | g)`) is one per-row weighted `iid`
      block per coefficient; a correlated slope (`(1 + x | g)`) one multivariate
      free-Sigma `miid` block. A slope term reports an `[n_groups x n_coefs]`
      BLUP matrix, a per-coefficient `sigma`, and (correlated) a `cor` matrix, all
      marginalized over the grid.
  BLUPs are in `fit$re` (a per-term list keyed by arm or `"<arm>:<var>"`) and via
  `ranef()`. `predict()` gains `type = "detection"` and sums each term's offset
  (weighting a slope by its covariate column in `newdata`), shrinking unseen /
  held-out levels to the population mean. Several crossed terms or a correlated
  slope grow the outer grid, so set `control$integration = "ccd"`; the free-Sigma
  grid uses a compact default, widened with `control$re.logchol.grid.p` /
  `re.logchol.grid.pos`. `simulate_occu_cover()` gains a `re_det` argument for
  crossed / nested / slope detection random-effect truth. As with the
  occupancy-arm RE the grid-integrated variances carry the binary / small-cluster
  inner-Laplace attenuation (a lower bound); the BLUPs recover the per-group
  structure. Needs `tulpa >= 0.0.39`.

## 0.0.47 (2026-06-17)

* `occu_cover()` and `occu_multiscale_cover()` joint nested-Laplace fits
  (`method = "nested_laplace"`) now default the outer Pareto-k accuracy
  diagnostic OFF (`control$diagnose.k = FALSE`), matching `occu_joint_coupled()`.
  Profiling traced the joint-fit runtime to this diagnostic: it re-solves the
  inner Laplace `k.samples` (200) times on the full areal field, each a
  super-linear sparse factorization, and accounted for 84-90% of wall time
  across field sizes -- the binding limit on per-species fits at fine spatial
  resolution (gcol33/tulpaObs#101). The diagnostic only reports the k-hat value;
  the fitted coefficients, SEs, and spatial field are byte-identical with it on
  or off. Re-enable it with `control$diagnose.k = TRUE` (and size the importance
  batch with `control$k.samples`).

## 0.0.46 (2026-06-17)

* `tobs_data(type = "cover")` gains `cover.floor` (default `0`): a cover value at
  or below the floor is stored as `NA` rather than as a positive observation,
  because the cover hurdle's positive arm is positive-only and a cover of `0` is
  an absence handled by the occurrence arm. This stops a `0` padded across
  unsampled cells (instead of left `NA`) from entering the positive arm as a
  fabricated zero, which, spread over a grid, flattens the spatial field. The
  conversion is reported with a one-line message; `cover.floor = -Inf` keeps
  every value verbatim.

## 0.0.45 (2026-06-17)

* Formula-native cross-arm coupling for `occu_cover()`. The cover (positive) arm
  carries a scaled copy of the occurrence arm's spatial effect, selected
  structurally with a constructor and no required name:
  `copy(spatial(), alpha = grid(c(0.25, 0.5, 1.0, 1.5)))`. The selector is
  type-carrying and reorder-stable: `copy(spatial(cell_idx))` disambiguates by
  grouping variable when several spatial effects are present, and an integer
  position is rejected. A per-component amplitude is
  `copy(spatial(), terms = list(intercept = grid(g0), time.sc = grid(g1)))`,
  keyed by the field's own block names and required to address every block. The
  coupling coefficient `alpha` (= sigma_pos / sigma_occ, the INLA `copy=`
  analogue) is marginalised over the grid; a scalar fixes it. Coupling is
  formula-native and explicit: a block with no `copy()` is decoupled (the field
  rides occupancy only), and decoupling is structural rather than a magic
  `alpha` of 0. `engine = "joint"` replaces `engine = "joint_coupled"`, and
  `occurrence =` / `positive =` read symmetrically with `detection =`
  (`formula =` stays as a deprecated alias). The fit is byte-identical to the old
  control-driven path (max abs difference 0).

* Full-model field-folded information criteria for `occu_cover()`. `tobs_waic()`,
  `tobs_dic()` and `tobs_cpo()` now fold the spatial intercept and trend fields and
  the per-visit detection process into the pointwise log-likelihood, so WAIC / DIC /
  CPO / LOO are numerically comparable to the INLA and spOccupancy criteria.
  `tobs_cpo()` additionally returns a LOO-PIT (`$pit`) for calibration checking and
  a per-observation LOO Pareto-k (`$failure`). Cross-checked against a brute-force
  pointwise evaluation (max abs difference < 1e-8).

* build: the `tulpa` dependency floor is raised to `tulpa (>= 0.0.38)` and the
  `Remotes` install reference to `gcol33/tulpa@v0.0.38`.

## 0.0.44 (2026-06-16)

* build: the `tulpa` dependency floor is raised to `tulpa (>= 0.0.37)` and the
  `Remotes` install reference to `gcol33/tulpa@v0.0.37`, version-matching tulpaObs
  to the current tulpa release.

## 0.0.43 (2026-06-16)

* build: the `tulpaMesh` dependency floor is raised to `tulpaMesh (>= 0.1.3)` and
  the stale `Remotes` install reference is corrected to `gcol33/tulpaMesh@v0.1.3`
  (the current tulpaMesh release).

## 0.0.42 (2026-06-16)

* build: the `tulpa` dependency floor is raised to `tulpa (>= 0.0.36)`, locking
  tulpaObs to the matching tulpa release (the two are ABI-coupled via
  `LinkingTo: tulpa`). The `Remotes` install reference stays at
  `gcol33/tulpa@v0.0.36`.

## 0.0.41 (2026-06-16)

* build: the `Remotes` install reference for `tulpa` is updated to
  `gcol33/tulpa@v0.0.36` (the current tulpa release). The `Imports` floor stays at
  `tulpa (>= 0.0.34)`, the minimum API tulpaObs requires.

## 0.0.40 (2026-06-16)

* `predict()` on an `occu_cover()` joint fit now propagates the visit-level
  positive-arm covariates from `newdata` into the conditional cover, instead of
  holding them at the reference (gcol33/tulpaObs#95). The positive arm splits a
  site design (intercept) from its visit-level covariate design at fit time; the
  predict handler rebuilt only the site design and zero-padded the visit columns,
  so a positive covariate supplied in `newdata` (e.g. the time axis of a
  `type = "change"` map) never entered the cover linear predictor and
  `delta_cover_cond` came out flat at zero. The handler now rebuilds the
  visit-level positive design from `newdata` with the same builder and column
  order as the fit. The model retains its visit-level formulas
  (`formulas$pos_visit`, `formulas$det_visit`) to support this. The occupancy /
  occurrence arm and the `cover()` hurdle were unaffected.

## 0.0.38 (2026-06-16)

* `tobs_data()` now preserves factor / character visit-level detection
  covariates as categorical. A column named in `det.covs` that is a factor or
  character is reshaped into a tagged character site x visit matrix (carrying its
  level set) instead of being coerced to numeric, and the detection / positive
  visit design expands it to k - 1 dummies for a k-level factor, with the first
  level (sorted-unique first value for a character column) the reference. Numeric
  `det.covs` follow the existing double-matrix path unchanged. The visit design
  builder now adds k - 1 contrast coding for any visit-level factor (previously a
  no-intercept visit formula expanded a factor to full k-dummy coding, collinear
  with the site-level intercept).

## 0.0.37 (2026-06-15)

* `occu_cover()` spatial NUTS (`method = "nuts"` with a `car_proper()` field on
  the occupancy arm, gcol33/tulpaObs#74) errored against `tulpa (>= 0.0.34)`. The
  fixed-hyper warm start passed the single-block `copy=` argument to
  `tulpa_nested_laplace_joint()`, which the single-block joint path no longer
  accepts; the copy coefficient is now declared on the cover arm as
  `field_coef = list(name = "alpha", grid = ...)`, the same convention the
  nested-Laplace `joint_coupled` path uses. Recovery, 95% interval coverage, and
  the calibration of the sampled coefficient SDs to the nested-Laplace SEs are
  restored.
* The "not yet supported" errors now name the supported route: a temporal term on
  `abun()` points to `dyn_abun()`; random effects with visit-level detection
  covariates on `abun()` name the two ways to proceed; and the multi-block
  nested-Laplace random-effect path lists the supported models (`iid`, `ar1`,
  `rw1`, `rw2`).
* `simulate_occu()` documents its actual return value (`y`, `data`, `truth`); the
  `coords` element it never produced is removed from the help.

## 0.0.36 (2026-06-15)

* Five cover / community families graduate from `status = "experimental"` to
  `"working"` after parameter-recovery and CI-coverage validation against
  simulated truth across seeds (gcol33/tulpaObs#96-100):
  * `occu_cover()` (#96): 95% Wald CI coverage holds near nominal on every path
    -- non-spatial `laplace` / `nuts` and the shared-field `nested_laplace`
    (`joint_coupled`) engine -- for both the beta and lognormal positive arms
    (measured pooled coverage 0.92-0.96). The recovery gates move from the 0.80
    experimental floor to the 0.85 working floor (pooled), and the beta arm and
    the shared-field paths gain explicit coverage gates.
  * `occu_multiscale_cover()` (#97): the four-arm recovery (cell psi, plot theta,
    visit p, cover) is validated on the `nested_laplace` and non-spatial
    `laplace` engines (pooled coverage ~0.95), and the coupled SVC-trend field
    recovers its shape. The availability / detection identifiability reduction is
    now surfaced: a fit with no within-plot replication (single releves) emits a
    message that theta and p collapse to the identified product theta * p, and
    the reduction is tested (psi and theta * p recover, the levels separately do
    not).
  * `ms_occu_cover()` (#98): community-mean 95% CIs are gated at the 0.85 working
    floor (measured ~0.92). The community-variance AGHQ debias cap
    (`re.aghq.maxdim`, default 4) is documented as a hard scope limit -- the
    tensor AGHQ is exponential in the total RE dimension -- with the EM variance
    above the cap explicitly tested as a lower bound. The reduced-rank
    spatial-factor (JSDM) path's loading / association recovery and per-species
    map calibration remain validated.
  * `ms_dyn_occu()` (#99): community-mean coverage gated at 0.85 (measured ~0.98;
    the shared colonisation / extinction dynamics cover at ~1.0) and the
    per-species first-season occupancy / detection variance components recover
    the realised spread.
  * `ms_int_occu()` (#100): the shared occupancy mean and the per-source
    detection components recover across seeds and more than one source; the
    community-mean coverage gate is tightened (measured ~0.89). For both
    `ms_dyn_occu()` and `ms_int_occu()`, a NUTS sampler and an areal-field path
    are a deliberate follow-up, not part of the working surface.

## 0.0.35 (2026-06-15)

* feat(tobs): `tobs()` gains a `by = "<species_col>"` argument for per-species
  batched fitting. Given a long / plot-level `data` frame and a species column,
  `tobs()` splits by species, builds each species' response onto one shared
  site x visit grid (via `tobs_data()`), and routes the per-species responses
  through the batched-independent driver, returning a `tobs_batch`. Scoped to
  `occu_cover()` and `cover()`. The split only reorganises the input: each
  species' fit is identical to the hand-built multi-response batch and to an
  independent single-species fit (equivalence tested to 1e-10).
* `tobs_batch_fit()` renamed to `tobs_get()`: it extracts one species' fit from a
  `tobs_batch`, and the old name read like a fitter.

* fix(cover): arm-specific spatial fields (single-arm `spatial(~ ... || node,
  to = "presence" / "positive")`, `method = "nested_laplace"`) are now projected
  at `predict()` time. The fields were estimated correctly (sigma > 0 on each
  arm) but every per-cell prediction came back flat, because the per-arm field
  block stored the fit-time per-observation node map and the accumulator skipped
  it whenever its length did not match the design (always true at predict, where
  the design has one row per cell). Arm-specific blocks now mirror the
  shared-field path: a block carries only its per-arm amplitude (membership) and,
  for a trend field, its covariate column name; the node map and per-cell weight
  are supplied by the consumer through `.tobs_joint_arm_eta` -- `predict()` passes
  the newdata cell map and column, the pointwise-log-likelihood consumer rebuilds
  the per-observation map and weight from `armspec_blocks`
  (`.tobs_armspec_obs_units` / `.tobs_armspec_obs_wfun`). The log-likelihood /
  WAIC / PPC consumer is numerically unchanged; the shared-copy and standalone
  paths are unaffected (gcol33/tulpaObs#95). Also fixes the spurious
  "this fit has 1 time-varying (trend) field(s); pass time_col" error raised by
  `predict()` on an intercept-only arm-specific fit: the positional "blocks 2.."
  trend convention no longer applies to arm-specific fits, whose blocks each
  carry their own weight column name.

## 0.0.34 (2026-06-12)

* perf(joint-coupled): the all-undetected occupancy-mixture detection (p, p)
  cross-Hessian in `occu()` / `occu_cover()` joint_coupled is no longer
  materialised as a dense V x V block per site (V = visits at that site). The
  block is analytically the rank-1 `a p p^T` (the cell density depends on the
  visits only through the scalar `P0 = prod_v (1 - p_v)`), so `nodet_mixture_block`
  / `occu_nodet_block` (`src/occu_coupling_shared.h`) now emit the per-site
  `(a, p)` into the engine's rank-1 self-cross descriptor (needs tulpa >= 0.0.34)
  and fold the rank-1 diagonal into the stored detection diagonal. Occupancy data
  is sparse, so nearly every site hits this branch; the per-iteration
  Gauss-Newton scatter -- the profiled inner-solve bottleneck at EVA scale
  (~99.9% of inner-solve wall time) -- drops from O(sum_s V_s^2) to O(sum_s V_s)
  (gcol33/tulpaObs#94; measured ~17x faster scatter at 32 visits/site on an
  occu_cover joint_coupled fit, the gap widening with visit count). The fit is
  numerically unchanged: a deterministic
  occu_cover joint_coupled fit (lognormal + beta arms) matches the former dense
  path to <= 1.4e-14 on coefficients, SDs, the spatial field, and logLik. Shared
  by the `occu_only`, `occu_cover`, and `occu_cover_latent` cell-coupling specs;
  `occu_multiscale_cover` keeps the dense path (its no-detection (p, p) block is a
  nested mixture, not a single rank-1).

## 0.0.33 (2026-06-12)

* perf(joint-coupled): the post-grid per-cell inner-covariance extraction in
  `occu()` / `occu_cover()` / `occu_multiscale_cover()` joint_coupled summaries
  no longer runs a serial-R `solve(Qk, E)` over the full betas + field latent
  per outer-grid cell. `.joint_inner_vcov_block()` now calls the engine's
  parallel selected-inversion primitive (`tulpa:::cpp_joint_inner_vcov_blocks`,
  needs tulpa >= 0.0.33): the betas block and betas x field cross are solved
  directly, the field marginal variances come from one Takahashi pass, and the
  cells run concurrently over `n.threads.outer` (gcol33/tulpaObs#93,
  gcol33/tulpa#112, #113). The SD summary is numerically unchanged; the betas
  block and field marginal variances match the former dense path to machine
  precision. `fit$joint_vcov` keeps its betas covariance and field marginal
  variances; its field x field off-diagonal now carries only the between-cell
  mode-dispersion term (the within-cell field cross-covariance, read by neither
  the summary nor the `Q_k`-direct `predict()` draws, is no longer formed).
* fix(cover): the single-field `cover(engine = "nested_laplace")` hurdle path is
  migrated off the `copy=` argument the joint single-block fitter dropped in the
  `(sigma, alpha)` reparam (it now hard-errors). The positive arm declares
  `field_coef = list(name = "alpha", grid = alpha_grid)` and the top-level `copy`
  is removed, mirroring `occu_cover()`'s single-block path; the multi-block /
  MCAR / trend branches (which still route through the copy-taking multi-block
  dispatch) are unchanged. Restores the single-field `cover()` nested-Laplace
  fits that errored at fit time (`test-cover-hurdle-nested-joint.R`).

## 0.0.32 (2026-06-12)

* fix(occu): the C++ NUTS occupancy likelihoods now read every observed visit
  when a missing visit precedes a valid one (#92). The single / dynamic /
  integrated occupancy kernels looped over `n_visits` (the count of valid visits)
  while indexing a `max_visits`-strided response that stores missing visits in
  place as a `-1` sentinel, so an interleaved or leading `NA` terminated the loop
  early and silently dropped the trailing valid visits, corrupting the
  likelihood, gradient, and posterior. The loops now stride the full
  `max_visits` dimension and skip the sentinels via the existing guard;
  `n_visits` is kept only for the "no surveys" early-out. Only `method = "nuts"`
  with non-trailing missing visits was affected (`method = "laplace"` always
  iterated the full dimension). `test-occu-interleaved-na.R` asserts that
  leading-vs-trailing `NA` encodings of identical data give identical fits under
  both methods (single and dynamic occupancy).
* refactor: the linear-predictor stability clamp is now a single
  `.tobs_clamp_eta()` helper over one `.TOBS_ETA_BOUND` constant, replacing ~20
  identical local `cl <- function(e) pmin(pmax(e, -30), 30)` definitions and the
  inline copies across the package; the C++ community-field kernel gains the
  matching `clamp_eta()` / `kEtaBound` (#89).
* refactor: the two-term no-detection log-likelihood in the occu_cover marginal
  is now the shared `.tobs_logsumexp2()` (a max-shifted `log1p` form), replacing
  the byte-identical block copied across the three occu_cover paths (#90).
* refactor: the visit-level design-matrix builder is now the shared
  `.tobs_build_visit_X()`. `occu()`, `abun()`, and `removal()` previously inlined
  their own builders and kept the visit `(Intercept)` column, which duplicated
  the site-level detection intercept; all four families now drop it through the
  one helper, so an intercept-bearing visit formula no longer makes the stacked
  detection design collinear (#91).

## 0.0.31 (2026-06-12)

* feat(cover): a `by = "factor"` argument on a cover spatial bar replicates the
  areal field across the factor's levels -- the graph becomes the block-diagonal
  `I_L (x) Q` (L disjoint copies) with each observation offset into its level's
  copy, sharing one precision across levels (one outer-grid axis). It composes with
  all three bar forms: shared (`||`), correlated (`|`, MCAR), and arm-specific, via
  `tulpa::tulpa_bar_field_replicate()`. Fixes the shared-`||` path, where the
  replication updated the adjacency but left the cached `n_spatial` at the base node
  count, so the replicated index was rejected as out of range; `n_spatial` now
  tracks the replicated graph, and the intentional L-component disconnection no
  longer emits the generic connectivity warning. `test-cover-spatial-bar-by.R`
  covers all three bar forms plus per-level field recovery.
* fix(dyn_occu): `logLik()` (and `AIC()` / `BIC()` / `glance()`) now return a
  finite value for `dyn_occu` fits (#87). The EM+Laplace packer left the
  log-likelihood unpopulated; `build_laplace_fit` now evaluates the exact
  HMM-forward marginal at the fixed-effect mode, and the dynamic exact-marginal
  refine moves the EM mode onto `colext`'s MLE. `tidy()` / `glance()` are
  re-exported so they resolve after `library(tulpaObs)`.
* fix(methods): a single unified convergence record across all families (#88).
  `cover()` stored its verdict at `fit$converged` while every other family uses
  `fit$convergence = list(converged, n_iter)`, so a mixed-family QC pass read `NA`
  for one location. `cover_fit` now carries the same record, and `convergence()` /
  `converged()` are exported as the documented accessors that normalise both
  layouts.

## 0.0.30 (2026-06-11)

* fix(dyn_occu): the dynamic-occupancy EM now uses an exact forward-backward
  (Baum-Welch) E-step. The colonization and extinction sufficient statistics are
  the smoothed pairwise joints `P(z_{t-1}=0, z_t=1 | y)` / `P(z_{t-1}=1, z_t=0 | y)`,
  and psi1 / detection use the smoothed marginals `P(z_t | y)`. Previously the
  E-step used only forward-**filtered** occupancy and a marginal-**product**
  approximation `(1 - w_{t-1}) w_t` for the transition events, which converged to a
  biased fixed point about 3.4 log-likelihood below `unmarked::colext` on the same
  data (inflated colonization / extinction). The fit now matches `colext`'s MLE to
  within the EM pseudo-count discretisation; a `colext` coefficient-equivalence
  gate is added to `test-dyn-int-occu-recovery.R` (#86). The shared single-species
  dynamic E-step also feeds the nested-Laplace and simplified-Laplace dynamic
  paths; the community dynamic family (`ms_dyn_occu`) has a separate E-step and is
  unchanged.

## 0.0.29 (2026-06-10)

* test(refimpl): `removal_laplace` and `distance_laplace` now gate head-to-head
  against `unmarked::multinomPois` and `unmarked::distsamp` (coefficients to
  ~5e-3, byte-identical log-likelihoods), extending the N-mixture-vs-`pcount`
  gold standard to the removal and distance families; plus a CI-runnable
  community-mean recovery smoke gate for `ms_abun` (gcol33/tulpaObs#83).
* test(recovery): multi-seed point recovery for `dyn_occu` (psi1, gamma, epsilon,
  p) with an independent R forward-recursion anchor for the dynamic-occupancy HMM
  marginal, and multi-seed recovery for single-source `int_occu` (gcol33/tulpaObs#84).
  NOTE: both families' deterministic Laplace standard errors are overconfident
  (the occupancy-intercept SE is ~an order of magnitude too small), so the gates
  assert point recovery; the SE-calibration gap is a separate kernel issue.
* docs: `DESCRIPTION` gains `URL` / `BugReports`; internal issue tokens stripped
  from rendered help; two non-ASCII characters removed from code comments; and
  the family front doors `abun()`, `cover()`, `distance()`, `ms_occu()`, and
  `occu_cover()` gain runnable `\donttest{}` fit-and-summary examples
  (gcol33/tulpaObs#85).

* feat(dyn_abun): a grouped random intercept on the **detection (`p`) arm** now
  fits on both engines (#82), alongside the existing initial-abundance (`lambda`)
  arm RE. Put the bar on the detection formula, e.g.
  `tobs(~ x, detection = ~ (1 | site), family = dyn_abun(), y = y)`. Unlike the
  initial-abundance arm -- where the predictor enters only the season-1 initial
  distribution, so the data-conditional weights are precomputed once and each
  quadrature node is an O(K) dot -- the detection predictor enters every season's
  observation pmf, so each AGHQ node re-evaluates the full exact HMM-forward
  marginal through a closed-form second-order `eta_p` forward-mode pass
  (`compute_dyn_abun_p_curv` / `cpp_dyn_abun_p_loglik`); NUTS adds a non-centered
  `p`-arm offset routed through the kernel's existing detection gradient. One
  grouping factor, on `lambda` OR `p` (a random effect on both arms in one fit is
  rejected; the AGHQ path integrates one arm at a time); survival / recruitment
  never carry random effects. Poisson and negative-binomial initial abundance.
  `ranef()` / `coef()` surface the detection RE as `sigma_p<t>_*` (AGHQ) and
  `log_sigma_p_*` (NUTS), with AGHQ debias on `sigma_p`.

## 0.0.28 (2026-06-10)

* fix(occu): the standalone `occu()` varying-coefficient (SVC) spatial bar
  `spatial(~ 1 + x || cell, graph = adj)` now fits through the joint direct-grid
  engine, single-arm (occupancy + detection, no cover arm), instead of the EM
  fixed-point nested-Laplace path (#81). The EM path oscillated and did not
  converge on real EVA-scale occupancy data -- a large-amplitude field, sparse
  detection, and a rich detection model drove the M-step to bounce and the fit
  was truncated at the iteration cap with `converged = FALSE`. The joint engine
  integrates the field hyperparameters on a direct outer grid with no
  fixed-point iteration, so it cannot oscillate; it is the same engine
  `occu_cover()` uses, with the cover arm removed. The occupancy mixture runs
  through a new `occu_only` cell-coupling spec that reuses the occupancy /
  detection derivatives of the `occu_cover` specs (single source of truth). The
  reroute is scoped to the SVC case: a plain single intercept field, a correlated
  MCAR bar, and temporal / random-effect structure stay on the EM path. `occu()`
  still reports `sigma` / `sigma_trend` (the field SDs marginalized over the
  grid), `spatial_field` / `trend_field`, and `predict(type = "occupancy" |
  "change")` now reads the shared areal field at each cell. Requires
  tulpa >= 0.0.30.

## 0.0.27 (2026-06-10)

* feat(formula): the varying-coefficient spatial bar `spatial(~ 1 + x || cell,
  graph = adj)` (a cell-indexed spatial intercept field plus a spatial trend
  field weighted by a per-site covariate) now fits on a standalone `occu()`
  nested-Laplace model, not only on `occu_cover()` (#67). The two-term spelling
  `icar(graph, group_var = "cell") + icar(graph, weight = x, group_var = "cell")`
  fits the same structure. This is the occupancy-only analogue of the
  `occu_cover()` coupled trend, with no cover arm and no coupling `alpha` -- the
  apples-to-apples match for an occupancy-only spatially-varying-coefficient
  model. Each field becomes one areal latent block on the existing multi-block
  nested-Laplace path: the intercept field is a plain icar block, and the trend
  field rides the same graph with a per-site `svc_weight` so its contribution is
  `weight_i * z[cell_i]`. The areal blocks are cell-indexed (one field node per
  graph cell, many sites per cell via `group_var`), so the field count can be far
  smaller than the site count. `occu()` reports `sigma` for the intercept field
  and `sigma_trend` for the trend field (the field SD marginalized over the outer
  grid), alongside `spatial_field` and `trend_field`. Requires tulpa >= 0.0.30
  (the single-arm driver gained the per-observation `svc_weight`). The correlated
  `|` bar (a free-Sigma MCAR field) stays on `occu_cover()` and errors with a
  pointer on `occu()`. Recovery-tested in `test-occu-spatial-svc-recovery.R`.
* The guard that previously rejected a weighted areal term on `occu()` now points
  to both the `occu_cover()` joint path and the standalone `occu()` nested-Laplace
  path; the term is rejected only on the engines that cannot carry it (the NUTS
  sampler and the single-Laplace path).

## 0.0.26 (2026-06-09)

* feat(formula): correlated (`|`) free-Sigma MCAR spatial coefficient fields on
  `occu_cover()` (#63). `spatial(~ 1 + x | cell, graph = adj)` on the occupancy
  formula declares the intercept and x-slope Besag fields as a separable MCAR with
  a free 2x2 cross-covariance `Sigma` (the within-arm covariance among the fields,
  integrated over the outer mode-centred CCD in log-Cholesky coordinates), then
  copies the whole correlated field onto the cover arm with one estimated
  amplitude `alpha`. Reports `sigma_mcar` / `rho_mcar` / `alpha_mcar`, marginalized
  over the grid. The independent `||` spelling (#61) is unchanged. Scoped to
  `icar`, the standard (non-latent) cover path, and at least two coefficient
  fields; correlated `|` does not compose with a per-group occupancy RE or the
  v2/v3 escape engines, which error with a pointer. This required a tulpa engine
  change (the coupled per-cell scatter now handles `INDEXED_MULTI` blocks);
  requires tulpa >= 0.0.29. Recovery-tested in `test-occu-cover-spatial-mcar.R`.
* A correlated `|` bar in an `occu_cover()` / `occu_multiscale_cover()` formula
  previously fell through to the independent `||` expansion **silently** -- the
  free cross-covariance was dropped with no warning. It now routes to the MCAR
  field (`occu_cover()`) or errors with a pointer (`occu_multiscale_cover()`),
  never a silent wrong model.

## 0.0.25 (2026-06-09)

* feat(formula): single-term varying-coefficient spatial bar in `cover()` /
  `occu_cover()` (#61). `spatial(~ 1 + time || cell, graph = adj, to =
  c("presence", "positive"))` is a compact spelling of the existing two
  weighted-areal-term coupled trend: the intercept column is the unweighted
  shared field, each covariate column a weight-scaled coefficient field, both on
  the bar node index, presence-anchored and copied to the positive arm with an
  estimated coupling `alpha`. It desugars to exactly the two-term form (`~ time +
  icar(graph, group_var = "cell") + icar(graph, weight = time, group_var =
  "cell")`), so the two spellings give the same fit; the bare `spatial()` /
  `weight =` forms are unchanged. `to =` validates against the family arm set,
  is order-free (presence anchor regardless of order), and defaults to both
  arms. Built on `tulpa::tulpa_is_spatial_bar()` / `tulpa::tulpa_bar_field_specs()`
  (tulpa#93). The correlated `|` bar (#64) and the arm-specific single-arm `||`
  bar (#65) have their own entries below.
* `cover()` / `occu_cover()` now emit an informative message when a bare `|` / `||`
  formula bar groups by the same factor an areal term uses as its graph-node
  `group_var` (#62): the bar is a random effect, not a spatial field, so the
  message points to `spatial(~ ... || cell, graph = adj)` (or the two-term
  `icar(graph, group_var) + icar(graph, weight, group_var)` form) for a spatial
  field. RE bars still fit as random effects; the message is suppressible and
  silent when the bar's factor is unrelated to any spatial term.
* feat(api): single-vector-response families accept the response on the top
  formula left-hand side, dropping `y =` (#66). `cover()` is the first such
  family: `tobs(cover.flat ~ predictors, data = dat, family = cover())` is
  equivalent to the one-sided `~ predictors` form with `y = dat$cover.flat`.
  A new `response` property on `obs_family()` (`"vector"` vs the default
  `"matrix"`) declares eligibility; matrix / array / list response families
  (`occu()`, `abun()`, the `ms_*` families, ...) keep `y =` and a two-sided
  formula for those errors. Supplying the response on both the LHS and `y =`
  errors.
* feat(formula): a correlated spatial bar (single `|`) on the cover hurdle wires
  the within-arm coefficient fields as a separable MCAR with a free 2x2 `Sigma`
  (#64). `spatial(~ 1 + time | cell, graph = adj)` makes the intercept and slope
  Besag fields correlated (free cross-covariance, integrated over the outer CCD
  grid in log-Cholesky coords), then copies the whole correlated field onto the
  positive arm with one estimated amplitude `alpha`. The fit reports `sigma_mcar`
  (per-field SDs), `rho_mcar` (cross-correlation), and `alpha_mcar`. `||`
  (independent fields) is unchanged. `nested_laplace` engine, intrinsic-CAR
  (icar) only; the simplified-Laplace correction over a correlated MCAR field is
  not yet wired (records `sla_status`). Requires tulpa >= 0.0.28.
* feat(formula): an INDEPENDENT (`||`) spatial bar with a single-arm `to =` fits
  an arm-specific separate latent field -- on that arm only, with its own
  precision and NO cross-arm copy (#65). `spatial(~ 1 + w || cell, graph = adj,
  to = "positive")` (or `"presence"`) places a non-copied areal field on the named
  arm; two separate single-arm calls give independent per-arm fields with no
  coupling between them. This is the free counterpart to the shared, presence-
  anchored, copied `to = c("presence", "positive")` field (#61): there one field
  is copied across arms with an estimated `alpha`, here each arm carries its own
  field, with that field's precision integrated on the outer nested-Laplace grid.
  No engine change -- a per-arm `spatial_idx = 0` sentinel makes the other arm's
  rows skip the block (the engine's `l_b > 0` scatter guard). `nested_laplace`
  engine; intrinsic-CAR / proper-CAR (`icar` / `car` / `car_proper`) only (the
  bym2 phi+theta mix is deferred). The fit reports `sigma_armspecific` (per-field
  SDs). Arm-specific fields do not compose with a shared field, a correlated `|`
  bar, a weighted trend, or `temporal()` / `re()` in the same formula (those
  couple the arms), and at most one field targets each arm. A single-arm
  correlated `|` bar stays copy-only and errors.
* The cover hurdle's two arms are now labelled `presence` (the `y > 0` Bernoulli
  arm) and `positive` (the `y | y > 0` arm) consistently across `summary()`,
  `print()`, and the `to =` argument (formula label == output label),
  replacing the earlier `occurrence` / `cover` headings.

## 0.0.24 (2026-06-08)

* feat(spatial): opt-in mode-centred central-composite design (CCD) for the outer
  field-hyperparameter integration of the in-package spatial / community fitters
  (#60). `control$integration = "ccd"` mode-finds the field hyperparameters
  (`tau`, `rho`, `sigma`, `range`) and places a CCD at the marginal-likelihood
  mode, scaled by the outer posterior covariance, reusing the engine's exported
  CCD primitives (`tulpa::ccd_grid()` / `ccd_to_theta()` / `ccd_weights()`) and
  surfacing the outer PSIS Pareto-k (`fit$spatial_pareto_k`). It declines to the
  fixed tensor grid when the outer curvature is ill-conditioned (a weakly-
  identified axis) and for a single positive hyperparameter, where the 1D grid is
  already cheap. The default `control$integration = "grid"` keeps the fixed tensor
  grid: each outer node is a full inner Laplace/EM solve, so the mode-find adds
  cost without a node-count saving on these already-coarse grids, and the CCD is
  most useful when a multi-axis hyperparameter posterior is well identified.
  Wired across the areal-BFGS families (`distance()`, `dyn_abun()`, `fp_occu()`),
  the community N-mixture Newton areal path, and the SPDE community path.
* Require `tulpa (>= 0.0.25)` and update the Remotes pin (the CCD primitives and
  the mode-centred-CCD machinery ship in tulpa 0.0.25).

## 0.0.23 (2026-06-08)

* **Breaking:** `control$trend` is removed. A spatially-varying trend is model
  structure, so it is now declared in the formula as a second, weighted areal
  term on the same graph as the intercept field --
  `icar(graph, weight = time, group_var)` (equivalently
  `spatial(graph, model = "icar", weight = time, group_var)`). The weighted
  term's contribution to each arm's predictor is `weight_i * z[cell_i]`, coupled
  onto the cover arm with its own scale (`fit$alpha_trend` / `fit$sigma_trend`)
  integrated over the outer grid (`control$alpha.grid.trend`, defaulting to
  `control$alpha.grid`). Requires `method = "nested_laplace"`; a leftover
  `control$trend` now errors with a migration pointer (#59).
* Require `tulpa (>= 0.0.18)` and update the Remotes pin. The committed s2z
  log-determinant guard test relies on the engine fix shipped in tulpa 0.0.18.
* feat(spatial): areal fields (ICAR, proper-CAR, BYM2) are wired on the
  abundance / occupancy arm across `removal()`, `distance()`, `dyn_abun()`,
  `fp_occu()`, `abun()`, and `occu_cover()` (#51), with the spatial fitters
  unified onto a shared field-spec areal driver.
* feat(re): Laplace AGHQ grouped random effects on one arm across `removal()`,
  `distance()`, `dyn_abun()`, and `fp_occu()` (occupancy and detection arms)
  (#51), plus NUTS sampling of a single intercept random effect on each arm
  including `abun()` (#51).
* feat(occu_cover): a shared spatial field with a per-group random intercept on
  the occupancy arm (#56).
* feat(occu_multiscale_cover / ms_occu_cover): `fitted()` and `predict()`
  (#53 part 1), a non-spatial Laplace path (#53 part 2), spatially-varying trend
  fields (#53 part 3), and AGHQ debias of the community variance components
  (#56).
* feat(ms_int_occu): partial / overlapping per-source site maps (#57).
* feat(dyn_abun): negative-binomial initial abundance (#52).
* feat(cover): fixed-effect priors thread through the nested-Laplace cover fit
  (#54).
* refactor(nuts): extract a shared single-arm-vector NUTS target oracle.

## 0.0.22 (2026-06-07)

* build: `tulpaMesh` moves from Suggests to Imports (the `spde()` term depends on
  it for mesh construction). `tulpaMesh::fem_matrices()` is re-exported so the
  mesh-assembly entry point is reachable directly from tulpaObs.
* build: drop the precompiled-header mechanism; each translation unit parses
  RcppEigen directly.
* docs(vignettes): correctness pass against the current API, and new
  documentation of the cover-hurdle row reductions (`control$aggregate.occ` and
  `control$aggregate.pos`, both default ON and byte-identical to the full
  per-plot fit). Corrected the spatial-occupancy WAIC interpretation (the areal
  field is not folded into the WAIC score), the `ms_occu` `ranef()` description
  (it returns the per-species deviations), and the `abun()` backend list (laplace,
  nested_laplace, nuts).

## 0.0.21 (2026-06-07)

* feat(cover): `control$aggregate.pos` now defaults ON for the beta positive
  arm (tulpaObs#49). The grouped sufficient-statistic collapse is byte-identical
  to the full per-plot beta arm on the single-block, coupled-trend and
  multi-block paths (`test-cover-hurdle-aggregate-pos.R`), with a multi-seed
  parameter-recovery suite behind the both-arms-aggregated default
  (`test-cover-hurdle-aggregate-recovery.R`); set `control$aggregate.pos = FALSE`
  for the full per-plot arm. An explicit `aggregate.pos = TRUE` still errors on a
  non-beta positive arm; the default leaves a non-beta arm untouched.

## 0.0.20 (2026-06-07)

* Require `tulpa (>= 0.0.16)` and update the Remotes pin so a fresh install
  resolves the grouped beta sufficient-statistic engine that `aggregate.pos`
  (0.0.19) depends on; the prior `>= 0.0.13` pin pre-dated that interface.

## 0.0.19 (2026-06-07)

* fix(occu_cover): the joint-coupled nested-Laplace parameter-surface
  covariance now carries the exact beta-hyperparameter cross-covariance and the
  full hyper-hyper covariance via the law of total covariance (the
  hyperparameters are the grid coordinates, so the within-cell term is zero and
  the cross-covariance is purely between-grid). Previously the hyperparameter
  block was diagonal, under-propagating the covariance of any derived quantity
  mixing a regression coefficient with `sigma` / `alpha`. Predicted occupancy
  and cover were unaffected (functions of the betas only). (tulpaObs#46)

* feat(cover): `control$aggregate.occ` (exact Binomial sufficient-statistic
  reduction of the cover-hurdle occurrence arm) now defaults to `TRUE`, backed
  by a multi-seed parameter-recovery suite on simulated beta-trend data
  (`test-cover-hurdle-aggregate-recovery.R`): the aggregated fit recovers truth
  with nominal 95% CI coverage on both arms' coefficients and the beta
  precision, and is byte-identical to the full per-plot fit. Set
  `aggregate.occ = FALSE` for the full occurrence arm. (tulpaObs#48)

* feat(cover): `control$aggregate.pos` (opt-in, default `FALSE`) adds the exact
  grouped-beta sufficient-statistic reduction of the positive (cover) arm. Plots
  sharing the positive design row and every per-observation latent component are
  collapsed to one row carrying `(n, sum log y, sum log(1 - y))`; tulpa's
  built-in beta likelihood reads those sufficient statistics, leaving the
  log-likelihood, gradient and Fisher Hessian pointwise unchanged. The fit is
  byte-identical to the full per-plot beta arm on the single-block and
  coupled-trend paths, alone and combined with `aggregate.occ`. Beta only (a
  lognormal positive arm errors with a pointer). (tulpaObs#49)

* docs(ms_occu_cover): the community joint fit now flags that its community
  VARIANCE components carry Laplace small-cluster attenuation (the community
  MEANS do not), via `print()`, a machine-readable
  `fit$ms_community$var_attenuation` marker, and `?ms_occu_cover`, so the
  reported between-species spread is not read as unbiased. (tulpaObs#47)

## 0.0.18

* test(ms-abun): correct two `ms_abun()` NUTS recovery assertions. The
  `fitted()$lambda` dimension check used the wrong site count (`40` instead of
  the fixture's `30`), and the per-species detection-coefficient recovery bar
  (`0.80`) was tighter than the arm actually recovers -- the realized
  correlation is `0.72` (abundance coefficients recover at `0.97`, with zero
  divergences), so the bar is now `0.65`. No change to the sampler or the
  model.

## 0.0.17

* fix(check): clears the `R CMD check --as-cran` ERRORs and WARNINGs (now only
  the GitHub-ecosystem CRAN-incoming WARNING remains). The namespace now
  imports every `stats`/`utils`/`methods` generic and function it uses
  (`nobs`, `simulate`, `update`, `model.matrix`, `glm`, `optim`, ...), which
  fixes the namespace-load failure that blocked the whole test suite; the two
  vignettes that failed to build are fixed (stale `summary()` column names;
  `predict()` returns a `tobs_prediction`, so the point estimate is `$mean`);
  non-ASCII characters and lost-brace Rd math are removed; the `tobs_test_*`
  goodness-of-fit helpers are documented; `model.matrix`-style globals and the
  `tulpaObs:::` self-reference are cleaned up.
* fix(build): drop the debug flags `-D_GLIBCXX_ASSERTIONS -g` from
  `Makevars.win` -- they triggered a spurious GCC 14 `-Warray-bounds` warning in
  `std::string` and bloated the shared object; the precompiled header now also
  rebuilds when `Makevars.win` changes.
* The progress reporter is on by default, so tests asserting silence pass
  `control$progress = FALSE`. Requires tulpa (>= 0.0.13).

## 0.0.16

* **feat(ms-abun): non-centered parameterization for the multi-species
  N-mixture NUTS.** The per-species block now holds standard-normal `z_s` and
  reconstructs the deviation per arm as `b_{s,arm} = C_arm z_{s,arm}` (`C_arm`
  the log-Cholesky factor of `Sigma_arm`). The community covariance leaves the
  `b`-prior (`z ~ N(0, I)`) and enters only the data term through `b = C z`,
  breaking the centered `b`/`Sigma` funnel that saturated the NUTS treedepth.
* **feat(progress): ETA reporting across every fitting loop**
  (gcol33/tulpaObs#43). Wires tulpa 0.0.12's unified progress reporter into all
  fitters, both channels ON by default -- a console bar plus a heartbeat file
  (written whenever `control$progress.file` is set, the channel that survives a
  detached run). The `cover()` / `occu_cover()` outer-grid paths flip progress
  ON by default (no longer tied to `verbose`); set `control$progress = FALSE`
  to silence the console bar.
* Require tulpa (>= 0.0.12) / `gcol33/tulpa@v0.0.12` for the shared progress
  reporter.

## 0.0.15

* Require tulpa (>= 0.0.10) / `gcol33/tulpa@v0.0.10`, which carries the fix for
  the threaded outer-grid nested-Laplace data race behind gcol33/tulpaObs#42
  (the coupled cover-arm dispersion was read lock-free across `n.threads.outer`
  threads). The threaded beta `cover()` / `occu_cover()` EVA-scale fits are now
  reproducible and crash-free; see the 0.0.14 entry and
  `dev_notes/issue42_root_cause.md`.

## 0.0.14

* **Fix: data race in threaded outer-grid nested-Laplace beta fits**
  (gcol33/tulpaObs#42). `cover(positive = "beta")` and
  `occu_cover(positive = "beta")` fits with `n.threads.outer > 1` could
  intermittently crash (native memory corruption) or hang at MOTIVATE/EVA scale.
  Root cause: in tulpa's threaded sparse joint outer-grid driver, the coupled
  (cover) arm's per-cell dispersion -- the beta precision on the `phi.grid.pos`
  axis -- was read lock-free from the shared `arms` during the inner Newton solve
  while a concurrent grid cell's `prep_at_grid` rewrote it under the phi-sync
  critical; every non-coupled arm already read a thread-local snapshot, but the
  coupled arm did not. Fixed in tulpa (`nested_laplace_joint_multi.{h,cpp}`) by
  snapshotting the coupled arms' dispersion under that critical and reading the
  per-thread snapshot in the coupled scatter / log-lik. Verified: a 220-region
  BYM2 beta cover fit is now identical serial vs `n.threads.outer = 6` to ~1e-10
  and finishes cleanly. The fix is in the tulpa dependency (root cause is there,
  compiled into `tulpaObs.dll` via the header-only joint driver); see
  `dev_notes/issue42_root_cause.md`.
* **Open-population (Dail-Madsen) N-mixture family `dyn_abun()`**
  (gcol33/tulpaObs#37). Latent abundance evolves across primary seasons:
  `N_1 ~ Poisson(lambda)`; for `t >= 2`, `N_t = Binomial(N_{t-1}, omega) +
  Poisson(gamma)`; observed via `Binomial(N_t, p)` over secondary visits. Unlike
  the static `abun()`, the latent abundance sequence is not closed form -- it is
  summed out by an exact HMM forward recursion over the abundance states, with
  analytic gradients from forward-mode differentiation of the scaled forward
  algorithm. Direct maximum-likelihood / Laplace fit (analytic-gradient BFGS,
  observed-information covariance) and a NUTS path over the same marginal
  (`method = "nuts"`, WAIC / LOO from the draws). Four site-level arms: initial
  abundance `lambda` (`formula`), detection `p` (`detection`), apparent survival
  `omega` (`omega_formula`), recruitment `gamma` (`gamma_formula`). The response
  is a 3D array `[n_sites x max_visits x n_seasons]`. `simulate_dyn_abun()`, full
  S3. Recovery / 95% coverage / NUTS recovery + WAIC, plus a correctness anchor
  (the C++ forward log-lik against an independent R forward recursion, exact to
  1e-9), an FD gradient check, and a C++ <-> R oracle cross-check in
  `test-dyn_abun.R`. Poisson initial abundance + constant recruitment this round;
  negative binomial, season-varying dynamics, and spatial / RE not yet wired.
  Per-site math `src/dyn_abun_kernel.h`; NUTS via the shared `src/nuts_engine.h`.
* **Multistate false-positive occupancy family `fp_occu()`** (gcol33/tulpaObs#40).
  The Miller et al. (2011) confirmed-detection design: each visit yields a state
  `y in {0, 1, 2}` (no detection / ambiguous detection / certain detection), with
  certain detections (state 2) only possible at occupied sites, which makes the
  model robustly identifiable. Four site-level logit arms -- occupancy `psi`
  (`formula`), true detection `p11` (`detection`), false-positive `p10`
  (`fp_formula`), certain-classification `b` (`b_formula`). The latent occupancy
  marginalises in closed form (two states); the Laplace fit maximises the exact
  marginal with an analytic gradient (BFGS) and an observed-information covariance
  (the inverse of the negative finite-difference Jacobian of the analytic
  gradient at the mode), and a NUTS path (`method = "nuts"`) samples the same
  marginal (WAIC / LOO from the draws). `simulate_fp_occu()`, full S3. Recovery /
  95% coverage / false-positive-arm covariate / NUTS recovery + WAIC, a
  correctness anchor (the two-state marginal against a direct computation, and the
  certain-detection identity), an FD gradient check, and a C++ <-> R oracle
  cross-check in `test-fp_occu.R`. Per-site math `src/fp_occu_kernel.h`; NUTS via
  the shared `src/nuts_engine.h` driver.
* **Binned distance-sampling family `distance()`** (gcol33/tulpaObs#38). Latent
  `N ~ Poisson(lambda)` (or negative binomial) in a covered region, observed
  through a half-normal or hazard-rate detection function over distance bins. The
  per-bin detected counts are multinomial over `(bin 1, ..., bin B, undetected)`
  with cell probabilities `pi_b = integral_bin g(x; sigma) f(x) dx` (line- or
  point-transect distance density `f`), integrated by Gauss-Legendre quadrature;
  the latent `N` is summed out in closed form (truncation `K_max`). Direct
  Laplace fit (`method = "laplace"`, Poisson or negbin, half-normal or
  hazard-rate with an estimated scalar shape) and a NUTS path over the same
  marginal (`method = "nuts"`, WAIC / LOO from the draws). The abundance formula
  is `tobs()`'s `formula`; the site-level `log sigma` model is `detection`; the
  response is an `n_sites x n_bins` count matrix; the bin edges and transect
  geometry travel with the family (`distance(cutpoints =, transect =)`).
  `simulate_distance()`, full S3. Recovery / 95% coverage / hazard-shape / NB
  dispersion / point-transect / NUTS recovery + WAIC, a closed-form correctness
  anchor (the Poisson distance marginal equals independent per-bin Poissons by
  thinning), and an analytic-observed-information vs finite-difference-Hessian
  check (the Louis curvature, including the second-derivative bin quadrature) in
  `test-distance.R`.
* Internal: the distance arm reuses the shared count-marginal core
  (`accumulate_count_moments` / `fill_nb_dispersion`, `src/nmix_kernel.h`) for
  the abundance / NB-dispersion math; the detection arm (site-level `log sigma`,
  optional scalar hazard shape, bin integrals + first/second eta-derivatives by
  quadrature) is `src/distance_quad.h` / `src/distance_kernel.h`. The tulpa NUTS
  engine plumbing is factored into a shared driver (`src/nuts_engine.h`) now used
  by both the count-marginal families and `distance()`. Byte-identical for the
  existing families (full `test-abun.R` / `test-abun-re.R` / `test-removal.R`
  suites unchanged).

## 0.0.13

* **Removal-sampling family `removal()`** (gcol33/tulpaObs#39). Sequential
  depletion: latent `N ~ Poisson(lambda)` (or negative binomial) observed
  through `K` ordered removal passes, where pass `k` removes
  `Binomial(N - sum_{l<k} y_l, p_k)` of the individuals still present. The
  depleting-binomial product is the multinomial-removal likelihood, and the
  latent `N` is summed out in closed form (truncation `K_max`), so the fit is a
  direct Laplace approximation (`method = "laplace"`, Poisson or negbin) with a
  NUTS path over the same marginal (`method = "nuts"`, WAIC / LOO from the
  draws). `simulate_removal()`, full S3 (`fitted`/`predict`/`simulate`/
  `residuals`/`coef`/`vcov`/`confint`/`logLik`). Recovery / 95% coverage / NB
  dispersion / NUTS recovery, plus a closed-form correctness anchor (the Poisson
  removal marginal equals independent Poissons) in `test-removal.R`.
* Internal: the N-mixture per-site moment / negative-binomial dispersion math is
  factored into shared helpers (`accumulate_count_moments`, `fill_nb_dispersion`
  in `src/nmix_kernel.h`), and the non-spatial count-marginal Laplace driver and
  NUTS machinery are now shared headers (`src/marginal_count_laplace.h`,
  `src/marginal_count_nuts.h`) instantiated by both `abun()` and `removal()` --
  one source of truth for the count-marginal fit. Byte-identical for `abun()`
  (full `test-abun.R` / `test-abun-re.R` recovery suites unchanged).

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
