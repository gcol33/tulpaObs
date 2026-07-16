# CLAUDE.md

Guidance for Claude Code in this repo. Caveman speak: terse, telegraphic. File
paths + function names exact. (Compacted; if a detail is missing, read the source
file named.)

## tulpaObs — hierarchical latent-state observation models on tulpa

Bayesian occupancy / abundance / distance / removal / cover. Built on
[`tulpa`](https://github.com/gcol33/tulpa) engine. R pkg, C++17 backend
(Rcpp/RcppEigen). Needs sibling `../tulpa` checkout (`LinkingTo: tulpa`).

**Public API:** `tobs()` + family ctors: `occu()`, `dyn_occu()`, `ms_occu()`,
`ms_dyn_occu()`, `ms_int_occu()`, `int_occu()`, `jsdm()`, `abun()`, `ms_abun()`,
`dyn_abun()`, `count()`, `ms_count()`, `distance()`, `ms_distance()`, `removal()`,
`fp_occu()`,
`cover()`, `occu_cover()`, `ms_occu_cover()`, `occu_multiscale_cover()`. S3 classes all `tobs_*`
(`tobs_fit/model/family/spatial/temporal/re/svc/latent/priors_spec`).

**Structured terms live in formula** (lme4/mgcv/INLA style), NOT `tobs()` args.
Registry `R/formula_terms.R` maps name -> ctor: spatial `icar() bym2() car()
car_proper() gp() multiscale_gp() spde()`; `re()`; `temporal()`; `svc()`;
`latent()`; `copy("id")` (share one realization state+detection). Ctors
`.tobs_term_*`, parser `.tobs_parse_formula`/`.tobs_bind_formulas` internal. No
exported `tobs_icar()`, no `spatial=`/`temporal=`/`re=` args. Term's process =
which formula it sits in; fitter derives `shared=c(occ,det)` via
`.tobs_structures_from_model()`.

**lme4 bars = sugar over `re()`.** `.tobs_desugar_bars()` (`R/formula_parse.R`)
rewrites bars to `re()` on AST before `terms()`. Forms (#10): `(1|g)`, `(x|g)`,
`(1+x+z|g)` (multi-slope, `cbind()` stack), `(x||g)` (uncorr), `(0+x|g)`
(slope-only), crossed/nested `(1|g:h)` `(1|g/h)`. Multi-slope + nested pure R-side
(`build_re_spec()` width->`n_coefs`); slope-only needs tulpa `re_has_intercept`
(ABI 22).

**RE both engines** (#11). NUTS fits all forms. Default `engine="laplace"` fits
iid intercept RE, uncorr slopes, AND corr slopes (`(1+x|g)`) on occ OR det arm of
single-season via variance-component EM (`R/em_laplace_re.R`,
`.tobs_em_laplace_re()`): splits RE into occ/det arm by `shared`
(`.tobs_re_split_arms()`), wraps tulpa `tulpa_laplace()` in occupancy missing-data
EM, EM/REML cov update `Sigma_k <- mean_g[b_g b_g' + Cov(b_g|y)]`. `Cov(b_g|y)`
from `tulpa_laplace(return_re_cov=TRUE)$cov_blocks`: occ arm scales by `M`
(pseudo-binomial inflation), det arm weighted binomial (`M=1`). Corr keeps full
`Sigma`; uncorr projected to diag each M-step. Gates (error, point to NUTS) via
`.validate_re_laplace()`: RE across BOTH arms, RE+spatial, RE+visit-level det, RE
on non-single families.

Raw EM variance components carry Laplace small-cluster bias for binary; FE SEs do
not. Default `re.aghq=TRUE` (`control=list(n.quad=)`) debiases via adaptive GH on
exact per-group marginal (single grouping factor, one arm, dim<=3; else EM
fallback). Engine `tulpa::tulpa_re_aghq()`; `R/re_aghq.R` (`.tobs_re_aghq()`)
wraps, supplies occ/det site marginal as `make_site`. Measured occ arm: n=8 sigma
bias ~18%(EM)->~4%(AGHQ); det arm ~70% attenuation->~1%, 88-96% coverage. Det-arm
RE params named for det process (`sigma_p<t>_*`, `re_p<t>_*`). Default LKJ
(`re.lkj=1.5`) per corr block regularizes weak correlation off +-1 (SDs untouched;
`re.lkj=1` disables). RE naming + BLUP recon both engines in `R/re_effects.R`
(`ranef()`/`coef()` in `R/methods.R`); corr off-diag `cor_<g>_<ci>_<cj>`.

Legacy internal: `.tobs_build_model()`, `.tobs_fit_model()`, `.tobs_laplace()`.
Family roster: `PLAN_tulpaObs.md`; dispatcher `R/obs_families.R` + `R/tobs.R`.

## Common commands

```r
Rcpp::compileAttributes()                  # after src/ or R/ doc edits
devtools::document()
devtools::install("../tulpa", quick = TRUE)# after tulpa ABI bump
devtools::load_all()                       # iterative dev
devtools::check(args = "--no-manual")      # full check (skip manual: dev non-ASCII)
devtools::test()                           # all tests
devtools::test(filter = "sla-cover-joint") # one file -> test-sla-cover-joint.R
testthat::test_file("tests/testthat/test-occu.R")
testthat::test_file("tests/testthat/test-occu.R", desc = "single fit recovers truth")
```

Convention: probes/repros/notes in `dev_notes/` (`_` prefix = runner, `probe_*`
= diagnostic, `repro_*` = minimal reproducer for upstream tulpa bugs).

## Testing: smoke-first (DEFAULT for iteration)

Full suite fits 15-20 models per seed x many seeds + NUTS/spatial recovery -> HOURS.
Do NOT run on every edit. Ladder:

1. **Iterating** -> only the test file(s) you touched:
   `testthat::test_file("tests/testthat/test-occu.R")` or
   `devtools::test(filter = "occu")`. Seconds.
2. **Whole-suite smoke** (plumbing/dispatch/closed-form, no fits) ->
   `Sys.setenv(TULPAOBS_FAST = "1"); devtools::test()`. **~22s**, ~1330 assertions;
   ~215 heavy fitting/recovery/NUTS blocks report as skips (never silently dropped).
   `skip_if_fast()` gates every fitting block; no-op when env var unset.
3. **Full recovery suite** (all seeds, NUTS, spatial) -> ONLY before committing to
   main, before release, or when asked. `Sys.unsetenv("TULPAOBS_FAST");
   devtools::test()`. Hours; uses `Config/testthat/parallel`.

Adding a slow test: pair `skip_if_fast()` + `skip_on_cran()` at top of any
multi-seed fit / NUTS block (helper `tests/testthat/helper-speed.R`). C++ recompiles
ccache-backed; only a killed/partial build needs `pkgbuild::clean_dll()`.

## Architecture

One C++ entry `cpp_occu_fit` for NUTS. Laplace via tulpa EM+Laplace
(`tulpa::tulpa_em_laplace`). Default `method="laplace"` = penalized EM, Gaussian
marginals, no post-EM correction; FE prior attached per M-step block as
`beta_prior`. MI/Gibbs opt-in: `method="laplace_gibbs"`/`"laplace_mi"`; same prior
threads into correction refits (tulpa#27), penalised unless `priors=FALSE`.

**nested-Laplace** (`em_nested_laplace.R`, `method="nested_laplace"`): non-conjugate
hyperpriors via multi-block latent prior. `.tobs_em_nested_laplace()` wraps
`.tobs_laplace(latent_prior=)`: multi-block prior from formula `spatial`/`temporal`/
`re` attached to state M-step block (`.tobs_laplace_nested()`), routed through
`tulpa_nested_laplace()`. Wired: single-season, integrated, community, dynamic occ.
Cover hurdle joint uses `tulpa_nested_laplace_joint()`.

**INLA-style NA-response prediction + calibrated CIs** (single-season): all-NA
detection sites (`.tobs_heldout_sites()`) field-interpolated; `predict(type="state")`
returns per-site psi posterior (`psi`, `psi_lower`, `psi_upper`, 95%) marginalised over
hyperparam grid. Calibrated by one exact-marginal refine
(`.tobs_occu_state_marginal_fit()`): integrate out z, each site Bernoulli on
D=1{>=1 detection}, mean q*sigma(eta), q=1-(1-p)^J. Fits tulpa `family="bernoulli"` w/
per-obs prob scale `det_prob=q` (scaled-Bernoulli): converged Hessian = marginal
curvature, `fitted_eta_var` = calibrated per-cell predictive variance. Per-row eta
posterior = Gaussian mixture over cells (`.nested_psi_mean()`,
`.nested_psi_quantiles()`). Coverage ~1.0, cor ~0.88, MAE ~0.11 on 10x10 icar/bym2
(`dev_notes/probe_nested_ci_coverage.R`). Old tulpa w/o `fitted_eta_var` -> NA intervals.

**Simplified-Laplace skew correction** (`simplified_laplace.R`, `sla_*` files):
orthogonal post-fit marginal refine, `approx="simplified_laplace"` (`*_sla` methods).
Computed for single/dynamic/integrated occ + cover hurdle; no-ops to Gaussian
(records `sla_status`) for jsdm.

**Backend coverage enforced centrally**: `.tobs_family_methods` in `R/tobs.R` =
single source of truth for which `method` each family supports; `tobs()` errors with
pointer, no silent downgrade. `nested_laplace` = occu/int_occu/dyn_occu + cover;
`*_sla` on nested = occu + cover only; cover hurdle has no
NUTS/`laplace_gibbs`/`laplace_mi`. `abun` = laplace + nuts (non-spatial) +
nested_laplace (areal). `ms_abun` = laplace + nested_laplace (shared areal field) +
nuts (non-spatial #14, and shared fixed-hyper proper-CAR field #73).
`occu_multiscale_cover` = nested_laplace ONLY. Community occupancy
`ms_occu`/`ms_dyn_occu`/`ms_int_occu` = laplace ONLY (shared community Laplace-EM,
`R/community_em.R`; per-species coef RE, per-arm community covariance).

Observation families (`removal`/`distance`/`fp_occu`/`dyn_abun`) = laplace + nuts,
non-spatial Pois/NB. Each: closed-form (or exact HMM-forward, dyn_abun) marginal over
latent N, analytic gradients, in-tree FullGradFn driving tulpa NUTS engine (shared
`src/nuts_engine.h`), draws -> WAIC/LOO. All filed observation-family issues shipped.

- **`removal`** (#39): sequential depletion, pass k sees `N - sum_{l<k} y_l` trials =
  depleting-binomial = multinomial removal; latent N to `K_max`. Shares count-marginal
  Laplace (`src/marginal_count_laplace.h`) + NUTS (`src/marginal_count_nuts.h`) w/ abun;
  per-site `src/removal_kernel.h` over shared
  `accumulate_count_moments`/`fill_nb_dispersion` (`src/nmix_kernel.h`). Areal
  icar/car_proper/bym2 on abundance arm via `nested_laplace` (#51):
  `compute_removal_site` returns same `NMixSiteResult` moments as Royle kernel -> reuses
  templated count-spatial driver (`src/nmix_count_spatial_driver.h`, byte-identical to
  nmix) via `cpp_nested_laplace_removal_*` (`src/removal_spatial.cpp`) + R packers
  (`R/removal_spatial.R`, `.tobs_fit_removal_spatial`); one unit/site, K_max=per-site
  total. spde/temporal + NUTS+spatial gated.
- **`distance`** (#38): latent N in covered region, per-bin counts multinomial over
  `(bin 1..B, undetected)`, `pi_b = int_bin g(x;sigma) f(x) dx` (half-normal/hazard
  key, line/point transect density `f`); bin integrals + 1st/2nd eta-derivs by
  Gauss-Legendre quad. Reuses count-marginal core (`src/nmix_kernel.h`); det arm =
  site-level `log sigma` + optional scalar hazard shape (`src/distance_quad.h`/
  `src/distance_kernel.h`), own Laplace driver (`src/distance_laplace.cpp`) + NUTS
  (`src/distance_nuts.cpp`). K_max default `3*max(rowSums)+100`. Laplace grouped RE on
  abundance arm (half-normal key, one grouping factor, dim<=3) via shared count AGHQ
  (`DistanceGroupedOracle` over `CountGroupedOracle`; theta = count layout
  `[beta_lambda|beta_sigma|log_r?]`); hazard-key + det-arm RE gated. Areal
  icar/car_proper/bym2 on abundance arm via `nested_laplace` (#51,
  `R/distance_spatial.R` over shared areal-BFGS driver `R/areal_bfgs.R`): distance
  marginal = bin-multinomial (NOT N-mixture, so not the count-spatial driver) but
  exposes analytic per-site gradient (`cpp_distance_site_sweep` over
  `compute_distance_site`), so BFGS + FD-Hessian recovers observed info. Half-normal
  only; Pois + NB; hazard-spatial/temporal + NUTS+spatial gated.
- **`fp_occu`** (#40): Miller 2011 multistate false-positive occupancy, `y in {0,1,2}`
  (none/ambiguous/certain); certain detections only at occupied sites identify it.
  Latent z summed (2-state); 4 logit arms psi/p11/p10/b (`fp_formula`/`b_formula`,
  default `~1`). Laplace = analytic-grad BFGS, vcov = inv of -FD-Jacobian of analytic
  grad (`src/fp_occu_kernel.h`, gradient only). NUTS `src/fp_occu_nuts.cpp`. Laplace
  grouped RE on psi OR p11 arm (one factor, dim<=3) via pure-R `make_site` AGHQ
  (`.tobs_fp_occu_re_aghq`, no native oracle, branches on arm): holding other arms
  fixed makes 2-state marginal `psi*A+(1-psi)*B` a fn of one scalar offset/site (psi
  shifts mixture weight; p11 shifts occupied emission `A`, uniform over a site's
  visits), both closed-form d1/d2. RE both arms at once rejected; p10/b never carry
  structured terms; NUTS samples psi-arm intercept RE only. Det-arm RE params named
  `sigma_p<t>_*` (p11 process). Areal icar/car_proper/bym2 on psi arm via
  `nested_laplace` (#51, `R/fp_occu_spatial.R` over shared areal-BFGS driver
  `R/areal_bfgs.R`, `.tobs_areal_bfgs_fit`, shared w/ dyn_abun): BFGS over two-state
  marginal (`cpp_fp_occu_total_log_lik` analytic grad) + CAR prior, FD-Hessian
  observed info; one unit/site. Occupancy fields more weakly identified than count
  (one binary site/node). temporal + NUTS+spatial gated.
- **`dyn_abun`** (#37): Dail-Madsen open N-mixture, `N_1~Pois/NB(lambda)`,
  `N_t=Binom(N_{t-1},omega_{t-1})+Pois(gamma_{t-1})`, `Binom(N_t,p)` obs. Latent N
  sequence summed by exact HMM forward over 0..K_max; analytic grad by forward-mode
  diff (`src/dyn_abun_kernel.h`). 4 arms lambda/p/omega/gamma (`omega_formula`/
  `gamma_formula`, default `~1`). NUTS `src/dyn_abun_nuts.cpp`. K_max default
  `max(count)+40` (forward ~cubic in K). Pois OR negbin init (`mixture="negbin"`,
  trailing `log_r`). **Season-varying omega/gamma (#80)**: the transition from
  season t-1 to t uses interval-(t-1) vital rates, so a covariate on
  `omega_formula`/`gamma_formula` carried as an `[n_sites x (T-1)]` matrix column of
  `data` drives the dynamics. Kernel takes interval-indexed `eta_omega`/`eta_gamma`
  (length T-1) with per-interval forward-mode gradients (direction iv born at its
  own transition, propagated through later ones); scalar overload broadcasts ->
  byte-identical constant-rate path. Binder `.tobs_dyn_abun_arm_design` unrolls to
  long-form `[(site x interval) x p]` only when a covariate is a `[n_sites x (T-1)]`
  matrix column; both laplace + NUTS scatter the per-interval score through the
  long-form design. `simulate_dyn_abun(beta_omega=, beta_gamma=)` season-varying
  truth. NUTS+temporal still gated. Laplace grouped RE on initial-abundance arm (one
  factor, dim<=3): RE shifts only eta_lambda -> per-site marginal factorises
  `L(eta_lambda)=sum_n1 pi_n1(eta_lambda) c(n1)`, conditional `c(n1)=P(all data|N_1=n1)`
  (O(K^2 T) HMM BACKWARD pass `compute_dyn_abun_init_weights`) precomputed ONCE per
  make_site (`cpp_dyn_abun_init_weights_mat`); each AGHQ node = O(K) dot
  (`cpp_dyn_abun_init_loglik`). Pure-R `make_site` AGHQ (`.tobs_dyn_abun_re_aghq`, no
  native oracle), Pois+NB; omega/gamma never structured. **Detection-arm RE (#82)**:
  `.tobs_dyn_abun_re_aghq` branches on arm; a p-arm RE shifts eta_p, which enters
  EVERY season's obs pmf, so `c(n1)` cannot be precomputed -- each AGHQ node
  re-evaluates the full O(K^2 T) forward marginal via a closed-form SECOND-ORDER
  eta_p forward-mode pass (`compute_dyn_abun_p_curv` -> `cpp_dyn_abun_p_loglik`:
  transition/initial are p-free, so `d obs`/`d2 obs` are the only source terms,
  propagated + renormalised per season; logL/d1/d2 FD-validated). NUTS p-arm RE
  (`src/dyn_abun_nuts.cpp` `re_arm` 0=lambda|1=p): non-centered offset to eta_p,
  grad via `grad_eta_p`; both-arm RE rejected (AGHQ integrates one arm). Det-arm RE
  params named `sigma_p<t>_*` / `log_sigma_p_*`. Shared pmf helpers
  (`da_obs_season_pmf`/`da_recruit_pmf`/`da_binom_pmf_row`) back the forward
  gradient kernel, the backward `c` pass, AND the eta_p second-order pass. Areal icar/car_proper/bym2 on initial-abundance
  arm via `nested_laplace` (#51, `R/dyn_abun_spatial.R` over shared areal-BFGS driver
  `R/areal_bfgs.R`): BFGS over exact forward-HMM marginal (`cpp_dyn_abun_total_log_lik`
  analytic grad) + CAR prior, FD-Hessian observed-info Laplace integrated over
  tau[,rho]; one unit/site, Pois+NB; temporal + NUTS+spatial gated.

### N-mixture abundance (`abun()`)

Royle 2004: `N_i ~ Pois(lambda_i)`, `y_ij|N_i ~ Binom(N_i, p_ij)`. NO EM — N
marginalises closed-form (sum to `K_max`), tulpa fits marginal directly, analytic
gradients + observed-Fisher. Wiring (`R/abun.R`): binder -> `model_type="nmix"` (site
`X_lambda`, long-form `y`/`site_idx`/`X_p`, drop NA visits), `.tobs_fit_nmix()` ->
in-tree fitter (`R/nmix_laplace*.R`, `src/nmix_*.{cpp,h}`), wrapped `build_nmix_fit()`.
Pseudo-draws MVN from JOINT (lambda,p) cov. Log link: `compute_intercepts()` reads
per-process `pi$link` (logit default, log for lambda). `simulate_abun()` +
`test-abun.R`.

- **Non-spatial** (`laplace`): `nmix_laplace()`. vcov = marginal observed-Fisher inv
  (full joint lambda/p block).
- **Areal** (`nested_laplace`): icar/bym2/car_proper on abundance ->
  `.tobs_fit_nmix_spatial()` -> `nmix_laplace_{icar,bym2,car_proper}` (one unit/site).
  Cov grid-integrated (law of total cov): kernels return per-grid `cov_blocks`, wrapper
  `V = sum_k w_k[cov_k + (m_k-mbar)(m_k-mbar)']` (`.nmix_grid_vcov()`). Rank-deficient
  intrinsic fields use sum-to-zero constraint penalty in `nmix_spatial_beta_cov()`/
  `nmix_beta_cov_bym2()`.
- **Negbin** (`abun(mixture="negbin")`): kernel P/NB threads both paths. NB adds `r`
  (`Var=lambda+lambda^2/r`). Non-spatial: `log_r` jointly estimated, trailing vcov
  coord (`nmix_dispersion`, delta `r_sd`). Spatial: `r` grid-integrated
  (`r_mean`/`r_sd` in `nmix_hyper$r`). Matches `unmarked::pcount(mixture="NB")`.
- **Grouped RE** (`abun()` + `(1|g)`/random slope, either arm; #13): site-level RE on
  abundance OR detection. Grouping = non-species; RE = subset of coefs on ONE arm.
  `.tobs_fit_nmix_re()` warm-starts no-RE betas, refines via `.tobs_nmix_re_aghq()`
  (thin `make_site` over `nmix_site_marginal()`) -> `tulpa::tulpa_re_aghq()`. Det arm:
  per-site offset uniform over visits. NB threads `log_r`. Native `NMixGroupedOracle`
  (`src/nmix_re_oracle.{h,cpp}`). Gates: RE+spatial, RE+visit-det, RE both arms.
  `test-abun-re.R`.
- **NUTS** (`method="nuts"`, #41; `R/abun_nuts.R`, `src/abun_nuts.cpp`): in-tree C++
  FullGradFn over closed-form marginal, byte-exact vs R oracle, warm-start + Laplace
  metric. NUTS + areal (#51, `.tobs_fit_abun_nuts_spatial`): FIXED-HYPER non-centered
  PROPER-CAR field on abundance — field precision tau Q(rho) fixed at nested-Laplace
  estimate, whitened raw ~ N(0,I), z = Linv %*% raw, sampled via OPTIONAL field block
  in shared count-NUTS header (`marginal_count_nuts.h`: n_field_units/field_map/
  field_Linv; non-spatial + removal NUTS byte-exact unchanged). 0 divergences, beta SDs
  calibrate to nested-Laplace SEs (tulpa#87 pattern). car_proper only — intrinsic
  icar's flat field-mean needs sum-to-zero reparam, so icar/bym2 NUTS+spatial gated.

### Community / multispecies N-mixture (`ms_abun()`)

spAbundance `msNMix`: per-species nmix w/ Gaussian community hyperpriors,
`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`, `beta_p_s ~ N(mu_p, Sigma_p)`.
`N_{s,i}` integrates closed-form per species-site (same kernel as `abun()`);
per-species deviations `b_s` = RE, community priors pin `(mu_lambda, mu_p)` fixed.
`y` = 3D array `[sites x max_visits x species]` or named list of count matrices;
`species=` required.

Wiring `R/ms_abun.R`: binder -> `model_type="ms_nmix"`, `.tobs_ms_nmix_longform()`
flattens 3D -> stacked `(y, site_idx, species_idx, X_p)`, `.tobs_fit_ms_nmix()` ->
in-tree C++ Laplace-EM `nmix_laplace_re()`. Builds native
`tulpaObs::NMixCommunityOracle` (subclass `tulpa::REGroupOracle`) via
`cpp_nmix_community_oracle()`, drives via in-tree `cpp_nmix_community_em()` (`n_quad=1`)
OR `tulpa::tulpa_re_aghq(oracle=)` for `n_quad>1` variance debias. Fit: per-species
coef mode-find (complete-data Fisher), closed-form EM cov M-step
`Sigma_k=mean_s[b_s b_s'+Cov(b_s|y)]`, FE SEs from marginal observed-info Schur
complement of b-block (Louis 1982). Per-site kernel `src/nmix_kernel.h` caches lgamma
(`compute_nmix_site_cached`); `compute_nmix_site()` Poisson delegates ->
single-species/spatial byte-identical.

`coef()` = community means; `ranef()` = per-species deviations; `vcov()`/`confint()` =
community-mean cov; `fitted()`/`simulate()` = per-species. `simulate_ms_abun()` +
`test-ms-abun.R` (community-mean recovery, 95% coverage 20 seeds, per-species coef,
S3).

- **Negbin** (`ms_abun(mixture="negbin")`): per-species dispersion RE
  `log_r_s ~ N(mu_log_r, sigma_log_r)` (native oracle widens per-species RE vector w/
  trailing `log_r_s`, community `mu_log_r` joins means), joint_grad / AGHQ
  (`n_quad=5` default). `coef()`/`vcov()` report `mu_log_r`; `ms_dispersion` carries
  `r`/`r_s`/`sigma_log_r`.
- **NUTS** (`method="nuts"`, #14; `R/ms_abun_nuts.R`, `src/ms_abun_nuts.cpp`): samples
  EXACT joint posterior over closed-form Royle marginal -> calibrated intervals +
  per-(species,site) WAIC/LOO (`.tobs_ploglik_ms_nmix`). NON-CENTERED: per-species block
  carries whitened `z_s ~ N(0,I)`, `b_{s,arm}=C_arm z_{s,arm}` (`C_arm`=log-Cholesky), so
  covariance enters only data term — breaks centered b<->Sigma funnel (322s->98s w/
  reparam + OpenMP at S=6). In-tree C++ FullGradFn parallelises per-species loop (OpenMP,
  deterministic reduction -> byte-exact vs R oracle `.tobs_ms_abun_nuts_logpost`); chol
  grad `chol_data_grad_noncentered` (`community_chol.h`, shared w/ #67). Warm-started at
  Laplace-EM mode + diagonal metric. Pois + NB, non-spatial. NUTS + areal (#73,
  car_proper only, Pois): the same per-species community sampler + a SHARED fixed-hyper
  non-centered proper-CAR field on abundance via an OPTIONAL field block in the eval
  (`n_field_units`/`field_map`/`field_Linv`; `f=Linv raw`, `raw~N(0,I)`, field precision
  `tau Q(rho)` fixed at the #12 nested-Laplace estimate). Field is a shared `eta_lambda`
  offset; its score sums `grad_eta_lambda` over species per unit -> `Linv' g_f - raw`.
  Field-off path byte-identical to non-spatial; full grad FD-validated.
  `.tobs_fit_ms_abun_nuts_spatial` warms from `nmix_community_laplace_car_proper` (#12).
  `test-ms-abun-nuts.R`. Design note `dev_notes/design_73.md`.

#### Areal-spatial community N-mixture (`ms_abun()` + shared field; sfMsNMix)

Shared ICAR/BYM2/proper-CAR field on abundance via `nested_laplace` (#12): term on
abundance -> `.tobs_fit_ms_nmix_spatial()` -> in-tree nested Laplace-EM
`nmix_community_laplace_{icar,bym2,car_proper}()` (`R/nmix_laplace_re_spatial.R`) over
`cpp_nmix_community_spatial_*` (`src/nmix_community_spatial.cpp`). Model
`log lambda_{s,i} = X_lambda_i.(mu_lambda + b_lambda_s) + f_{u(i)}`, one `f` shared. Top
block `(mu_lambda, mu_p, f)` shares single-species spatial layout
(`field_start = p_lambda+p_p`), so `nmix_spatial_kernel*.h` apply. Per outer grid point
EM iterates joint `(mu,f,{b_s})` mode-find (block-elim Newton, `b_s` Schur-folded) +
closed-form `Sigma` M-step; R grid-integrates. Community means FLAT (no ridge). NB `r`
grid-integrated. `fit$spatial_field` = posterior-mean field. `test-ms-abun-spatial.R`.
NUTS + shared field: #73 wires the #14 community sampler + an optional fixed-hyper
non-centered car_proper field block on the abundance arm (Poisson; 0 divergences,
field cor ~0.97); icar/bym2 + temporal/RE NUTS gated to n-L.

### Boundary: tulpaObs vs tulpa

- **tulpaObs owns**: family likelihoods (`src/*_likelihood.h`, `src/nmix_*.{cpp,h}`),
  E-step weights, M-step encoding, family S3/diagnostics. RE AGHQ: occ/det per-site
  marginal (`make_site` in `R/re_aghq.R`), native `NMixGroupedOracle`
  (`src/nmix_re_oracle.{h,cpp}`, `R/nmix_re_aghq.R`), `NMixCommunityOracle`
  (`src/nmix_community_oracle.{h,cpp}`); both exported as `XPtr<tulpa::REGroupOracle>`
  via `tulpa::tulpa_re_aghq(oracle=)`.
- **tulpa owns**: EM engine, MI/Gibbs, Rubin pooling, AGHQ engine (`tulpa_re_aghq()`:
  quad grid / log-Cholesky / LKJ / marginal Hessian), generic S3/diagnostics, NUTS/HMC.

NUTS crash for component w/ correct `populate_*` here = bug in tulpa
`hmc_sampler.cpp`, not here. File against `gcol33/tulpa`.

### Key design rules

- **Never pass `Rcpp::Nullable<T>` to header helpers** — MinGW crashes. Unwrap in
  `.cpp`, pass concrete types into headers.
- **Composition over registry** — families x spatial/temporal/re/svc/latent orthogonal,
  no per-combination branches.
- **tulpaObs defines likelihoods, tulpa handles structure.**
- **Dotted arg names (spOccupancy/base-R), never underscores.** Prefer single word
  (`visits`, `priors`, `control`); compound separated `.` not `_`: `n.iter`,
  `n.warmup`, `n.chains`, `n.thin`, `n.threads`, `adapt.delta`, `max.treedepth`,
  `max.iter`, `sigma.beta`. Governs `tobs()` args + every `control=list(...)` key
  (splatted into `.tobs_fit_model()`). Internal helper formals may stay underscore;
  anything user types = dotted.

## What works (tested)

Detail in Architecture above + per-family detail sections below. `n-L` = nested_laplace.

| Feature | Laplace | NUTS | Notes |
|---|---|---|---|
| Single-season occupancy | Yes | Yes | parity w/ inlaocc |
| Dynamic (HMM) | Yes | Yes | colonization/extinction. **Season-varying gamma/epsilon (#124)**: a covariate carried as a `[n_sites x (T-1)]` matrix column of `data` drives interval-indexed colonization/extinction (`colonization = ~ cov`), the dyn_abun #80 recipe ported to the colext forward. Detected via `.tobs_interval_arm_design` (shared w/ dyn_abun); the arm design unrolls long-form `[(site x interval) x p]`. E-step forward-backward uses per-interval transition matrices `A**[i, iv]` (constant rates recycle -> byte-identical to pre-#124); M-step encodes the transition arm as a WEIGHTED LOGISTIC (response = transition prob given origin state `xi/n_mat`, weight = origin-state prob `n_mat`) NOT the per-site M=1000 pseudo-binomial (which, applied per interval, makes each row near-separable so the inner Newton returns ~0 slope). The exact-marginal refine handles SV via a pure-R season-varying HMM-forward marginal `.tobs_dyn_occu_marginal_nlp_sv` (`R/dyn_occu_marginal.R`; the cpp `cpp_occu_dynamic_ploglik` reads one rate/site) -- escapes the EM local optimum AND calibrates SEs. `simulate_dyn_occu(beta_gamma=, beta_epsilon=)`. Laplace/gibbs/mi only; NUTS+SV gated (C++ forward per-site). 20-seed recovery + coverage (`test-dyn-occ.R`) |
| Community single-season (`ms_occu`) | Yes | Yes | per-arm community RE, shared community Laplace-EM (`R/community_em.R`, `R/ms_occu.R`); msPGOcc. NUTS non-spatial samples community means/deviations/covariances jointly (#69, `R/ms_occu_nuts.R`); shared areal field icar/bym2/car_proper on occ arm via n-L (#75, `R/ms_occu_spatial.R`, sfMsNMix analogue). NUTS+field -> n-L. **SVC (`svcMsPGOcc`, #118)**: a varying-coefficient bar `spatial(~ 1 + w \|\| cell, graph)` on the occ arm -> intercept + SVC field(s) via BLOCK COORDINATE ascent (`R/ms_occu_field.R`, `.ms_occu_field_solve` = two-state-marginal occupancy field Newton + ICAR prior; community occ EM w/ field as psi offset). Both fields recover ~0.96. icar only; plain intercept field stays on the C++ path (no regression). `test-ms-occu-field.R` |
| Community dynamic (`ms_dyn_occu`) | Yes | — | per-species psi1/p RE + shared gamma/eps; HMM-forward; `R/ms_dyn_occu.R`. **Shared areal field on psi1 (stMsPGOcc, #123)**: `~ 1 + icar(graph=adj)` under `nested_laplace` -> `R/ms_dyn_occu_spatial.R`. KEY MATH: psi1 sets ONLY the initial HMM mixing weight, so the per-(species,site) marginal is LINEAR in psi1 -- `L=(1-psi1)A+psi1 B` w/ A/B = HMM-forward likelihood cond on season-1 state 0/1 -- IDENTICAL mixture to the single-season `ms_occu` field oracle, so the Louis score/curv (`score=r-psi1`, `curv=psi1(1-psi1)-r(1-r)`, `r=psi1 B/L`) carry over verbatim (only A/B differ). Routed through the shared block-coordinate driver (`R/community_latent.R`, `.tobs_community_latent_ascent`): community EM (field=psi1 offset) alternated w/ areal Newton. Kernels VECTORIZED over sites (`.ms_dyn_occu_condAB_vec`/`.ms_dyn_occu_fb_vec`) + ANALYTIC Fisher-identity gradient (fwd-bwd smoothed w1/pairwise xi -> grad_psi1=X'(w1-psi1), grad_p=X'(sum_t w_t(ndet-nvalid p)), grad_gamma=X'(col_y-gamma col_n)) so the community EM skips its O(U^2) FD Hessian (FD-validated 1e-9; a per-site-loop version was ~10x too slow). icar only; field cor ~0.94, community-mean coverage >=0.85. `simulate_ms_dyn_occu(field=)`. NUTS/bym2/car follow-ups |
| Community integrated (`ms_int_occu`) | Yes | — | per-species psi + per-source det RE; multi-source two-state marginal; `R/ms_int_occu.R` |
| Integrated multi-source | Yes | Yes | shared psi |
| JSDM (`jsdm`) | Yes | Yes | `jsdm()` = the COMMUNITY GLMM on observed presence/absence (#121): per-species coefs + Gaussian community covariance, NO detection/latent state = `ms_count()` w/ logit link, so it SHARES that binder (`.tobs_build_ms_count(response="bernoulli")`, model_type `ms_count`), community EM, latent driver, NUTS target + S3. `latent(n)` = lfJSDM; + shared field = sfJSDM; icar/car_proper/bym2/spde field via n-L. NUTS = exact joint community posterior over the Bernoulli response (`MSC_BERN` in `src/ms_count_nuts.cpp`, byte-exact vs R oracle). laplace_sla/gibbs/mi were single-block routes, dropped |
| Cover hurdle (joint) | Yes | — | `family_cover_hurdle.R`, `sla_cover_*`. positive = beta/lognormal/lognormal_trunc/ordinal/`beta_oi`. `beta_oi` (#108) = one-inflated Beta: ceiling (cover=1) plots = a point mass (constant pi = ceiling share, binomial SE), interior Beta on (0,1); encode splits `is_pos` to interior, `enc$oi` carries pi, decode reports `pi_one`, predict conditional cover = `pi + (1-pi)*mu` (`.tobs_cover_mu`). `control$aggregate.occ` (ON, #48) collapses occ arm to Binomial suff-stat; `control$aggregate.pos` (ON beta arm, #49) collapses beta pos arm to grouped suff-stat (tulpa `slog_y`/`slog_1my`), errors on non-beta. Both byte-identical to per-plot |
| Cover hurdle spatial coef fields (`\|\|` / `\|`) | n-L | — | `spatial(~ 1 + w \|\| node, graph, to=)` independent (#61, two coupled ICAR blocks, per-field alpha) OR `\| ` correlated (#64, one separable-MCAR block sharing free Sigma). `\|` both-arm `to=c("presence","positive")` = copied to pos arm w/ one alpha (#64); `\|` single-arm `to="presence"`/`"positive"` = free-Sigma field on that arm alone, NO copy (#109, 0-sentinel `spatial_idx` on the other arm via `mc$to`, `copy=NULL`, `alpha_mcar`=NA). `\|` -> `.cover_build_mcar_spec`/`.fit_cover_hurdle_joint_mcar` (tulpa `type="mcar"` block, copy only when both-arm); reports `sigma_mcar`/`rho_mcar`/`alpha_mcar`. SLA on `\|` no-op. icar only |
| Cover hurdle arm-specific fields (single-arm `to`) | n-L | — | `spatial(~ 1 + w \|\| cell, graph, to="positive")` (or `"presence"`); separate single-arm calls = independent per-arm fields, NO cross-arm copy (#65). NO engine change: per-arm `spatial_idx=0` makes the other arm's rows skip the block (tulpa `l_b>0` scatter guard), own precision grid-integrated. `.tobs_armspecific_bar_fields` (formula_terms.R) -> `enc$armspec` -> `.fit_cover_hurdle_joint_armspecific` (non-copied per-arm blocks, no `copy=`). `armspec_blocks` carries per-block arm/slot/type; `.tobs_joint_draws_cover_armspecific` scatters each block onto its arm only (amp 0 on other). icar/car/car_proper AND bym2 (#107): a bym2 block is the non-copied length-2 (phi ICAR + iid theta) block, paired (sigma,rho) grid; the draw projection reconstructs the rho-mixed unit field `z = sqrt(rho)*sf*phi + sqrt(1-rho)*theta` so predict/WAIC see the full mix. `\|\|` only (`\|` arm-specific undefined: copy-only). No mix w/ shared field/trend/temporal/re; one field per arm. SLA no-op |
| Joint occu + cover | Yes | Yes | `occu_cover()` — see below. NUTS non-spatial (in-tree FullGradFn over exact two-state marginal, beta/lognormal; `R/occu_cover_nuts.R`, `src/occu_cover_nuts.cpp`) AND spatial fixed-hyper coupled proper-CAR field (#74, `car_proper()` only; icar/bym2 NUTS gated -> n-L). Sampled-field (estimated-variance) route = `ms_occu_cover()` factor |
| occu + cover + areal field + per-group RE | n-L | — | `occu_cover()` + icar/bym2 + `re(g)`/`(1\|g)` on psi; one iid RE block (#56); `sigma_re` + BLUPs; intercept RE only |
| occu + cover independent cover-arm field (single-arm `to="positive"`) | n-L | — | `occu_cover()` + `spatial(~ 1 + w \|\| cell, graph, to="positive")` on occurrence formula (#110): NON-copied ICAR block(s) on the cover arm alone, decoupled from the occupancy field's alpha copy. Composes w/ the shared occupancy field (psi + `delta_cover_exp` keep it) so `delta_cover_cond` varies instead of collapsing when `alpha->0`. Parse: `.occu_cover_spatial_fields` splits single-arm `to="positive"` bar -> `spatial_info$pos_armspec` via `.tobs_armspecific_bar_fields`. Fitter appends non-copied ICAR blocks (`spatial_idx` psi=0/p=0/pos=cell, `tau_grid`, svc_weight for trend) after occ fields, before RE blocks; `ctx$field_specs` labels each field block shared-vs-pos. Postprocess partitions occ (1..n_occ_fields) vs pos blocks; reports `sigma_pos_field`/`sigma_pos_field_<col>` from `b<k>.tau`; surfaces `fit$pos_field`/`pos_field_table(s)`. Draw substrate (`.tobs_joint_draws_occu_cover`) reads `field_specs`: pos block amp_occ=0, amp_pos=1/sqrt(tau); predict/WAIC add it to `field_pos` automatically. Predict weight override skips pos blocks. Per-visit cover (`cover_aggregate="none"`); icar only (bym2/car->icar); NOT w/ MCAR `\|`, latent cover RE, or batch. Static intercept field weakly ID'd vs alpha copy; time-weighted trend cleanly ID'd. `test-occu-cover-pos-field.R` |
| occu + cover independent detection-arm field (single-arm `to="detection"`) | n-L | — | `occu_cover()` + `spatial(~ 0 + time \|\| cell, graph, to="detection")` on the `detection` formula (or lifted `to="detection"`) (tulpa#140): spatially-structured detection prob. Same non-copied arm-specific block machinery as the cover-arm field (#110) -- builder `arm_field_blocks(af, "p")` sets slot 2 (`spatial_idx` psi=0/p=cell/pos=0). Enters via detection arm `field_coef=1` (set when `det_armspec` present, `.occu_cover_build_joint_arms(det_field=)`) so the block scatters onto the p rows; shared occ field kept off detection by its `spatial_idx=0` sentinel -- identical to the #102 detection-RE mechanism, so NO tulpa engine change. Reports `sigma_p_field`/`sigma_p_field_<col>`. icar only; per-visit cover. `test-occu-cover-pos-field.R` |
| occu + cover obs-arm RE (detection / pos) | n-L | — | `occu_cover()` + RE on `detection=`/`positive=`; per-visit grouping. Intercepts: `(1\|g)`, crossed `(1\|g)+(1\|h)`, nested `(1\|g/h)` = N iid blocks (#102 single, #103 crossed/nested). Slopes (#103, tulpa>=0.0.39): uncorr `(x\|\|g)`/`(0+x\|g)` = weighted iid per coef; corr `(1+x\|g)` = miid free-Sigma block. `sigma_re_p`/`sigma_re_pos` (+`_<coef>`, `cor_re_p_*`) + BLUPs |
| occu + cover + detection / cover-arm RE | n-L | — | `occu_cover()` + `(1\|g)`/`re(g)` on `detection=` or `positive=` (#102); per-VISIT grouping iid block on that arm, composes w/ psi field; `sigma_re_p`/`sigma_re_pos` + BLUPs; `fit$re` per-arm list; `predict(type="detection")` adds BLUP offset, unseen->pop mean. Detection arm `field_coef=1` (not 0) so the iid block scatters, field still skipped via `spatial_idx=0` sentinel. pos RE needs `cover_aggregate="none"`. Intercept only; slope/correlated/non-spatial/NUTS gated |
| Community joint occu + cover | Yes | — | `ms_occu_cover()` — see below |
| Spatial-factor community occu+cover (JSDM) | Yes | Yes | `ms_occu_cover()` + icar/car_proper/bym2 shared field, per-species loadings (tulpa#67). Laplace-EM (`R/ms_occu_cover_spatial.R`) + NUTS (`src/ms_occu_cover_spatial_nuts.cpp`). Cover-arm factor, `tobs_associations()`, per-species `predict()` maps |
| Multiscale occu + cover | n-L | — | `occu_multiscale_cover()` — 3-level cell/plot/visit + cover; spatial joint only |
| Count / relative-abundance GLMM | Yes | — | `count(response=)` (spAbundance `abund`): GLMM on the observed response directly, NO detection, NO latent state -- the abundance analogue of `jsdm()`. Pois/negbin (log link) + gaussian (identity) + **binomial (logit, `k`-of-`n`, `trials=`; spOccupancy `svcPGBinom`, #125)**. One tulpa GLMM block (`build_count_callbacks`, `R/laplace_callbacks.R`); negbin size / gaussian variance by an outer dispersion loop in `.dispatch_count` (tulpa_laplace takes fixed phi), reported `fit$count_dispersion`; binomial has NO dispersion (variance pinned by `n`) so it takes the Poisson no-loop path + carries per-site `model$n_trials`. `simulate_count()`, `.tobs_ploglik_count` (WAIC), `count_methods.R` (fitted = expected successes `n*p`, predict newdata = per-trial prob, residuals binomial deviance/pearson). Community (msAbund) binomial = `ms_count("binomial")`; NUTS pending. Non-spatial binomial matches `glm()` MLE to ~1e-3. `test-count.R` (3x 20-seed recovery + binomial trials>1 & trials=1) |
| Community spatial-factor (`ms_count`) | n-L | — | `ms_count()` + `icar()` + `latent(n)` (spatial-factor `sfMsAbund`; #117): shared areal field AND latent factors on ONE formula, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}+sum_q lambda_{s,q} eta_{q,i}`. One block-coordinate loop runs BOTH the field update and the factor update; when both present the factor loadings are CENTERED across species (`sum_s lambda_{s,q}=0`) so the field owns the shared spatial mean, factors own between-species residual co-occ -- both recover (field cor ~0.99, residual cor ~0.99). UNIFIED fitter `.tobs_fit_ms_count_latent` (field-only/factor-only/both = one source of truth; old fitters = thin wrappers). `test-ms-count-factor.R` |
| Community latent factor (`ms_count`) | Yes | — | `ms_count()` + `latent(n)` (spAbundance `lfMsAbund`; #117): residual species co-occurrence via Q per-site latent factors + per-species loadings, `log mu_{s,i}=X_i(mu+b_s)+sum_q lambda_{s,q} eta_{q,i}`. BLOCK COORDINATE via the shared driver (`R/community_latent.R`): community EM w/ factor offset <-> factor update (alternating Newton eta\|lambda, lambda\|eta + unit-variance anchor). loadings/factors identified up to ROTATION -> recoverable target is residual cov `Sigma_res=lambda lambda'` (`fit$ms_factor$residual_cov`/`residual_cor`/`loadings`/`factors`); residual cor recovers ~0.95. `fitted()`/WAIC factor-aware via `model$count_factor_offset`. Poisson; not composed w/ areal field yet. `test-ms-count-factor.R` |
| Community SVC (`ms_count`) | n-L | — | `ms_count()` + `spatial(~ 1 + w \|\| cell, graph)` (spAbundance `svcMsAbund`; #117/#118): intercept field + one shared spatially-varying-coefficient field per covariate, `log mu_{s,i}=X_i(mu+b_s)+f0+w_i*f1`. Same block-coordinate scheme; the field solve generalizes to K covariate-weighted ICAR fields solved jointly (K x K sparse block system, per-field tau). `fit$spatial_field`=intercept field, `fit$trend_field(s)`=SVC field(s). Both recover cor ~0.98. icar only. `test-ms-count-spatial.R` |
| Community count + areal (`ms_count`) | n-L | — | `ms_count()` + `icar()` (spAbundance `sfMsAbund`; #117): ONE shared areal field across species, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}`. Fit by BLOCK COORDINATE ASCENT (`R/ms_count_spatial.R`) -- community Laplace-EM w/ the field as a per-site OFFSET (captured in sp_ll closure, NO community_em.R change) alternated w/ a self-contained Poisson-ICAR field update (analytic sparse Newton `.ms_count_field_solve` + closed-form tau M-step). PURE R, no C++ (sidesteps the community EM's FD-Hessian which doesn't scale to an O(n_sites) field). Field informed by all species/site -> recovers cleanly (cor ~0.98); `fitted()`/WAIC field-aware via `model$count_field_offset`. Poisson + icar ONLY (overdispersed community count NOT identified vs per-site field; bym2/car/group_var/bar follow-ups). `test-ms-count-spatial.R` (community-mean + pooled coverage + field recovery, 20 seeds) |
| Community count (`ms_count`) | Yes | Yes | `ms_count(response=)` (spAbundance `msAbund`; #117): community relative-abundance GLMM -- per-species GLMM on the observed response w/ Gaussian community hyperpriors on the coefficients, NO detection/latent state. The community analogue of `count()`, abundance analogue of `ms_occu()`. Pois/negbin (per-species dispersion RE = second arm)/gaussian (per-species resid var, outer loop)/**binomial (logit `k`-of-`n`, `trials=`; community `svcPGBinom`, #125 -- laplace non-spatial only, NUTS + shared field/latent gated w/ pointer; community-mean intercept carries a small O(1/n_species) first-order-Laplace bias, absent at `trials=1`, slope unbiased)**. Reuses shared community Laplace-EM (`.tobs_community_em`, `R/community_em.R`) -- PURE R, no C++ -- w/ count `sp_ll`/`sp_grad` (`R/ms_count.R`). `y` = `[n_sites x n_species]` or named list. `coef`=community means, `ranef`=per-species deviations, `simulate_ms_count()`, `.tobs_ploglik_ms_count` (WAIC). Gaussian/Pois community means unbiased; negbin slope carries mild first-order-Laplace attenuation (~10%, documented). **NUTS (Pois/negbin/gaussian, #117)**: samples the exact joint posterior via a family-aware in-tree C++ FullGradFn (`src/ms_count_nuts.cpp`, `R/ms_count_nuts.R`), the reduced `ms_abun` NUTS (no detection/latent-N); non-centered `b_{s,arm}=C_arm z_{s,arm}`; byte-exact vs the R oracle (1e-7), warm-started at the Laplace-EM mode, 0 divergences, NUTS==Laplace. negbin adds a per-species dispersion RE `log_r_s~N(mu_log_r,sigma_log_r^2)` as a SECOND community arm (block-diag chol, mirrors `ms_abun`); gaussian adds S FREE per-species `log_phi_s` (no community prior, weakly-informative N(pooled log-var, 2^2), matching the Laplace outer loop) -- NUTS removes the negbin Laplace attenuation, `fit$ms_dispersion` carries `r_s`/`variance`. Missing (`NA`) site x species entries are dropped from the per-(species,site) data sum (NA-in-`y` IS the mask; both C++ + R oracle skip identically -> byte-exact w/ NA present, `N` counts observed), matching the Laplace path's per-species `valid` subsets. `test-ms-count.R` + `test-ms-count-nuts.R` |
| Count + areal | n-L | — | `count()` + `icar()`/`car_proper()` (spAbund; #117): a plain areal field on the abundance formula. Response is observed (no latent state), so the fit is ONE `tulpa_nested_laplace()` call over the count block with the areal field as its GMRF prior -- NOT the occupancy EM. Grid-integrated FE (law of total cov via per-cell `keep_grid_hessians`) + field/sigma via shared `.tobs_nested_attach_field_summary`; dedicated `.tobs_fit_count_spatial` (`R/count_spatial.R`), reuses `.tobs_to_multi_block_prior`. `fitted()` field-aware in-sample. **Poisson OR binomial** -- Poisson has no dispersion; binomial's variance is pinned by `n` so it is identified against a per-node field (spOccupancy `svcPGBinom`, #125; works at `trials=1` too -- one Bernoulli/node is weak but not confounded). negbin size / gaussian residual variance are NOT jointly identified with the field under the fixed-phi loop (degenerate: size->Inf, resid var->0), gated w/ pointer; bym2 + improper car() gated to icar/car_proper. `test-count-spatial.R` (20-seed FE coverage + field recovery, car_proper, gates, + binomial trials>1 & Bernoulli field recovery) |
| N-mixture + areal | n-L | — | `abun()`+icar/bym2/car_proper; Pois/NB (r grid-int); grid-int cov (constrained intercept) |
| Community N-mixture | Yes | Yes | `ms_abun()` (msNMix); per-species coef RE, in-tree Laplace-EM (`nmix_laplace_re`) -> `NMixCommunityOracle` AGHQ, Schur SEs; Pois + negbin. NUTS (#14) |
| Community N-mixture + areal | n-L | Yes | `ms_abun()`+icar/bym2/car_proper (sfMsNMix; #12); shared field + per-species RE; `nmix_community_spatial.cpp`; Pois/NB. NUTS (#73, car_proper only): the #14 non-centered community sampler + a SHARED fixed-hyper non-centered proper-CAR field on abundance (tau Q(rho) fixed at the #12 nested-Laplace estimate, raw ~ N(0,I), f=Linv raw; optional field block in `src/ms_abun_nuts.cpp`, FD-validated, field-off byte-identical to #14). 0 divergences, field cor ~0.97; Pois only. icar/bym2 NUTS gated to n-L |
| Community N-mixture latent factor (`ms_abun`) | Yes | — | `ms_abun()` + `latent(n)` (spAbundance `lfMsNMix`; #117): residual species co-occurrence on the ABUNDANCE arm via Q per-site factors + per-species loadings. The latent N still marginalises closed-form per species-site, so the whole latent structure sits on `eta_lambda` and the family reduces to ONE working oracle over the Royle marginal (`R/ms_abun_latent.R`) driving the shared block-coordinate engine: `score = grad_eta_lambda`, `curv = info_eta_lambda - var_N * score_wt_lambda^2` = the Louis (1982) (1,1) block (abundance curvature, detection arm profiled out), which `nmix_site_marginal()` ALREADY exposes -- no new kernel. Oracle FD-validated (score 2.6e-9; curv diff = FD truncation; cor 1.000000 both). Residual cor recovers ~0.99 at N=120/S=12 (counts are information-rich -- cleanly better than occupancy's ~0.9). Poisson only (a negbin size is a second per-site dispersion, not identified vs a per-site latent). `test-ms-abun-factor.R` |
| Community N-mixture spatial-factor (`ms_abun`) | n-L | — | `ms_abun()` + `icar()`/`car_proper()`/`bym2()`/`spde()` + `latent(n)` (#117): shared field AND factors on the abundance arm; the centred loadings separate them (field cor ~0.98, residual cor ~0.94). Same block-coordinate driver. A plain field with NO factors KEEPS the dedicated C++ path (#12) -- faster + already recovery-tested, deliberately not replaced |
| N-mixture + grouped RE | Yes | — | `abun()`+`(1\|g)`/`(x\|g)` either arm (#13); non-species grouping; Pois/NB; `NMixGroupedOracle`. Gated: RE+spatial, RE+visit-det, RE both arms |
| Community distance (`ms_distance`) | Yes | — | `ms_distance(key=, transect=, cutpoints=)` (spAbundance `msDS`; #117): per-species binned distance sampling w/ Gaussian community hyperpriors on the abundance (`log lambda`) + detection-scale (`log sigma`) coefs. `y` = `[n_sites x n_bins x n_species]` or named list. NO new C++ — the latent N marginalises closed-form per species-site, and `cpp_distance_site_sweep` already returns `log_lik`/`grad_lam`/`info_lam`/`var_N`/`swl`, so the shared community Laplace-EM (`R/ms_distance.R`) reads its per-species score straight off it. Hazard-key log-shape = a community `global` (shared across species) via the EM's `global` slot. Community means UNBIASED over 10 seeds (z of bias 1.07/1.68/-0.66), 95% Wald coverage 1.0/0.8/0.9; single-species `distance()` control agrees (bias 0.02/-0.01). NB: a single seed shows a ~0.18 paired lambda-down/sigma-up shift that looks like bias and is NOT — it is one draw on the lambda/sigma ridge. Poisson only (negbin size not yet a per-species RE); NUTS not wired. `simulate_ms_distance()` draws through `cpp_simulate_distance` (a separate R-side quadrature would simulate from a pi the model is not fit against). `test-ms-distance.R` |
| Community distance latent factor / spatial-factor | Yes | — | `ms_distance()` + `latent(n)` (`lfMsDS`) and + `icar()`/`car_proper()`/`bym2()`/`spde()` (`sfMsDS`), #117: the SAME shared block-coordinate driver, whose working oracle is the SAME Louis formula as `ms_abun` — `score = grad_lam`, `curv = info_lam - var_N * swl^2` (both are count marginals with `B_i = diag(info) - var_N v v'`; only the kernel supplying the pieces differs). Oracle FD-validated (score 1.6e-9, curv = FD truncation, cor 1.000000). Poisson only. Slower than `ms_abun`'s latent path: the distance kernel exposes the per-site Louis pieces but no assembled block, so this family does not yet pass `sp_info` and the EM finite-differences its per-species Hessian |
| Removal (Pois/NB) | Yes | Yes | `removal()` (#39) — see Architecture. NUTS single intercept RE; Laplace grouped RE one arm; areal icar/car_proper/bym2 abundance arm; areal+temporal AND temporal-only AR1/RW1/RW2/iid on abundance arm via shared areal-BFGS (#78/#114); NUTS+areal car_proper field on abundance arm (#72) |
| Distance (Pois/NB) | Yes | Yes | `distance(key=, transect=, cutpoints=)` (#38); `formula`=log lambda, `detection`=log sigma, `y`=`n_sites x n_bins` — see Architecture. NUTS single abundance intercept RE; Laplace grouped RE abundance arm (half-normal); areal field (half-normal + hazard key); areal+temporal AND temporal-only on abundance arm (#78/#114); DETECTION-arm areal field on log sigma (`detection=~icar()`, spatially-varying detection scale, half-normal, #114); NUTS+areal car_proper field on abundance arm (half-normal, #72) |
| False-positive occupancy | Yes | Yes | `fp_occu()` (#40) — see Architecture. NUTS single psi intercept RE; Laplace grouped RE psi OR p11; areal psi arm; areal+temporal AND temporal-only on psi arm (#78/#114); NUTS+areal car_proper field on psi arm (#72) |
| Open N-mixture (Dail-Madsen) | Yes | Yes | `dyn_abun()` (#37); y 3D `[n_sites x J x T]` — see Architecture. Pois/NB init; season-varying omega/gamma via `[n_sites x (T-1)]` covariate, interval-indexed forward kernel, all backends (#80). NUTS single intercept RE on init-abundance OR detection arm; Laplace grouped RE on init-abundance OR detection arm (#82, p-arm = per-node full-forward second-order eta_p pass); areal init arm; areal+temporal AND temporal-only on init-abundance arm (#78/#114); NUTS+areal car_proper field on init-abundance arm (#72); NUTS+temporal-only fixed-hyper ar1/rw1/rw2/iid field on init-abundance arm (#114, same field block as areal, `field_map`=period, temporal whitened loading; areal+temporal under NUTS gated to n-L) |
| Spatial ICAR/BYM2/NNGP | — | Yes | |
| Spatial + dynamic | — | Yes | |
| Nested-Laplace (areal) | n-L | — | icar/bym2/car (+temporal/iid) on occu/int_occu/dyn_occu. RECOVERY-TESTED on icar for all three (`test-{occu,int-occu,dyn-occu}-areal-recovery.R`) + SVC bar (`test-dyn-occu-svc-recovery.R`). State arm encodes at **M = 1** for ANY nested latent block (`nested_block <- !is.null(latent_prior)` in all three `build_*_callbacks`): a state row = ONE binary occupancy obs, so the old M=1000 (areal) / M=4 (spde) overstated its info M-fold, swamped the field prior, read between-cell noise as field, inflated the slope via `sqrt(1+0.346 sigma^2)`. Monotone in M, no arm regresses at M=1. `use_louis` (laplace_helpers.R) covers single/dynamic/integrated — state score `x_i(z_i-psi_i)` is family-generic; without it the state block falls to `.se_from_laplace_fit()` -> **NA SEs** (dynamic + integrated both shipped that way). dynamic needs `weights[,1]` (season-1 col), integrated's is already per-site. **A null-field fixture CANNOT test this path**: grid = `b1.tau` 9 log cells [0.3,30] (sigma [1.83,0.18]), so truth sigma=0 is OFF-grid and the marginal pins on the last cell (weight 0.60) — "shrunk" and "pinned" read identical. Test with an INTERIOR field + assert peak cell off both boundaries. `tau` = ICAR CONDITIONAL precision, NOT 1/sd^2: marginal field sd 1.0 -> sigma~0.40, peak cell 6/9 — assert on field sd/cor, never on `sigma`. `dev_notes/finding_dyn_nested_laplace_field.md` |
| NA-response prediction | n-L | — | `predict(type="state")`: all-NA single-season sites field-interpolated, calibrated 95% from exact-marginal bernoulli pass (coverage ~1.0) |
| Formula RE (intercept) | Yes | Yes | `(1\|g)`; variance-component EM, occ OR det arm (#11) |
| Formula RE (uncorr slope) | Yes | Yes | `(x\|\|g)`, `(0+x\|g)`, `(1+x\|\|g)`; either arm |
| Formula RE (corr slope) | Yes | Yes | `(1+x\|g)`; EM M-step consumes tulpa `cov_blocks`; either arm (#11) |
| Formula RE on detection | Yes | Yes | `detection=~(1\|g)`; separate det-arm RE block, AGHQ branches on arm |
| All S3 methods | Yes | Yes | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, `$.tobs_fit` |
| Diagnostics | Yes | Yes | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatial_range |
| Simulation | Yes | Yes | `simulate_occu/ms_occu/dyn_occu/int_occu/dyn_ms_occu/int_ms_occu/cover/cover_joint` |
| Spatial prediction | — | Yes | `tobs_predict_spatial` (IDW on field) |
| Components | Yes | Yes | `tobs_re`, `tobs_temporal`, `tobs_svc`, `tobs_latent`, `tobs_community_re`, `tobs_areal` |

### `occu_cover()` detail

Joint occu-detection + cover hurdle. **Non-spatial** (`laplace`, `R/occu_cover.R`):
cell psi + per-visit detection + per-visit cover (beta/lognormal) on exact two-state
marginal. **Non-spatial NUTS** (`nuts`, `R/occu_cover_nuts.R`, `src/occu_cover_nuts.cpp`):
samples that exact two-state coefficient marginal (no field/RE) via in-tree C++
FullGradFn over shared `run_tulpa_nuts`, warm-started at Laplace mode w/ diagonal
Laplace metric. Reuses `nodet_mixture_block` + `Lognormal/BetaPositive` (cover density
+ eta-grad + new `grad_logdisp` dispersion score, single-source in
`occu_coupling_shared.h`); coef grad = design-sandwiched eta-grad. Layout
`c(beta_psi, beta_p, beta_pos, log_disp)`, weak N(0,sigma.beta^2) coef + broad
N(0,sigma.logdisp^2=25) log-disp priors. Byte-exact R oracle
`.tobs_occu_cover_nuts_logpost` vs `cpp_occu_cover_nuts_joint_logpost`; calibrated
WAIC/LOO + Rhat/ESS (shared `.tobs_nuts_rhat_ess` in `nuts_chains.R`). 0 divergences,
NUTS==Laplace mode, recovery + 95% coverage (`test-occu-cover-nuts.R`). **Spatial NUTS**
(`nuts` + `car_proper()` on psi, #74, `R/occu_cover_nuts.R::.tobs_fit_occu_cover_nuts_spatial`):
FIXED-HYPER non-centered coupled PROPER-CAR field — psi-arm field `f` (one/cell) enters psi
linearly, copied to pos arm w/ scaling `alpha`; field precision `tau Q(rho)` + `alpha` FIXED
at the nested-Laplace joint proper-CAR estimate (`.tobs_occu_cover_nuts_carproper_warm`:
one warm `tulpa_nested_laplace_joint(type="car_proper")` fit), whitened `raw ~ N(0,I)` (`f =
Linv %*% raw`) sampled jointly. Param vector `c(beta_psi, beta_p, beta_pos, log_disp,
raw_field)`. Reuses the abun#51 field-block recipe (tulpa#87): the optional field block in
`src/occu_cover_nuts.cpp` (`n_field_units`/`field_map`/`field_Linv`/`field_alpha`), byte-exact
vs R oracle's field branch; field-off path byte-identical to non-spatial NUTS. car_proper only
(full-rank precision -> well-conditioned geometry); icar/bym2 + SVC/trend/temporal/RE gated ->
n-L. group_var maps sites>cells. 0 divergences, field cor ~0.78, 95% slope coverage ~0.92, beta
SD calibrates to n-L SEs (`test-occu-cover-spatial-nuts.R`). predict() needs the joint object
(non-spatial laplace AND nuts both error w/ pointer); sampled-field (estimated-variance) route =
`ms_occu_cover()` factor (tulpa#67). **Spatial
default** (`nested_laplace`, `R/occu_cover_joint.R`):
`joint` engine via `tulpa_nested_laplace_joint()` w/
`occu_cover_{lognormal,beta}` cell-coupling spec (tulpa#32) — 3-arm joint
nested-Laplace, outer-grid over `(sigma, alpha)`, per-cell occupancy mixture
closed-form derivs drive inner Newton. ~150-300x faster than v3_nested at N=100,
completes N=200+. 18-seed lognormal + beta recovery
(`test-occu-cover-joint-coupled.R`); status `"working"` (#96). Shared-field occ
SLOPE Wald CI mildly anti-conservative small-N (pooled ~0.94; NUTS non-spatial
calibrated). Outer Pareto-k diagnostic (`control$diagnose.k`) defaults OFF
(#101): 84-98% of joint-fit wall time (re-solves the inner Laplace `k.samples`=200x
on the full field vs the grid's ~30-70). Per-phase profiling
(`dev_notes/_profile_pareto_k.R`) corrected the cause: the binding per-solve cost
is the per-Newton-iter Hessian/grad SCATTER (beta arm digamma/trigamma fill, 73-83%),
NOT the factorize (flat ~0.5ms, 8-12%, not super-linear <=1100 cells). tulpa#118 sped
the diagnostic 2-4x (Shamanskii reuse `.K_DIAG_REFRESH` -> grad-only scatter; loosened
inner tol `.K_DIAG_TOL`=1e-4; near-neighbour batch order) with the k-hat byte-stable
(externally validated == `loo::psis`/`posterior::pareto_khat` on real EVA ratios; knobs
`tulpa.kdiag.{refresh,tol,reorder,capture}`). Still OFF by default: reports k-hat only,
fit byte-identical on/off, opt-in; matches `occu_joint` /
`occu_multiscale_cover`. `control$diagnose.k = TRUE` re-enables. Same default flip on
`occu_multiscale_cover_joint.R`.

**Cover-arm intercept prior (#32)**: on shared-field path cover intercept confounds
w/ field level over detected cells. `.occu_cover_coupled_arm_priors()` hands pos arm
`cover_priors()` weakly-informative intercept prior **by default** (not engine flat
ridge); else cover intercept floats to huge SD, `predict()` conditional cover blows up
via Jensen. `priors = FALSE`/`"none"` disables all three arms.

**Cell-aggregated cover (`cover_aggregate`, #33)**: per-visit cover gives cover arm one
row/visit, so a cell w/ many detected plots drives field more than its single occupancy
obs. `cover_aggregate = "mean"` (default spatial) / `"median"` collapses cover arm to
ONE row/occupancy unit; `"none"` keeps per-visit. ONLY on spatial `joint`.
Needs cell-level positive design (from `data`); visit-level `positive` keeps per-visit.
C++ compile-time `Aggregated` flag on `OccuCoverCoupling`
(`src/cell_coupling_occu_cover.h`), registered `occu_cover_{lognormal,beta}_agg`; R
`.occu_cover_build_joint_arms(cover_aggregate=)`. `test-occu-cover-coupling.R`,
`test-occu-cover-aggregate.R`.

**Latent cover-per-unit (`cover_aggregate = "latent"`)**: principled mean/median
alternative — cover arm carries per-unit cover RE `u_i ~ N(0, sigma_u^2)` shared across
the unit's detected visits, integrated out per unit (keeps every visit). Unit-level
predictor -> per-unit marginal `log M_i` SCALAR in one eta, reuses one-row-per-unit
layout. Lognormal = closed form (`src/occu_cover_latent.h::LognormalLatent`); beta =
adaptive GH over `u_i` reusing `BetaPositive` (`BetaLatent`, `control$n.quad` default
15). Within-unit dispersion pre-fit + held FIXED; `sigma_u` rides pos arm `phi_grid`
(`control$sigma.u.grid`), reported `phi_pos`. Stateful spec
(`OccuCoverLatentCoupling<PosLatent>`, `src/cell_coupling_occu_cover_latent.h`),
(re)registered per fit via `cpp_register_occu_cover_{lognormal,beta}_latent_coupling()`.
Shared det-branch in `occu_det_psi_p_block`/`occu_nodet_block`
(`src/occu_coupling_shared.h`). Same gates as aggregation. `test-occu-cover-latent.R`.

**Coupled SVC/trend fields** (#15): extra shared areal fields = WEIGHTED areal terms in
psi formula — `icar(graph=adj, weight=year)` couples spatially-varying coef on `year`
atop unweighted intercept field. N fields compose (each own outer-grid `alpha`) via
tulpa multi-block copy. Intercept field = `fit$spatial_field`; weighted fields
`fit$trend_field`/`fit$trend_fields`, scale `alpha_trend` (`control$alpha.grid.trend`).
p arm excluded via `field_coef=0` (NOT `svc_weight=0`). Resolved by
`.occu_cover_spatial_fields()`; off-path errors via `.tobs_reject_weighted_spatial()`.
`test-occu-cover-trend.R`.

**Escape hatches**: `control$engine="v3_nested"` (pure-R outer-BFGS,
`R/occu_cover_nested.R`, lognormal only), `"v2_joint"` (v2 joint Laplace).

**`group_var` (sites > cells)**: `group_var="<col>"` on icar/bym2 maps each site ->
field node, so `n_sites` > `n_cells`. Field length `n_cells` while psi/p/cover run over
`n_sites`; per-arm `spatial_idx` (field node) + `cell_obs_map` (occupancy unit)
decouple. Layout: site = cell x time-period. R-side only (`.dispatch_occu_cover`,
`.occu_cover_build_joint_arms`); joint only. `test-occu-cover-group-var.R`.

**Per-group RE on shared-field path (#56)**: one random intercept on psi —
`re(g)`/`(1|g)` — alongside field joins joint fit as one `iid` prior block (consumer of
tulpa#86). Variance integrates on outer grid (reported `sigma_re`); BLUPs in `fit$re` +
`ranef()`. `.occu_cover_spatial_fields` extracts `tobs_re`, fitter appends iid block
(`obs_idx` = group on psi rows, 0 on p/pos), `.occu_cover_jc_postprocess` extracts
`sigma_re` + BLUPs. Scope: ONE random intercept; slope/correlated/RE-without-field
rejected; v2/v3 hatches carry no RE. Per-arm community variances don't scale here (joint
engine grid-integrates every variance component), so COMMUNITY spatial occu_cover =
reduced-rank Laplace-EM `ms_occu_cover()` + `icar()` (tulpa#67, below), NOT this path.
`re.sigma.grid` knob. `test-occu-cover-field-re.R`. `fit$re` is a per-TERM flat
list keyed by arm (lone term) or `"<arm>:<var>"` (crossed); psi entry `fit$re$psi`.

**Per-group RE on detection / cover arm (#102 intercept, #103 crossed/nested)**:
random intercept(s) on `detection=`/`positive=` (`(1|g)`/`re(g)`), per-VISIT
grouping (one code per (site,visit)), composes w/ the required psi field on the
nested_laplace joint path. `.occu_cover_obs_re_parse` (occu_cover.R) strips ALL re
terms off the obs formula BEFORE copy-extraction + design build (rejects other
structured terms; copy/re allowed), returns `$terms` (LIST of specs, crossed/nested)
+ `$has_slope`; `.occu_cover_obs_re_design` resolves each term's per-(site,visit)
codes site-major from `data` (site-level, broadcast) or `visits` (visit-level) via
`.occu_cover_obs_flat_eval`, levels from observed visits only, builds slope `Z`
(intercept + covariate cols) for slope terms, attaches per-term LIST
`model$re_det`/`re_pos`. **Slopes (#103, tulpa>=0.0.39)**: NO gate -- the
DESCRIPTION floor enforces the engine. Builder (`add_re_term`, arms builder passes
keep-subset `Z` via `re_det_terms`/`re_pos_terms`) per term: intercept -> one
scalar iid block; uncorr slope (`!correlated`) -> one weighted iid block per coef
(tulpa `svc_weight` = `Z[,c]`, intercept col all-ones = scalar iid); corr slope
(`(1+x|g)`) -> one `miid` block (tulpa#114: `Q=I` mcar, `n_fields`=n_coefs,
`field_weight`=Z cols, free Sigma log-Cholesky). **Slope covariate STANDARDIZED**
to unit SD in `.occu_cover_obs_re_design` (`coef_scales`, intercept scale 1) so the
fixed Sigma grid is scale-invariant; BLUP/sigma + the predict draws (joint_substrate)
back-transformed `/scale` to natural units (cor scale-free). miid grid =
`.occu_cover_miid_logchol_grid` (p=2 principled compact: SYMMETRIC rho incl 0 +
strong +/-, log-spaced SD) to stay under the engine's 2048 outer-grid cap, knob
`control$re.logchol.grid.p`/`.pos`. `re_descs` = ONE desc per TERM
(block_start/n_blocks span). **Key fix (#102)**: det arm's `field_coef`=1 when
`model$re_det` present so the iid block scatters; field skipped by `spatial_idx=0`.
Postprocess: per term gather its blocks' latent cols (uncorr = n_coefs iid blocks;
corr = one miid, latent coef-major `(c-1)*ng+g`) -> `[n_groups x n_coefs]` BLUP
(centred per coef); sigma per coef from `b<P>.sigma` (uncorr) or marginalized from
`b<P>.L<ij>` log-Cholesky (corr, reuses `.cover_mcar_logchol_to_L` + `put_derived`)
-> per-coef `sigma` + `cor` matrix on `re_terms`. names `.occu_cover_re_sigma_names`
(base `sigma_re`/`_p`/`_pos`, `_<var>` when >1 term/arm) + `_<coef>` per coef +
`cor_re_*_<ci>_<cj>`. `fit$re` = flat list (arm/var/sigma/cor/blup
[vec if n_coefs==1 else matrix]/blup_sd/n_coefs/coef_names/covariate_names/
correlated/levels/latent_idx). predict: `.tobs_joint_draws` draws latents
(coef-major), `.occu_cover_re_offset` SUMS terms on the arm; slope coef weighted by
`newdata[[coef_name]]` (intercept=1), unseen -> 0 = pop mean. Gates: obs RE needs
`nested_laplace`; pos RE needs `cover_aggregate="none"`; not composed w/ correlated
MCAR / latent cover / batch. Corr slope / crossed grow the grid ->
`control$integration="ccd"`. sigma carries binary small-cluster attenuation (lower
bound); BLUPs recover (cor>0.5). knobs `re.sigma.grid.p`/`.pos`,
`re.logchol.grid.p`/`.pos`. `test-occu-cover-obs-re.R`.

### `ms_occu_cover()` detail

Community `occu_cover()` (`R/ms_occu_cover.R`); per-species coef RE w/ Gaussian
community covariances across psi/p/cover arms, shared dispersion. Latent z integrates
closed-form per species-cell (reuses `.occu_cover_site_ll`); per-species deviations by
in-tree pure-R Laplace-EM — arrowhead joint Newton (RE Schur-folded, analytic grads
`.occu_cover_eta_grad`) + closed-form community-cov M-step. Community-mean SEs =
marginal observed info (Louis 1982). Beta+lognormal. Non-spatial Laplace only
(structured term any arm errors+pointer). Recovery + 20-seed coverage
(`test-ms-occu-cover.R`); status `"working"` (#98). Community VARIANCE carries Laplace
small-cluster attenuation (means do not); AGHQ-debiased BY DEFAULT below
`re.aghq.maxdim` (4), above the cap EM variance = tested lower bound (tensor AGHQ
exp in total RE dim). Flagged via `print.tobs_fit` +
`fit$ms_community$var_attenuation` marker + `?ms_occu_cover` (#47). NUTS/negbin/
dispersion RE pending.

### `occu_multiscale_cover()` detail

Three-level occupancy + cover hurdle (#29; `R/occu_multiscale_cover.R`, fitter
`R/occu_multiscale_cover_joint.R`). For data where "visits" are spatially
distinct PLOTS aggregated into `(cell, period)`, not temporal revisits (EVA/MOTIVATE
vegetation; Nichols 2008, Mordecai 2011). `occu_cover()` treats plots as detection
replicates -> conflates within-cell prevalence into detection (Kendall & White 2009);
this family adds explicit middle level:

```
z_c        ~ Bernoulli(psi_c)        # cell/range occupancy
a_cj|z=1   ~ Bernoulli(theta_cj)     # plot availability/use
y_cjv|a=1  ~ Bernoulli(p_cjv)        # detection
cover|y=1  ~ f_pos(eta_pos, disp)    # hurdle (beta/lognormal)
```

Both z (cells) + a (plots) marginalize closed-form (two states each) -> exact joint
marginal LL, reuses occu_cover nested-Laplace cell-coupling machinery.

**Inputs**: `y`/`y_pos` = `[n_plots x max_visits]`. State `formula` = cell-level psi,
MUST carry areal field naming per-plot cell col: `icar(graph=adj, group_var="cell")`.
`availability=~...` = plot-level theta (default `~1`); `detection` = per-visit p;
`positive=~...` = cover. **Engine** (`method="nested_laplace"` only): 4-arm
generalization of occu_cover joint via
`tulpa_nested_laplace_joint(cell_coupling="occu_multiscale_cover_*")`. Field coupling:
psi `field_coef=1`; theta/p `0`; pos `list(name="alpha")`. Cell spec
(`src/cell_coupling_occu_multiscale_cover.{cpp,h}`, per-fit) reuses occu_cover helpers
(`src/occu_coupling_shared.h`).

**Identifiability**: theta + p separate only with replication WITHIN a plot. Single
releves -> identifies psi (cell) + product theta*p, reduces to occu_cover; surfaced
via `message()` in `.dispatch_occu_multiscale_cover` on single-releve data (#97).
**Scope** (`"working"`, #97): three engines — `nested_laplace` (shared + SVC-trend
coupled field), `laplace` + `nuts` (both non-spatial: iid cells, field fixed at 0;
cell-declaring areal term supplies plot->cell map only, graph ignored). NUTS rejects
a coupled SVC/trend field (single cell-declaring term only).
`test-occu-multiscale-cover-{recovery,coupling,nuts}.R`. `simulate_occu_multiscale_cover()`.

## NUTS coverage status

`temporal`, multi-term `re`, `latent` smoke-tested 2026-05-20
(`dev_notes/probe_blocked_nuts.R`, single-season occ) -> return `tobs_fit` w/o crash,
but gradient correctness / calibration / convergence NOT verified. Treat as "not
blocked", not "validated". Family NUTS paths (#37/#38/#39/#40/#41/#14/#67) ARE
recovery-tested.

**SVC = two distinct flavors, do NOT conflate (#118):**
- **Areal spatially-varying coefficient** (the `svcPGOcc` analogue): a WEIGHTED
  areal bar `spatial(~ 1 + w || cell, graph)` on `occu()` via `nested_laplace`,
  rerouted through the joint direct-grid engine (`.tobs_fit_occu_joint`, #81).
  RECOVERY-TESTED: `test-occu-spatial-svc-recovery.R` +
  `test-occu-svc-joint-recovery.R` recover the known intercept + trend surfaces
  (`fit$spatial_field`/`fit$trend_field`, `cor > 0.75`), the SD hyperparameters,
  and the coefficients. This arm arrives as a `spatial` term, not `svc`.
- **Continuous NNGP `svc()` term** (`svc(lon, lat, indices=)`): wired to
  single-season `occu()` NUTS ONLY (`populate_svc`, `src/populate_helpers.h` ->
  `data.svc_data`). The eta-assembly + NNGP prior + gradient live in the compiled
  UPSTREAM tulpa engine (not tulpaObs). The fitted surface IS now exposed as
  `fit$svc_field` (an n_obs vector / n_obs x n_svc matrix of posterior means,
  per-draw surface on `attr(., "draws")`), sliced by position from the layout the
  engine reports as `fit$svc_layout`; the block is named (`svc_w[i,j]`,
  `log_sigma2_svc[j]`, `log_phi_svc[j]`) instead of falling through to
  `param[k]`. Requires `prior_range = c(r0, alpha)` (PC prior on the range,
  `P(range < r0) = alpha`); tulpa ships the anchors unset and refuses without
  them, and neither package defaults them.
  **RECOVERY-VALIDATED** (#119, tulpa 0.0.82): `test-occu-svc-nngp-recovery.R`
  measures seeds 1/2/3/11 at N=150, J=6, p=0.6, truth phi=0.25 / sigma=1.3 ->
  divergences 72-83% -> **0 on every seed**, phi ~4 -> 0.14-0.23, sigma
  1.06-2.31. Two upstream causes: the Uniform-behind-a-wall range prior
  (gcol33/tulpa#144) and the SVC marginal-SD half-Cauchy being improper on its
  sampled coordinate (nothing bounded sigma above); fixing the SD prior is what
  pulled phi onto truth, the two being ends of the GP ridge. Surface cor did NOT
  move (0.73 mean vs 0.66) -- it is information-bounded at these settings, not
  sampler-bounded, so the calibration test asserts divergences + phi + sigma and
  deliberately does NOT assert on cor. tulpa cannot make this measurement:
  `svc()` is a tulpaObs term and `cpp_tulpa_fit_generic` is a plain LM, so no
  tulpa-side fit reaches the NNGP SVC path -- which is how #144 survived there.
  On any OTHER family, or `occu()` under laplace/nested_laplace, `svc()` ERRORS
  with a pointer to the areal bar (`.tobs_fit_model` guard, #118) rather than
  silently dropping. `test-svc-guard.R`.

## Performance

N=200, single-season (J=3), 2026-05-24: `laplace` (prior-aware penalized EM, default)
**~0.2s**; `laplace_gibbs` ~1.7s; `nuts` ~13s (historical). inlaocc 0.7s, spOccupancy
0.9s. Penalized EM trades speed to break psi-p identifiability ridge at small J;
Gibbs/MI add Rubin-pooled correction.

## Progress + ETA (all backends, #43)

Every fitting loop reports progress bar + ETA. ONE config surface: `tobs()` sets scoped
option `tulpa.nl_progress` from `control$progress[.every/.throttle/.file]`
(`.tobs_progress_opt`, tobs.R). Two channels, both ON by default: Rcout console bar
(`progress`, set `control$progress = FALSE` to silence; NOT tied to `verbose`) + a
heartbeat file (`progress.file`, only signal surviving detached Start-Process/nohup
stdout buffer). File wire format = one overwritten line `"<done> <total> <elapsed_s>
<eta_s>"`. Use `control[["progress"]]` (exact), never `control$progress` — `$`
prefix-matches `progress.file`.

Backend -> reporter:
- outer-grid (nested-Laplace areal, cover/occu_cover/multiscale joint, nmix spatial) ->
  C++ `tulpa_progress::GridProgress` (unit "cells").
- NUTS, ALL families -> `GridProgress` via active-pointer `g_active_grid_progress`,
  ticked once/iter under omp critical. Console auto-suppressed in across-chain parallel
  region (file is channel); byte-exact preserved (tick touches only clock/counter/file,
  never RNG). `make_nuts_progress` reads option (unit "iter").
- EM-Laplace (occu/dyn/int/jsdm), community EM (ms_occu/ms_dyn/ms_int, ms_occu_cover),
  RE-EM (em_laplace_re.R), fp_occu/dyn_abun optim -> tulpa R loop via
  `tulpa:::.tulpa_iter_progress` (R/progress_iter.R).
- count-marginal Laplace (abun/removal/distance) + community N-mixture EM (ms_abun,
  cpp_nmix_community_em) -> C++ `make_grid_progress_from_option` (nmix_progress.h).

ETA = upper bound to max_iter, finalised by `finish()` on early convergence. Test:
`test-progress-all-variants.R` (+ `test-occu-cover-progress.R`).

## File organization

Most files named inline above; non-obvious ones:

```
R/
  tobs.R / obs_families.R / occu.R   — dispatcher+print; family ctors; .tobs_build_model()
  community_em.R            — shared community Laplace-EM .tobs_community_em() (ms_occu/ms_dyn_occu/ms_int_occu/ms_count/jsdm/ms_abun-latent/ms_distance). Three OPTIONAL args, each defaulting to the previous behaviour byte-identically: `sp_info` (analytic per-species observed info; the FD fallback costs 2(P+G) marginal sweeps per species per Newton step -- supply it when the kernel exposes the Louis block), `init_b`/`init_Sigma` (warm start; a block-coordinate caller re-enters the EM once per outer pass, so a cold restart each time dominates)
  ms_{occu,dyn_occu,int_occu}.R      — community single/dynamic/integrated
  abun.R / abun_nuts.R               — nmix family + non-spatial NUTS (#41)
  count_spatial.R / count_methods.R  — areal count fitter .tobs_fit_count_spatial (#117); fitted/predict/residuals
  ms_count.R                — community count / relative-abundance GLMM (msAbund, #117); .tobs_fit_ms_count over shared community_em.R
  ms_count_nuts.R / src/ms_count_nuts.cpp — community count NUTS (msAbund NUTS, #117); in-tree C++ FullGradFn (reduced ms_abun_nuts: no detection/latent-N), R oracle .tobs_ms_count_nuts_logpost, warm-start from Laplace-EM
  ms_count_spatial.R        — community count + shared areal field (sfMsAbund) + SVC bar (svcMsAbund, #117/#118); block coordinate ascent (community EM offset <-> multi-field Poisson-ICAR), pure R
  community_latent.R        — the factor Newton (`.tobs_latent_factor_update`) is a RAW Newton (no line search, unlike `.tobs_community_em` which backtracks); `safe_step()` rejects a non-finite/singular step + holds the previous iterate (byte-identical when finite). A latent-count marginal (nmix/distance) can return non-finite curvature far from the mode; the NaN used to surface only later as a non-finite `sd()` in the rescale. Found by a 6-SEED loop, never by a single fit. SHARED latent-structure engine for EVERY community family (#119/#120/#121): one block-coordinate ascent (community EM w/ the latent as an offset <-> field / factor updates) + the areal Newton, the factor update, bym2/car_proper/spde hyper grids. A family supplies ONE callback `working(eta) -> list(score, curv)` (per-(site,species) score+curvature wrt an additive offset on the structured arm): Poisson `(y-mu, mu)`, occupancy two-state, Bernoulli `(y-psi, psi(1-psi))`. Field solve is `t(A) diag(w) A + tau Q`, so the site->node map slot takes an areal group_var incidence OR an spde barycentric projector unchanged. Adding a family to every latent route = one callback, not a new fitter
  ms_occu_field.R           — community occupancy SVC (svcMsPGOcc, #118); block coordinate ascent (community occ EM psi offset <-> two-state-marginal occupancy field solve), intercept + SVC field(s), pure R; plain intercept -> C++ ms_occu_spatial.R
  ms_abun.R / ms_abun_nuts.R         — community nmix + NUTS (#14)
  ms_abun_latent.R          — community nmix + latent() factors (lfMsNMix) / + shared field (spatial-factor); the ONLY new piece is the working oracle over the Royle marginal: score=grad_eta_lambda, curv=info_eta_lambda-var_N*score_wt_lambda^2 (the Louis (1982) (1,1) block = abundance curvature w/ the detection arm profiled out), which nmix_site_marginal() already exposes. Supplies `sp_info` (the design-sandwiched per-site Louis block) so the community EM skips its FD Hessian: 387s -> see table. Plain field w/o factors KEEPS the C++ #12 path
  ms_distance.R             — community binned distance sampling (msDS, #117) + latent() factors (lfMsDS) + shared field (sfMsDS). NO new C++: cpp_distance_site_sweep already returns log_lik/grad_lam/info_lam/var_N/swl, so the community EM reads its per-species score from it and the driver oracle is the SAME Louis formula as ms_abun (curv = info_lam - var_N*swl^2). Hazard-key log-shape = a community `global` (shared across species). simulate_ms_distance() draws through cpp_simulate_distance (the kernel the likelihood integrates against) -- a separate R-side quadrature simulates from a pi the model is not fit against and biases recovery. Poisson only; NUTS not wired
  nmix_laplace{,_re,_re_spatial,_spatial}.R — non-spatial / community / sfMsNMix / areal fitters
  nmix_re_aghq.R / nmix_site_marginal.R — grouped RE -> NMixGroupedOracle; per-site AGHQ callback
  occu_fit.R / occu_priors.R / laplace.R — .tobs_fit_model(); occu_priors()+beta_prior; .tobs_laplace()+EM cbs
  em_nested_laplace.R / simplified_laplace.R / sla_*.R — nested-Laplace EM; SLA wrapper + paths
  family_cover_hurdle.R    — .dispatch_cover() (two-Laplace hurdle), large
  occu_cover.R             — joint occu-det+cover wiring + .occu_cover_eta_from_par()
  occu_cover_nuts.R        — non-spatial occu_cover NUTS: R oracle + .tobs_fit_occu_cover_nuts
  nuts_chains.R            — multi-chain pooling + shared split-Rhat/bulk-ESS (.tobs_nuts_rhat_ess)
  ms_occu_cover.R / ms_occu_cover_spatial{,_nuts}.R — community joint; spatial-factor JSDM (tulpa#67) Laplace-EM + NUTS
  occu_multiscale_cover{,_joint}.R — 3-level occu+cover (#29) + 4-arm fitter
  occu_cover_spatial.R / occu_cover_nested.R — v2_joint / v3_nested escape hatches
  formula_terms.R / formula_parse.R  — term registry+ctors; AST parser
  inputs.R                 — single source of truth for response/site/visit input: .tobs_check_site_count() (site-count cross-check every binder used to hand-roll), .tobs_input_dims()/fit$dims canonical totals, .tobs_unpack_frame() (tobs_data -> data/y/visits)
  spatial.R / methods.R / diagnostics.R / data.R / within_between.R — precompute; S3; diags; data+sims; decomposition
  RcppExports.R            — generated, do not edit
src/
  occu_fit.cpp / populate_helpers.h  — unified C++ entry; populate_spatial/temporal/re/svc/latent
  occ_*.h / dyn_occ_*.h / integrated_occ_*.h — single/dynamic(HMC fwd)/integrated
  cell_coupling_occu_cover.h / cell_coupling_occu_multiscale_cover.{cpp,h} / occu_coupling_shared.h — coupling specs + shared helpers
  ms_occu_cover_spatial_nuts.cpp / abun_nuts.cpp / ms_abun_nuts.cpp / occu_cover_nuts.cpp — NUTS (#67/#41/#14; occu_cover non-spatial)
  nuts_engine.h            — shared run_tulpa_nuts driver for the in-tree FullGradFn targets
  community_chol.h         — shared log-Cholesky helpers (#14 non-centered, #67 centered)
  RcppExports.cpp          — generated, do not edit
  Makevars.win             — CXX_STD=CXX17, OpenMP, -Wa,-mbig-obj (large-obj MinGW)
tests/testthat/ — test files;  vignettes/ — cover-hurdle/spatial-spde/vs-inla Rmd;  dev_notes/ — probes/repros
```

## Roxygen / Rd

After exported-doc edits or `@export` changes:

```r
Rcpp::compileAttributes()   # only if src/ changed
devtools::document()
```

`devtools::document()` regenerating Rd w/ non-Latin-1 Unicode (`n⁴`, `≤`, `→`) breaks
CRAN PDF manual. Stick ASCII / Latin-1 Supplement; see global CLAUDE.md "ASCII-Only in
Roxygen/Rd" for safe/unsafe map.
