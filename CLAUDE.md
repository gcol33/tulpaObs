# CLAUDE.md

Guidance for Claude Code in this repo. Caveman speak: terse, telegraphic. File
paths + function names exact. (Compacted from a much longer original; if a
detail you need is missing, read the source file named.)

## tulpaObs — hierarchical latent-state observation models on tulpa

Bayesian occupancy / abundance / distance / removal / cover. Built on
[`tulpa`](https://github.com/gcol33/tulpa) engine. R pkg, C++17 backend
(Rcpp/RcppEigen). Needs sibling `../tulpa` checkout (`LinkingTo: tulpa`).

**Public API:** `tobs()` + family ctors: `occu()`, `dyn_occu()`, `ms_occu()`,
`ms_dyn_occu()`, `ms_int_occu()`, `int_occu()`, `jsdm()`, `abun()`, `ms_abun()`,
`dyn_abun()`, `distance()`, `removal()`, `fp_occu()`, `cover()`, `occu_cover()`,
`ms_occu_cover()`, `occu_multiscale_cover()`. S3
classes all `tobs_*` (`tobs_fit/model/family/spatial/temporal/re/svc/latent/priors_spec`).

**Structured terms live in formula** (lme4/mgcv/INLA style), NOT `tobs()` args.
Registry `R/formula_terms.R` maps name -> ctor: spatial `icar() bym2() car()
car_proper() gp() multiscale_gp() spde()`; `re()`; `temporal()`; `svc()`;
`latent()`; `copy("id")` (share one realization across state + detection).
Ctors `.tobs_term_*` + parser `.tobs_parse_formula`/`.tobs_bind_formulas`
internal. No exported `tobs_icar()` etc, no `spatial=`/`temporal=`/`re=` args.
Term's process = which process formula it sits in; fitter derives `shared =
c(occ, det)` via `.tobs_structures_from_model()`.

**lme4 bars = sugar over `re()`.** `.tobs_desugar_bars()` (`R/formula_parse.R`)
rewrites bars to `re()` on AST before `terms()`. One parser, one term type.
Forms (tulpaObs#10): `(1|g)`, `(x|g)`, `(1+x+z|g)` (multi-slope, `cbind()`
stack), `(x||g)` (uncorr), `(0+x|g)` (slope-only), crossed/nested `(1|g:h)`
`(1|g/h)`. Multi-slope + nested pure R-side (`build_re_spec()` width->`n_coefs`);
slope-only needs tulpa `re_has_intercept` flag (ABI 22).

**RE both engines** (tulpaObs#11). NUTS fits all forms. Default
`engine="laplace"` fits iid intercept RE, uncorr slopes, AND corr slopes
(`(1+x|g)`) on occ OR det arm of single-season via variance-component EM
(`R/em_laplace_re.R`, `.tobs_em_laplace_re()`): splits RE into occ/det arm by
`shared` (`.tobs_re_split_arms()`), wraps tulpa `tulpa_laplace()` in occupancy
missing-data EM, EM/REML cov update `Sigma_k <- mean_g[b_g b_g' + Cov(b_g|y)]`.
`Cov(b_g|y)` from `tulpa_laplace(return_re_cov=TRUE)$cov_blocks`: occ arm scales
by `M` (pseudo-binomial inflation), det arm weighted binomial (`M=1`). Corr term
keeps full `Sigma`; uncorr projected to diag each M-step. Gates (error, point to
NUTS) via `.validate_re_laplace()`: RE shared across BOTH arms, RE+spatial,
RE+visit-level det, RE on non-single families.

Raw EM variance components (sigma, cor) carry Laplace small-cluster bias for
binary; fixed-effect SEs do not. Default `re.aghq=TRUE` (`control=list(n.quad=)`)
debiases variance components after EM via adaptive Gauss-Hermite quad on the exact
per-group marginal (single grouping factor, one arm, RE dim<=3; else falls back to
EM). Engine `tulpa::tulpa_re_aghq()`; `R/re_aghq.R` (`.tobs_re_aghq()`) thin
wrapper supplying occ/det site marginal as `make_site`. Measured occ arm:
per-group-n=8 sigma bias ~18% (EM)->~4% (AGHQ); det arm ~70% attenuation->~1%,
88-96% coverage. Det-arm RE params named for det process (`sigma_p<t>_*`,
`re_p<t>_*`). Default LKJ (`re.lkj=1.5`) per corr block regularizes weak
correlation off +-1 (SDs untouched; `re.lkj=1` disables). RE naming + BLUP recon
both engines in `R/re_effects.R` (`ranef()`/`coef()` in `R/methods.R`); corr
off-diag `cor_<g>_<ci>_<cj>`.

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

The full suite fits 15-20 models per seed across many seeds plus NUTS/spatial
recovery -> it takes HOURS. Do NOT run it on every edit. Use this ladder:

1. **While iterating** -> run only the test file(s) covering the code you
   touched: `testthat::test_file("tests/testthat/test-occu.R")` or
   `devtools::test(filter = "occu")`. Seconds.
2. **Whole-suite smoke** (plumbing/dispatch/closed-form, no model fits) -> set
   the fast tier: `Sys.setenv(TULPAOBS_FAST = "1"); devtools::test()`. **~22s**,
   ~1330 assertions, the ~215 heavy fitting/recovery/NUTS blocks report as skips
   (never silently dropped). `skip_if_fast()` gates every model-fitting block;
   it is a no-op when the env var is unset, so the full run is unchanged.
3. **Full recovery suite** (all seeds, NUTS, spatial) -> ONLY before committing
   to main, before a release, or when explicitly asked. `Sys.unsetenv(
   "TULPAOBS_FAST"); devtools::test()`. Hours; uses `Config/testthat/parallel`.

Adding a slow test: pair `skip_if_fast()` with `skip_on_cran()` at the top of
any block whose cost is a multi-seed fit or a NUTS sample (helper
`tests/testthat/helper-speed.R`). C++ recompiles are ccache-backed, so a clean
edit-rebuild is fast; only a killed/partial build needs `pkgbuild::clean_dll()`.

## Architecture

One C++ entry `cpp_occu_fit` for NUTS. Laplace via tulpa EM+Laplace
(`tulpa::tulpa_em_laplace`). Default `method="laplace"` = penalized EM, Gaussian
marginals, no post-EM correction; fixed-effect prior attached per M-step block
as `beta_prior`. MI/Gibbs opt-in: `method="laplace_gibbs"`/`"laplace_mi"`; same
prior threads into correction refits (tulpa#27), penalised like laplace unless
`priors=FALSE`.

**nested-Laplace** engine (`em_nested_laplace.R`, `method="nested_laplace"`):
non-conjugate hyperpriors via multi-block latent prior. `.tobs_em_nested_laplace()`
thin wrapper over `.tobs_laplace(latent_prior=)`: multi-block prior from formula
`spatial`/`temporal`/`re` attached to state M-step block (`.tobs_laplace_nested()`),
routed through `tulpa_nested_laplace()`. Wired: single-season, integrated,
community, dynamic occ. Cover hurdle joint uses `tulpa_nested_laplace_joint()`.

**INLA-style NA-response prediction + calibrated CIs** (single-season): all-NA
detection sites (`.tobs_heldout_sites()`) interpolated by field;
`predict(type="state")` returns per-site psi posterior (`psi`, `psi_lower`,
`psi_upper`, 95%) marginalised over hyperparam grid. Calibrated by one
exact-marginal refine pass (`.tobs_occu_state_marginal_fit()`): integrate out z,
each site Bernoulli on D=1{>=1 detection}, mean q*sigma(eta), q=1-(1-p)^J. Fits
tulpa `family="bernoulli"` w/ per-obs prob scale `det_prob=q` (scaled-Bernoulli):
converged Hessian = marginal curvature, `fitted_eta_var` = calibrated per-cell
predictive variance. Per-row eta posterior = Gaussian mixture over cells
(`.nested_psi_mean()`, `.nested_psi_quantiles()`). Coverage ~1.0, cor ~0.88,
MAE ~0.11 on 10x10 icar/bym2 (`dev_notes/probe_nested_ci_coverage.R`). Old tulpa
w/o `fitted_eta_var` -> `psi` only, NA intervals.

**Simplified-Laplace skew correction** (`simplified_laplace.R`, `sla_*` files):
orthogonal post-fit marginal refine, `approx="simplified_laplace"` (`*_sla`
methods). Computed for single/dynamic/integrated occ + cover hurdle; no-ops to
Gaussian (records `sla_status`) for jsdm.

**Backend coverage enforced centrally**: `.tobs_family_methods` in `R/tobs.R` =
single source of truth for which `method` each family supports; `tobs()` errors
with pointer, no silent downgrade. `nested_laplace` = occu/int_occu/dyn_occu +
cover; `*_sla` on nested = occu + cover only; cover hurdle has no
NUTS/`laplace_gibbs`/`laplace_mi`. `abun` = laplace + nuts (non-spatial) +
nested_laplace (areal). `ms_abun` = laplace + nested_laplace (shared areal field).
`occu_multiscale_cover` = nested_laplace ONLY. Community occupancy
`ms_occu`/`ms_dyn_occu`/`ms_int_occu` = laplace ONLY (shared community Laplace-EM,
`R/community_em.R`; per-species coef RE, per-arm community covariance).

Observation families below (`removal`/`distance`/`fp_occu`/`dyn_abun`) = laplace
+ nuts, non-spatial Pois/NB. Each: closed-form (or exact HMM-forward, dyn_abun)
marginal over latent N, analytic gradients, in-tree FullGradFn driving tulpa's
NUTS engine (shared `src/nuts_engine.h`), draws -> WAIC/LOO; spatial/RE pending.

- `removal` (tulpaObs#39): sequential depletion, pass k sees `N - sum_{l<k} y_l`
  trials = depleting-binomial product = multinomial removal; latent N summed to
  `K_max`. SHARES the count-marginal Laplace driver (`src/marginal_count_laplace.h`)
  + NUTS (`src/marginal_count_nuts.h`) with `abun`; per-site `src/removal_kernel.h`
  over shared `accumulate_count_moments`/`fill_nb_dispersion` (`src/nmix_kernel.h`).
  Areal spatial (icar/car_proper) on the abundance arm via `nested_laplace`
  (tulpaObs#51): the removal per-site marginal (`compute_removal_site`) returns the
  same `NMixSiteResult` moments as the Royle kernel, so it reuses the family-
  agnostic nested-Laplace inner Newton + outer-grid orchestration
  (`src/nmix_count_spatial_driver.h`, templated on the site kernel; extracted from
  `nmix_spatial.cpp`, byte-identical for nmix) via `cpp_nested_laplace_removal_{icar,
  car_proper}` (`src/removal_spatial.cpp`) and the shared R packer
  `.count_spatial_pack_common` / `.count_spatial_pack_bym2_common`
  (`R/removal_spatial.R`, `.tobs_fit_removal_spatial`); one spatial unit per site,
  K_max = per-site removal total. icar / car_proper / bym2 (the templated bym2
  driver `run_count_nested_laplace_bym2` is shared with nmix too); spde/temporal +
  NUTS+spatial gated.
- `distance` (tulpaObs#38): latent N in covered region, per-bin detected counts
  multinomial over `(bin 1..B, undetected)`, `pi_b = int_bin g(x;sigma) f(x) dx`
  (half-normal / hazard-rate key, line/point transect density `f`), bin integrals
  + 1st/2nd eta-derivs by Gauss-Legendre quad. REUSES count-marginal core
  (`src/nmix_kernel.h`); det arm = site-level `log sigma` + optional scalar hazard
  shape (`src/distance_quad.h`/`src/distance_kernel.h`), own Laplace driver
  (`src/distance_laplace.cpp`) + NUTS (`src/distance_nuts.cpp`). K_max default
  `3*max(rowSums)+100`. Laplace grouped RE on the abundance arm (half-normal key,
  one grouping factor, dim<=3) via the shared count AGHQ path
  (`DistanceGroupedOracle` over `CountGroupedOracle`; one sigma row/site, so the
  half-normal theta is the count layout `[beta_lambda|beta_sigma|log_r?]`);
  hazard-key + detection-arm RE gated. Areal icar()/car_proper()/bym2() field on
  the abundance arm via `nested_laplace` (tulpaObs#51, `R/distance_spatial.R` over
  the shared areal-BFGS driver `R/areal_bfgs.R`): the distance marginal is the
  bin-multinomial (not the binomial N-mixture, so not the count-spatial driver),
  but it exposes the analytic per-site gradient (`cpp_distance_site_sweep` over
  `compute_distance_site`), so the driver's BFGS + FD-Hessian recovers the
  documented observed info `diag(info_lam,info_sig)-var_N v v'`,
  `v=(score_wt_lambda,vN_sigma)`, `vN_sigma=-p_sigma/(1-p)`. Half-normal key only;
  Pois + NB (log_r jointly estimated, as non-spatial); hazard-key spatial /
  temporal + NUTS+spatial gated.
- `fp_occu` (tulpaObs#40): Miller 2011 multistate false-positive occupancy, `y in
  {0,1,2}` (none/ambiguous/certain), certain detections only at occupied sites
  identify it. Latent z summed (2-state); 4 site-level logit arms psi/p11/p10/b
  (`fp_formula`/`b_formula`, default `~1`). Laplace = analytic-grad BFGS, vcov =
  inv of -FD-Jacobian of analytic grad (`src/fp_occu_kernel.h`, gradient only).
  NUTS `src/fp_occu_nuts.cpp`. Laplace grouped RE on the psi OR p11 arm (one
  grouping factor, dim<=3) via the pure-R `make_site` AGHQ path
  (`.tobs_fp_occu_re_aghq`, no native oracle, branches on arm): holding the other
  arms fixed makes the 2-state marginal `psi*A+(1-psi)*B` a function of one scalar
  offset per site -- psi shifts the mixture weight, p11 shifts the occupied-state
  emission `A` (uniform across a site's visits), both closed-form d1/d2. RE on
  BOTH arms at once rejected (one arm per AGHQ pass); p10/b never carry structured
  terms; NUTS samples the psi-arm intercept RE only. Det-arm RE params named
  `sigma_p<t>_*` for the p11 process. Areal icar()/car_proper() field on the
  occupancy (psi) arm via `nested_laplace` (tulpaObs#51, `R/fp_occu_spatial.R`):
  the shared areal-BFGS driver (`R/areal_bfgs.R`, `.tobs_areal_bfgs_fit`, shared
  with dyn_abun) -- BFGS over the two-state marginal (`cpp_fp_occu_total_log_lik`
  analytic gradient) + CAR prior, FD-Hessian observed info; one unit/site.
  Occupancy fields are more weakly identified than count fields (one binary site
  per node); icar/car_proper/bym2 (the shared areal-BFGS field spec
  `.areal_field_{car,bym2}`); temporal + NUTS+spatial gated.
- `dyn_abun` (tulpaObs#37): Dail-Madsen open N-mixture, `N_1~Pois(lambda)`,
  `N_t=Binom(N_{t-1},omega)+Pois(gamma)`, `Binom(N_t,p)` obs. Latent N sequence
  summed by exact HMM forward over states 0..K_max; analytic gradient by
  forward-mode diff (`src/dyn_abun_kernel.h`). 4 arms lambda/p/omega/gamma
  (`omega_formula`/`gamma_formula`, default `~1`). NUTS `src/dyn_abun_nuts.cpp`.
  K_max default `max(count)+40` (forward cost ~cubic in K). Poisson init +
  constant recruitment v1 (negbin/season-varying pending). Laplace grouped RE on
  the initial-abundance arm (one grouping factor, dim<=3): the RE shifts only
  eta_lambda, which enters only the season-1 initial distribution, so the per-site
  marginal factorises as `L(eta_lambda)=sum_n1 pi_n1(eta_lambda) c(n1)` with the
  conditional likelihood `c(n1)=P(all data|N_1=n1)` (the O(K^2 T) HMM BACKWARD
  pass, `compute_dyn_abun_init_weights`) precomputed ONCE per make_site
  (`cpp_dyn_abun_init_weights_mat`); each AGHQ node is then an O(K) dot product
  (`cpp_dyn_abun_init_loglik`). Pure-R `make_site` AGHQ path
  (`.tobs_dyn_abun_re_aghq`, no native oracle), Pois+NB; detection-arm RE gated
  (not yet wired), omega/gamma never carry structured terms; spatial pending.
  Shared pmf helpers (`da_obs_season_pmf`/`da_recruit_pmf`/`da_binom_pmf_row`)
  back both the forward gradient kernel and the backward `c` pass. Areal
  icar()/car_proper() field on the initial-abundance arm via `nested_laplace`
  (tulpaObs#51, `R/dyn_abun_spatial.R` over the shared areal-BFGS driver
  `R/areal_bfgs.R`, `.tobs_areal_bfgs_fit`, shared with fp_occu): BFGS over the
  exact forward-HMM marginal (`cpp_dyn_abun_total_log_lik` analytic gradient) +
  the CAR prior, FD-Hessian observed-info Laplace marginal integrated over
  tau[,rho]; one unit/site, Pois+NB; icar/car_proper/bym2 (shared areal-BFGS field
  spec); temporal + NUTS+spatial gated.

All filed observation-family issues shipped; no planned-status families remain.

### N-mixture abundance (`abun()`)

Royle 2004: `N_i ~ Pois(lambda_i)`, `y_ij|N_i ~ Binom(N_i, p_ij)`. NO EM — N
marginalises closed-form (sum to `K_max`), tulpa fits marginal directly, analytic
gradients + observed-Fisher. Family wiring (`R/abun.R`): binder ->
`model_type="nmix"` (site `X_lambda`, long-form `y`/`site_idx`/`X_p`, drop NA
visits), `.tobs_fit_nmix()` -> in-tree fitter (`R/nmix_laplace*.R`,
`src/nmix_*.{cpp,h}`), wrapped `build_nmix_fit()`. Pseudo-draws MVN from JOINT
(lambda,p) cov. Abundance log link: `compute_intercepts()` reads per-process
`pi$link` (default logit, log for lambda). `simulate_abun()` + `test-abun.R`.

- **Non-spatial** (`laplace`): `nmix_laplace()`. vcov = marginal observed-Fisher
  inverse (full joint lambda/p block).
- **Areal spatial** (`nested_laplace`): `icar()`/`bym2()`/`car_proper()` on
  abundance formula -> `.tobs_fit_nmix_spatial()` ->
  `nmix_laplace_{icar,bym2,car_proper}` (one unit/site). Cov grid-integrated (law
  of total cov): kernels return per-grid `cov_blocks`, wrapper `V = sum_k w_k[cov_k
  + (m_k-mbar)(m_k-mbar)']` (`.nmix_grid_vcov()`). Rank-deficient intrinsic fields
  use a sum-to-zero constraint penalty in
  `nmix_spatial_beta_cov()`/`nmix_beta_cov_bym2()` (intercept variance reflects
  constraint, not field-mean confounding).
- **Negbin** (`abun(mixture="negbin")`): kernel P/NB threads both paths. NB adds
  overdispersion `r` (`Var=lambda+lambda^2/r`). Non-spatial: `log_r` jointly
  estimated, trailing vcov coord (`nmix_dispersion` w/ delta `r_sd`). Spatial: `r`
  integrated over outer grid, `r_mean`/`r_sd` in `nmix_hyper$r`. Matches
  `unmarked::pcount(mixture="NB")`; `test-abun.R`.
- **Random effects** (`abun(...) + (1|g)` / random slope, either arm;
  tulpaObs#13): site-level grouped RE on abundance OR detection of single-species
  nmix. Grouping = non-species; RE = subset of coefs on ONE arm.
  `.tobs_fit_nmix_re()` warm-starts no-RE Laplace betas, refines via
  `.tobs_nmix_re_aghq()` (thin `make_site` over `nmix_site_marginal()`) ->
  `tulpa::tulpa_re_aghq()` adaptive GH. Det arm: per-site RE offset uniform over
  visits -> marginal = fn of one scalar shift/site. NB threads `log_r`. Gates
  (error+pointer): RE+areal spatial, RE+visit-level det, RE both arms. Native
  `NMixGroupedOracle` (`src/nmix_re_oracle.{h,cpp}`). `test-abun-re.R`.
- **NUTS** (`method="nuts"`, tulpaObs#41; `R/abun_nuts.R`, `src/abun_nuts.cpp`):
  in-tree C++ FullGradFn over the closed-form marginal, byte-exact vs R oracle,
  warm-start + Laplace metric, draws -> WAIC/LOO; non-spatial Pois/NB.

### Community / multispecies N-mixture (`ms_abun()`)

spAbundance `msNMix`: per-species nmix w/ Gaussian community hyperpriors,
`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`, `beta_p_s ~ N(mu_p, Sigma_p)`.
`N_{s,i}` integrates closed-form per species-site (same kernel as `abun()`);
per-species deviations `b_s` = RE, community priors pin `(mu_lambda, mu_p)` fixed.
`y` = 3D array `[sites x max_visits x species]` or named list of count matrices;
`species=` required.

Family wiring `R/ms_abun.R`: binder -> `model_type="ms_nmix"`,
`.tobs_ms_nmix_longform()` flattens 3D -> stacked `(y, site_idx, species_idx,
X_p)`, `.tobs_fit_ms_nmix()` -> in-tree C++ Laplace-EM `nmix_laplace_re()`. Builds
native `tulpaObs::NMixCommunityOracle` (subclass `tulpa::REGroupOracle`) via
`cpp_nmix_community_oracle()`, drives via in-tree `cpp_nmix_community_em()`
(`n_quad=1`) OR `tulpa::tulpa_re_aghq(oracle=)` for `n_quad>1` variance debias.
Fit: per-species coef mode-find (complete-data Fisher), closed-form EM cov M-step
`Sigma_k=mean_s[b_s b_s'+Cov(b_s|y)]`, fixed-effect SEs from marginal observed-info
Schur complement of b-block (Louis 1982). Per-site kernel `src/nmix_kernel.h`
caches `lgamma` (`compute_nmix_site_cached`); `compute_nmix_site()` Poisson
delegates -> single-species/spatial byte-identical.

`coef()` = community means; `ranef()` = per-species deviations; `vcov()`/`confint()`
= community-mean cov; `fitted()`/`simulate()` = per-species lambda/p/counts.
`simulate_ms_abun()` + `test-ms-abun.R` cover community-mean recovery, 95% coverage
(20 seeds), per-species coef recovery, S3.

- **Poisson + negbin.** `ms_abun(mixture="negbin")`: per-species dispersion RE
  `log_r_s ~ N(mu_log_r, sigma_log_r)` (native oracle widens per-species RE vector
  w/ trailing `log_r_s` coord, community `mu_log_r` joins means), joint_grad / AGHQ
  (`n_quad=5` default). `coef()`/`vcov()` report `mu_log_r`; `ms_dispersion` carries
  `r`/`r_s`/`sigma_log_r`. `test-ms-abun.R` covers NB recovery.
- **NUTS** (`method="nuts"`, tulpaObs#14; `R/ms_abun_nuts.R`, `src/ms_abun_nuts.cpp`):
  samples the EXACT joint posterior (means, per-species deviations, covariances)
  over the closed-form Royle marginal -> calibrated intervals + per-(species,site)
  WAIC/LOO (`.tobs_ploglik_ms_nmix`). NON-CENTERED: per-species block carries
  whitened `z_s ~ N(0,I)`, `b_{s,arm}=C_arm z_{s,arm}` (`C_arm`=log-Cholesky), so
  covariance enters only the data term -- breaks the centered b<->Sigma funnel that
  saturated treedepth (322s->98s with reparam + OpenMP at S=6). In-tree C++
  FullGradFn parallelises the per-species loop (OpenMP, deterministic reduction ->
  byte-exact vs R oracle `.tobs_ms_abun_nuts_logpost`); chol gradient
  `chol_data_grad_noncentered` (`community_chol.h`, shared w/ #67). Data-driven
  K_max. Warm-started at Laplace-EM mode + diagonal metric. Pois + NB, non-spatial.
  `test-ms-abun-nuts.R`.

#### Areal-spatial community N-mixture (`ms_abun()` + shared field; sfMsNMix)

Shared ICAR/BYM2/proper-CAR field on abundance via `nested_laplace`
(tulpaObs#12): term on abundance formula -> `.tobs_fit_ms_nmix_spatial()` ->
in-tree nested Laplace-EM `nmix_community_laplace_{icar,bym2,car_proper}()`
(`R/nmix_laplace_re_spatial.R`) over `cpp_nmix_community_spatial_*`
(`src/nmix_community_spatial.cpp`). Model `log lambda_{s,i} = X_lambda_i.(mu_lambda
+ b_lambda_s) + f_{u(i)}`, one `f` shared across species. Top block
`(mu_lambda, mu_p, f)` shares single-species spatial layout (`field_start =
p_lambda+p_p`), so `nmix_spatial_kernel*.h` apply. Per outer grid point EM
iterates joint `(mu,f,{b_s})` mode-find (block-elim Newton, `b_s` Schur-folded) +
closed-form `Sigma` M-step; R grid-integrates means/field/Sigma. Community means
FLAT (no ridge). NB `r` integrated over outer grid (`ms_dispersion`/`ms_hyper`).
`fit$spatial_field` = posterior-mean field. `test-ms-abun-spatial.R`. NUTS pending.

### Boundary: tulpaObs vs tulpa

- **tulpaObs owns**: family likelihoods (`src/*_likelihood.h`,
  `src/nmix_*.{cpp,h}`), E-step weights, M-step encoding, family S3/diagnostics.
  RE AGHQ: occ/det per-site marginal (`make_site` in `R/re_aghq.R`), native
  `NMixGroupedOracle` (`src/nmix_re_oracle.{h,cpp}`, `R/nmix_re_aghq.R`),
  `NMixCommunityOracle` (`src/nmix_community_oracle.{h,cpp}`); both exported as
  `XPtr<tulpa::REGroupOracle>` via `tulpa::tulpa_re_aghq(oracle=)`.
- **tulpa owns**: EM engine, MI/Gibbs, Rubin pooling, AGHQ engine
  (`tulpa_re_aghq()`: quad grid / log-Cholesky / LKJ / marginal Hessian), generic
  S3/diagnostics, NUTS/HMC.

NUTS crash for component w/ correct `populate_*` here = bug in tulpa
`hmc_sampler.cpp`, not here. File against `gcol33/tulpa`.

### Key design rules

- **Never pass `Rcpp::Nullable<T>` to header helpers** — MinGW crashes. Unwrap in
  `.cpp`, pass concrete types into headers.
- **Composition over registry** — families x spatial/temporal/re/svc/latent
  orthogonal, no per-combination branches.
- **tulpaObs defines likelihoods, tulpa handles structure.**
- **Dotted arg names (spOccupancy/base-R), never underscores.** Prefer single
  word (`visits`, `priors`, `control`); compound separated `.` not `_`:
  `n.iter`, `n.warmup`, `n.chains`, `n.thin`, `n.threads`, `adapt.delta`,
  `max.treedepth`, `max.iter`, `sigma.beta`. Governs `tobs()` args + every
  `control=list(...)` key (splatted into `.tobs_fit_model()`). Internal helper
  formals may stay underscore; anything user types = dotted.

## What works (tested)

| Feature | Laplace | NUTS | Notes |
|---|---|---|---|
| Single-season occupancy | Yes | Yes | parity w/ inlaocc |
| Dynamic (HMM) | Yes | Yes | colonization/extinction |
| Community single-season (`ms_occu`) | Yes | — | independent per-arm community RE via shared community Laplace-EM (`R/community_em.R`, `R/ms_occu.R`); spOccupancy `msPGOcc`. Legacy generic-engine community path removed (tulpaObs#30); NUTS/spatial deferred |
| Community dynamic (`ms_dyn_occu`) | Yes | — | per-species psi1/p RE + shared gamma/eps; HMM-forward marginal; `R/ms_dyn_occu.R`; recovery `test-ms-dyn-occu.R` |
| Community integrated (`ms_int_occu`) | Yes | — | per-species psi + per-source detection RE; multi-source two-state marginal (analytic grad); `R/ms_int_occu.R`; recovery `test-ms-int-occu.R` |
| Integrated multi-source | Yes | Yes | shared psi |
| JSDM | Yes | — | no detection |
| Cover hurdle (joint) | Yes | — | `family_cover_hurdle.R`, `sla_cover_*`. `control$aggregate.occ` (ON, #48) collapses occurrence arm to exact Binomial suff-stat; `control$aggregate.pos` (ON for beta arm, #49) collapses beta positive arm to grouped `(n, sum log y, sum log(1-y))` (tulpa beta spec `slog_y`/`slog_1my`); explicit `aggregate.pos=TRUE` errors on non-beta arm. Both byte-identical to per-plot; `.cover_aggregate_occ`/`.cover_aggregate_pos`, scattered via `.cover_arm_keys_from_blocks`/`.cover_scatter_arm_keys` |
| Joint occu + cover | Yes | — | `occu_cover()` — see below |
| Community joint occu + cover | Yes | — | `ms_occu_cover()` — see below |
| Spatial-factor community occu + cover (JSDM) | Yes | Yes | `ms_occu_cover()` + `icar()`/`car_proper()`/`bym2()` shared field, per-species loadings (tulpa#67). Laplace-EM (`R/ms_occu_cover_spatial.R`) + NUTS (`src/ms_occu_cover_spatial_nuts.cpp`). Cover-arm factor, auto-K, `tobs_associations()`, per-species `predict()` maps |
| Multiscale occu + cover | n-L | — | `occu_multiscale_cover()` — 3-level cell/plot/visit + cover; spatial joint only — see below |
| N-mixture (Pois/NB) | Yes | Yes | `abun(mixture=)`; in-tree `nmix_laplace`, joint-vcov draws, calibrated CIs. NB jointly-est `log_r`. NUTS (#41, `src/abun_nuts.cpp`). `test-abun.R`. See below |
| N-mixture + areal spatial | n-L | — | `abun()`+icar/bym2/car_proper, `nested_laplace`; Pois/NB (r grid-int); grid-int cov (constrained intercept) |
| Community N-mixture | Yes | Yes | `ms_abun()` (msNMix); per-species coef RE, in-tree C++ Laplace-EM (`nmix_laplace_re`) -> `NMixCommunityOracle` via AGHQ, Schur SEs; Pois + negbin. NUTS (#14, `src/ms_abun_nuts.cpp`). `test-ms-abun.R`/`-nuts.R`. See below |
| Community N-mixture + areal spatial | n-L | — | `ms_abun()`+icar/bym2/car_proper, `nested_laplace` (sfMsNMix; #12); shared field on log lambda + per-species RE; nested Laplace-EM (`nmix_community_spatial.cpp`); Pois/NB; `test-ms-abun-spatial.R`. NUTS pending |
| N-mixture + grouped RE | Yes | — | `abun()`+`(1\|g)`/`(x\|g)` either arm (tulpaObs#13); non-species grouping; Pois/NB; AGHQ via `NMixGroupedOracle`. Gated: RE+spatial, RE+visit-det, RE both arms |
| Removal sampling (Pois/NB) | Yes | Yes | `removal()` (#39); `R/removal{,_nuts,_spatial}.R`. See Architecture. `test-removal.R`. NUTS samples a single intercept RE (abundance OR detection arm, #51). Laplace fits a site-level grouped RE on one arm via the shared count-model AGHQ path (`RemovalGroupedOracle`). Areal icar()/car_proper()/bym2() field on the abundance arm via `nested_laplace` (#51), reusing the templated count-spatial driver (`nmix_count_spatial_driver.h`); spde/temporal + NUTS+spatial gated |
| Distance sampling (Pois/NB) | Yes | Yes | `distance(key=, transect=, cutpoints=)` (#38); `formula`=log lambda, `detection`=log sigma, `y`=`n_sites x n_bins`. See Architecture. `test-distance.R`. NUTS samples a single abundance-arm intercept RE (#51). Laplace fits a site-level grouped RE on the abundance arm (half-normal key, dim<=3, one grouping factor) via the shared count-model AGHQ path (`DistanceGroupedOracle` over `CountGroupedOracle`); hazard-key/detection-arm RE gated. Areal icar()/car_proper()/bym2() field on the abundance arm via `nested_laplace` (#51, `R/distance_spatial.R` over the shared `R/areal_bfgs.R` driver, `cpp_distance_site_sweep`); hazard-spatial + NUTS+spatial gated |
| False-positive occupancy (multistate) | Yes | Yes | `fp_occu()` (#40); `R/fp_occu{,_nuts}.R`. See Architecture. `test-fp_occu.R`. NUTS samples a single occupancy (psi)-arm intercept RE (#51). Laplace fits a site-level grouped RE on the psi OR p11 (detection) arm via the pure-R `make_site` AGHQ path (branches on arm); both-arms-at-once rejected. Areal icar()/car_proper() field on the psi arm via `nested_laplace` (#51, `R/fp_occu_spatial.R` over the shared `R/areal_bfgs.R` driver, icar/car_proper/bym2); NUTS+spatial gated |
| Open N-mixture (Dail-Madsen) | Yes | Yes | `dyn_abun()` (#37); y is 3D `[n_sites x J x T]`. See Architecture. `test-dyn_abun.R`. NUTS samples a single initial-abundance intercept RE (#51). Laplace fits a site-level grouped RE on the initial-abundance arm (one grouping factor, dim<=3) via the backward-`c` precompute + `make_site` AGHQ path; detection-arm RE gated. Areal icar()/car_proper() field on the initial-abundance arm via `nested_laplace` (#51, `R/dyn_abun_spatial.R` over the shared `R/areal_bfgs.R` driver, icar/car_proper/bym2, BFGS + FD-Hessian); NUTS+spatial gated |
| Spatial ICAR/BYM2/NNGP | — | Yes | |
| Spatial + dynamic | — | Yes | |
| Nested-Laplace (areal) | n-L | — | `nested_laplace`: icar/bym2/car (+temporal/iid) on occu/int_occu/dyn_occu |
| NA-response prediction | n-L | — | `predict(type="state")`: all-NA single-season sites field-interpolated, calibrated 95% `psi_lower`/`psi_upper` from exact-marginal bernoulli pass (coverage ~1.0) |
| Formula RE (intercept) | Yes | Yes | `(1\|g)`; variance-component EM, occ OR det arm (tulpaObs#11) |
| Formula RE (uncorr slope) | Yes | Yes | `(x\|\|g)`, `(0+x\|g)`, `(1+x\|\|g)`; either arm |
| Formula RE (corr slope) | Yes | Yes | `(1+x\|g)`; EM M-step consumes tulpa `cov_blocks`; either arm (tulpaObs#11) |
| Formula RE on detection | Yes | Yes | `detection=~(1\|g)`; separate det-arm RE block, AGHQ branches on arm |
| All S3 methods | Yes | Yes | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, `$.tobs_fit` |
| Diagnostics | Yes | Yes | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatialRange |
| Simulation | Yes | Yes | `simulate_occu/ms_occu/dyn_occu/int_occu/dyn_ms_occu/int_ms_occu/cover/cover_joint` |
| Spatial prediction | — | Yes | `tobs_predict_spatial` (IDW on field) |
| Components | Yes | Yes | `tobs_re`, `tobs_temporal`, `tobs_svc`, `tobs_latent`, `tobs_community_re`, `tobs_areal` |

### `occu_cover()` detail

Joint occu-detection + cover hurdle. **Non-spatial** (`laplace`,
`R/occu_cover.R`): cell psi + per-visit detection + per-visit cover (beta or
lognormal) on exact two-state marginal. **Spatial default** (`nested_laplace`,
`R/occu_cover_joint_coupled.R`): `joint_coupled` engine via
`tulpa_nested_laplace_joint()` w/ `occu_cover_{lognormal,beta}` cell-coupling
spec (tulpa#32) — 3-arm joint nested-Laplace, outer-grid over `(sigma, alpha)`,
per-cell occupancy mixture closed-form derivs drive inner Newton. ~150-300x faster
than v3_nested at N=100, completes N=200+. Recovery: 10-seed lognormal + beta
(`test-occu-cover-joint-coupled.R`); status `"experimental"`.

**Cover-arm intercept prior (tulpaObs#32)**: on the shared-field path the cover
intercept confounds with the field level over detected cells.
`.occu_cover_coupled_arm_priors()` hands the pos arm the `cover_priors()`
weakly-informative intercept prior **by default** (not the engine's flat ridge);
else the cover intercept floats to huge SD and `predict()`'s conditional cover
blows up via Jensen. `priors = FALSE`/`"none"` disables all three arms.

**Cell-aggregated cover (`cover_aggregate`, tulpaObs#33)**: per-visit cover gives
the cover arm one row per visit, so a cell with many detected plots drives the
field more than its single occupancy obs. `cover_aggregate = "mean"` (default on
spatial path) / `"median"` collapses cover arm to ONE row per occupancy unit;
`"none"` keeps per-visit. Wired ONLY on spatial `joint_coupled` path. Needs a
cell-level positive design (from `data`); visit-level `positive` keeps per-visit
arm. C++ compile-time `Aggregated` flag on `OccuCoverCoupling`
(`src/cell_coupling_occu_cover.h`), registered `occu_cover_{lognormal,beta}_agg`;
R `.occu_cover_build_joint_coupled_arms(cover_aggregate=)`.
`test-occu-cover-coupling.R`, `test-occu-cover-aggregate.R`.

**Latent cover-per-unit (`cover_aggregate = "latent"`)**: principled mean/median
alternative — cover arm carries a per-unit cover RE `u_i ~ N(0, sigma_u^2)` shared
across the unit's detected visits, integrated out per unit (keeps every visit).
Unit-level predictor -> per-unit marginal `log M_i` SCALAR in one eta, reusing the
one-row-per-unit layout, no within-arm Hessian coupling. Lognormal = closed form
(`src/occu_cover_latent.h::LognormalLatent`); beta = adaptive GH over `u_i` reusing
`BetaPositive` (`BetaLatent`, `control$n.quad` default 15). Within-unit dispersion
pre-fit + held FIXED; `sigma_u` rides pos arm `phi_grid` (`control$sigma.u.grid`),
reported `phi_pos`. Stateful spec (`OccuCoverLatentCoupling<PosLatent>`,
`src/cell_coupling_occu_cover_latent.h`), (re)registered per fit via
`cpp_register_occu_cover_{lognormal,beta}_latent_coupling()`. Shared det-branch in
`occu_det_psi_p_block`/`occu_nodet_block` (`src/occu_coupling_shared.h`). Same
gates as aggregation. `test-occu-cover-latent.R`.

**Coupled SVC/trend fields** (tulpaObs#15): extra shared areal fields = WEIGHTED
areal terms in psi formula — `icar(graph=adj, weight=year)` couples a
spatially-varying coef on `year` atop unweighted intercept field. N fields compose
(each own outer-grid `alpha`) via tulpa multi-block copy. Intercept field =
`fit$spatial_field`; weighted fields `fit$trend_field`/`fit$trend_fields`, own
scale `alpha_trend` (`control$alpha.grid.trend`). p arm excluded via `field_coef=0`
(NOT `svc_weight=0`). Resolved by `.occu_cover_spatial_fields()`;
`control=list(trend=list(weight="<col>"))` = back-compat alias. Off-path errors via
`.tobs_reject_weighted_spatial()`. `test-occu-cover-trend.R`.

**Escape hatches**: `control$engine="v3_nested"` (pure-R outer-BFGS,
`R/occu_cover_nested.R`, lognormal only), `"v2_joint"` (v2 joint Laplace).

**`group_var` (sites > cells)**: `group_var="<col>"` on `icar()`/`bym2()` maps each
site -> field node, so `n_sites` > `n_cells`. Field length `n_cells` while
psi/p/cover run over `n_sites`; per-arm `spatial_idx` (field node) + `cell_obs_map`
(occupancy unit) decouple. Motivating layout: site = cell x time-period. R-side only
(`.dispatch_occu_cover`, `.occu_cover_build_joint_coupled_arms`); joint_coupled only.
`test-occu-cover-group-var.R`.

### `ms_occu_cover()` detail

Community version of `occu_cover()` (`R/ms_occu_cover.R`); per-species coef RE w/
Gaussian community covariances across psi/p/cover arms, shared dispersion. Latent
z integrates closed-form per species-cell (reuses `.occu_cover_site_ll`);
per-species deviations by in-tree pure-R Laplace-EM — arrowhead joint Newton (RE
Schur-folded, analytic grads `.occu_cover_eta_grad`) + closed-form community-cov
M-step. Community-mean SEs = marginal observed info (Louis 1982). Beta+lognormal.
Non-spatial Laplace only (structured term any arm errors+pointer). Recovery +
15-seed coverage (`test-ms-occu-cover.R`); status `"experimental"`. Community
VARIANCE carries Laplace small-cluster attenuation (means do not); flagged via
`print.tobs_fit` + `fit$ms_community$var_attenuation` marker + `?ms_occu_cover`
(tulpaObs#47). NUTS/negbin/dispersion RE/AGHQ debias pending.

### `occu_multiscale_cover()` detail

Three-level occupancy + cover hurdle (tulpaObs#29; `R/occu_multiscale_cover.R`,
fitter `R/occu_multiscale_cover_joint_coupled.R`). For data where "visits" are
spatially distinct PLOTS aggregated into `(cell, period)`, not temporal revisits
(EVA/MOTIVATE vegetation; Nichols 2008, Mordecai 2011). `occu_cover()` treats
plots as detection replicates -> conflates within-cell prevalence into detection
(Kendall & White 2009); this family adds explicit middle level:

```
z_c        ~ Bernoulli(psi_c)        # cell/range occupancy
a_cj|z=1   ~ Bernoulli(theta_cj)     # plot availability/use
y_cjv|a=1  ~ Bernoulli(p_cjv)        # detection
cover|y=1  ~ f_pos(eta_pos, disp)    # hurdle (beta/lognormal)
```

Both z (cells) + a (plots) marginalize closed-form (two states each) -> exact
joint marginal LL, reuses occu_cover nested-Laplace cell-coupling machinery.

**Inputs**: `y`/`y_pos` = `[n_plots x max_visits]` matrices. State `formula` =
cell-level psi, MUST carry the areal field naming the per-plot cell col:
`icar(graph=adj, group_var="cell")`. `availability=~...` = plot-level theta
(default `~1`); `detection` = per-visit p; `positive=~...` = cover.

**Engine** (`method="nested_laplace"` only): 4-arm generalization of occu_cover
joint_coupled via `tulpa_nested_laplace_joint(cell_coupling="occu_multiscale_cover_*")`.
Field coupling: psi `field_coef=1`; theta/p `0`; pos `list(name="alpha")`. Cell
spec (`src/cell_coupling_occu_multiscale_cover.{cpp,h}`, per-fit) reuses occu_cover
helpers (`src/occu_coupling_shared.h`).

**Identifiability**: theta + p separate only with replication WITHIN a plot.
Single releves -> identifies psi (cell) + product theta*p, reduces to occu_cover.

**Scope** (status `"experimental"`): spatial joint nested-Laplace only (SVC trend
+ non-spatial not wired; `.dispatch_occu_multiscale_cover` rejects laplace +
non-spatial state formula). `test-occu-multiscale-cover-recovery.R`, `-coupling.R`.
`simulate_occu_multiscale_cover()`.

## NUTS coverage status

`temporal`, multi-term `re`, `svc`, `latent` smoke-tested 2026-05-20
(`dev_notes/probe_blocked_nuts.R`, single-season occ) -> return `tobs_fit` w/o
crash, but gradient correctness / calibration / convergence NOT verified. Treat as
"not blocked", not "validated". Family NUTS paths (#37/#38/#39/#40/#41/#14/#67) ARE
recovery-tested.

## Performance

N=200, single-season (J=3), 2026-05-24:
- `laplace` (prior-aware penalized EM, default): **~0.2s**;
  `laplace_gibbs` ~1.7s; `nuts` ~13s (historical). inlaocc 0.7s, spOccupancy 0.9s.

Penalized EM trades speed to break the psi-p identifiability ridge at small J;
Gibbs/MI add Rubin-pooled correction.

## Progress + ETA (all backends, gcol33/tulpaObs#43)

Every fitting loop reports progress bar + ETA. ONE config surface: `tobs()` sets
scoped option `tulpa.nl_progress` from `control$progress[.every/.throttle/.file]`
(`.tobs_progress_opt`, tobs.R). Two channels, both ON by default: Rcout console bar
(`progress`, set `control$progress = FALSE` to silence; NOT tied to `verbose`) and
a heartbeat file (`progress.file`, the only signal surviving a detached
Start-Process/nohup stdout buffer). File wire format = one overwritten line
`"<done> <total> <elapsed_s> <eta_s>"`. Use `control[["progress"]]` (exact), never
`control$progress` -- `$` prefix-matches `progress.file`.

Backend -> reporter:
- outer-grid (nested-Laplace areal, cover/occu_cover/multiscale joint, nmix
  spatial) -> C++ `tulpa_progress::GridProgress` (unit "cells").
- NUTS, ALL families (`run_hmc_chain_cpp` + parallel core) -> `GridProgress` via
  active-pointer `g_active_grid_progress`, ticked once/iter under omp critical.
  Console auto-suppressed in across-chain parallel region (file is the channel);
  byte-exact preserved (tick touches only clock/counter/file, never RNG).
  `make_nuts_progress` reads the option (unit "iter").
- EM-Laplace (occu/dyn/int/jsdm) -> tulpa `tulpa_em_laplace` R loop via
  `tulpa:::.tulpa_iter_progress` (R/progress_iter.R).
- community EM (ms_occu/ms_dyn/ms_int, ms_occu_cover) + RE-EM (em_laplace_re.R) +
  fp_occu/dyn_abun optim -> the same R reporter.
- count-marginal Laplace (abun/removal/distance) + community N-mixture EM
  (ms_abun, cpp_nmix_community_em) -> C++ `make_grid_progress_from_option`
  (nmix_progress.h, unit "iter").

ETA = upper bound to max_iter, finalised by `finish()` on early convergence.
Test: `test-progress-all-variants.R` (+ `test-occu-cover-progress.R`).

## File organization

```
R/
  tobs.R                    — tobs() dispatcher + print.tobs_fit
  obs_families.R            — family ctors (occu, dyn_occu, int_occu, jsdm, …, cover)
  occu.R                    — .tobs_build_model() (single/dynamic/integrated/jsdm/nmix)
  community_em.R            — shared community Laplace-EM .tobs_community_em() (arrowhead Newton + per-arm cov M-step + marginal info); drives ms_occu/ms_dyn_occu/ms_int_occu
  ms_occu.R                 — community single-season (msPGOcc): build, fit (reuses ms_int kernel), S3, .tobs_richness_ms_occu, ms_occu()
  ms_dyn_occu.R             — community dynamic: HMM-forward marginal, fit (psi1/p RE + gamma/eps global), S3, simulate, ms_dyn_occu()
  ms_int_occu.R             — community integrated: multi-source two-state marginal + analytic grad, fit, S3, simulate, ms_int_occu()
  abun.R                    — nmix family: build_abun, fit_nmix, build_nmix_fit, S3, simulate_abun
  abun_nuts.R               — non-spatial nmix NUTS (#41): .tobs_abun_nuts_logpost, .tobs_fit_abun_nuts
  ms_abun.R                 — community nmix: build_ms_abun, ms_nmix_longform, fit_ms_nmix (-> nmix_laplace_re), S3, simulate_ms_abun
  ms_abun_nuts.R            — community nmix NUTS (#14): .tobs_ms_abun_nuts_logpost, layout/marginal/metric helpers, .tobs_fit_ms_abun_nuts
  nmix_laplace.R            — in-tree non-spatial nmix (Royle 2004) Laplace (Pois+NB)
  nmix_laplace_re.R         — in-tree community nmix (msNMix), .nmix_re_oracle()
  nmix_laplace_re_spatial.R — spatial community nmix (sfMsNMix): _icar/_bym2/_car_proper over cpp_nmix_community_spatial_*
  nmix_re_aghq.R            — single-species grouped RE: .tobs_nmix_re_aghq() -> NMixGroupedOracle
  nmix_laplace_spatial.R    — nmix_laplace_icar/_bym2/_car_proper areal fitters
  nmix_site_marginal.R      — per-site marginal as composable AGHQ callback
  occu_fit.R                — .tobs_fit_model() (Laplace default, NUTS fallback; nmix -> fit_nmix)
  occu_priors.R             — occu_priors() + per-block beta_prior plumbing
  laplace.R                 — .tobs_laplace() + EM callbacks per model type
  em_nested_laplace.R       — nested Laplace EM for non-conjugate hyperpriors
  simplified_laplace.R      — SLA wrapper for sla_*.R
  sla_dyn_occu.R / sla_int_occu.R / sla_cover_hurdle.R / sla_cover_hurdle_joint.R — SLA paths
  family_cover_hurdle.R     — .dispatch_cover() (two-Laplace hurdle), large
  occu_cover.R              — joint occu-det + cover hurdle: wiring + shared .occu_cover_eta_from_par()
  ms_occu_cover.R           — community joint: build_ms_occu_cover, per-species kernel, Laplace-EM fit, S3, simulate
  ms_occu_cover_spatial.R   — reduced-rank spatial-factor community occu_cover (JSDM, tulpa#67): Laplace-EM, constrained-loading param, associations/maps, S3
  ms_occu_cover_spatial_nuts.R — R NUTS target (.ms_ocs_joint_logpost), chol helpers, .tobs_fit_ms_occu_cover_spatial_nuts
  occu_multiscale_cover.R   — 3-level occu+cover (#29): builder, dispatcher, simulate_occu_multiscale_cover
  occu_multiscale_cover_joint_coupled.R — 4-arm joint nested-Laplace fitter
  occu_cover_spatial.R      — v2 joint-Laplace (escape hatch control$engine="v2_joint")
  occu_cover_nested.R       — v3 nested-Laplace (control$engine="v3_nested", lognormal only)
  sim_cover_hurdle.R        — cover hurdle simulators (incl joint)
  formula_terms.R           — term registry + ctors (.tobs_term_*), print, .tobs_term_to_tulpa_spatial
  formula_parse.R           — AST parser: parse_formula/parse_processes/resolve_terms/bind_formulas
  spatial.R                 — precompute (adjacency_to_csr, compute_bym2_scale, compute_nngp_neighbors)
  methods.R                 — S3 on tobs_fit, $.tobs_fit, predict_spatial, checkIdentifiability
  diagnostics.R             — tobs_waic/ppc/test_*/pit_residuals
  data.R                    — tobs_format/data, summary/plot.tobs_data, simulators
  within_between.R          — within_between() decomposition
  RcppExports.R             — generated, do not edit
src/
  occu_fit.cpp                — unified C++ entry
  populate_helpers.h          — populate_spatial/temporal/re/svc/latent
  occ_data.h / occ_likelihood.h — single-season occupancy
  dyn_occ_data.h / dyn_occ_likelihood.h — HMC forward algorithm
  integrated_occ_data.h / integrated_occ_likelihood.h — integrated multi-source
  jsdm_likelihood.h           — JSDM (Bernoulli, no detection)
  cell_coupling_occu_cover.h  — occu_cover joint cell-coupling spec
  cell_coupling_occu_multiscale_cover.{cpp,h} — 4-arm multiscale cell-coupling spec
  occu_coupling_shared.h      — shared coupling helpers (CSR/field-demean/inner-vcov/rmvn); reused by spatial-factor NUTS marginal (nodet_mixture_block, Lognormal/BetaPositive)
  ms_occu_cover_spatial_nuts.cpp — spatial-factor community occu_cover NUTS (#67): marginal LL + full joint log-post gradient (FullGradFn), field R(h) layer, cpp_ms_ocs_nuts
  abun_nuts.cpp               — non-spatial nmix NUTS (#41): abun_nuts_eval + FullGradFn; byte-exact vs R oracle
  ms_abun_nuts.cpp            — community nmix NUTS (#14): ms_abun_nuts_eval (non-centered b=Cz log-post + gradient, OpenMP species-parallel, lgamma cache), cpp_ms_abun_nuts; byte-exact vs R oracle
  community_chol.h            — shared C++ log-Cholesky helpers for community NUTS targets (#14 non-centered, #67 centered)
  RcppExports.cpp             — generated, do not edit
  Makevars.win                — CXX_STD=CXX17, OpenMP, -Wa,-mbig-obj (large-obj MinGW)
tests/testthat/  — test files
vignettes/       — cover-hurdle-*.Rmd, occupancy-spatial-spde.Rmd, occupancy-vs-inla.Rmd
dev_notes/       — probes, repros, notes (kept in-tree)
```

## Roxygen / Rd

After exported-doc edits or `@export` changes:

```r
Rcpp::compileAttributes()   # only if src/ changed
devtools::document()
```

`devtools::document()` regenerating Rd w/ non-Latin-1 Unicode (`n⁴`, `≤`, `→`)
breaks CRAN PDF manual. Stick ASCII / Latin-1 Supplement; see global CLAUDE.md
"ASCII-Only in Roxygen/Rd" for safe/unsafe map.
