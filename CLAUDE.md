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
missing-data EM (arm RE mode fed to psi/p at E-step), EM/REML cov update
`Sigma_k <- mean_g[b_g b_g' + Cov(b_g|y)]`. `Cov(b_g|y)` from
`tulpa_laplace(return_re_cov=TRUE)$cov_blocks`: occ arm scales by `M`
(pseudo-binomial inflation), det arm genuine weighted binomial (`M=1`). Corr
term keeps full `Sigma`; uncorr projected to diag each M-step. Gates (error,
point to NUTS, not silent drop) via `.validate_re_laplace()`: RE shared across
BOTH arms, RE+spatial, RE+visit-level det, RE on non-single families.

Raw EM variance components (sigma, cor) carry Laplace small-cluster bias for
binary (glmer nAGQ=1 regime, attenuated at small per-group n); fixed-effect SEs
(natural scale) do not. Default `re.aghq=TRUE` (`control=list(n.quad=)`)
debiases variance components after EM via adaptive Gauss-Hermite quad on exact
per-group marginal (nAGQ>1 fix; single grouping factor, one arm, RE dim<=3;
falls back to EM else, incl. split-across-arms). Engine
`tulpa::tulpa_re_aghq()` (quad grid, per-group mode, log-Cholesky cov, LKJ
penalty, joint opt, marginal Hessian); `R/re_aghq.R` (`.tobs_re_aghq()`) thin
wrapper supplying occ/det site marginal as `make_site` callback. Marginal
branches on arm: occ-arm moves `psi`, det-arm moves `p` (binomial-in-p, FD
verified). Measured occ arm: per-group-n=8 sigma bias ~18% (EM)->~4% (AGHQ),
matches NUTS; corr sigmas ~1%. Det arm: EM attenuates sigma ~70% (only occupied
sites inform p)->AGHQ ~1%, 88-96% fixed-effect CI coverage. Det-arm RE params
named for det process (`sigma_p<t>_*`, `re_p<t>_*`). Default LKJ
(`re.lkj=1.5`) penalty per corr block regularizes weak correlation off +-1
boundary toward 0 (SDs untouched; `re.lkj=1` disables); NUTS for full posterior.
RE naming + per-group BLUP recon both engines in `R/re_effects.R`
(`ranef()`/`coef()` in `R/methods.R`); corr off-diag reported `cor_<g>_<ci>_<cj>`.

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
non-conjugate hyperpriors via multi-block latent prior.
`.tobs_em_nested_laplace()` thin wrapper over `.tobs_laplace(latent_prior=)`:
multi-block prior from formula `spatial`/`temporal`/`re` attached to state
("occ") M-step block (`.tobs_laplace_nested()`), routed through
`tulpa_nested_laplace()`. Wired: single-season, integrated, community, dynamic
occ (shared `build_*_callbacks`; `.tobs_state_block_dims()` maps state row ->
spatial unit, so community shares site-level field across species). Cover hurdle
joint uses `tulpa_nested_laplace_joint()` separately.

**INLA-style NA-response prediction + calibrated CIs** (single-season): all-NA
detection sites (`.tobs_heldout_sites()`) interpolated by field;
`predict(type="state")` returns per-site psi posterior (`psi`, `psi_lower`,
`psi_upper`, 95%) marginalised over hyperparam grid. Calibrated by one
exact-marginal refine pass (`.tobs_occu_state_marginal_fit()`): integrate out z,
each site Bernoulli on D=1{>=1 detection}, mean q*sigma(eta), q=1-(1-p)^J
(held-out -> q=0 -> dropped, field-interpolated). Fits tulpa
`family="bernoulli"` with per-obs prob scale `det_prob=q` (scaled-Bernoulli
family): converged Hessian = marginal curvature, `fitted_eta_var` = calibrated
per-cell predictive variance, grid weights natural scale (no M-inflation).
EM M-inflated pseudo-binomial kept only for mode/det estimate (reading
variance/weights off it under-covers + collapses grid). Per-row eta posterior =
Gaussian mixture over cells (`.nested_psi_mean()` GH mean,
`.nested_psi_quantiles()` mixture-CDF quantiles). Measured held-out coverage
~1.0 (conservative), cor ~0.88, MAE ~0.11 on 10x10 icar/bym2
(`dev_notes/probe_nested_ci_coverage.R`). Other model types keep EM occ fit (NA
mapping not wired). Old tulpa w/o `fitted_eta_var` -> `psi` only, NA intervals.

**Simplified-Laplace skew correction** (`simplified_laplace.R`, `sla_*` files):
orthogonal post-fit marginal refine, `approx="simplified_laplace"` (`*_sla`
methods). Computed for single/dynamic/integrated occ + cover hurdle; no-ops to
Gaussian (records `sla_status`) for jsdm.

**Backend coverage enforced centrally**: `.tobs_family_methods` in `R/tobs.R` =
single source of truth for which `method` each family supports; `tobs()` errors
with pointer, no silent downgrade. `nested_laplace` = occu/int_occu/dyn_occu +
cover; `*_sla` on nested = occu + cover only; cover hurdle has no
NUTS / `laplace_gibbs` / `laplace_mi`. `abun` = laplace (non-spatial) +
nested_laplace (areal). `ms_abun` = laplace + nested_laplace (shared areal field
on abundance arm). `occu_multiscale_cover` = nested_laplace ONLY (spatial joint).
Community occupancy families `ms_occu`/`ms_dyn_occu`/`ms_int_occu` = laplace ONLY
(shared community Laplace-EM, `R/community_em.R`; per-species coef RE, per-arm
community covariance). Planned (error via `.stop_planned_family()`): `dyn_abun`,
`distance`, `removal`, `fp_occu`.

### N-mixture abundance (`abun()`)

Royle 2004: `N_i ~ Pois(lambda_i)`, `y_ij|N_i ~ Binom(N_i, p_ij)`. NO EM — N
marginalises closed-form (sum to `K_max`), tulpa fits marginal directly,
analytic gradients + observed-Fisher. tulpaObs owns family wiring (`R/abun.R`):
binder -> `model_type="nmix"` (site `X_lambda`, long-form `y`/`site_idx`/`X_p`,
drop NA visits), `.tobs_fit_nmix()` -> in-tree fitter (`R/nmix_laplace*.R`,
`src/nmix_*.{cpp,h}`), wrapped `build_nmix_fit()`. Pseudo-draws MVN from JOINT
(lambda,p) cov -> cross-arm cov propagates. Abundance log link:
`compute_intercepts()` reads per-process `pi$link` (default logit, log for
lambda). `simulate_abun()` + `test-abun.R` cover point recovery + 95% coverage.

- **Non-spatial** (`laplace`): `nmix_laplace()`. vcov = marginal observed-Fisher
  inverse (full joint lambda/p block).
- **Areal spatial** (`nested_laplace`): `icar()`/`bym2()`/`car_proper()` on
  abundance formula -> `.tobs_fit_nmix_spatial()` ->
  `nmix_laplace_{icar,bym2,car_proper}` (one unit/site). Cov grid-integrated
  (law of total cov over hyperparam grid): kernels return per-grid `cov_blocks`
  (beta-block of grid mode `H^{-1}`); wrapper `V = sum_k w_k[cov_k +
  (m_k-mbar)(m_k-mbar)']` (`.nmix_grid_vcov()`). Rank-deficient intrinsic fields
  (ICAR, BYM2 structured): within-grid block under sum-to-zero constraint (big
  `(sum field)^2` penalty in `nmix_spatial_beta_cov()`/`nmix_beta_cov_bym2()`)
  so intercept variance reflects constraint, not flat intercept/field-mean
  confounding.
- **Negbin** (`abun(mixture="negbin")`): `mixture` poisson/negbin -> kernel
  P/NB, threads both paths. NB adds overdispersion `r` (`Var=lambda+lambda^2/r`).
  Non-spatial: `log_r` jointly estimated, trailing vcov coord; `build_nmix_fit()`
  carries it (means/vcov/draws gain `log_r` col, `nmix_dispersion` w/ delta `r_sd`);
  autoscale leaves trailing coord untouched; coef/confint/vcov report 2 arms only.
  Spatial: `r` integrated over outer grid (tau/rho/sigma), no `log_r` coord,
  summarized `r_mean`/`r_sd` in `nmix_hyper$r`. `simulate()` draws NB. Matches
  `unmarked::pcount(mixture="NB")`; `test-abun.R` covers recovery/dispersion/
  coverage/spatial-NB.
- **Random effects** (`abun(...) + (1|g)` / random slope, either arm;
  tulpaObs#13): site-level grouped RE on abundance OR detection of single-species
  nmix. Grouping = non-species (station, observer-per-site, cluster); RE = subset
  of coefs on ONE arm; `theta` = plain fixed-effect vector (no community-mean,
  wrong shape for community oracle; see `R/nmix_re_aghq.R` header).
  `.tobs_fit_nmix_re()` warm-starts no-RE Laplace betas, refines via
  `.tobs_nmix_re_aghq()` (thin `make_site` over `nmix_site_marginal()`). Engine
  `tulpa::tulpa_re_aghq(make_site=)`: adaptive GH (`n_quad=1` Laplace, default
  from `control$n.quad`; `>1` debiases small-cluster sigma). Det arm: per-site RE
  offset uniform over visits -> marginal = fn of one scalar shift/site (the
  per-row separability `make_site` needs). Visit-level p RE (observers vary
  WITHIN site) unreachable from public API (binder doesn't parse structured terms
  in `det_visit_formula`); per-site latent-N coupling would break factorization
  anyway -> needs different engine. NB threads `log_r` trailing theta. Gates
  (error+pointer): RE+areal spatial, RE+visit-level det, RE both arms.
  `test-abun-re.R` covers recovery/coverage/S3/gates. Native
  `NMixGroupedOracle` (`src/nmix_re_oracle.{h,cpp}`): per-group iterates only that
  group's sites, stays in C++ no round trip. ~tens of s -> ~2s at N=100/J=4/10
  groups (`dev_notes/probe_nmix_re_oracle_speedup.R`); 27-test suite ~110s.
- **Pending**: no NUTS (no nmix HMC likelihood).

### Community / multispecies N-mixture (`ms_abun()`)

spAbundance `msNMix`: per-species nmix w/ Gaussian community hyperpriors on
per-species abundance + det coefs, `beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`,
`beta_p_s ~ N(mu_p, Sigma_p)`. `N_{s,i}` integrates closed-form per species-site
(same kernel as `abun()`); per-species deviations `b_s` = random effects,
community priors pin `(mu_lambda, mu_p)` fixed (no sum-to-zero, unlike spatial).
`y` = 3D array `[sites x max_visits x species]` or named list of count matrices;
`species=` required.

Family wiring `R/ms_abun.R`: binder -> `model_type="ms_nmix"`,
`.tobs_ms_nmix_longform()` flattens 3D -> stacked `(y, site_idx, species_idx,
X_p)`, `.tobs_fit_ms_nmix()` -> in-tree C++ Laplace-EM `nmix_laplace_re()` (fast
path). Builds native `tulpaObs::NMixCommunityOracle` (subclass
`tulpa::REGroupOracle`, `<tulpa/aghq_oracle.h>`) via
`cpp_nmix_community_oracle()`, then drives via in-tree `cpp_nmix_community_em()`
(default `n_quad=1`) OR `tulpa::tulpa_re_aghq(oracle=)` for `n_quad>1` variance
debias. Fit: per-species coef mode-find (complete-data Fisher, PSD), closed-form
EM cov M-step `Sigma_k=mean_s[b_s b_s'+Cov(b_s|y)]`, fixed-effect SEs from
marginal observed-info Schur complement of b-block (`Var[N|y]` rank-1 cross-arm
coupling; Louis 1982). Per-site kernel `src/nmix_kernel.h` caches eta-independent
`lgamma` (`compute_nmix_site_cached`); single-shot `compute_nmix_site()` Poisson
delegates -> single-species/spatial byte-identical (regression passes).

`coef()` = community means; `ranef()` = per-species deviations (long:
species/arm/term/estimate); `vcov()`/`confint()` = community-mean cov;
`fitted()`/`simulate()` = per-species lambda/p/counts. `simulate_ms_abun()` +
`test-ms-abun.R` cover community-mean recovery, 95% coverage (20 seeds),
per-species coef recovery, S3.

- **Poisson only.** Native oracle already carries global NB `log_r` + kernel
  branches on `is.finite(r)`; just needs `mixture` plumbed.

#### Areal-spatial community N-mixture (`ms_abun()` + shared field; sfMsNMix)

Shared ICAR/BYM2/proper-CAR field on abundance via `nested_laplace`
(tulpaObs#12): term on abundance formula -> `.tobs_fit_ms_nmix_spatial()` ->
in-tree nested Laplace-EM `nmix_community_laplace_{icar,bym2,car_proper}()`
(`R/nmix_laplace_re_spatial.R`) over `cpp_nmix_community_spatial_*`
(`src/nmix_community_spatial.cpp`). Model `log lambda_{s,i} = X_lambda_i.(mu_lambda
+ b_lambda_s) + f_{u(i)}`, one `f` shared across species. tulpa AGHQ/nested-Laplace
don't combine per-group RE block + shared field on one predictor -> dedicated
tulpaObs fit: top block `(mu_lambda, mu_p, f)` same layout as single-species
spatial state (`field_start = p_lambda+p_p`), so `nmix_spatial_kernel*.h` helpers
apply. Per outer grid point (tau[,rho]/sigma,rho[,r]) EM iterates joint
`(mu,f,{b_s})` mode-find (block-elim Newton on complete-data-Fisher arrowhead,
`b_s` Schur-folded) + closed-form `Sigma` M-step; at convergence accumulate
grid-point Laplace log-marginal + b/field-folded community-mean cov `vcov_mu`,
then R grid-integrates (law of total cov for means; weighted mean for
field/Sigma/BLUPs). Community means kept FLAT (no ridge): abundance intercept +
field constant mode exactly flat dir; ridge would drive intercept->0, field
absorbs level (then sum-to-zero deletes). NB `r` integrated over outer grid
(field-agnostic) -> hyperparameter (`ms_dispersion`/`ms_hyper`), no `log_r`
coord. `fit$spatial_field` = posterior-mean field. `test-ms-abun-spatial.R`:
community-mean + field-shape recovery (icar/bym2/car_proper, Pois+NB), 95%
coverage, S3, `spAbundance::sfMsNMix()` interop SMOKE (plumbing only). Numerical
head-to-head = offline manual bench
(`dev_notes/probe_ms_abun_spatial_vs_spabundance.R`), not test (converged
sfMsNMix chain runs hours). NUTS pending.

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
| Cover hurdle (joint) | Yes | — | `family_cover_hurdle.R`, `sla_cover_*` |
| Joint occu + cover | Yes | — | `occu_cover()` — see below |
| Community joint occu + cover | Yes | — | `ms_occu_cover()` — see below |
| Spatial-factor community occu + cover (JSDM) | Yes | Yes | `ms_occu_cover()` + `icar()`/`car_proper()`/`bym2()` shared field, per-species loadings (gcol33/tulpa#67). Laplace-EM (`R/ms_occu_cover_spatial.R`) + NUTS (`method="nuts"`, in-tree C++ FullGradFn `src/ms_occu_cover_spatial_nuts.cpp` driving tulpa's sampler; byte-exact vs the R target). Cover-arm factor, auto-K (Laplace only), `tobs_associations()`, per-species `predict()` maps |
| Multiscale occu + cover | n-L | — | `occu_multiscale_cover()` — 3-level cell/plot/visit + cover; spatial joint only — see below |
| N-mixture (Pois/NB) | Yes | — | `abun(mixture=)`; in-tree `nmix_laplace`, joint-vcov draws, calibrated CIs (`test-abun.R`). NB jointly-est `log_r`. NUTS pending |
| N-mixture + areal spatial | n-L | — | `abun()`+icar/bym2/car_proper, `nested_laplace`; Pois/NB (r grid-int); grid-int cov (constrained intercept) |
| Community N-mixture | Yes | — | `ms_abun()` (msNMix); per-species coef RE, in-tree C++ Laplace-EM (`nmix_laplace_re`) -> native `NMixCommunityOracle` via AGHQ, Schur SEs; Pois only; recovery+20-seed (`test-ms-abun.R`). NUTS/negbin pending |
| Community N-mixture + areal spatial | n-L | — | `ms_abun()`+icar/bym2/car_proper, `nested_laplace` (sfMsNMix; tulpaObs#12); shared field on log lambda + per-species RE; nested Laplace-EM (`nmix_community_spatial.cpp`), joint `(mu,f,{b_s})` mode-find, field/Sigma grid-int; Pois/NB; `test-ms-abun-spatial.R`. NUTS pending |
| N-mixture + grouped RE | Yes | — | `abun()`+`(1\|g)`/`(x\|g)` either arm (tulpaObs#13); non-species grouping; Pois/NB; AGHQ via `NMixGroupedOracle`. Gated: RE+spatial, RE+visit-det, RE both arms. ~2s/fit (N=100,J=4,10g,n.quad=5) |
| Spatial ICAR/BYM2/NNGP | — | Yes | |
| Spatial + dynamic | — | Yes | |
| Nested-Laplace (areal) | n-L | — | `nested_laplace`: icar/bym2/car (+temporal/iid) on occu/int_occu/dyn_occu |
| NA-response prediction | n-L | — | `predict(type="state")`: all-NA single-season sites field-interpolated, calibrated 95% `psi_lower`/`psi_upper` from exact-marginal bernoulli pass (coverage ~1.0) |
| Formula RE (intercept) | Yes | Yes | `(1\|g)`; variance-component EM, occ OR det arm (tulpaObs#11) |
| Formula RE (uncorr slope) | Yes | Yes | `(x\|\|g)`, `(0+x\|g)`, `(1+x\|\|g)`; either arm |
| Formula RE (corr slope) | Yes | Yes | `(1+x\|g)`; EM M-step consumes tulpa `cov_blocks`; either arm (tulpaObs#11) |
| Formula RE on detection | Yes | Yes | `detection=~(1\|g)`; separate det-arm RE block, AGHQ branches on arm |
| All S3 methods | Yes | Yes | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, `$.tobs_fit` |
| Diagnostics | Yes | Yes | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatialRange, temporalCorr |
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
per-cell occupancy mixture closed-form derivs drive inner Newton. Lognormal+beta
both. ~150-300x faster than v3_nested at N=100
(`dev_notes/probe_bench_v3_vs_jc.R`), completes N=200+ where v3 trips on
missing-value compare in outer BFGS. Recovery: 10-seed lognormal + 10-seed beta
(`test-occu-cover-joint-coupled.R`); status `"experimental"`.

**Cover-arm intercept prior (tulpaObs#32)**: on the shared-field path the cover
intercept confounds with the field level over detected cells (the cover arm sees
the field only where detected; the sum-to-zero field constraint pins only the
global field mean). `.occu_cover_coupled_arm_priors()` therefore hands the pos
arm the `cover_priors()` weakly-informative intercept prior **by default** (like
the load-bearing detection-arm prior), not the engine's flat 1e-4 ridge —
otherwise the cover intercept floats to a huge posterior SD (occupancy stays
tight, being regularised + observing every cell) and `predict()`'s conditional
cover blows up via Jensen. `priors = FALSE`/`"none"` disables all three arms.

**Cell-aggregated cover (`cover_aggregate`, tulpaObs#33)**: per-visit cover gives
the cover arm one row per *valid visit*, so a cell with many detected plots
drives the shared field far more than its single occupancy obs (field flattens).
`cover_aggregate = "mean"` (default on the spatial path) / `"median"` collapses
the cover arm to ONE row per occupancy unit (the per-site mean/median cover over
detected visits) so the two arms inform the field with comparable weight;
`"none"` keeps per-visit. Wired ONLY on the spatial `joint_coupled` path (v2/v3 +
non-spatial laplace reject explicit aggregation, default falls back to per-visit).
Aggregation needs a **cell-level** positive design (resolved from `data`, not
`visits`): a visit-level `positive` covariate keeps the per-visit arm (bare
default falls back, explicit request errors). C++: compile-time `Aggregated`
template flag on `OccuCoverCoupling` (`src/cell_coupling_occu_cover.h`) evaluates
the pos density once per cell at arm-2 row 0; registered as
`occu_cover_{lognormal,beta}_agg`. R: `.occu_cover_build_joint_coupled_arms(
cover_aggregate=)` builds the one-row-per-detected-site pos arm; the fitter picks
the `_agg` spec + pre-fits dispersion on the aggregated values. FD checks
(`test-occu-cover-coupling.R`), recovery + gates (`test-occu-cover-aggregate.R`).

**Latent cover-per-unit (`cover_aggregate = "latent"`)**: the principled version
of mean/median — instead of collapsing a unit's detected covers to one number,
the cover arm carries a per-unit cover RE `u_i ~ N(0, sigma_u^2)` shared across
the unit's detected visits and integrates it out per unit (keeping every detected
visit). Because the cover predictor is unit-level the per-unit marginal `log M_i`
is a SCALAR function of one eta, so it reuses the one-row-per-unit layout with no
within-arm Hessian coupling. Lognormal = closed form (compound-symmetry suff
stats `m,T1,T2` → `src/occu_cover_latent.h::LognormalLatent`); beta = adaptive
Gauss-Hermite over scalar `u_i` reusing `BetaPositive` (`BetaLatent`, GH nodes
from the engine's exported `<tulpa/gauss_hermite.h>`, `control$n.quad` default
15). Two dispersions: the within-unit one (`sigma_eps` / beta precision) is
pre-fit from the **within-unit** spread (NOT the total — that would swallow
`sigma_u`) and held FIXED in the spec; `sigma_u` rides the pos arm's `phi_grid`
axis (`control$sigma.u.grid`), integrated on the outer grid and reported as the
`phi_pos` hyperparameter (`joint_fit$theta_mean[["phi_pos"]]`). The spec is
**stateful** (`OccuCoverLatentCoupling<PosLatent>`,
`src/cell_coupling_occu_cover_latent.h`): it captures the per-unit cover data +
fixed dispersion at construction and is (re)registered per fit via
`cpp_register_occu_cover_{lognormal,beta}_latent_coupling()` (last-writer-wins; the
joint driver holds the resolved shared_ptr for the fit). Shared det-branch psi/p
+ nodet logic is factored into `occu_det_psi_p_block` / `occu_nodet_block`
(`src/occu_coupling_shared.h`), used by both the per-visit/agg template and the
latent spec. Same gates as aggregation. FD vs brute-force integration + recovery
in `test-occu-cover-latent.R`.

**Coupled SVC/trend fields** (tulpaObs#15): extra shared areal fields = WEIGHTED
areal terms in psi formula — `icar(graph=adj, weight=year)` couples a
spatially-varying coef on `year` atop unweighted intercept field
`icar(graph=adj)`. N fields compose (each own outer-grid `alpha`) via tulpa
multi-block copy (`prior=list(icar_intercept, icar_trend, ...)`, one `copy`/block
onto pos arm, per-block `svc_weight`). Intercept field (`fit$spatial_field`)
couples uniformly; each weighted field (`fit$trend_field`, `fit$trend_fields` for
>1) weighted by `weight_i` on psi (per-cell) + pos (per-visit broadcast), own
scale (`alpha_trend`, `control$alpha.grid.trend`). p (detection) arm excluded via
`field_coef=0` (per-arm scalar zeroes field every block — `svc_weight=0` NOT
equivalent: leaves p-intercept biased through cross-arm Hessian; verified
`dev_notes/_phase1_altwire.R`). Weighted term = formula-DSL of field list,
resolved `.occu_cover_spatial_fields()` (exactly 1 unweighted base + N weighted
SVC, one graph); `control=list(trend=list(weight="<col>"))` = back-compat alias
for single trend field (one way or other, not both — errors). Weighted areal off
joint occu_cover path errors via `.tobs_reject_weighted_spatial()`. Recovery:
10-seed 2-field lognormal + route-equiv (`test-occu-cover-trend.R`).

**Escape hatches**: `control$engine="v3_nested"` (pure-R outer-BFGS,
`R/occu_cover_nested.R`, lognormal only, `test-occu-cover-spatial.R`),
`control$engine="v2_joint"` (v2 joint Laplace, kept for (z,alpha,sigma) ridge
investigation).

**`group_var` (sites > cells)**: `group_var="<col>"` on `icar()`/`bym2()` maps
each site (y/data row, one occupancy state) -> field node, so `n_sites` >
`n_cells`. Field stays length `n_cells` (graph) while psi/p/cover run over
`n_sites`; per-arm `spatial_idx` (field node) + `cell_obs_map` (occupancy unit)
decouple — psi `spatial_idx=site_cell`, pos `spatial_idx=cell_of_visit`, both
`cell_obs_map` over sites; SVC trend weight per-site. Motivating layout: site =
cell x time-period (plots-in-cell-period = detection replicates) ->
detection-corrected occupancy trend on shared cell field. R-side only
(`.dispatch_occu_cover` resolves `site_cell`,
`.occu_cover_build_joint_coupled_arms` splits site/cell); joint_coupled only
(v2/v3 reject). Recovery: `test-occu-cover-group-var.R`.

### `ms_occu_cover()` detail

Community version of `occu_cover()` (`R/ms_occu_cover.R`); per-species coef RE w/
Gaussian community covariances across psi/p/cover arms, shared dispersion. Latent
z integrates closed-form per species-cell (reuses `.occu_cover_site_ll`);
per-species deviations integrated by in-tree pure-R Laplace-EM — arrowhead joint
Newton (per-species RE Schur-folded, analytic grads `.occu_cover_eta_grad`,
FD-of-gradient observed-info) + closed-form community-cov M-step
`Sigma_arm=mean_s[b_s b_s'+Cov(b_s|y)]`. Community-mean SEs = marginal observed
info (Louis 1982 Schur of b-block), natural scale. Beta+lognormal. Non-spatial
Laplace only; `nested_laplace` NOT offered (community analogue needs upstream
tulpa support combining per-group RE + shared field; structured term any arm
errors+pointer). Community-mean recovery + 15-seed coverage + per-species coef
recovery (`test-ms-occu-cover.R`); status `"experimental"`. Binary-detection
community VARIANCE carries Laplace small-cluster attenuation (EM n_quad=1, like
single-species RE); community MEANS do not. NUTS/negbin/per-species dispersion
RE/AGHQ variance debias pending.

### `occu_multiscale_cover()` detail

Three-level occupancy + cover hurdle (tulpaObs#29; `R/occu_multiscale_cover.R`,
fitter `R/occu_multiscale_cover_joint_coupled.R`). For data where "visits" are
spatially distinct PLOTS aggregated into `(cell, period)`, not temporal revisits
(EVA/MOTIVATE vegetation; Nichols 2008, Mordecai 2011). `occu_cover()` treats
plots as detection replicates of one occupancy state -> conflates within-cell
prevalence into detection (Kendall & White 2009); this family adds explicit
middle level:

```
z_c        ~ Bernoulli(psi_c)        # cell/range occupancy
a_cj|z=1   ~ Bernoulli(theta_cj)     # plot availability/use
y_cjv|a=1  ~ Bernoulli(p_cjv)        # detection
cover|y=1  ~ f_pos(eta_pos, disp)    # hurdle (beta/lognormal)
```

Both z (cells) + a (plots) marginalize closed-form (two states each) -> exact
joint marginal LL, reuses occu_cover nested-Laplace cell-coupling machinery.

**Inputs**: `y`/`y_pos` = `[n_plots x max_visits]` matrices (plots = rows /
availability units, visits = cols). State `formula` = cell-level psi, MUST carry
areal field naming per-plot cell col: `icar(graph=adj, group_var="cell")`.
`availability=~...` = plot-level theta (default `~1`); `detection` = per-visit p;
`positive=~...` = cover (default = detection formula). `y_pos` read only where
`y==1`.

**Engine** (`method="nested_laplace"` only): 4-arm generalization of occu_cover
joint_coupled via `tulpa_nested_laplace_joint(cell_coupling=
"occu_multiscale_cover_*")`. Arm field coupling: psi row/cell `field_coef=1`
(shared field); theta row/plot `field_coef=0`; p row/valid-visit `field_coef=0`;
pos row/valid-visit `field_coef=list(name="alpha")`. Single shared intercept
field, outer `(sigma, alpha)` grid. Cell-coupling spec
(`src/cell_coupling_occu_multiscale_cover.{cpp,h}`, registered per-fit carrying
this fit's per-cell plot structure) writes closed-form derivs; visit rows
compacted to valid + laid plot-major within cell. Reuses occu_cover CSR /
field-demean / inner-vcov / rmvn helpers (shared `src/occu_coupling_shared.h`).

**Identifiability**: theta + p separate only with replication WITHIN a plot.
Single releves -> plain fit identifies psi (cell) + product theta*p (plot),
reduces to occu_cover w/ `p:=theta*p`. Within-plot temporal replication
(plot resurvey later period) makes 3rd level estimable.

**Scope** (status `"experimental"`): spatial joint nested-Laplace only. SVC
trend fields + non-spatial Laplace not wired. Dispatcher
`.dispatch_occu_multiscale_cover` rejects laplace + non-spatial state formula.
Recovery: 8 fixed-effect coefs within 0.25 mean, 95% Wald coverage >=0.80,
field cor(z_hat, f_true)>0.70, lognormal + beta arms
(`test-occu-multiscale-cover-recovery.R`); spec FD-derivative + 2-level-reduction
checks (`test-occu-multiscale-cover-coupling.R`). `simulate_occu_multiscale_cover()`.

## NUTS coverage status

Prior CLAUDE.md flagged `temporal`, multi-term `re`, `svc`, `latent` as NUTS
crashers. Smoke-tested 2026-05-20 (`dev_notes/probe_blocked_nuts.R`, N=40, 50
iter/25 warmup, single-season occ) — return `tobs_fit` w/o crash. C++ population
not segfaulting. (The community NUTS path the original probe used for `latent`
was removed with the legacy community engine; tulpaObs#30.)

NOT verified: gradient correctness, posterior calibration, production-iter
convergence, stacking (spatial+temporal+multi-RE). No `tests/testthat/` exercise
these under NUTS — only Laplace/nested-Laplace covered. Treat as "not blocked",
not "validated".

## Performance

N=200, single-season (J=3), 2026-05-24:
- `laplace` (prior-aware penalized EM, default): **~0.2s**
- `laplace_gibbs` (EM+Gibbs, n.gibbs=10): ~1.7s
- `nuts`: ~13s (historical)
- inlaocc 0.7s; spOccupancy MCMC 0.9s (historical ref)

Old "0.01s" predates prior-aware default; penalized EM trades speed to break
psi-p identifiability ridge at small J. Gibbs/MI add Rubin-pooled correction.

## File organization

```
R/
  tobs.R                    — tobs() dispatcher + print.tobs_fit
  obs_families.R            — family ctors (occu, dyn_occu, int_occu, jsdm, …, cover)
  occu.R                    — .tobs_build_model() (single/dynamic/integrated/jsdm/nmix)
  community_em.R            — shared community Laplace-EM engine .tobs_community_em() (arrowhead Newton + per-arm covariance M-step + marginal info); drives ms_occu/ms_dyn_occu/ms_int_occu
  ms_occu.R                 — community single-season occupancy (msPGOcc): build, fit (reuses ms_int single-source kernel), S3, .tobs_richness_ms_occu, ms_occu() ctor
  ms_dyn_occu.R             — community dynamic occupancy: build, HMM-forward marginal, fit (psi1/p RE + gamma/eps global), S3, simulate, ms_dyn_occu() ctor
  ms_int_occu.R             — community integrated occupancy: build, multi-source two-state marginal + analytic grad, fit, S3, simulate, ms_int_occu() ctor
  abun.R                    — nmix family: build_abun, fit_nmix, build_nmix_fit, S3, simulate_abun
  ms_abun.R                 — community nmix: build_ms_abun, ms_nmix_longform, fit_ms_nmix (-> nmix_laplace_re), build_ms_nmix_fit, S3, simulate_ms_abun
  nmix_laplace.R            — in-tree non-spatial nmix (Royle 2004) Laplace (Pois+NB)
  nmix_laplace_re.R         — in-tree community nmix (msNMix), .nmix_re_oracle()
  nmix_laplace_re_spatial.R — spatial community nmix (sfMsNMix): _icar/_bym2/_car_proper over cpp_nmix_community_spatial_*, grid-int means/field/Sigma
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
  ms_occu_cover.R           — community joint: build_ms_occu_cover, per-species kernel, Laplace-EM fit, build fit, S3, simulate
  ms_occu_cover_spatial.R   — reduced-rank spatial-factor community occu_cover (JSDM, tulpa#67): simulate, field-structure abstraction (icar/car/bym2), penalised joint LL+grad, constrained-loading param, Laplace-EM fitter, build fit, associations/maps, S3
  ms_occu_cover_spatial_nuts.R — R NUTS target (.ms_ocs_joint_logpost = the C++ oracle), chol/hyperprior helpers, .ms_ocs_nuts_spec marshalling, .tobs_fit_ms_occu_cover_spatial_nuts (method="nuts" runner)
  occu_multiscale_cover.R   — 3-level occu+cover (tulpaObs#29): builder, dispatcher, simulate_occu_multiscale_cover
  occu_multiscale_cover_joint_coupled.R — 4-arm joint nested-Laplace fitter (cell/plot/visit + cover)
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
  occu_coupling_shared.h      — shared coupling helpers (CSR/field-demean/inner-vcov/rmvn); reused by the spatial-factor NUTS marginal (nodet_mixture_block, Lognormal/BetaPositive)
  ms_occu_cover_spatial_nuts.cpp — spatial-factor community occu_cover NUTS: marginal LL + full joint log-posterior gradient (FullGradFn), field R(h) layer (icar/car/bym2), cpp_ms_ocs_nuts driving tulpa's sampler (tulpa#67)
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
