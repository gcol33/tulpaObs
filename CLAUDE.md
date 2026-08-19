# CLAUDE.md

Guidance for Claude Code in this repo. Caveman speak: terse, telegraphic. File
paths + function names exact. (Compacted; if a detail is missing, read the source
file named.)

**This file has a hard 150k-char budget and is near it.** Keep RULES + CONTRACTS
here; put the MEASUREMENTS they rest on (fixtures, seeds, wall times, per-seed
recovery numbers, cost ratios) in `NOTES_measurements.md`, committed alongside and
Rbuildignored. Adding a paragraph of numbers here costs someone else's context.

**Adding a family / prior / backend / diagnostic -> read `DESIGN.md` FIRST.** This
file is the ROSTER (what each family supports); `DESIGN.md` is the EXTENSION POINTS
(which single registry to edit, which duplication is deliberate). Short by design;
read it in full. Do not restate the roster there or the extension points here.

## tulpaObs — hierarchical latent-state observation models on tulpa

Bayesian occupancy / abundance / distance / removal / cover. Built on
[`tulpa`](https://github.com/gcol33/tulpa) engine. R pkg, C++17 backend
(Rcpp/RcppEigen). Needs sibling `../tulpa` checkout (`LinkingTo: tulpa`).

**Public API:** `tobs()` + family ctors: `occu()`, `dyn_occu()`, `ms_occu()`,
`ms_dyn_occu()`, `ms_int_occu()`, `int_occu()`, `jsdm()`, `abun()`, `ms_abun()`,
`dyn_abun()`, `count()`, `ms_count()`, `distance()`, `ms_distance()`, `removal()`,
`fp_occu()`, `royle_nichols()`, `occu_ttd()`, `occu_multi()`, `double_observer()`,
`dyn_int_occu()`, `t_occu()`, `gdistremoval()`, `distsamp_open()`,
`occu_categorical()`,
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
   `Sys.setenv(TULPAOBS_FAST = "1"); devtools::test()`. Baseline 2026-07-20
   SERIAL (`TESTTHAT_PARALLEL=false`, `test_dir()` after `load_all()`, tulpa
   0.0.93): 3322 assertions, 629 skips, **0 fail, 0 error**. Earlier 5 errors =
   #148 (tulpa fed its own nested-Laplace validator a length-1 `re_idx`), fixed
   tulpa 0.0.93; assertion count rose b/c those blocks now run. Wall time NOT
   measured (box shared) -- time it yourself, do NOT trust old ~153s / ~22s
   figures. `skip_if_fast()` gates every fitting block (628 sites, 166 files);
   no-op when env var unset.
3. **Full recovery suite** (all seeds, NUTS, spatial) -> CI cron
   (`full-recovery.yaml`), or on request. NOT a release gate, NOT pre-commit:
   ~9h serial, one file never terminates -> local whole-tier run finishes only
   by luck. `Sys.unsetenv("TULPAOBS_FAST")`, then parallel note below --
   `devtools::test()` cannot use `Config/testthat/parallel` here.

**Release gate = what the diff changed, NOT the whole tier.** Tarball
`R CMD check` + whole-suite smoke vs INSTALLED package + recovery files for the
families the diff touches. Calibration evidence (every multi-seed recovery /
coverage loop) rests on `full-recovery.yaml`'s last green run + the version it
ran against; a release states which run it inherits. Policy:
`tests/testthat/helper-speed.R`. As of 2026-07-28 that workflow completed ZERO
times (one run 2026-07-25, cancelled at 350-min cap) -> calibration evidence
ABSENT, not stale; any claim resting on it = unverified.

**`devtools::test()` can NEVER run this suite in parallel (#151, won't-fix
upstream).** It hardcodes `load_package = "source"` per worker
(`devtools:::load_package_for_testing()`, no override) -> N `callr` workers each
`pkgload::load_all()` the same tree. Harmless pure-R; here `src/tulpaObs.dll`
~165MB, so N workers recompile into one `src/` and race the build artifacts.
Reproduced in isolated copy: cold `src/` + 2 concurrent `load_all()` corrupted
DLL registration in <1 min (`Error in getDLLRegisteredRoutines.DLLInfo(dll, ...)
: must specify DLL via a "DLLInfo" object`), wrapped by testthat
`cli_abort(..., parent = msg$error)` -- real error in that `parent`, NOT the
printed message (`rlang::cnd_message(e, inherit = TRUE)` / `e$parent$message`).
Silent 10-min hang w/ zero output = same race, different landing: partial DLL
load blocks in the loader w/o emitting its startup handshake -> `queue$poll(Inf)`
waits forever.

Safe recipe = what `.github/scripts/run-tests.R` does, works locally too:
install ONCE, then call testthat directly (never via `devtools::test()`), every
worker loading the built package instead of recompiling:

```r
devtools::install(quick = TRUE)
testthat::test_dir("tests/testthat", package = "tulpaObs",
                   load_package = "installed")
```

Verified: `load_package = "installed"` runs multi-file parallel in ~2-3s, warm
or cold `src/` -- each worker's `library(tulpaObs)` = one read of one built DLL,
safe at any worker count. CI runs `smoke.yaml` + `full-recovery.yaml` w/
`TESTTHAT_PARALLEL: true` for this reason.

Slow test -> pair `skip_if_fast()` + `skip_on_cran()` atop any multi-seed fit /
NUTS block (`tests/testthat/helper-speed.R`). C++ recompiles ccache-backed; only
killed/partial build needs `pkgbuild::clean_dll()`.

### CI (#149)

`.github/workflows/`:

- `R-CMD-check.yaml` -- push/PR + weekly. Checks a **built tarball**, not
  `load_all()`, so a missing NAMESPACE export surfaces (#147 shipped `abun()`
  unexported precisely because `load_all()` resolves internals regardless).
  `--no-manual` (dev non-ASCII in Rd). ubuntu on push; ubuntu+windows+macOS on
  weekly cron + `workflow_dispatch`. Serial; parallel here goes via R's own
  `test_check()`, which already hardcodes `load_package = "installed"` -> not
  subject to #151. Left alone: this job = tarball correctness, not speed.
- `smoke.yaml` -- push/PR, tier 2 (`TULPAOBS_FAST=1`) vs INSTALLED package.
  Catches #148-class breakage the day it lands. `TESTTHAT_PARALLEL: true`
  (#151, safe here).
- `full-recovery.yaml` -- weekly cron + dispatch, tier 3, `NOT_CRAN=true` +
  `TULPAOBS_REQUIRE_SPDE=1`. Hours; carries the calibration evidence.
  `TESTTHAT_PARALLEL: true` -- the tier #151 was filed to unblock.

Both test workflows -> `.github/scripts/run-tests.R` (one runner; logs which
tier ran -- smoke vs broken-full otherwise report similar counts). Every job
first runs `.github/scripts/check-engine-pin.R`: fails if DESCRIPTION `Imports`
floor, `Remotes` tag, and installed `tulpa`/`tulpaMesh` disagree -> #150 skew
cannot reopen silently.

`R-CMD-check` also sets `TULPAOBS_FAST=1`: check runs `tests/` itself, ungated
that = tier 3, blew the job cap. Green baseline 2026-07-20: check `Status: OK`
(0/0/0), smoke 3320 assertions / 630 skips / 0 fail on ubuntu.

**Linux runner != local box; two test classes feel it.** Both surfaced on the
first CI run, neither reproduces on Windows:

- *Exact float equality across two BLAS call shapes.* Batched ploglik
  (`R/diagnostics.R`) builds eta w/ one `[S x p] x [p x n]` GEMM; R oracle
  builds it per draw as GEMV. Different routine -> different accumulation order
  -> eta differs in the last ULP, and a marginal summing ~100 terms through
  `exp`/`lgamma` carries it (max relative 6.5e-16). Reference BLAS collapses
  both to the same naive ordering -> equality held on Windows ALONE. Assert a
  tolerance AND report the magnitude: bare `expect_true(a == b)` cannot separate
  1 ULP from a real bug, and the two want opposite responses. Thread-count
  invariance IS exact (same GEMM, same kernel) -> stays `expect_identical`.
- *Single-seed recovery under an aggregate tolerance.* `expect_equal(coefs,
  truth, tolerance = t)` scores the WHOLE vector -> a one-sided bias in one
  coefficient sits under `t` on one platform, over it on another. That was #153
  (lfJSDM slope 1.44x truth, high in 15/16 seeds, biased on BOTH platforms,
  tolerance hiding it). Test flips -> measure the bias over many seeds BEFORE
  touching the tolerance; widening it deletes the signal.

  Two things make that tolerance looser than it reads; take out BOTH before a
  community recovery number means anything:

  - *Estimand = wrong constant.* Community simulator draws per-species coefs
    around a POPULATION mean (`rnorm(S, 0.2, 0.4)`) -> the seed's realized mean
    sits `beta_sd/sqrt(S)` off that constant, so scoring vs the constant spends
    most of the budget on draw noise. Score vs `colMeans(bs)` (seed's own
    realized mean) -> pure estimator budget (jsdm fixture: honest tolerance
    0.35 -> 0.10/0.12). `.jsdmc_sim` returns `beta_real`; most other community
    simulators still return only the nominal constant (#155).
  - *`tolerance` silently switches scale.* `all.equal.numeric` = relative while
    target exceeds tolerance, absolute below -> one `tolerance = 0.35` meant
    +-0.28 on a 0.8 slope, +-0.35 on a 0.2 intercept. Assert absolutely
    (`expect_lt(abs(est - truth), tol)`) when the budget is meant to be uniform.

  One-sided shift = property of the MEAN deviation over seeds, NOT of one fit
  -> assert in a multi-seed loop. Single-seed assertion = gross-regression
  guard, budget it as one.

Suggests are NOT all installed (`_R_CHECK_FORCE_SUGGESTS_: false`); INLA is a
large non-CRAN install and every INLA vignette chunk is already `eval =
have_inla`, every INLA/unmarked/spAbundance test `skip_if_not_installed()`.

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
hyperparam grid. Calibrated by ONE exact-marginal refine
(`.tobs_occu_state_marginal_fit()`): integrate out z, each site Bernoulli on
D=1{>=1 detection}, mean q*sigma(eta), q=1-(1-p)^J. Fits tulpa `family="bernoulli"` w/
per-obs prob scale `det_prob=q` (scaled-Bernoulli) -> converged Hessian = marginal
curvature, `fitted_eta_var` = calibrated per-cell predictive variance. Per-row eta
posterior = Gaussian mixture over cells (`.nested_psi_mean()`,
`.nested_psi_quantiles()`). Coverage ~1.0, cor ~0.88, MAE ~0.11 on 10x10 icar/bym2.
Old tulpa w/o `fitted_eta_var` -> NA intervals.

**Simplified-Laplace skew correction** (`simplified_laplace.R`, `sla_*` files):
orthogonal post-fit marginal refine, `approx="simplified_laplace"` (`*_sla` methods).
Computed for single/dynamic/integrated occ + cover hurdle; no-ops to Gaussian
(records `sla_status`) for jsdm.

**Backend coverage enforced centrally**: `.tobs_family_methods` in `R/tobs_helpers.R`
(NOT `R/tobs.R`, which only calls `.tobs_validate_family_method()`) = single source of
truth for which `method` each family supports; `tobs()` errors w/ pointer, no silent
downgrade. READ THAT OBJECT before trusting any support claim here -- roster below
drifts. `nested_laplace` = occu/int_occu/dyn_occu + cover; `*_sla` on nested = occu +
cover only; cover hurdle has NO `laplace_gibbs`/`laplace_mi` but DOES have nuts
(`R/cover_nuts.R`, `src/cover_nuts.cpp`). `abun` = laplace + nuts (non-spatial) +
nested_laplace (areal). `ms_abun` = laplace + nested_laplace (shared areal field) +
nuts (non-spatial #14, + shared fixed-hyper proper-CAR field #73).
`occu_multiscale_cover` = laplace + nested_laplace + nuts (last two non-spatial;
`R/occu_multiscale_cover_nuts.R`, `src/occu_multiscale_cover_nuts.cpp`).
`occu_categorical` = laplace ONLY. Community occupancy
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
  total. **Detection-arm areal field (#114)**: `detection=~icar()/car_proper()/bym2()`
  loads a spatially-varying capture logit on eta_p via the shared areal-BFGS driver
  (`.tobs_fit_removal_spatial_bfgs(det_arm=TRUE)`); removal's detection design is
  per-PASS, so the per-pass grad_eta_p is `rowsum`-summed to a per-site field
  gradient + the per-site offset expanded via site_idx. `fit$spatial_field_arm`
  labels it. NUTS carries the field on the abundance arm only (detection-arm ->
  nested_laplace). spde + NUTS+detection-arm gated.
- **`distance`** (#38): latent N in covered region, per-bin counts multinomial over
  `(bin 1..B, undetected)`, `pi_b = int_bin g(x;sigma) f(x) dx` (half-normal/hazard
  key, line/point transect density `f`); bin integrals + 1st/2nd eta-derivs by
  Gauss-Legendre quad. Reuses count-marginal core (`src/nmix_kernel.h`); det arm =
  site-level `log sigma` + optional scalar hazard shape (`src/distance_quad.h`/
  `src/distance_kernel.h`), own Laplace driver (`src/distance_laplace.cpp`) + NUTS
  (`src/distance_nuts.cpp`). K_max default `3*max(rowSums)+100`. Laplace grouped RE on
  abundance arm (one grouping factor, dim<=3) via shared count AGHQ
  (`DistanceGroupedOracle` over `CountGroupedOracle`; theta = count layout
  `[beta_lambda|beta_sigma|log_r?]`). **hazard-key grouped RE (#114)**: the scalar
  log-shape eta_b is not a per-site design column, so rather than a second global
  theta slot it is PROFILED over the AGHQ log-marginal (`.tobs_distance_re_aghq`
  optimises eta_b in an outer loop, each candidate a full AGHQ fit at a FIXED shape
  the `DistanceGroupedOracle` carries via `key_code`/`eta_b`); the profile log-shape
  row/col is inserted into the AGHQ vcov (`.tobs_distance_insert_shape_vcov`, profile
  SE, off-diag 0). Pois + NB; det-arm RE gated (couples a site's bins through the
  latent N). Areal icar/car_proper/bym2 on abundance arm via `nested_laplace` (#51,
  `R/distance_spatial.R` over shared areal-BFGS driver `R/areal_bfgs.R`): distance
  marginal = bin-multinomial (NOT N-mixture, so not the count-spatial driver) but
  exposes analytic per-site gradient (`cpp_distance_site_sweep` over
  `compute_distance_site`), so BFGS + FD-Hessian recovers observed info. **DETECTION-arm
  areal field (#114)**: `detection=~icar()` loads on the per-site detection scale
  eta_sigma (spatially-varying detection scale, `det_arm` in `.tobs_fit_distance_spatial`
  routes offset->eta_sigma, grad_eta=sw$grad_sig); works under BOTH the half-normal key
  AND the hazard key (the field on eta_sigma, the log-shape eta_b threaded as a global
  in the same eval). Pois + NB; NUTS+spatial gated.
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
  (one binary site/node). **Detection-arm areal field (#114)**: `detection=~icar()`
  loads a spatially-varying true-positive logit on the per-site eta_p11 (per-site
  detection design, so no aggregation -- `det_arm` routes offset->eta_p11,
  grad_eta=grad_eta_p11); p10/b never structured. NUTS carries the field on the psi
  arm only (detection-arm -> nested_laplace).
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
  tau[,rho]; one unit/site, Pois+NB. **Detection-arm areal field (#114)**:
  `detection=~icar()` loads a spatially-varying detection logit on the per-site eta_p
  (per-site detection design applied across every season, so no aggregation --
  `det_arm` routes offset->eta_p, grad_eta=grad_eta_p); omega/gamma never structured.
  NUTS carries the field on the initial-abundance arm only (detection-arm ->
  nested_laplace).

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
- **Latent-N truncation is PER SITE, and guarded.** Sum runs on `[max(y_i), K]`:
  lower end always the site's own maximum, only the ceiling shared -> a shared
  ceiling makes every site pay for the largest count ANYWHERE in the data.
  Cost linear in the state count -> one heavy-tailed species-site drags every other
  site's evaluation up by that ratio. `K_max = NULL` (default) caps each site at
  `max(y_i) + 100` (`.nmix_truncation`, `R/nmix_site_marginal.R`) = SAME headroom
  the site holding the global maximum already had; an explicit `K_max` keeps its
  documented meaning (hard global truncation) and is NEVER capped, since a caller
  raising it is compensating for the one regime headroom cannot see. Threaded as a
  `headroom` argument (`< 0` = no cap) through `nmix_precompute_site` -- the single
  place `K_hi` is set, so every caller inherits it. Do NOT give `compute_nmix_site`
  a headroom argument: its arity IS the `CountKernelFn` function-pointer contract
  the shared count-NUTS / count-Laplace drivers take, and a defaulted 7th parameter
  silently breaks that conversion (`abun_nuts.cpp`). Capping callers build the cache
  themselves instead.
  **The cap is verified, NOT assumed.** Exact where the posterior over N decays
  inside its window, wrong where it does not; no fixed headroom survives `p -> 0`
  (never-detected count ~`Poisson(lambda (1-p)^J)`). Guard is NOT the boundary mass
  `nmix_laplace()` warns on -- that bounds the error where the fit STOPPED, and a
  fixture passed it cleanly while sitting well below the uncapped optimum, b/c the
  optimiser's PATH ran through the truncated region. What IS tested: the answer is
  also stationary under the shared ceiling -- score at the fitted coefficients under
  both truncations must agree to `.NMIX_SCORE_TOL` (1e-4). Separates by eight orders
  of magnitude where boundary mass separated by nothing. Failing that -> window
  widens 4x, fit redone, escalating to uncapped -> a guarded fit is NEVER worse than
  the shared-ceiling fit, at worst it costs the extra fits. Verified bit-identical
  (0.000e+00 under the uncapped likelihood) on ordinary counts (guard never fires),
  an unidentified lambda/p ridge (escalates to uncapped), and an identified
  high-abundance fixture (settles at 400). Wired into `nmix_laplace()`,
  `nmix_laplace_re()` + the `ms_abun() + latent()` path (scores at the predictor the
  fit ran on, field + factor offsets included -- a latent surface can push a site
  above what its own counts suggest = the direction that exhausts a window). NB /
  zero-inflated keep the shared ceiling: the check runs a Poisson marginal, would
  understate a heavier tail. Dedicated C++ areal community path (#12) + the removal /
  distance kernels untouched + uncapped -- a capped oracle beside an uncapped field
  solve in ONE fit is worse than neither.
  **The guard is live, not just wired.** `.nmix_community_score_gap()` sits behind a
  `tryCatch(error = NA)` and `is.finite(NA)` is FALSE -> an ERRORING guard and a
  PASSING guard both leave the fit untouched and look identical from outside; a
  capped fit matching its uncapped twin is NOT evidence the check ran. Measured live
  on a fitted `ms_abun() + latent()`: finite on every call, gap rises further when
  the fit's own factor offset is included -- design rationale showing up in the
  measurement, since the latent surface lifts exactly the sites a tight window
  starves. **What it buys** is confined to heavy-tailed seeds and is NOTHING
  where counts are not heavy-tailed, so this is a tier-3 fix, NOT a way to
  shrink `test-ms-abun-factor.R` -- its lever is `factor.starts`. Numbers in
  `NOTES_measurements.md`.
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
- **Zero-inflated** (`abun(mixture="zip"/"zinb")`, #116): structural-zero mixture
  `L_i = omega*1{all y_i=0} + (1-omega)*L_royle_i`. PURE-R additive layer over
  `nmix_site_marginal()` (no kernel change; plain P/NB paths byte-identical),
  `.tobs_fit_nmix_zip` (`R/nmix_zip.R`) L-BFGS-B over theta
  `[beta_lambda|beta_p|logit_omega|(log_r)]`, vcov = inv numeric Hessian. omega
  intercept-only; `logit_omega` a model coord (in `means`/`vcov`/`sds`, surfaced
  by `summary` like `log_r`, NOT the per-process `coef()` list); `fit$zi_omega`
  = omega, `fit$zero_inflated`. Box-bounds only the pathological corners
  (`logit_omega`, `log_r`) so a degenerate seed can't diverge on the ZINB
  zero-source ridge; betas unbounded (no interior bias). ZIP recovers cleanly
  (betas/omega, 95% slope coverage >=0.85, 15 seeds); ZINB weakly identified
  (structural omega vs NB overdispersion both absorb zeros), recovers in a
  higher-count regime (median over 12 seeds). Gated at `.dispatch_abun`
  (non-spatial laplace only; nuts/nested_laplace/structured-term -> pointer, no
  silent drop). `simulate_abun(mixture=, omega=)`. `test-abun-zip.R`. NUTS +
  areal + zero-inflation covariate design = follow-ups (the marginal is the
  additive layer they share).
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
  **Latent-N ceiling PER (species, site) (#233)**: warm-mode-lambda tail
  quantile, an EXCURSION margin re-checked every fit; kernel arity untouched.
  Contract + numbers in `NOTES_measurements.md`.

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
field cor ~0.97); icar/bym2 also sample + centre via the #71 sum-to-zero
shared-field loading (#113); temporal/RE NUTS still gated to n-L.

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

### Diagnostic doors = S3 methods, NEVER a second name

ONE verb per diagnostic, owned by whoever owns the CONCEPT; tulpaObs registers a
`tobs_fit` method + re-exports the generic (`R/reexports.R`), so `library(tulpaObs)`
alone reaches every door and no session masks anything. Roster:
`waic()`/`loo()` = **loo**'s generics (so `loo_compare()` reads a `tobs_fit`;
`loo.tobs_fit` returns a real `psis_loo` via `.tobs_loo_one()`, shared w/
`tobs_stack()`); `sbc()` `pit_residuals()` `test_uniformity()` `test_dispersion()`
`test_outliers()` `test_zero_inflation()` `dic()` `cpo()` `check_model()` =
**tulpa**'s; `ppc()` = tulpaObs's own (no owner elsewhere). Old `tobs_*` spellings
DELETED, not aliased.

NEVER define a generic another package owns (`waic`, `loo`, `pp_check`) --
attaching both masks it and their `waic(matrix)` breaks. Method's first formal must
MATCH the generic's (`loo::waic(x, ...)` -> `waic.tobs_fit(x, ...)`; tulpa's take
`object`).

A method must honour the generic's CONTRACT. `ppc()` did NOT fold into
`tulpa::pp_check()` for that reason: pp_check DRAWS, ppc returns a Bayesian
p-value. Same topic, different return -> different verb.

`check_model.tobs_fit()` (`R/diagnostics.R`, panel `.tobs_check_panel`) DID fold,
by meeting the contract: prints the roll-up AND draws the panel (`plot=FALSE`
prints only). tulpa's `.default` cannot run here at all -- it wants `fitted()` /
`residuals()` as numeric VECTORS and a latent-state fit has `list(psi,p,z)` /
`list(occ,det)`, so it dies resolving the response. Method reads the family's own
doors instead (`pit_residuals`/`ppc`/`test_dispersion`/`moran_i`) and plots ONLY
what the report already computed -> nothing simulated twice (hence
`test_dispersion()` returning `$sim`).

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
| Single-season occupancy | Yes | Yes | parity w/ inlaocc. **`method="pg_gibbs"`** (#126, spOccupancy `PGOcc`): a REAL Polya-Gamma Gibbs chain over the exact posterior via PG data augmentation (Polson-Scott-Windle, tulpa `cpp_rpg`), NOT the stochastic-EM `laplace_gibbs`. z\|. Bernoulli, then omega~PG conjugate-Gaussian beta updates both arms (`R/occu_pg_gibbs.R`, `.tobs_fit_occu_pg_gibbs`). Posterior MATCHES the Laplace observed-Fisher fit (means <1 SE, SDs <20%); Rhat~1.00, ESS 700-900; nominal 95% coverage (15 seeds). Wired via method route table + `.tobs_control_allow(pg_gibbs="sampler")`. **spPGOcc** (`occu() + icar()` + pg_gibbs, `.tobs_fit_occu_pg_gibbs_spatial`): intrinsic ICAR field on the psi logit, jointly updated w/ beta as a GMRF (dense `(p+n)x(p+n)` precision `t(W)Omega W + blkdiag(B0inv, tau Q)`, W=[X\|I]); conjugate Gamma tau; field sum-to-zero'd each sweep, level moved to the intercept (else eta leaks). Field cor ~0.74, intercept/detection recover, Rhat~1.00. NB: site-level occupancy covariate spatially-confounds w/ the saturated 1-node-per-site field (known property, spOccupancy too) -> recovery target = field+intercept+detection. icar only (bym2/car follow-ups). `test-occu-pg-gibbs.R` + `test-occu-pg-gibbs-spatial.R` |
| Dynamic (HMM) | Yes | Yes | colonization/extinction. **Season-varying gamma/epsilon (#124)**: a covariate carried as a `[n_sites x (T-1)]` matrix column of `data` drives interval-indexed colonization/extinction (`colonization = ~ cov`), the dyn_abun #80 recipe ported to the colext forward. Detected via `.tobs_interval_arm_design` (shared w/ dyn_abun); the arm design unrolls long-form `[(site x interval) x p]`. E-step forward-backward uses per-interval transition matrices `A**[i, iv]` (constant rates recycle -> byte-identical to pre-#124); M-step encodes the transition arm as a WEIGHTED LOGISTIC (response = transition prob given origin state `xi/n_mat`, weight = origin-state prob `n_mat`) NOT the per-site M=1000 pseudo-binomial (which, applied per interval, makes each row near-separable so the inner Newton returns ~0 slope). The exact-marginal refine handles SV via a pure-R season-varying HMM-forward marginal `.tobs_dyn_occu_marginal_nlp_sv` (`R/dyn_occu_marginal.R`; the cpp `cpp_occu_dynamic_ploglik` reads one rate/site) -- escapes the EM local optimum AND calibrates SEs. `simulate_dyn_occu(beta_gamma=, beta_epsilon=)`. **Season-varying DETECTION (#124)**: a `[n_sites x T]` matrix column of `data` on the `detection` formula (`detection=~det_cov`) drives per-SEASON detection `logit p_it`. Shares the period-agnostic unroller (`.tobs_interval_arm_design` -> `.tobs_period_arm_design`, w/ a season (T) wrapper `.tobs_season_arm_design`; single source of truth): E-step emission reads per-`(site,season)` `p_mat[i,t]`, M-step (+ `hard_encode`) encodes one detection row per `(site,season)`, marginal refine + gate carry `det_season_varying`. `.tobs_dynamic_smoothed_z` (`fitted()$z`) indexes detection per season + transitions per interval (also FIXED the smoothed state for season-varying colext). `simulate_dyn_occu(beta_det_season=)`. Laplace/gibbs/mi only; NUTS+SV gated (C++ forward per-site). 20-seed recovery + coverage (`test-dyn-occ.R`) |
| Multi-season AR1 year effect (`t_occu`) | — | — | `t_occu()` (#124, spOccupancy `tPGOcc`): NOT colext -- per-`(site,season)` Bernoulli `z_{i,t}~Bern(psi_{i,t})`, `logit psi = X beta + eta_t` w/ a SHARED AR1 year effect `eta_t = rho eta_{t-1} + w_t` on the logit, NO colonization/extinction. Given the year effects the seasons factorise -> EXACT Polya-Gamma Gibbs (`method="pg_gibbs"`, the engine spOccupancy uses; `R/t_occu.R`, `.tobs_fit_t_occu_pg_gibbs`): draw z (Bernoulli), then joint `(beta_occ, eta)` as ONE GMRF update -- dense precision `blkdiag(t(Xf)Omega Xf + B0inv, Z'Omega Z + Q/sigma^2)` w/ stationary AR1 precision `Q` (`.t_occu_ar1_Q`, tridiag) as the year-effect prior (as spPGOcc's ICAR precision is for a spatial field), then `beta_p`, then AR1 hypers (`sigma^2` conjugate Inverse-Gamma, `rho` on a grid over the AR1 log-density). Year effect sum-to-zero'd each sweep, level moved to the intercept; AR1 hypers update on RAW `eta_raw` BEFORE centering (else centering distorts the AR1 increments -> rho collapses). `fit$temporal_field` = posterior-mean year effect. Recovers year-effect surface (cor ~0.97), occ/det coefs, `sigma`, `Rhat`~1.00 (20-seed `test-t-occu.R`). AR1 correlation `rho` WEAKLY IDENTIFIED at few seasons -- even given true year effects `rho_hat` climbs ~0.08 (T=8) -> 0.43 (T=20) -> 0.59 (T=200) for truth 0.6, a property of the short AR1 series NOT the sampler -> reported, not asserted tightly. `y` = 3D `[n_sites x n_seasons x max_visits]` or list of per-season matrices; `simulate_t_occu()`. v1 = site-level occ/det covariates, `pg_gibbs` only (season-varying covariates, areal psi field, NUTS = follow-ups) |
| Community single-season (`ms_occu`) | Yes | Yes | per-arm community RE, shared community Laplace-EM (`R/community_em.R`, `R/ms_occu.R`); msPGOcc. **`method="pg_gibbs"`** (#115/#126, spOccupancy `msPGOcc`): hierarchical Polya-Gamma Gibbs (`R/ms_occu_pg_gibbs.R`) -- per-species PG-augmented conjugate coef updates + conjugate community mean + Inverse-Gamma (near-Jeffreys) community variance draws (diagonal per-arm covariance), tulpa `cpp_rpg`. Gives a CALIBRATED community-VARIANCE posterior: Laplace-EM `sd_psi`/`sd_p` attenuate (documented lower bound), Gibbs recovers them (sd_psi 0.667/0.417 vs truth 0.6/0.4 vs Laplace 0.595/0.342; sd = posterior MEDIAN, robust to the variance skew). Community means + per-species coefs (cor ~0.9) recover; Rhat~1.00. `test-ms-occu-pg-gibbs.R`. NUTS non-spatial samples community means/deviations/covariances jointly (#69, `R/ms_occu_nuts.R`); shared areal field icar/bym2/car_proper on occ arm via n-L (#75, `R/ms_occu_spatial.R`, sfMsNMix analogue). NUTS+field -> n-L. **SVC (`svcMsPGOcc`, #118)**: varying-coefficient bar `spatial(~ 1 + w \|\| cell, graph)` on the occ arm -> intercept + SVC field(s) via BLOCK COORDINATE ascent (`R/ms_occu_field.R`, `.ms_occu_field_solve` = two-state-marginal occupancy field Newton + ICAR prior; community occ EM w/ field as psi offset). Both fields recover ~0.96. icar only; plain intercept field stays on the C++ path (no regression). `test-ms-occu-field.R` |
| Community dynamic (`ms_dyn_occu`) | Yes | Yes | per-species psi1/p RE + shared gamma/eps; HMM-forward; `R/ms_dyn_occu.R`. **NUTS (`method="nuts"`, #115, `R/ms_dyn_occu_nuts.R`, `src/ms_dyn_occu_nuts.cpp`)**: non-spatial community sampler over the exact per-(species,site) HMM-forward marginal via an in-tree C++ FullGradFn -- samples community means, per-species psi1/p deviations, the two independent per-arm community covariances, AND the shared gamma/eps globals jointly. Non-centered `b_{s,arm}=C_arm z_{s,arm}` (dynamic analogue of the `ms_occu` #69 target; data term is the HMM-forward marginal `.ms_dyn_occu_fwd_ll_vec` + two extra SHARED transition arms carrying no RE; per-arm eta scores from `.ms_dyn_occu_fb_vec`, shared verbatim w/ the stMsPGOcc field fitter). Byte-exact vs the R oracle; warm-started at the Laplace-EM mode; 0 divergences; de-attenuates the community variance the EM under-reports. Non-spatial only (a structured term -> nested_laplace). `test-ms-dyn-occu-nuts.R`. **`method="pg_gibbs"`** (#115/#126, spOccupancy `tMsPGOcc`, `R/ms_dyn_occu_pg_gibbs.R`): the msPGOcc community PG machinery + a 2-state HMM forward-filter backward-sample latent step (per-species psi1/p RE w/ community hyperpriors; SHARED gamma/eps from the aggregated 0->/1-> transitions across all species). Calibrated community-variance posterior (Laplace-EM attenuates); shared gamma/eps recover tightly (informed by all species), sd_psi1 recovers, Rhat~1.00. Constant transitions, site-level detection. `test-ms-dyn-occu-pg-gibbs.R`. **Shared areal field on psi1 (stMsPGOcc, #123)**: `~ 1 + icar(graph=adj)` under `nested_laplace` -> `R/ms_dyn_occu_spatial.R`. KEY MATH: psi1 sets ONLY the initial HMM mixing weight, so the per-(species,site) marginal is LINEAR in psi1 -- `L=(1-psi1)A+psi1 B` w/ A/B = HMM-forward likelihood cond on season-1 state 0/1 -- IDENTICAL mixture to the single-season `ms_occu` field oracle, so the Louis score/curv (`score=r-psi1`, `curv=psi1(1-psi1)-r(1-r)`, `r=psi1 B/L`) carry over verbatim (only A/B differ). Routed through the shared block-coordinate driver (`R/community_latent.R`, `.tobs_community_latent_ascent`): community EM (field=psi1 offset) alternated w/ areal Newton. Kernels VECTORIZED over sites (`.ms_dyn_occu_condAB_vec`/`.ms_dyn_occu_fb_vec`) + ANALYTIC Fisher-identity gradient (fwd-bwd smoothed w1/pairwise xi -> grad_psi1=X'(w1-psi1), grad_p=X'(sum_t w_t(ndet-nvalid p)), grad_gamma=X'(col_y-gamma col_n)) so the community EM skips its O(U^2) FD Hessian (a per-site-loop version was ~10x too slow). icar only; field cor ~0.94, community-mean coverage >=0.85. `simulate_ms_dyn_occu(field=)`. **svcTMsPGOcc (#123)**: a weighted areal bar `spatial(~ 1 + w || cell, graph)` adds a shared covariate-weighted field alongside the intercept field -- NO new code, the psi1 oracle already returns per-site/species score+curv and the K-field weighted-ICAR block solve is the SAME `community_latent.R` machinery as svcMsAbund, so it flows through the existing dynamic-spatial dispatch unchanged. Both fields recover (cor ~0.90/0.89). `simulate_ms_dyn_occu(trend=)`. NUTS/bym2/car follow-ups |
| Community integrated (`ms_int_occu`) | Yes | Yes | per-species psi + per-source det RE; multi-source two-state marginal; `R/ms_int_occu.R`. **NUTS (`method="nuts"`, #115, `R/ms_int_occu_nuts.R`, `src/ms_int_occu_nuts.cpp`)**: multi-source generalisation of the `ms_occu`/`ms_dyn_occu` community samplers -- samples community means, per-species occupancy/per-source detection deviations, and the D+1 independent per-arm community covariances jointly over the exact multi-source two-state per-(species,site) marginal via in-tree C++ FullGradFn (per-arm non-centered `b=Cz`, NO shared globals), warm-started at the Laplace-EM mode; reuses shared `.ms_ocs_*` epilogue (#128). Byte-exact vs R oracle; 0 divergences; de-attenuates the community variance the Laplace-EM under-reports. Non-spatial only. `test-ms-int-occu-nuts.R`. **`method="pg_gibbs"`** (#115/#126, `R/ms_int_occu_pg_gibbs.R`): msPGOcc generalized to D detection arms -- per species draw the single latent z (occupied if any source detects, else Bernoulli on the pooled occupied-undetected mass), PG-conjugate `beta_psi_s` + D per-source `beta_pd_s` (each at that species' occupied covered sites), conjugate community mean + IG variance per arm. Calibrated community-VARIANCE posterior (sd_psi ~0.53 vs truth 0.5; Laplace-EM attenuates), community means recover, Rhat<1.1. `ms_community` layout matches the Laplace fit (`Sigma_/sd_/coef_/blup_<arm>`). `test-ms-int-occu-pg-gibbs.R` |
| Integrated multi-source | Yes | Yes | shared psi |
| Multi-season integrated occupancy | Yes | — | `dyn_int_occu()` (#122, spOccupancy `tIntPGOcc`): product of dynamic occupancy (multi-season HMM: psi1/gamma/eps) + integrated occupancy (per-season emission pooling S detection sources). Pure-R two-state HMM forward w/ multi-source pooled emission (`R/dyn_int_occu.R`, the `dyn_occu_marginal.R` forward generalised), optim BFGS + observed-Fisher vcov, no new C++. `y`=list of S `[sites x visits x seasons]` arrays; `colonization=~`/`extinction=~` required (as dyn_occu); shared per-source detection design (`p_<src>` arms). Recovers psi1/gamma/eps + per-source detection (20 seeds + pooled 95% coverage). **Partial season overlap**: source absent at a (site,season) -> NA -> contributes nothing to that season's emission (nvalid=0); a (site,season) unobserved by EVERY source is marginalised (e0=e1=1) by the forward -> staggered surveys NA-padded to the common `[n_sites x max_visits x T]` grid (`simulate_dyn_int_occu(source_seasons=list(1:4, 3:6))`, 20-seed recovery). Two boundary ANCHORS pin the reduction: one source (others all-NA) reproduces `dyn_occu()`; one season of data (`T=2`, season 2 all-NA) reproduces `int_occu()` (both to ~0.001). `fit$means` names `psi1_*`/`gamma_*`/`eps_*`/`p_<src>_*`. **stIntPGOcc**: `~ 1 + icar(graph)` on the occupancy formula loads a shared areal field on psi1 via the areal-BFGS driver (`.tobs_fit_dyn_int_occu_spatial`); field gradient = Fisher-identity psi1 score `w1 - psi1` (`.dio_fb`, FD-validated; also sped the non-spatial fit); field cor ~0.8. **svcTIntPGOcc (#122)**: a `spatial(~ 1 + w || cell, graph)` bar adds a covariate-weighted (SVC) field alongside the intercept field -- the areal-BFGS driver already takes a LIST of field blocks and scatters `w1 - psi1` to each, so the weighted block = ICAR field w/ a `w`-weighted loading (`.areal_field_car(weight=)`, byte-identical unweighted). `fit$spatial_field`=intercept, `fit$trend_field(s)`=SVC surface; both recover by cor (~0.85 / ~0.72). `simulate_dyn_int_occu(field=, trend=)`. icar only. v1 = constant transitions, site-level detection; season-varying rates (the #124 recipe) + bym2/car_proper + NUTS = follow-ups. Full S3 + WAIC. `test-dyn-int-occu.R` + `test-dyn-int-occu-areal-recovery.R` |
| JSDM (`jsdm`) | Yes | Yes | `jsdm()` = COMMUNITY GLMM on observed presence/absence (#121): per-species coefs + Gaussian community covariance, NO detection/latent state = `ms_count()` w/ logit link -> SHARES that binder (`.tobs_build_ms_count(response="bernoulli")`, model_type `ms_count`), community EM, latent driver, NUTS target + S3. `latent(n)` = lfJSDM; + shared field = sfJSDM; icar/car_proper/bym2/spde field via n-L. NUTS = exact joint community posterior over the Bernoulli response (`MSC_BERN` in `src/ms_count_nuts.cpp`, byte-exact vs R oracle). **`method="pg_gibbs"`** (#126): hierarchical Polya-Gamma Gibbs -- community logistic GLMM has NO latent state -> pure per-species conjugate coef update + community mean + near-Jeffreys Inverse-Gamma community variance (`R/ms_count_pg_gibbs.R`, `.tobs_fit_ms_count_pg_gibbs`, tulpa `cpp_rpg`), shared w/ `ms_count("binomial")` (n=trials, jsdm n=1). Recovers the community variance the Laplace-EM attenuates (sd_mu 0.674/0.439 vs truth 0.7/0.5 vs Laplace 0.635/0.401); community means + per-species cor (~0.9) recover; Rhat~1.00. Non-spatial (bernoulli/binomial only; poisson/negbin/gaussian rejected w/ pointer). `test-ms-count-pg-gibbs.R`. laplace_sla/gibbs/mi were single-block routes, dropped |
| Cover hurdle (joint) | Yes | Yes | `family_cover_hurdle.R`, `sla_cover_*`, NUTS `R/cover_nuts.R` + `src/cover_nuts.cpp` (no `laplace_gibbs`/`laplace_mi`). positive = beta/lognormal/lognormal_trunc/ordinal/`beta_oi`. `beta_oi` (#108) = one-inflated Beta: ceiling (cover=1) plots = a point mass (constant pi = ceiling share, binomial SE), interior Beta on (0,1); encode splits `is_pos` to interior, `enc$oi` carries pi, decode reports `pi_one`, predict conditional cover = `pi + (1-pi)*mu` (`.tobs_cover_mu`). `control$aggregate.occ` (ON, #48) collapses occ arm to Binomial suff-stat; `control$aggregate.pos` (ON beta arm, #49) collapses beta pos arm to grouped suff-stat (tulpa `slog_y`/`slog_1my`), errors on non-beta. Both byte-identical to per-plot |
| Cover hurdle spatial coef fields (`\|\|` / `\|`) | n-L | — | `spatial(~ 1 + w \|\| node, graph, to=)` independent (#61, two coupled ICAR blocks, per-field alpha) OR `\| ` correlated (#64, one separable-MCAR block sharing free Sigma). `\|` both-arm `to=c("presence","positive")` = copied to pos arm w/ one alpha (#64); `\|` single-arm `to="presence"`/`"positive"` = free-Sigma field on that arm alone, NO copy (#109, 0-sentinel `spatial_idx` on the other arm via `mc$to`, `copy=NULL`, `alpha_mcar`=NA). `\|` -> `.cover_build_mcar_spec`/`.fit_cover_hurdle_joint_mcar` (tulpa `type="mcar"` block, copy only when both-arm); reports `sigma_mcar`/`rho_mcar`/`alpha_mcar`. SLA on `\|` no-op. icar only |
| Cover hurdle arm-specific fields (single-arm `to`) | n-L | — | `spatial(~ 1 + w \|\| cell, graph, to="positive")` (or `"presence"`); separate single-arm calls = independent per-arm fields, NO cross-arm copy (#65). NO engine change: per-arm `spatial_idx=0` makes the other arm's rows skip the block (tulpa `l_b>0` scatter guard), own precision grid-integrated. `.tobs_armspecific_bar_fields` (formula_terms.R) -> `enc$armspec` -> `.fit_cover_hurdle_joint_armspecific` (non-copied per-arm blocks, no `copy=`). `armspec_blocks` carries per-block arm/slot/type; `.tobs_joint_draws_cover_armspecific` scatters each block onto its arm only (amp 0 on other). icar/car/car_proper AND bym2 (#107): a bym2 block is the non-copied length-2 (phi ICAR + iid theta) block, paired (sigma,rho) grid; the draw projection reconstructs the rho-mixed unit field `z = sqrt(rho)*sf*phi + sqrt(1-rho)*theta` so predict/WAIC see the full mix. `\|\|` only (`\|` arm-specific undefined: copy-only). No mix w/ shared field/trend/temporal/re; one field per arm. SLA no-op |
| Joint occu + cover | Yes | Yes | `occu_cover()` — see below. NUTS non-spatial (in-tree FullGradFn over exact two-state marginal, beta/lognormal; `R/occu_cover_nuts.R`, `src/occu_cover_nuts.cpp`) AND spatial coupled areal field (#74, car_proper recovers; icar/bym2 sample + centre via the #71 sum-to-zero coupled field, #113) with its HYPERS SAMPLED (#204: field SD, bym2/car mixing rho, copy amplitude alpha — fixed basis, no per-step decomposition; `fit$nuts$sampled_hyper`/`$fixed_hyper` report per fit). Community sampled-field (per-species loadings) route = `ms_occu_cover()` factor |
| occu + cover + areal field + per-group RE | n-L | — | `occu_cover()` + icar/bym2 + `re(g)`/`(1\|g)` on psi; one iid RE block (#56); `sigma_re` + BLUPs; intercept RE only |
| occu + cover independent cover-arm field (single-arm `to="positive"`) | n-L | — | `occu_cover()` + `spatial(~ 1 + w \|\| cell, graph, to="positive")` on occurrence formula (#110): NON-copied ICAR block(s) on the cover arm ALONE, decoupled from the occupancy field's alpha copy. Composes w/ the shared occupancy field (psi + `delta_cover_exp` keep it) -> `delta_cover_cond` varies instead of collapsing when `alpha->0`. Parse: `.occu_cover_spatial_fields` splits the single-arm `to="positive"` bar -> `spatial_info$pos_armspec` via `.tobs_armspecific_bar_fields`. Fitter appends non-copied ICAR blocks (`spatial_idx` psi=0/p=0/pos=cell, `tau_grid`, svc_weight for trend) after occ fields, before RE blocks; `ctx$field_specs` labels each block shared-vs-pos. Postprocess partitions occ (1..n_occ_fields) vs pos blocks; reports `sigma_pos_field`/`sigma_pos_field_<col>` from `b<k>.tau`; surfaces `fit$pos_field`/`pos_field_table(s)`. Draw substrate (`.tobs_joint_draws_occu_cover`) reads `field_specs`: pos block amp_occ=0, amp_pos=1/sqrt(tau) -> predict/WAIC add it to `field_pos` automatically. Predict weight override skips pos blocks. Per-visit cover (`cover_aggregate="none"`); icar only (bym2/car->icar); NOT w/ MCAR `\|`, latent cover RE, or batch. Static intercept field weakly ID'd vs alpha copy; time-weighted trend cleanly ID'd. `test-occu-cover-pos-field.R` |
| occu + cover independent detection-arm field (single-arm `to="detection"`) | n-L | — | `occu_cover()` + `spatial(~ 0 + time \|\| cell, graph, to="detection")` on the `detection` formula (or lifted `to="detection"`) (tulpa#140): spatially-structured detection prob. Same non-copied arm-specific block machinery as the cover-arm field (#110) -- builder `arm_field_blocks(af, "p")` sets slot 2 (`spatial_idx` psi=0/p=cell/pos=0). Enters via detection arm `field_coef=1` (set when `det_armspec` present, `.occu_cover_build_joint_arms(det_field=)`) so the block scatters onto the p rows; shared occ field kept off detection by its `spatial_idx=0` sentinel -- identical to the #102 detection-RE mechanism, so NO tulpa engine change. Reports `sigma_p_field`/`sigma_p_field_<col>`. icar only; per-visit cover. `test-occu-cover-pos-field.R` |
| occu + cover obs-arm RE (detection / pos) | n-L | Yes | `occu_cover()` + RE on `detection=`/`positive=`; per-visit grouping. Intercepts: `(1\|g)`, crossed `(1\|g)+(1\|h)`, nested `(1\|g/h)` = N iid blocks (#102 single, #103 crossed/nested). Slopes (#103, tulpa>=0.0.39): uncorr `(x\|\|g)`/`(0+x\|g)` = weighted iid per coef; corr `(1+x\|g)` = miid free-Sigma block. `sigma_re_p`/`sigma_re_pos` (+`_<coef>`, `cor_re_p_*`) + BLUPs. **NUTS = INTERCEPTS, non-spatial (#205)**: each grouping = one shared `ReBlock` (`src/nuts_re_block.h`, reused not re-derived) trailing the field block, non-centered `b_g=sigma_re*z_g`, group SD SAMPLED (prior `N(0,1.5^2)` on `log sigma_re`, the observation-family width -- written ON the sampled coordinate, so NO change-of-variables term; asserted by reducing the target to the prior at `z=0`). Row-indexed, NOT site-indexed: one `re_group` code per `(site,visit)` row, 0 = padded/unseen level (engine's scatter sentinel). Designs from the SAME `.occu_cover_obs_re_parse`/`_design` as n-L (`.occu_cover_attach_obs_re`, one dispatch-level call site). `fit$re` keys + `ranef()` match n-L; adds `sigma_median`/`sigma_draws`, `fit$nuts$re_sigma_rhat`. `fit$draws`/`vcov` stay the coefficient block (WAIC reads it -> scores RE at 0, as the sampled field already does). Off-block byte-identity `expect_identical`. Gated -> n-L: slopes, RE+areal field. `test-occu-cover-nuts-re.R` |
| occu + cover + detection / cover-arm RE | n-L | — | `occu_cover()` + `(1\|g)`/`re(g)` on `detection=` or `positive=` (#102); per-VISIT grouping iid block on that arm, composes w/ psi field; `sigma_re_p`/`sigma_re_pos` + BLUPs; `fit$re` per-arm list; `predict(type="detection")` adds BLUP offset, unseen->pop mean. Detection arm `field_coef=1` (not 0) so the iid block scatters, field still skipped via `spatial_idx=0` sentinel. pos RE needs `cover_aggregate="none"`. Intercept only; slope/correlated/non-spatial/NUTS gated |
| Community joint occu + cover | Yes | Yes | `ms_occu_cover()` — see below. positive = beta/lognormal/gaussian (`gaussian` #127, delta-normal `mu=eta`). **NUTS (`method="nuts"`, non-spatial, #115 B7, `R/ms_occu_cover_nuts.R`, `src/ms_occu_cover_nuts.cpp`)**: joint-cover analogue of the `ms_occu`/`ms_int_occu` community samplers — three non-centered per-species arms (occ/p/pos), each w/ a log-Cholesky community covariance + ONE shared community log-dispersion scalar, over the per-(species,cell) two-state occu_cover marginal via in-tree C++ FullGradFn (WRAPS existing `occu_coupling_shared.h` cover-density kernels + `community_chol.h`, reuses `.ms_ocs_*` #128), warm-started at the Laplace-EM mode. Byte-exact vs R oracle; 0 divergences; REMOVES the Laplace community-variance attenuation (`var_attenuation` caveat) — sampled per-arm SD de-attenuates toward truth where the EM collapses. **Per-species dispersion RE** via `control=list(dispersion.re=TRUE)` (#115 B7): shared cover log-dispersion becomes a 4th 1-D community arm `log_disp_s=mu_ld+sigma_ld*z_ld_s` (ms_abun log_r_s analogue), byte-exact vs the RE oracle, shared-disp path byte-identical (extracted per-cell `msoc_cell_sweep`); `fit$ms_dispersion$sigma_log_disp`/`log_disp_species`. Non-spatial lognormal/beta/gaussian; negbin N/A (continuous cover arm). `test-ms-occu-cover-nuts.R`. Shared-field spatial-factor variant (a field on the occ arm) stays lognormal-only |
| Spatial-factor community occu+cover (JSDM) | Yes | Yes | `ms_occu_cover()` + icar/car_proper/bym2 shared field, per-species loadings (tulpa#67). Laplace-EM (`R/ms_occu_cover_spatial.R`) + NUTS (`src/ms_occu_cover_spatial_nuts.cpp`). Cover-arm factor, `tobs_associations()`, per-species `predict()` maps |
| Multiscale occu + cover | Yes | Yes | `occu_multiscale_cover()` — 3-level cell/plot/visit + cover. THREE engines: `nested_laplace` (shared + SVC-trend coupled field) and `laplace` + `nuts` (both non-spatial, iid cells, field fixed at 0; `R/occu_multiscale_cover_nuts.R`, `src/occu_multiscale_cover_nuts.cpp`). Spatial = joint only. positive = beta/lognormal/gaussian (`gaussian` #127, threaded through all three engines via the shared `pos_*` code dispatch + `OccuMultiscaleCoverGaussianCoupling`) |
| Count / relative-abundance GLMM | Yes | — | `count(response=)` (spAbundance `abund`): GLMM on the observed response directly, NO detection, NO latent state = abundance analogue of `jsdm()`. Pois/negbin (log link) + gaussian (identity) + **binomial (logit, `k`-of-`n`, `trials=`; spOccupancy `svcPGBinom`, #125)**. ONE tulpa GLMM block (`build_count_callbacks`, `R/laplace_callbacks.R`); negbin size / gaussian variance by an outer dispersion loop in `.dispatch_count` (tulpa_laplace takes fixed phi), reported `fit$count_dispersion`; binomial has NO dispersion (variance pinned by `n`) -> takes the Poisson no-loop path + carries per-site `model$n_trials`. `simulate_count()`, `.tobs_ploglik_count` (WAIC), `count_methods.R` (fitted = expected successes `n*p`, predict newdata = per-trial prob, residuals binomial deviance/pearson). Community (msAbund) binomial = `ms_count("binomial")`; NUTS pending. Non-spatial binomial matches `glm()` MLE to ~1e-3. `test-count.R` (3x 20-seed recovery + binomial trials>1 & trials=1) |
| Community spatial-factor (`ms_count`) | n-L | — | `ms_count()` + `icar()` + `latent(n)` (spatial-factor `sfMsAbund`; #117): shared areal field AND latent factors on ONE formula, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}+sum_q lambda_{s,q} eta_{q,i}`. ONE block-coordinate loop runs BOTH field + factor updates; when both present the factor loadings are CENTERED across species (`sum_s lambda_{s,q}=0`) -> field owns the shared spatial mean, factors own between-species residual co-occ; both recover (field cor ~0.99, residual cor ~0.99). UNIFIED fitter `.tobs_fit_ms_count_latent` (field-only/factor-only/both = one source of truth; old fitters = thin wrappers). `test-ms-count-factor.R` |
| Community latent factor (`ms_count`) | Yes | — | `ms_count()` + `latent(n)` (spAbundance `lfMsAbund`; #117): residual species co-occurrence via Q per-site latent factors + per-species loadings, `log mu_{s,i}=X_i(mu+b_s)+sum_q lambda_{s,q} eta_{q,i}`. BLOCK COORDINATE via the shared driver (`R/community_latent.R`): community EM w/ factor offset <-> factor update (alternating Newton eta\|lambda, lambda\|eta + unit-variance anchor). Loadings/factors identified up to ROTATION -> recoverable target = residual cov `Sigma_res=lambda lambda'` (`fit$ms_factor$residual_cov`/`residual_cor`/`loadings`/`factors`); residual cor recovers ~0.95. `fitted()`/WAIC factor-aware via `model$count_factor_offset`. Poisson; not composed w/ areal field yet. `test-ms-count-factor.R` |
| Community SVC (`ms_count`) | n-L | — | `ms_count()` + `spatial(~ 1 + w \|\| cell, graph)` (spAbundance `svcMsAbund`; #117/#118): intercept field + one shared spatially-varying-coefficient field per covariate, `log mu_{s,i}=X_i(mu+b_s)+f0+w_i*f1`. Same block-coordinate scheme; the field solve generalizes to K covariate-weighted ICAR fields solved jointly (K x K sparse block system, per-field tau). `fit$spatial_field`=intercept field, `fit$trend_field(s)`=SVC field(s). Both recover cor ~0.98. icar only. `test-ms-count-spatial.R` |
| Community count + areal (`ms_count`) | n-L | — | `ms_count()` + `icar()` (spAbundance `sfMsAbund`; #117): ONE shared areal field across species, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}`. Fit by BLOCK COORDINATE ASCENT (`R/ms_count_spatial.R`) -- community Laplace-EM w/ the field as a per-site OFFSET (captured in sp_ll closure, NO community_em.R change), alternated w/ a self-contained Poisson-ICAR field update (analytic sparse Newton `.ms_count_field_solve` + closed-form tau M-step). PURE R, no C++ (sidesteps the community EM FD-Hessian, which does not scale to an O(n_sites) field). Field informed by all species/site -> recovers cleanly (cor ~0.98); `fitted()`/WAIC field-aware via `model$count_field_offset`. Poisson + icar ONLY (overdispersed community count NOT identified vs per-site field; bym2/car/group_var/bar follow-ups). `test-ms-count-spatial.R` (community-mean + pooled coverage + field recovery, 20 seeds) |
| Community count (`ms_count`) | Yes | Yes | `ms_count(response=)` (spAbundance `msAbund`; #117): community relative-abundance GLMM -- per-species GLMM on the observed response w/ Gaussian community hyperpriors on the coefficients, NO detection/latent state = community analogue of `count()`, abundance analogue of `ms_occu()`. Pois/negbin (per-species dispersion RE = second arm)/gaussian (per-species resid var, outer loop)/**binomial (logit `k`-of-`n`, `trials=`; community `svcPGBinom`, #125 -- laplace non-spatial ONLY, NUTS + shared field/latent gated w/ pointer; community-mean intercept carries a small O(1/n_species) first-order-Laplace bias, absent at `trials=1`, slope unbiased)**. Reuses shared community Laplace-EM (`.tobs_community_em`, `R/community_em.R`) -- PURE R, no C++ -- w/ count `sp_ll`/`sp_grad` (`R/ms_count.R`). `y` = `[n_sites x n_species]` or named list. `coef`=community means, `ranef`=per-species deviations, `simulate_ms_count()`, `.tobs_ploglik_ms_count` (WAIC). Gaussian/Pois community means unbiased; negbin slope carries mild first-order-Laplace attenuation (~10%, documented). **NUTS (Pois/negbin/gaussian, #117)**: samples the exact joint posterior via a family-aware in-tree C++ FullGradFn (`src/ms_count_nuts.cpp`, `R/ms_count_nuts.R`) = reduced `ms_abun` NUTS (no detection/latent-N); non-centered `b_{s,arm}=C_arm z_{s,arm}`; byte-exact vs R oracle, warm-started at the Laplace-EM mode, 0 divergences, NUTS==Laplace. negbin adds a per-species dispersion RE `log_r_s~N(mu_log_r,sigma_log_r^2)` as a SECOND community arm (block-diag chol, mirrors `ms_abun`); gaussian adds S FREE per-species `log_phi_s` (no community prior, weakly-informative N(pooled log-var, 2^2), matching the Laplace outer loop) -- NUTS removes the negbin Laplace attenuation, `fit$ms_dispersion` carries `r_s`/`variance`. Missing (`NA`) site x species entries are dropped from the per-(species,site) data sum (NA-in-`y` IS the mask; both C++ + R oracle skip identically -> byte-exact w/ NA present), matching the Laplace path's per-species `valid` subsets. `test-ms-count.R` + `test-ms-count-nuts.R` |
| Count + areal / SPDE / GP | n-L | — | `count()` + `icar()`/`car_proper()`/`bym2()`/`spde()`/`gp()` (spAbund; #117, bym2+spde #116, gp #116): plain areal OR continuous-mesh/NNGP field on the abundance formula. **bym2 (#116)**: generic nested-Laplace field summary `.tobs_nested_attach_field_summary` (R/laplace.R) gained a bym2 branch reading the `b<n>.sigma`/`b<n>.rho` grid axes + the (phi\|theta) mode slices to reconstruct `z=sqrt(rho/scale)*phi+sqrt(1-rho)*theta` (Riebler 2016), grid-averaged; field cor ~0.93. SAME summary drives occu() etc. -> a bym2 areal field there is reconstructed too (was dropped/NULL). icar/car_proper byte-identical. **spde (#116)**: `spde(lon,lat,...)` fits the same nested-Laplace path -- the Matern field lives on the mesh nodes (`fit$spatial_field`, length `n_mesh`), reconstructed by the icar-style else-branch of the field summary, and the barycentric projector `fit$spatial$tulpa_spec$A` (`n_sites x n_mesh`) maps it to sites. `.count_spatial_field_offset` (`R/count_spatial.R`) gained an A-projection branch so `fitted()`/`predict()` add the projected per-site field; areal fits carry no projector (`tulpa_spec$A` NULL) so that branch is skipped, byte-identical. Field cor ~0.93, slope recovers (`test-count-spatial.R`, gated `skip_if_no_tulpamesh()`). **gp (#116)**: `gp(lon,lat,prior_range=c(r0,alpha))` routes to tulpa's single-block `nngp` nested-Laplace kernel via `.tobs_fit_count_gp` (over `tulpa::spatial_gp`+`prior_from_spec`): the GP marginal variance + range are integrated on the kernel's OWN outer grid and the field Schur-folded out, so the fit reports grid-integrated FE + the GP hyper posterior (`fit$gp_hyper`) but NO per-cell map (`fit$spatial_field=NULL`; use spde() for a reconstructed field). `multiscale_gp()` (two-scale) NOT hosted by nested-Laplace -> errors w/ pointer to gp()/spde(). count dispatch gate `ftypes` includes `gp`. Response observed (no latent state) -> fit = ONE `tulpa_nested_laplace()` call over the count block w/ the areal field as its GMRF prior, NOT the occupancy EM. Grid-integrated FE (law of total cov via per-cell `keep_grid_hessians`) + field/sigma via shared `.tobs_nested_attach_field_summary`; dedicated `.tobs_fit_count_spatial` (`R/count_spatial.R`), reuses `.tobs_to_multi_block_prior`. `fitted()` field-aware in-sample. **Poisson OR binomial** -- Poisson has no dispersion; binomial variance pinned by `n` -> identified against a per-node field (spOccupancy `svcPGBinom`, #125; works at `trials=1` too -- one Bernoulli/node weak but NOT confounded). negbin size / gaussian residual variance NOT jointly identified w/ the field under the fixed-phi loop (degenerate: size->Inf, resid var->0) -> gated w/ pointer; bym2 + improper car() gated to icar/car_proper. `test-count-spatial.R` (20-seed FE coverage + field recovery, car_proper, gates, + binomial trials>1 & Bernoulli field recovery) |
| N-mixture (`abun`) | Yes | Yes | `abun()` = laplace + nuts (non-spatial, #41) + nested_laplace (areal). Royle 2004; see the detail section below |
| N-mixture + areal | n-L | Yes | `abun()`+icar/bym2/car_proper; Pois/NB (r grid-int); grid-int cov (constrained intercept). NUTS+areal = car_proper ONLY (#51, fixed-hyper non-centered field); icar/bym2 NUTS+spatial gated |
| Community N-mixture | Yes | Yes | `ms_abun()` (msNMix); per-species coef RE, in-tree Laplace-EM (`nmix_laplace_re`) -> `NMixCommunityOracle` AGHQ, Schur SEs; Pois + negbin. NUTS (#14) |
| Community N-mixture + areal | n-L | Yes | `ms_abun()`+icar/bym2/car_proper (sfMsNMix; #12); shared field + per-species RE; `nmix_community_spatial.cpp`; Pois/NB. NUTS (#73, car_proper only) = #14 non-centered community sampler + a SHARED fixed-hyper non-centered proper-CAR field on abundance (tau Q(rho) fixed at the #12 nested-Laplace estimate, raw ~ N(0,I), f=Linv raw; optional field block in `src/ms_abun_nuts.cpp`, FD-validated, field-off byte-identical to #14). 0 divergences, field cor ~0.97; Pois only. icar/bym2 also sample + centre via the #71 sum-to-zero shared-field loading (#113) |
| Community N-mixture latent factor (`ms_abun`) | Yes | — | `ms_abun()` + `latent(n)` (spAbundance `lfMsNMix`; #117): residual species co-occurrence on the ABUNDANCE arm via Q per-site factors + per-species loadings. Latent N still marginalises closed-form per species-site -> whole latent structure sits on `eta_lambda`, family reduces to ONE working oracle over the Royle marginal (`R/ms_abun_latent.R`) driving the shared block-coordinate engine: `score = grad_eta_lambda`, `curv = info_eta_lambda - var_N * score_wt_lambda^2` = Louis (1982) (1,1) block (abundance curvature, detection arm profiled out), which `nmix_site_marginal()` ALREADY exposes -> no new kernel. Oracle FD-validated. Residual cor recovers ~0.99 at N=120/S=12 (counts information-rich, cleanly better than occupancy's ~0.9). Poisson only (a negbin size = second per-site dispersion, not identified vs a per-site latent). `test-ms-abun-factor.R` |
| Community N-mixture spatial-factor (`ms_abun`) | n-L | — | `ms_abun()` + `icar()`/`car_proper()`/`bym2()`/`spde()` + `latent(n)` (#117): shared field AND factors on the abundance arm; the centred loadings separate them (field cor ~0.98, residual cor ~0.94). Same block-coordinate driver. A plain field with NO factors KEEPS the dedicated C++ path (#12) -- faster + already recovery-tested, deliberately not replaced |
| N-mixture + grouped RE | Yes | — | `abun()`+`(1\|g)`/`(x\|g)` either arm (#13); non-species grouping; Pois/NB; `NMixGroupedOracle`. Gated: RE+spatial, RE+visit-det, RE both arms |
| Community distance (`ms_distance`) | Yes | — | `ms_distance(key=, transect=, cutpoints=)` (spAbundance `msDS`; #117): per-species binned distance sampling w/ Gaussian community hyperpriors on the abundance (`log lambda`) + detection-scale (`log sigma`) coefs. `y` = `[n_sites x n_bins x n_species]` or named list. NO new C++ — latent N marginalises closed-form per species-site, and `cpp_distance_site_sweep` already returns `log_lik`/`grad_lam`/`info_lam`/`var_N`/`swl` -> shared community Laplace-EM (`R/ms_distance.R`) reads its per-species score straight off it. Hazard-key log-shape = community `global` (shared across species) via the EM `global` slot. Community means UNBIASED over 10 seeds, 95% Wald coverage 1.0/0.8/0.9; single-species `distance()` control agrees. NB: a single seed shows a ~0.18 paired lambda-down/sigma-up shift that LOOKS like bias and is NOT -- one draw on the lambda/sigma ridge. Poisson only (negbin size not yet a per-species RE); NUTS not wired. `simulate_ms_distance()` draws through `cpp_simulate_distance` (a separate R-side quadrature would simulate from a pi the model is not fit against). `test-ms-distance.R` |
| Community distance latent factor (`ms_distance`) | Yes | — | `ms_distance()` + `latent(n)` (`lfMsDS`, #117). Factors ALONE = plain block-coordinate Laplace-EM, so `method="laplace"`; `.dispatch_ms_distance` REJECTS `nested_laplace` without a field ("needs a shared field ... use method = \"laplace\""). Same shared driver + oracle as the spatial-factor row below |
| Community distance spatial-factor (`ms_distance`) | n-L | — | `ms_distance()` + `icar()`/`car_proper()`/`bym2()`/`spde()` (`sfMsDS`, #117), optionally + `latent(n)`. Shared field REQUIRES `method="nested_laplace"` (`.dispatch_ms_distance` errors under `laplace`) — same gate as sibling `ms_count`/`ms_abun` spatial-factor rows. Field on the ABUNDANCE arm only; a `spatial()`/`latent()` term on the detection formula errors. Working oracle = SAME Louis formula as `ms_abun` -- `score = grad_lam`, `curv = info_lam - var_N * swl^2` (both count marginals w/ `B_i = diag(info) - var_N v v'`; only the kernel supplying the pieces differs). Oracle FD-validated. Poisson only; temporal/re/svc not wired. **Per-species Hessian IS assembled now (#161)**: used to be left to the community EM finite difference at `2U` sweeps per species per Newton step, every sweep summing over the latent N -> `test-ms-distance.R` the most expensive file in the recovery tier (5.12h at tier 3 = 88% of `full-recovery.yaml`'s 350-min cap on its own). `.tobs_ms_distance_info_block()` builds the per-site Louis block `B_i = diag(info_lam, info_sig) - Var(N_i|y) v v'` from what `cpp_distance_site_sweep` already returns (`info_lam`, `info_sig_obs`, `var_N`, `swl`, `vN_sig`); distance has ONE unit per site -> three `crossprod`s against diagonal weights, no per-site loop (N-mixture needs one only b/c its visit count varies). Community EM **2.3x** faster, fit UNCHANGED to ~1e-13. **Sign inside `v` is NOT the sibling's**: `v = (-swl, -vN_sig)` -- distance kernel stores `vN_d` already negated (`-p_k/(1-p)`) where `nmix_site_marginal()`'s `v` carries `+p`. Invisible on the diagonal, wrong by ~2x on the cross block -> asserted against the finite difference, NOT reasoned across (`test-ms-distance-info.R`, which also asserts the cross block >10% of the diagonal so the agreement cannot be vacuous, and that the flipped convention is REJECTED). Half-normal only: under the hazard key `grad_b`/`info_b` come back already summed over sites + per-site detection cross terms are not exported -> shared log-shape global cannot be sandwiched, that key keeps the FD fallback |
| Removal (Pois/NB) | Yes | Yes | `removal()` (#39) — see Architecture. **continuous NNGP `svc()` varying coefficients on the abundance arm under laplace/nested_laplace, composing with the areal/temporal blocks (#144)**; NUTS single intercept RE; Laplace grouped RE one arm; areal icar/car_proper/bym2 abundance arm; **detection-arm areal field on the capture logit (`detection=~icar()`, #114)**; areal+temporal AND temporal-only AR1/RW1/RW2/iid on abundance arm via shared areal-BFGS (#78/#114); NUTS+areal icar/car_proper/bym2 field on abundance arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| Distance (Pois/NB) | Yes | Yes | `distance(key=, transect=, cutpoints=)` (#38); `formula`=log lambda, `detection`=log sigma, `y`=`n_sites x n_bins` — see Architecture. **continuous NNGP `svc()` varying coefficients on the abundance arm under laplace/nested_laplace (#144)**; NUTS single abundance intercept RE; Laplace grouped RE abundance arm (half-normal AND hazard key -- hazard log-shape profiled over the AGHQ log-marginal, #114); areal field (half-normal + hazard key); areal+temporal AND temporal-only on abundance arm (#78/#114); DETECTION-arm areal field on log sigma (`detection=~icar()`, spatially-varying detection scale, half-normal AND hazard key, #114); NUTS+areal icar/car_proper/bym2 field on abundance arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| False-positive occupancy | Yes | Yes | `fp_occu()` (#40) — see Architecture. **continuous NNGP `svc()` varying coefficients on the psi arm under laplace/nested_laplace (#144)**; NUTS single psi intercept RE; Laplace grouped RE psi OR p11; areal psi arm; **detection-arm areal field on the p11 logit (`detection=~icar()`, #114)**; areal+temporal AND temporal-only on psi arm (#78/#114); NUTS+areal icar/car_proper/bym2 field on psi arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| Presence + nominal class | Yes | — | `occu_categorical()` (#106, `R/occu_categorical.R`): categorical analogue of the `cover()` hurdle -- presence Bernoulli + an UNORDERED class (multinomial logit / softmax), not a continuous magnitude. Arms factorise exactly: `P(y=0)=1-psi`, `P(y=k)=psi*p_k`, `p=softmax(X beta_class)`. Multinomial math = FD-validated tulpa kernel `cpp_multinomial_logit_terms` (`src/multinomial_logit.h`); non-spatial Laplace fit = vectorised R Newton over the same closed forms. laplace ONLY (native multi-process LikelihoodSpec path for spatial fields / NUTS = documented follow-up). `simulate_occu_categorical`, S3 (`print`/`coef`/`predict`). `test-occu-categorical.R` |
| Royle-Nichols occupancy | Yes | — | `royle_nichols()` (#116, unmarked `occuRN`): abundance-induced detection heterogeneity, `p_site=1-(1-r)^N`, `N~Pois(lambda)`. Pure-R closed-form Poisson-sum marginal over per-site suff-stats (k,n), optim BFGS + observed-Fisher vcov (`R/royle_nichols.R`). Site-level detection; laplace only. `simulate_royle_nichols`, WAIC. `test-royle-nichols.R` |
| Time-to-detection occupancy | Yes | — | `occu_ttd()` (#116, unmarked `occuTTD`): survey records TIME to first detection; `z~Bern(psi)`, `t\|z=1~Exp(rate lambda)` censored at `surveyLength`. Pure-R two-state marginal (`.ttd_site_loglik` over suff-stats k/S_det/Tcens), optim BFGS + observed-Fisher vcov (`R/occu_ttd.R`); same recipe as royle_nichols w/ a continuous censored-exponential emission. `y`=`[N x J]` TTD matrix (val<Tmax=detection, >=Tmax=censored, NA=not surveyed). Site-level rate; exponential only (Weibull shape + visit-varying rate + areal follow-ups). `simulate_occu_ttd`, full S3 + WAIC. `test-occu-ttd.R` |
| Double-observer abundance | Yes | — | `double_observer(type=)` (#116, unmarked `multinomPois` double-observer pi): `N~Pois(lambda)`. **`type="independent"`** (default): two independent observers `p1`/`p2`; Poisson-multinomial thinning makes the 3 observable cells (obs1-only/obs2-only/both) INDEPENDENT Poissons `n_c~Pois(lambda*pi_c)` -> closed-form marginal, NO latent-N sum. `y`=`[N x 3]` cell counts. **`type="dependent"`** (removal-style, #116 follow-up): PRIMARY observer records what it detects, SECONDARY records only the primary's misses -> two cells `n1~Pois(lambda*ppri)`, `n2~Pois(lambda*(1-ppri)*psec)` (`.dobs_dep_site_loglik`). Single fixed primary = 2 cells for 3 params (lambda,p1,p2) = a ridge -> NOT identified; observer ROLE-SWAPPING (obs 1 primary at some sites, obs 2 at others, per-site `primary` in {1,2}) -> 4 cell means, identifies all three. `y`=`[N x 2]` (primary-detected, secondary-only) + `primary=`; single-observer `primary` warns. optim BFGS + observed-Fisher vcov (`R/double_observer.R`). `formula`=log lambda, `detection`=shared design w/ separate p1/p2 observer arms; `fit$means` names `lambda_*`/`p1_*`/`p2_*`. Both types recover lambda + both detection probs (dependent: 5-seed, cov ~0.9). `simulate_double_observer(type=)`, full S3 + WAIC. `test-double-observer.R` |
| Joint distance + removal | Yes | — | `gdistremoval()` (#116, unmarked `gdistremoval`, Amundson 2014): SINGLE-SEASON (NOT the open-population `distsampOpen` HMM). `N~Pois(lambda)`; detected birds cross-classified by a distance band AND a removal period. Total detected = binomial thinning of N, Pois closed under thinning -> `ysum~Pois(lambda*pdist*prem)` CLOSED-FORM (no latent-N sum); band + period counts = two conditional multinomials (`double_observer` Poisson-multinomial pattern). Half-normal band integrals CLOSED-FORM line + point (`.gdr_dist_cp`, == integrate to machine precision); depleting-removal `pi_k=r(1-r)^(k-1)` (`.gdr_rem_cp`). 3 site-level arms: log lambda (`formula`), log sigma (`detection`), logit r (`removal=~`). optim BFGS + observed-Fisher (`R/gdistremoval.R`, `.gdr_site_loglik`). `y`=`[N x Jdist]` band counts, `y_rem`=`[N x Jrem]` period counts (rowSums must match). `simulate_gdistremoval`, full S3 + WAIC. 20/20 recovery, bias<=0.011, pooled cov ~0.975 (`test-gdistremoval.R`). v1 = halfnorm key, line/point, Pois, constant r, availability phi FIXED 1 (single period does not ID phi); hazard key / NB-ZIP / multi-period phi arm = follow-ups (all closed-form under the same thinning). Closed-form validated == unmarked brute-force K-sum, diff 0 |
| Open-population distance sampling | Yes | — | `distsamp_open()` (#116, unmarked `distsampOpen`, Dail-Madsen + distance): open-population counterpart of `gdistremoval`. `N_1~Pois(lambda)`, `N_t=Binom(N_{t-1},omega)+Pois(gamma)`; each primary period the detected birds are distance-sampled into bins. Band allocation conditional on the period total -> factors OUT of the abundance HMM (the gdistremoval trick) -> marginal REUSES the validated `cpp_dyn_abun_total_log_lik` kernel fed **eta_p=logit(pdist)** (overall distance detection as the detection logit) + y=per-period detected totals (J=1) + per-period band multinomials. NO new HMM kernel. 4 site-level arms: log lambda (`formula`), log sigma (`detection`), logit omega (`omega=~`), log gamma (`gamma=~`). optim BFGS + ANALYTIC gradient (`.dso_grad`: `dyn_abun` returns grad_eta_*, sigma block chains through pdist via a distance-integral FD; ONE kernel call/eval; FD-validated) + observed-Fisher (`R/distsamp_open.R`). `y`=`[n_sites x n_bins x n_seasons]`. K_max default = DETECTION-CORRECTED (`max(ntot)/pdist + infl*sqrt + 10`, infl=4 Pois / 5 NB-ZI, NOT 3*max -- the forward is cubic in K). **NB/ZIP/ZINB init (`mixture=`, #116 follow-up)**: `negbin` threads the kernel's `use_nb`/`eta_logr` through `.dso_negll`/`.dso_grad` (trailing `log_r` coord, analytic `grad_eta_logr` FD-validated) + `.dso_unpack`/`.dso_site_loglik`; `zip`/`zinb` = pure-R additive layer `.tobs_fit_distsamp_open_zip` over the composed per-site marginal (HMM `log_lik_site` + band multinomials), `omega*1{all bands 0}+(1-omega)*L_open`, per-site eta/sigma grads weighted by structural-zero posterior w_i (band part=0 on all-zero sites) + ZI-logit score `(1-om)-w_i`, log_r via central FD (summed-only score); ZI coord `zi_logit`/`zi_omega` (NOT survival `omega`). `simulate_distsamp_open(mixture=,size=,zi=)`. lambda/sigma/zi recover on all mixtures; omega/gamma/NB-size on the weak ridge at short series. Full S3 + WAIC. 20/20 recovery (Pois), per-param cov all >=0.90 (`test-distsamp-open.R`). **Alternative dynamics (`dynamics=`, #116 follow-up)**: `constant` (default) / `notrend` / `trend` / `autoreg` / `ricker` / `gompertz`, the unmarked `distsampOpen` tp1..tp5. `.dso_dyn_meta(dynamics)` sets the active-arm layout + kernel code; `constant`/`notrend` keep the analytic-gradient path (notrend ties `gamma=(1-omega)*lambda`), the density-dependent forms (`trend`/`autoreg`/`ricker`/`gompertz`) use a value-only forward kernel `compute_dyn_abun_site_dyn` (`src/dyn_abun_kernel.h`, `int dynamics` 2..5; brute-force validated) exported as `cpp_dyn_abun_dynamics_log_lik`, driven by numeric-gradient BFGS (`.dso_dyn_negll`/`.tobs_fit_distsamp_open_dyn`, dynamics-aware binder returns only the active arms). `ricker`/`gompertz` estimate a carrying capacity `K` (log link, `omega=~` slot) + growth `r` (`gamma=~` slot); `trend` drops survival; `autoreg` gamma is per-capita. lambda/sigma recover on every dynamics (`test-distsamp-open.R` sec 8); K/r/omega/gamma on the short-series ridge, and the density-regulated forms grow toward K so bound `K_max` (forward cubic in K). `fitted()`/`predict()`/`simulate()`/WAIC dynamics-aware (`.dso_draw_dyn`). v1 = halfnorm, line/point, site-level arms; season-varying detection, NUTS = follow-ups |
| Multi-species co-occurrence occupancy | Yes | — | `occu_multi()` (#116, unmarked `occuMulti`, Rota 2016): joint state `z in {0,1}^S` from a log-linear model w/ first-order (per species) + second-order (per pair = the INTERACTION) natural params; conditional per-species detection. Pure-R marginal enumerating the `2^S` states (`.occu_multi_site_loglik`), optim BFGS + observed-Fisher vcov (`R/occu_multi.R`). `y`=list of S `[N x J]` matrices or 3D array; `species=`. Shared covariate design (separate coefs per natural param), site-level detection, laplace only (per-param formulas, higher order, visit-level det, areal = follow-ups). NB: log-linear natural-param SLOPES trade off (weakly ID'd individually); MARGINAL psi + the interaction recover cleanly -> recovery targets those. `fit$means` names `f_<sp>`/`f_<sp>_<sp>`/`p_<sp>`; `fitted()$psi`=marginal occupancy. `simulate_occu_multi`, full S3 + WAIC. `test-occu-multi.R` |
| Open N-mixture (Dail-Madsen) | Yes | Yes | `dyn_abun()` (#37); y 3D `[n_sites x J x T]` — see Architecture. Pois/NB init; **ZIP/ZINB (`mixture="zip"/"zinb"`, #116)**: structural-zero site never occupied in any season -> all its counts 0, per-site marginal `omega*1{all y=0}+(1-omega)*L_DailMadsen` (abun ZIP additive layer over the forward-HMM marginal); PURE-R `.tobs_fit_dyn_abun_zip`, analytic-gradient BFGS, FD-Jacobian observed-info vcov; ZI coord named `zi_logit` (NOT `omega_*` = SURVIVAL arm); `simulate_dyn_abun(zi=)`; non-spatial laplace intercept-only omega, field/RE/NUTS stay Pois/NB w/ pointer (`test-dyn-abun-zip.R`). season-varying omega/gamma via `[n_sites x (T-1)]` covariate, interval-indexed forward kernel, all backends (#80). **continuous NNGP `svc()` varying coefficients on the initial-abundance arm under laplace/nested_laplace (#144)**; NUTS single intercept RE on init-abundance OR detection arm; Laplace grouped RE on init-abundance OR detection arm (#82, p-arm = per-node full-forward second-order eta_p pass); areal init arm; **detection-arm areal field on the detection logit (`detection=~icar()`, #114)**; areal+temporal AND temporal-only on init-abundance arm (#78/#114); NUTS+areal icar/car_proper/bym2 field on init-abundance arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113); NUTS+temporal-only fixed-hyper ar1/rw1/rw2/iid field on init-abundance arm (#114, same field block as areal, `field_map`=period, temporal whitened loading; areal+temporal under NUTS gated to n-L) |
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
(`nuts` + an areal term on psi, #74/#113/#204,
`R/occu_cover_nuts.R::.tobs_fit_occu_cover_nuts_spatial`): non-centered coupled areal field
— psi-arm field `z` (one/cell) enters psi linearly, copied to pos arm w/ amplitude `alpha`.
Param vector `c(beta_psi, beta_p, beta_pos, log_disp, raw_field, u_sigma?, u_rho?, u_alpha?)`.

**Hypers SAMPLED, not pinned (#204)** — the point being that a fit conditioned on the outer
grid's own point estimate cannot serve as an independent reference for that grid. Every areal
kind's loading factors as a FIXED basis w/ hyper-dependent column weights, so a leapfrog step
costs a rescale, never a re-decomposition (`src/nuts_field_hyper.h`, R mirror `.ochf_*` in
`occu_cover_nuts.R`):
`z = sigma*(B1 %*% (s1(rho)*raw1) + s2(rho)*raw2)`.
icar `B1` = sum-to-zero eigen-loading of intrinsic `Q` (#71), `s1=1`, no `raw2`;
bym2 same `B1`, `s1=sqrt(rho/sf)`, `raw2` iid w/ `s2=sqrt(1-rho)` (Riebler);
car_proper `Q(rho)=D-rho W = D^{1/2}(I-rho Lambda)D^{1/2}` in the eigenbasis of the
symmetrically normalised adjacency `D^{-1/2} W D^{-1/2}=U Lambda U'` -> `B1 = D^{-1/2}U`
FIXED, `s1_j=(1-rho lambda_j)^{-1/2}`. That last one is why car_proper rho is NOT the O(n^3)
per-step Cholesky the issue scoped it as (measured dense `chol` at n=2025: 0.974 s/call ->
27.1 h per 1e5 leapfrog steps; the eigen route pays 9.7 s once and nothing per step).
Each sampled hyper rides `t = t_lo + (t_hi-t_lo)*plogis(u)`, `value = inv_link(t)` (log for
sigma/alpha, logit for rho) — bounded, so no wall, Jacobian in the target.
**Prior = flat in `t` over the WARM FIT'S OWN outer-grid span** (`fit$theta_grid` column
range), which is the measure the nested-Laplace grid integrates against (tulpa defaults to
equal weight per cell over log-spaced sigma/alpha + logit-spaced rho nodes;
`R/nested_laplace_joint_hyperpriors.R` shows the optional pc.prec/half_normal families are
NOT set by this path). So `control$sigma.grid`/`alpha.grid`/`rho.car.grid` move BOTH backends.
Flat prior + change of variables = normalised `log(e)+log(1-e)`, `e=plogis(u)`.
alpha's grid `0` atom is not HMC-representable -> bounds take the positive nodes only.
icar pins `rho=1` (intrinsic precision has no mixing param); an axis the grid pinned to one
node stays pinned. `fit$nuts$sampled_hyper` / `$fixed_hyper` (CHARACTER vectors, empty when
nothing pinned — deliberately a type change from the old `fixed_hyper=TRUE`, so a stale
`isTRUE()` read fails loudly) + `$fixed_hyper_values`; `fit$hyper_draws` cols
`sigma|rho|alpha|field_sd`. **`field_sd`** = geo-mean marginal SD the block implies at that
draw = the ONLY field-scale summary comparable across kinds (the three normalise their
precisions differently), so it is what a simulation truth is stated in: the simulator's `f`
carries geo-mean marginal variance 1 (Sorbye-Rue), truth = `sigma`. `control$fixed.hyper=TRUE`
restores #74 conditioning, byte-identical to the old loading (pinned = the degenerate
configuration of the same block, not a second path).
Validated: C++ == R oracle to ~1e-13, analytic == central FD to ~1e-9 incl. every hyper
coordinate, prior verified flat in `t` by holding `raw=0` (field vanishes -> the target IS the
prior; catches a missing Jacobian, which the gradient check cannot).
Was: `inv_metric` sized `n_cells` while bym2's `n_raw = 2n-1` -> the engine took a
short metric pointer w/o a length check; now sized `n_raw + n_hyper`.
**Sampled terms reach the criteria (#211).** `.tobs_occu_cover_components()` returns the
structured terms as OFFSETS beside the coefficient draws: per SITE `field_occ`/`field_pos`
(#204's sampled field off `fit$field_draws` + the per-draw `alpha`; the v3 route still reads
`field_table`), per VISIT `off_det`/`off_pos` (#205's sampled obs-arm RE, `re_draws` ->
`sigma_re*z` mapped through the view's `flat_idx`). All four diagnostics fold them in --
ploglik (WAIC/LOO/CPO), PPC, PIT/LOO-PIT -- b/c the per-visit offset enters the SHARED `Arms`
view (`src/occu_cover_ragged.h`, `eta_p_visit`/`eta_pos_visit`), so ONE change reaches all
three kernels and dense == compact by construction (#185). A 0-column matrix = "arm carries
none" -> null pointer -> the no-offset path byte-identical. Cell-aggregated cover scores one
cover row per detected UNIT, so a per-visit offset errors there w/ a pointer. `test-occu-cover-nuts-ic.R`.
Was: `fit$draws` = the coefficient block + `ncol(draws)` read positionally as `log_disp`, so
the RE draws never reached the scorer, and `.tobs_occu_cover_v3_field()` read only
`field_table` -> both terms scored at ZERO on a fit that carried them (elpd moved 31 nats on a
det-arm RE, 549 on a sampled field). The grid-integrated (`nested_laplace`) route reaches the
same criteria (#215): `.tobs_joint_draws()` returns its RE latents on `bundle$re` in the layout
the offset builder already reads, so the components builder feeds them through the same
per-visit path (elpd_waic moved 11.5 nats on a det-arm RE at sigma 1.1, 61.5 on a cover-arm
RE). An occupancy-arm RE (#56) is per SITE, so the dispatcher stores its group codes on the
model (`model$re_psi`, the occupancy counterpart of `model$re_det`/`re_pos`) and the components
builder adds the per-group draw to `field_occ` -- the per-site offset every kernel already reads
(elpd_waic moved 9.34 nats, elpd_loo 9.41, at sigma_re 1.57; a fit carrying none byte-identical).
`.occu_cover_spatial_fields()` also carries the term's `var` + factor `levels` through, so
`fit$re$psi`/`ranef()` label the grouping like the obs arms and `predict(newdata=)` matches it.
**A SECOND (SVC / trend) field samples too (#214).** The block is a LIST now
(`hyper_field_build_list()`, `src/nuts_field_hyper.h`): each field carries its own basis,
site->node map, per-site design WEIGHT (`field_weight`; absent = the intercept field) and
its own sampled `(sigma, rho, alpha)` -- two fields share no hyper. Site i loads
`sum_b w_b(i) z_b[cell(i)]` on psi and `sum_b alpha_b w_b(i) z_b[cell(i)]` on cover; the
three places that loading is written are `hyper_field_site_{value,offsets,score}()`, so the
eta assembly and the score cannot express it differently. Spec spelling `field_blocks`
(list); the top-level entries still read as one block, and a one-block fit is
`expect_identical` on lp+grad AND on a whole fit's means/sds/field/draws. Layout per block
`[raw, sampled hypers]`, blocks back to back, RE blocks after -- so `n.iter`-for-`n.iter`
the one-field vector is unchanged. **The warm fit is the multi-block coupled path**
(`multi = TRUE` arms + one `icar/bym2/car_proper` block per field w/ `svc_weight` + a copy
spec per block, `alpha.grid` / `alpha.grid.trend`), and it FORCES `integration = "grid"`:
above 3 axes the engine switches to a mode-centred CCD star whose column range is a design
radius, not an integrated span, and the sampler reads each axis's span as its flat prior's
support (measured: a CCD warm put alpha bounds at 1.49-4.70 for a requested 0.2-2 axis, and
the fit sampled alpha to 2.79). A DEFAULTED axis is thinned to 3 nodes over the SAME span
when a second field is present (the prior is defined by the span alone -> unchanged; the
tensor is a product over blocks). Reported: `fit$trend_field`/`trend_fields` (named by
weight column), `fit$trend_field_draws`, per-block suffixed hypers (`sigma_trend`,
`alpha_trend`, `field_sd_trend`, indexed when several) in `hyper_draws` / `sampled_hyper` /
`fixed_hyper`, and `fit$spatial$field_suffix`/`field_weights` so
`.tobs_occu_cover_sampled_field()` sums every block's loading into the criteria (#211 rule:
a block the scorer cannot see is scored at ZERO). Both surfaces recover w/ 0 divergences;
numbers (+ the CCD-vs-tensor measurement) in `NOTES_measurements.md`.
`test-occu-cover-nuts-svc.R`. Correlated `|` (one free-Sigma MCAR block), temporal + RE
still gated -> n-L. group_var maps sites>cells; predict() needs the joint object
(non-spatial laplace AND nuts both error w/ pointer); sampled-field (estimated-variance) route =
`ms_occu_cover()` factor (tulpa#67). **Spatial default** (`nested_laplace`,
`R/occu_cover_joint.R`): `joint` engine via `tulpa_nested_laplace_joint()` w/
`occu_cover_{lognormal,beta}` cell-coupling spec (tulpa#32) — 3-arm joint
nested-Laplace, outer-grid over `(sigma, alpha)`, per-cell occupancy mixture
closed-form derivs drive inner Newton. Much faster than v3_nested, completes at
sizes v3_nested does not. 18-seed lognormal + beta recovery
(`test-occu-cover-joint-coupled.R`); status `"working"` (#96). Shared-field occ
SLOPE Wald CI mildly anti-conservative small-N (pooled ~0.94; NUTS non-spatial
calibrated). Outer Pareto-k diagnostic (`control$diagnose.k`) defaults OFF
(#101): dominates joint-fit wall time (re-solves the inner Laplace on the full
field vs the grid's node count). tulpa#118 sped it up (Shamanskii reuse via
`.K_DIAG_REFRESH` -> grad-only scatter; loosened inner tol `.K_DIAG_TOL`=1e-4;
near-neighbour batch order), k-hat byte-stable (externally validated ==
`loo::psis`/`posterior::pareto_khat`). Numbers in `NOTES_measurements.md`. Still
OFF by default: reports k-hat only, fit byte-identical on/off, opt-in; matches
`occu_joint` / `occu_multiscale_cover`. `control$diagnose.k = TRUE` re-enables.
Same default flip on `occu_multiscale_cover_joint.R`.

**Cover-arm intercept prior (#32)**: on the shared-field path the cover intercept
confounds w/ field level over detected cells. `.occu_cover_coupled_arm_priors()` hands
the pos arm a `cover_priors()` weakly-informative intercept prior **by default** (NOT
the engine flat ridge); else the cover intercept floats to huge SD and `predict()`
conditional cover blows up via Jensen. `priors = FALSE`/`"none"` disables all three arms.

**Cell-aggregated cover (`cover_aggregate`, #33)**: per-visit cover gives the cover arm
one row/visit -> a cell w/ many detected plots drives the field more than its single
occupancy obs. `cover_aggregate = "mean"` (default spatial) / `"median"` collapses the
cover arm to ONE row/occupancy unit; `"none"` keeps per-visit. ONLY on spatial `joint`.
Needs cell-level positive design (from `data`); visit-level `positive` keeps per-visit.
C++ compile-time `Aggregated` flag on `OccuCoverCoupling`
(`src/cell_coupling_occu_cover.h`), registered `occu_cover_{lognormal,beta}_agg`; R
`.occu_cover_build_joint_arms(cover_aggregate=)`. `test-occu-cover-coupling.R`,
`test-occu-cover-aggregate.R`.

**Latent cover-per-unit (`cover_aggregate = "latent"`)**: principled mean/median
alternative — cover arm carries a per-unit cover RE `u_i ~ N(0, sigma_u^2)` shared across
the unit's detected visits, integrated out per unit (keeps EVERY visit). Unit-level
predictor -> per-unit marginal `log M_i` SCALAR in one eta, reuses the one-row-per-unit
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

**Compact (ragged) input**: `tobs_data(compact=TRUE)` (the DEFAULT under
`method="nested_laplace"`, `R/tobs.R:433`) returns a `tobs_ragged` carrier -- one row
per VALID visit in `order(site, visit)` -- instead of a padded `[n_sites x max_visits]`
grid -> memory O(observations), NO per-site visit cap. Binder
`.tobs_build_occu_cover_ragged` (`R/occu_cover.R`) sets `ragged=TRUE` +
`site_of_visit`/`y_det_visit`/`y_pos_visit`/V-row visit designs; a compact model carries
**NO `model$y` / `model$y_pos` / `model$valid`** (all NULL). Gated to the joint
nested-Laplace path + `cover_aggregate="none"`. `test-occu-cover-compact.R`.

**A per-visit diagnostic reads `.occu_cover_visit_view()`, NEVER `model$y`/`$valid`**
(#185). ONE length-V view for BOTH layouts: a compact fit's stored visit rows, or a
dense grid flattened by `.occu_cover_dense_ragged()` (site-major, visit-ascending = the
order the dense `rowSums` accumulates). Also derives `n_valid` + `any_det` -> per-site
detection summary has ONE definition. All three consumers -- pointwise loglik
(`cpp_occu_cover_ploglik_ragged`), PPC (`cpp_occu_cover_ppc`), PIT/LOO-PIT CDF limits
(`cpp_occu_cover_cdf_limits`) -- go through it, and the three kernels assemble per-draw
predictors from one shared `Arms` view (`src/occu_cover_ragged.h`). One kernel per
diagnostic => dense == compact TO THE BIT (0.000e+00 on elpd_waic, elpd_loo, LOO-PIT,
PIT residuals, PPC fit.y/fit.y.rep/bayesian.p). Reading the dense grid instead is what
made `cpo()`/`ppc()` error w/ "'x' must be an array of at least two
dimensions" on every compact fit while `waic()` worked. The
aggregated (mean/median/latent) PPC keeps the padded grid + `cpp_occu_cover_ppc_agg`
-- aggregation is dense-only by gate. Cover density gates on `detected AND finite`
everywhere (a detected visit may carry NA cover, missing-at-random); the PPC used to
score that NA and returned `fit.y = NA` for every draw.

**`group_var` (sites > cells)**: `group_var="<col>"` on icar/bym2 maps each site ->
field node, so `n_sites` > `n_cells`. Field length `n_cells` while psi/p/cover run over
`n_sites`; per-arm `spatial_idx` (field node) + `cell_obs_map` (occupancy unit)
decouple. Layout: site = cell x time-period. R-side only (`.dispatch_occu_cover`,
`.occu_cover_build_joint_arms`); joint only. `test-occu-cover-group-var.R`.

**Per-group RE on shared-field path (#56)**: ONE random intercept on psi —
`re(g)`/`(1|g)` — alongside the field joins the joint fit as one `iid` prior block
(consumer of tulpa#86). Variance integrates on the outer grid (reported `sigma_re`);
BLUPs in `fit$re` + `ranef()`. `.occu_cover_spatial_fields` extracts `tobs_re`, fitter
appends the iid block (`obs_idx` = group on psi rows, 0 on p/pos),
`.occu_cover_jc_postprocess` extracts `sigma_re` + BLUPs. Scope: ONE random intercept;
slope/correlated/RE-without-field rejected; v2/v3 hatches carry no RE. Per-arm community
variances do NOT scale here (joint engine grid-integrates every variance component) ->
COMMUNITY spatial occu_cover = reduced-rank Laplace-EM `ms_occu_cover()` + `icar()`
(tulpa#67, below), NOT this path. `re.sigma.grid` knob.
`test-occu-cover-field-re.R`. `fit$re` = per-TERM flat list keyed by arm (lone term) or
`"<arm>:<var>"` (crossed); psi entry `fit$re$psi`.

**Per-group RE on detection / cover arm (#102 intercept, #103 crossed/nested)**:
random intercept(s) on `detection=`/`positive=` (`(1|g)`/`re(g)`), per-VISIT
grouping (one code per (site,visit)), composes w/ the required psi field on the
nested_laplace joint path. `.occu_cover_obs_re_parse` (occu_cover.R) strips ALL re
terms off the obs formula BEFORE copy-extraction + design build (rejects other
structured terms; copy/re allowed), returns `$terms` (LIST of specs, crossed/nested)
+ `$has_slope`. `.occu_cover_obs_re_design` resolves each term's per-(site,visit)
codes site-major from `data` (site-level, broadcast) or `visits` (visit-level) via
`.occu_cover_obs_flat_eval`, levels from observed visits ONLY, builds slope `Z`
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
marginal observed info (Louis 1982). Beta, lognormal, AND identity-Gaussian
(`ms_occu_cover("gaussian")`, #127: delta-normal magnitude, `mu = eta`, shared residual
`sigma_pos`; recovers to lognormal parity). Non-spatial Laplace ONLY (structured term
on any arm errors + pointer). Recovery + 20-seed coverage (`test-ms-occu-cover.R`);
status `"working"` (#98). Community VARIANCE carries Laplace small-cluster attenuation
(means do NOT); AGHQ-debiased BY DEFAULT below `re.aghq.maxdim` (4), above the cap EM
variance = tested lower bound (tensor AGHQ exp in total RE dim). Flagged via
`print.tobs_fit` + `fit$ms_community$var_attenuation` marker + `?ms_occu_cover` (#47).
WAIC/DIC/CPO via `.tobs_ploglik_community_occu_cover` (`R/community_ploglik.R`, #116):
exact per-(species,cell) two-state occu_cover marginal (`.occu_cover_site_ll`, beta/
lognormal/gaussian) scored over community-mean pseudo-draws w/ per-species BLUP plugged
in; routed via `.tobs_ploglik_ms_community`. NUTS/negbin/dispersion RE pending.

### `occu_multiscale_cover()` detail

Three-level occupancy + cover hurdle (#29; `R/occu_multiscale_cover.R`, fitter
`R/occu_multiscale_cover_joint.R`). For data where "visits" are spatially distinct
PLOTS aggregated into `(cell, period)`, NOT temporal revisits (EVA/MOTIVATE
vegetation; Nichols 2008, Mordecai 2011). `occu_cover()` treats plots as detection
replicates -> conflates within-cell prevalence into detection (Kendall & White 2009);
this family adds an explicit middle level:

```
z_c        ~ Bernoulli(psi_c)        # cell/range occupancy
a_cj|z=1   ~ Bernoulli(theta_cj)     # plot availability/use
y_cjv|a=1  ~ Bernoulli(p_cjv)        # detection
cover|y=1  ~ f_pos(eta_pos, disp)    # hurdle (beta/lognormal/gaussian)
```

Both z (cells) + a (plots) marginalize closed-form (two states each) -> exact joint
marginal LL, reuses occu_cover nested-Laplace cell-coupling machinery.

**Inputs**: `y`/`y_pos` = `[n_plots x max_visits]`. State `formula` = cell-level psi,
MUST carry areal field naming per-plot cell col: `icar(graph=adj, group_var="cell")`.
`availability=~...` = plot-level theta (default `~1`); `detection` = per-visit p;
`positive=~...` = cover. **Engine** (`nested_laplace` = SPATIAL engine; `laplace` +
`nuts` also fit, both non-spatial -- see Scope): 4-arm generalization of occu_cover
joint via `tulpa_nested_laplace_joint(cell_coupling="occu_multiscale_cover_*")`. Field
coupling: psi `field_coef=1`; theta/p `0`; pos `list(name="alpha")`. Cell spec
(`src/cell_coupling_occu_multiscale_cover.{cpp,h}`, per-fit) reuses occu_cover helpers
(`src/occu_coupling_shared.h`).

**Identifiability**: theta + p separate ONLY w/ replication WITHIN a plot. Single
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
recovery-tested. `dyn_occu`/`int_occu` NUTS are recovery + CI-coverage tested
(#139, `test-family-nuts-coverage.R`; the `cpp_occu_fit` path reports NO
per-parameter `fit$sds`, so those blocks read the 95% CI off `fit$draws` via
`.nuts_ci_cover_draws`). **`n.iter` = POST-WARMUP samples kept per chain on ALL
NUTS paths; the total run per chain is `n.iter + n.warmup`** (every R sampler
call site passes `n_iter = n.iter + n.warmup` to the engine, which returns
`n_iter - n_warmup` draws). The occu/dyn_occu/int_occu + cover/occu_cover paths
used to pass `n_iter = n.iter` (n.iter = TOTAL) -> kept `n.iter - n.warmup` and
returned ZERO draws (NaN means) at `n.iter == n.warmup`; unified to the sampling
convention (bugfix, occu_fit.R/cover_nuts.R/occu_cover_nuts.R).

**Convergence diagnostics on EVERY NUTS path** (#174): ONE writer
`.tobs_nuts_attach_convergence()` (`R/nuts_chains.R`) fills the record
`summary.tobs_fit` + `print.tobs_fit` read -- `convergence$parameter/rhat/
ess_bulk/ess_tail` -- plus scalars `fit$max_rhat`/`fit$min_ess`. tulpa owns the
estimator (`tulpa::diagnostics()`: rank-normalized split-Rhat, bulk ESS,
5%/95% tail-indicator ESS, Vehtari 2021); tulpaObs only hands it per-chain
matrices + names. Wired at shared choke points: `.tobs_count_nuts_attach()`
(abun/removal/distance/fp_occu/dyn_abun), `.tobs_nuts_field_attach()` (their
fixed-hyper field paths), `.ms_ocs_finalize_nuts_fit()` (community samplers),
plus ms_count/jsdm, cover, occu_cover (both), occu_multiscale_cover, and the
spatial-factor community fit via `.ms_ocs_attach_spatial_convergence()`.
`convergence$parameter` MUST carry the names `summary()` puts on its rows, else
record present but unreadable -> writer takes `par_names` (fit's `fixed_names`)
+ `cols` (sampler coordinates those name). Most community builders report
moment-matched pseudo-draws around the posterior mean, NOT the chain -> their
record comes from the sampler's own `rc$chains` at `par_cols` (default `lay$mu`;
ms_dyn_occu adds `lay$global`, ms_occu_cover adds `lay$log_disp` / `lay$mu_ld`).
Computing it off `fit$draws` there = diagnosing a chain never run. `converged`
on a sampled fit = every reported split-Rhat < 1.01 (print's warn threshold),
NOT the warm-start optimiser flag. `test-nuts-convergence-contract.R` fits every
family advertising `nuts`, asserts the record resolves an Rhat for every
`summary()` row (199 assertions, ~173s); its name-set check fails when a new
NUTS family lands w/o a case. `.tobs_nuts_rhat_ess()` (`list(rhat, ess)`
accessor for the `fit$nuts` full-coordinate block + PG-Gibbs summariser) reads
the same table. PG-Gibbs still keeps rhat/ess at `fit$rhat`/`fit$ess` ONLY,
where summary/print do not look.

**Areal field on the occu NUTS path** (#142): `occu() + icar()/bym2()` under
`method="nuts"` EXPOSES its field. `occu_fit.cpp` emits `spatial_layout` (engine
ParamLayout) + names the columns (`spatial_field[i]`/`spatial_theta[i]`/
`log_tau_spatial`/`log_sigma_spatial`/`logit_rho_spatial`, no more `param[k]`);
`.tobs_areal_field()` (occu_fit.R) sets `fit$spatial_field` = centred per-cell
surface for icar/car_proper (level confounded w/ intercept -> centred). bym2
field = Riebler rho-mix of both blocks x graph scale factor -> columns named,
reconstruction left to the draws. icar field cor ~0.81
(`test-occu-areal-nuts-recovery.R`). car_proper NOT a wired occu-NUTS term
(errors at fit, pre-existing).

**SVC = two flavors behind ONE verb (#118, unified #146).** Both written under
`spatial()`; `model =` picks which. `spatial(~ 1 + w || cell, graph)` = areal,
`spatial(lon, lat, model = "svc", coefficients = )` = continuous. `svc()` stays
as direct ctor (like `icar()`). Continuous form selects columns BY NAME
(`coefficients = c("(Intercept)", "elev")`), matched against the arm design at
fit time by `.tobs_svc_columns()` (`R/occu_svc.R`, single source
of truth -- BOTH the Laplace field builder `.tobs_svc_field_blocks()` and the
NUTS packer in `R/occu_fit.R` call it, so the two backends cannot drift).
`indices =` (column positions) is the lower-level form, still supported; giving
both errors. Fits report `fit$svc_indices` (resolved positions) +
`fit$svc_coefficients` on both backends. A bar with `model = "svc"` errors with
a pointer to the continuous form. `test-svc-spatial-umbrella.R`. The two flavors
are still DIFFERENT ENGINES with different coverage -- do not conflate what they
fit, only how they are spelled:
- **Areal spatially-varying coefficient** (the `svcPGOcc` analogue): a WEIGHTED
  areal bar `spatial(~ 1 + w || cell, graph)` on `occu()` via `nested_laplace`,
  rerouted through the joint direct-grid engine (`.tobs_fit_occu_joint`, #81).
  RECOVERY-TESTED: `test-occu-spatial-svc-recovery.R` +
  `test-occu-svc-joint-recovery.R` recover the known intercept + trend surfaces
  (`fit$spatial_field`/`fit$trend_field`, `cor > 0.75`), the SD hyperparameters,
  and the coefficients. This arm arrives as a `spatial` term, not `svc`.
- **Continuous NNGP `svc()` term** (`spatial(lon, lat, model="svc",
  coefficients=)`, or direct `svc()`): wired to single-season `occu()` on NUTS
  (`populate_svc`, `src/populate_helpers.h` -> `data.svc_data`) AND on the
  deterministic backends (#143, below). On NUTS the eta-assembly + NNGP prior +
  gradient live in the compiled UPSTREAM tulpa engine, NOT tulpaObs. Surface
  exposed as `fit$svc_field` (n_obs vector / n_obs x n_svc matrix of posterior
  means, per-draw on `attr(., "draws")`), sliced by position from
  `fit$svc_layout`; block named (`svc_w[i,j]`, `log_sigma2_svc[j]`,
  `log_phi_svc[j]`), no longer `param[k]`. NEEDS `prior_range = c(r0, alpha)`
  (PC prior, `P(range < r0) = alpha`) -- tulpa ships anchors unset + refuses
  without them; neither package defaults them.
  **RECOVERY-VALIDATED** (#119, tulpa 0.0.82, `test-occu-svc-nngp-recovery.R`):
  divergences 72-83% -> **0 every seed**, phi onto truth. Two upstream causes:
  Uniform-behind-a-wall range prior (gcol33/tulpa#144) + SVC marginal-SD
  half-Cauchy improper on its sampled coordinate (nothing bounded sigma above);
  fixing the SD prior is what pulled phi on, the two = ends of the GP ridge.
  Surface cor did NOT move -- information-bounded at these settings, not
  sampler-bounded, so the calibration test asserts divergences + phi + sigma and
  deliberately NOT cor. tulpa CANNOT make this measurement: `svc()` = tulpaObs
  term, `cpp_tulpa_fit_generic` = plain LM, so no tulpa-side fit reaches the
  NNGP SVC path -- how #144 survived there. Numbers in `NOTES_measurements.md`.
  **Laplace backends (#143, `R/occu_svc.R`)**: `occu() + svc()` also fits under
  `method="laplace"` / `"nested_laplace"`. K surfaces = latent field blocks on
  the psi logit -> rides the SHARED areal-BFGS driver (`.tobs_areal_bfgs_fit`,
  `R/areal_bfgs.R`); two new pieces only: `.tobs_svc_nngp_field()` (continuous
  NNGP block w/ optional per-site design weight, continuous sibling of
  `.areal_field_car(weight=)`) + `.tobs_occu_svc_marginal()` (exact two-state
  occupancy marginal, Fisher-identity gradients `w-psi` / `w(y-p)`,
  FD-validated). Vecchia precision `Q=(I-A)'D^-1(I-A)` assembled in R
  (`.tobs_nngp_precision`) from the term's OWN neighbour structure w/ the
  compiled kernel's kernels/jitter/variance floor -> both backends integrate the
  SAME density, asserted == tulpa `cpp_test_svc_nngp_twins` to 1e-8. Hypers
  (sigma, phi) grid-integrated on both routes (`laplace` == `nested_laplace`
  here) -> `fit$svc_hyper`; surface -> `fit$svc_field` (NUTS naming). Surface cor
  matches the NUTS path on the same truth (information-bounded, NOT
  backend-bounded; numbers in `NOTES_measurements.md`). `fitted()` adds the
  surface in-sample via `model$occ_eta_offset`; `predict(newdata=)` does NOT
  krige to new locations (as on NUTS). Gated on occu(): detection-arm svc,
  spatial/temporal/re term alongside svc, `pg_gibbs` -- all error w/ pointer.
  `test-occu-svc-laplace-recovery.R` + `test-svc-guard.R`.
  **Observation families (#144, `laplace`/`nested_laplace`)**: `removal()`,
  `distance()`, `fp_occu()`, `dyn_abun()` carry `svc()` too, with NO family-specific
  code. Those four ALREADY ride `.tobs_areal_bfgs_fit`, and an svc surface IS just
  another latent block on the arm their `eval(theta, offset)` already exposes, so
  the whole wiring is `.tobs_svc_field_blocks()` (single source of truth for the
  term's validation + hyper grid) + `.tobs_build_field_spec(svc=, X_svc=)`
  appending one NNGP block per `indices` entry AFTER the areal / temporal blocks +
  `.tobs_attach_field_results(svc=, has_spatial=)` slicing the trailing blocks into
  `fit$svc_field`/`svc_hyper`/`svc_field_arm`. Composes with an areal and/or
  temporal field on the same arm. Surfaces load on the STATE arm only (log lambda /
  psi); a detection-arm areal field alongside svc errors (`.tobs_check_svc_arm()` --
  the driver exposes ONE `grad_eta`, so the surfaces would otherwise be fit against
  the detection arm). The N-mixture families (`abun`/`ms_abun`) do NOT get it: their
  areal path is the C++ count-spatial driver, not this one. NUTS still errors
  everywhere but single-season occu(). Surface cor in `NOTES_measurements.md`.
  `test-svc-families-recovery.R`.
  The driver also returns `res$eta_offset`, the marginalised per-observation offset
  the blocks jointly load; the family wrappers read it instead of re-deriving each
  block's site map, which also fixed the temporal-only fp_occu path (it indexed a
  length-`n_t` field by site via `res$field_mean[map]`).

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
  tobs.R / obs_families.R / occu.R   — top-level router+print; family ctors; .tobs_build_model()
  tobs_dispatch.R           — the per-family `.dispatch_<family>()` bodies (~1370 lines). tobs.R routes HERE; it does not hold the dispatch itself
  tobs_helpers.R            — `.tobs_family_methods` (the method-support single source of truth) + shared dispatch helpers. NOT tobs.R
  engine_defaults.R         — `.TOBS_ENGINE_DEFAULTS` / `.TOBS_FAMILY_DEFAULTS` (#183): per-engine SAMPLER defaults, ONE table. `.tobs_control_defaults(control, engine, family)` fills every knob left unset (resolve once per dispatcher branch); `.tobs_default(engine, knob, family)` reads one knob inline. Scope = sampler knobs ONLY -- `max.iter`/`tol` are per-ROUTE Laplace-EM values (`ms_occu_cover()` iterates 30 at 1e-3 on its own EM, warm-starts its sampler at 200/1e-4; `ms_occu()` plain areal C++ EM 100, its latent fitter 200) -> stay at their call sites. Log-link families (`ms_count`/`jsdm`/`ms_abun`) `sigma.beta=10` under nuts, logit-link occupancy families 5 (matches C++ model defaults) -> a `.TOBS_FAMILY_DEFAULTS` row, NOT a uniform "community NUTS" value. `test-engine-defaults.R` pins every resolved profile.
    **Every sampler fitter reads it now (#188)** -- table used to be read by the 14 community dispatch branches ONLY, while 29 fitters restated knobs in their own formals. Everything routed through `.tobs_fit_model()` (occu/dyn_occu/int_occu + abun/removal/distance/fp_occu/dyn_abun/count) had DEAD formals: that entry forwards explicit values, so `abun(method="nuts")` ran `.tobs_fit_model`'s `n.iter = 2000`, never its own `1000L` (measured 2050 iterations, 2000 kept). `cover_nuts`/`occu_cover_nuts`/`occu_multiscale_cover_nuts` have own dispatch -> THEIR formals live + diverged. Now every sampler knob = `NULL` formal + `.tobs_fill_sampler(environment(), engine, ...)` as first statement; `test-engine-defaults.R` asserts structurally that no fitter carries a literal + each calls the filler -> new sampler family cannot reopen it.
    Divergences DELETED (proven artifacts): `n.iter` 2000 -> 1000 on occu/cover/occu_cover NUTS (commit 8975470 changed `n.iter` TOTAL -> kept-post-warmup on exactly those paths, left the literal -> a default that always kept 2000-1000=1000 silently began keeping 2000; that commit's message records the behaviour change); occu `pg_gibbs` 2000/1000 -> shared 3000/1500 (was keeping 1000 where every `ms_*` sibling kept 1500, no recorded reason; `test-occu-pg-gibbs{,-spatial}.R` 29/29 green after).
    Divergences RECORDED (`.TOBS_SINGLE_SPECIES_NUTS` / `_LAPLACE`, via `.tobs_single_species_defaults(engine)`): `.tobs_fit_model()` entry keeps `sigma.beta = 10` (vs community 5, on laplace/nested_laplace too), `adapt.delta = 0.8` (vs 0.9), `seed = 42` (vs 1). Each = value the single-species recovery/coverage tests were calibrated against; override belongs to the ENTRY, not the nine families passing through -> one constant, not nine family rows.
    **`n.iter` means the OPPOSITE thing under `pg_gibbs`**: TOTAL sweeps, warmup taken out of it (kept = `length(seq.int(n.warmup+1, n.iter, by=n.thin))` = 1500); under `nuts` it is the KEPT count and the run = `n.iter + n.warmup`. Stated in the table's `pg_gibbs` block, pinned by a test recomputing the kept count from the profile.
    **`sd.load` (1.0) + `re.lkj` (1.5) = laplace rows (#189)**; `n.quad` NOT one number, deliberately -- `.TOBS_NQUAD_ROUTES` / `.tobs_n_quad(route)` enumerates seven marginals (`re_aghq` 9, `ms_nmix` 1, `ms_nmix_scalar` 2, `ms_occu_cover` 5, `cover_latent_beta` 15 vs `cover_latent_lognormal` 1 (closed form, needs none), `community_latent` 5). `?tobs` lists per route, no single default most routes do not use
  occu_categorical.R        — presence + nominal K-class hurdle (#106); Bernoulli presence arm + baseline-category multinomial logit class arm over tulpa `cpp_multinomial_logit_terms`
  occu_cover_dispatch.R     — formula-native cross-arm copy()/spatial-field DAG coupling dispatch for occu_cover()
  occu_joint.R              — standalone occu() SVC-spatial-bar nested-Laplace path (#81): occu_cover()'s joint direct-grid engine with the cover arm removed
  cover_hurdle_joint.R      — the joint nested-Laplace cover fit (lognormal / beta), incl. occurrence-arm suff-stat aggregation `.cover_aggregate_occ`. One of the largest files in the package
  cover_nuts.R / src/cover_nuts.cpp — standalone cover() NUTS
  joint_substrate.R         — shared joint-fit draw / Pareto-k substrate (.tobs_joint_fit, .tobs_joint_draws*) across the occu_cover / cover nested-Laplace joint routes.
    **Outer-grid placement promoted (#187)**: `.tobs_promote_outer_grid(jf)` lifts `outer_grid_placement` ("fixed"/"auto_recentered"), `_recenter_attempts`, `_prior_added`, `_recenter_declined` (reason a "fixed" placement stayed fixed, tulpa#293) to `tobs_fit` top level, spliced everywhere `.tobs_promote_pareto_k` is (occu_cover postprocess, occu_joint, occu_multiscale_cover_joint) + the cover decode. NOT gated on the grid having MOVED -- a declined recenter is exactly the case worth seeing; an inert one invisible across a whole batch is what filed it. `.tobs_glance_outer_grid(g, x)` adds the two columns in BOTH `glance.tobs_fit()` + `glance.tobs_multiarm_fit()`; latter terminal for `cover_fit` (class order cover_fit/tobs_multiarm_fit/tobs_fit) -> never reaches the former.
    **Defaulted grids declare themselves (#186, needs tulpa >= 0.0.132)**: engine auto-recenter decides axis PROVENANCE, not field presence -- moves an axis that is absent, `auto_grid()`-marked, or exactly equal to its own default; anything else = user pin. tulpaObs writes a grid on EVERY joint fit -> unmarked default reads as a pin, rescue goes inert. `.tobs_default_{sigma,alpha,bym2_rho}_grid()` now RETURN `tulpa::auto_grid(...)` (they are the layer that chose the values; a user grid never passes through them); `.tobs_mark_auto(x, auto)` re-applies the mark wherever a site reshapes a defaulted grid, since `sort()`/`[`/`c()`/`as.numeric()`/`expand.grid()` all drop the attribute -- cover arm-specific tau translation, occu_cover pos-arm tau, copy `alpha_grid`s, RE `sigma_grid`s, EM-path bym2/ar1 pairings. Verified end to end: defaulted axis reports `declined = "grid_not_collapsed"`, a `control$sigma.grid` one `"axis_pinned"`
  pg_gibbs_shared.R         — shared Polya-Gamma machinery (.tobs_pg_draw_beta, .tobs_pg_community_update, .tobs_pg_finalize_fit) behind EVERY pg_gibbs fitter
  stacking.R                — `tobs_stack()`, LOO-weighted predictive stacking across fitted tobs_fit objects of any family
  sbc.R                     — `sbc.tobs_fit()` (#207): posterior SBC, the `tobs_fit` method of `tulpa::sbc()`. Registry `.TOBS_SBC_REGISTRY` per family (simulate/refit/draws/+statistic); pooling, group_ids, arms, controls shared. Replicate on FRESH cells (block-diag graph) -> premise VERIFIED. occu_cover registered+verified; field SD not gated (#212). Numbers in NOTES_measurements.md
  aggregation_scan.R        — cell-size / NN-spacing / changepoint scan utilities
  community_em.R            — shared community Laplace-EM .tobs_community_em() (ms_occu/ms_dyn_occu/ms_int_occu/ms_count/jsdm/ms_abun-latent/ms_distance). Three OPTIONAL args, each defaulting to the previous behaviour byte-identically: `sp_info` (analytic per-species observed info; the FD fallback costs 2(P+G) marginal sweeps per species per Newton step -- supply it when the kernel exposes the Louis block), `init_b`/`init_Sigma` (warm start; a block-coordinate caller re-enters the EM once per outer pass, so a cold restart each time dominates)
  ms_{occu,dyn_occu,int_occu}.R      — community single/dynamic/integrated
  abun.R / abun_nuts.R               — nmix family + non-spatial NUTS (#41)
  count_spatial.R / count_methods.R  — areal count fitter .tobs_fit_count_spatial (#117); fitted/predict/residuals
  ms_count.R                — community count / relative-abundance GLMM (msAbund, #117); .tobs_fit_ms_count over shared community_em.R
  ms_count_nuts.R / src/ms_count_nuts.cpp — community count NUTS (msAbund NUTS, #117); in-tree C++ FullGradFn (reduced ms_abun_nuts: no detection/latent-N), R oracle .tobs_ms_count_nuts_logpost, warm-start from Laplace-EM
  ms_count_spatial.R        — community count + shared areal field (sfMsAbund) + SVC bar (svcMsAbund, #117/#118); block coordinate ascent (community EM offset <-> multi-field Poisson-ICAR), pure R
  community_latent.R        — SHARED latent-structure engine for EVERY community family (#119/#120/#121): one block-coordinate ascent (community EM w/ the latent as an offset <-> field / factor updates) + the areal Newton, the factor update, bym2/car_proper/spde hyper grids. A family supplies ONE callback `working(eta) -> list(score, curv)` (per-(site,species) score+curvature wrt an additive offset on the structured arm): Poisson `(y-mu, mu)`, occupancy two-state, Bernoulli `(y-psi, psi(1-psi))`. Field solve is `t(A) diag(w) A + tau Q`, so the site->node map slot takes an areal group_var incidence OR an spde barycentric projector unchanged. Adding a family to every latent route = one callback, not a new fitter. **The measured evidence behind every number below -- fixtures, seeds, wall times, per-family screens -- is in `NOTES_measurements.md`; the rules here are what it concluded.**

**Backtracking + guards.** Factor Newton (`.tobs_latent_factor_update`) backtracks: local `ascend()` halves the step until the penalized objective improves, holds previous iterate if never; `nstep()` ridge-bumps singular curvature. Non-finite guards are inline (`if (all(is.finite(Dz)))` / `(Dl)`), not a named helper; non-finite `working()` score/curv `break`s the pass. Field solve `.tobs_latent_field_solve` has its own local `safe_solve()` (ridge retry for a singular Hessian only, Newton update unconditional) -- do not confuse the two, there is no `safe_step()` anywhere in the repo.

**Loadings by MARGINAL likelihood, NOT the joint mode (#153 -> #156).** Factor update holds zeta at its joint mode -> `(zeta, lambda)` is a joint-likelihood estimate with `Ns*Q` incidental params growing with the sample = Neyman-Scott, inconsistent. Site factors' estimation error lands in the fitted co-occurrence and lambda absorbs it, over-fit growing with Q/S. Fix = `.tobs_latent_factor_mmle()`: EM on the SAME joint site marginal over all S*Q loadings (E-step = posterior `p(z_i|y_i)` off `.tobs_latent_joint_grid()`; M-step = per-species Qk-dim weighted Newton, backtracked on expected complete-data ll). Numbers in `NOTES_measurements.md`.

**ONE estimator, ONE state (#156).** `.tobs_latent_factor_update()` + `.tobs_latent_factor_scale()` run ONCE, outer pass 1, purely to INITIALIZE -- the marginal's lambda-gradient vanishes at lambda=0, and the 1-D bracket is a global magnitude search the local EM cannot do. Running the joint-mode update every pass alongside the MMLE diverges both ways, so it never repeats.

**Offset by SCORE-MATCHING, NOT `zeta t(lambda)` (#156).** `.tobs_latent_factor_offset()` solves `score(eta+off) = E_z[score(eta+lambda_s'z)]` per cell (scalar Newton) -- plug-in and integrated stationary conditions match for ANY family, no link-specific derivation; reduces to `lambda'zhat + v/2` on a Poisson log link. `fit$model[[offset_slot]]` reads THIS, not `zeta t(lambda)`.

**Block-coordinate callers MUST warm-start `init_b`/`init_Sigma`.** Omitting cold-restarts every per-species deviation and both community covariances on every outer pass.

**`max.outer`: factor path 150, field path 25 by default.** Field reaches `tol` and breaks early; factor does NOT (alternates with the coefficient block along a slow mode) -- 25 leaves real community-mean bias. `factor.outer` is per-family, set from that family's own measurement: ms_count/jsdm/ms_occu 150, everything else 25 until measured. Do NOT globalize -- cost is not transferable.

**`factor.starts` (multi-start width) dominates a latent-N fit, NOT `max.outer`.** #157's basin escape runs K candidate starting directions on the first factor pass, each a full loading-EM to convergence. Per family, measured from that family's OWN recovery suite (a random seed screen alone is not enough -- the committed regression test #157 was built for caught a value ms_count's own screen missed): ms_abun 1, ms_occu 1, ms_count/jsdm + ms_distance at the driver default 8. NEVER copy a family's value -- cost and benefit both scale with how expensive one oracle eval is.

**NEITHER summary screens a fit alone.** `residual_cor` is row-normalized -> blind to a magnitude regression; `mag_ratio` is rotation-invariant -> blind to a direction regression. Screen on both.

**`n.quad` NOT threaded from any caller** -- driver default 5 is what every community latent fit actually uses; `control$n.quad` silently does nothing here.

TRAP surviving all of the above: magnitude MUST come from the JOINT marginal, never a per-species one (identifiability ridge between sigma and coef scale). Assert on `sqrt(tr(Sigma_res))` (rotation-invariant), as `test-ms-count-factor.R` / `test-ms-occu-factor.R` do
  ms_occu_field.R           — community occupancy SVC (svcMsPGOcc, #118); block coordinate ascent (community occ EM psi offset <-> two-state-marginal occupancy field solve), intercept + SVC field(s), pure R; plain intercept -> C++ ms_occu_spatial.R
  ms_abun.R / ms_abun_nuts.R         — community nmix + NUTS (#14)
  ms_abun_latent.R          — community nmix + latent() factors (lfMsNMix) / + shared field (spatial-factor); the ONLY new piece is the working oracle over the Royle marginal: score=grad_eta_lambda, curv=info_eta_lambda-var_N*score_wt_lambda^2 (the Louis (1982) (1,1) block = abundance curvature w/ the detection arm profiled out), which nmix_site_marginal() already exposes. Supplies `sp_info` (the design-sandwiched per-site Louis block) so the community EM skips its FD Hessian: 387s -> see table. Plain field w/o factors KEEPS the C++ #12 path
  ms_distance.R             — community binned distance sampling (msDS, #117) + latent() factors (lfMsDS) + shared field (sfMsDS). NO new C++: cpp_distance_site_sweep already returns log_lik/grad_lam/info_lam/var_N/swl, so the community EM reads its per-species score from it and the driver oracle is the SAME Louis formula as ms_abun (curv = info_lam - var_N*swl^2). `.tobs_ms_distance_info_block()` assembles the full per-species observed information from the same sweep (#161) so the EM does not finite-difference it; sign inside `v` differs from the N-mixture's and is asserted, not inherited -- see the spatial-factor row. Hazard-key log-shape = a community `global` (shared across species), and the hazard key keeps the FD fallback (its per-site cross terms are not exported). simulate_ms_distance() draws through cpp_simulate_distance (the kernel the likelihood integrates against) -- a separate R-side quadrature simulates from a pi the model is not fit against and biases recovery. Poisson only; NUTS not wired
  nmix_laplace{,_re,_re_spatial,_spatial}.R — non-spatial / community / sfMsNMix / areal fitters
  nmix_re_aghq.R / nmix_site_marginal.R — grouped RE -> NMixGroupedOracle; per-site AGHQ callback
  occu_fit.R / occu_priors.R / laplace.R — .tobs_fit_model(); occu_priors()+beta_prior; .tobs_laplace()+EM cbs
  em_nested_laplace.R / simplified_laplace.R / sla_*.R — nested-Laplace EM; SLA wrapper + paths
  family_cover_hurdle.R    — .dispatch_cover() (two-Laplace hurdle), large
  occu_cover.R             — joint occu-det+cover wiring + .occu_cover_eta_from_par()
  occu_cover_nuts.R        — non-spatial occu_cover NUTS: R oracle + .tobs_fit_occu_cover_nuts
  nuts_chains.R            — multi-chain pooling + shared split-Rhat/bulk-ESS (.tobs_nuts_rhat_ess)
  ms_occu_cover.R / ms_occu_cover_spatial{,_nuts}.R — community joint; spatial-factor JSDM (tulpa#67) Laplace-EM + NUTS
  occu_multiscale_cover{,_joint,_nuts}.R — 3-level occu+cover (#29) + 4-arm joint fitter + non-spatial NUTS (src/occu_multiscale_cover_nuts.cpp)
  occu_cover_spatial.R / occu_cover_nested.R — v2_joint / v3_nested escape hatches
  formula_terms.R / formula_parse.R  — term registry+ctors; AST parser
  inputs.R                 — single source of truth for response/site/visit input: .tobs_check_site_count() (site-count cross-check every binder used to hand-roll), .tobs_input_dims()/fit$dims canonical totals, .tobs_unpack_frame() (tobs_data -> data/y/visits)
  spatial.R / methods.R / diagnostics.R / data.R / within_between.R — precompute; S3; diags; data+sims; decomposition
  RcppExports.R            — generated, do not edit
src/
  occu_fit.cpp / populate_helpers.h  — unified C++ entry; populate_spatial/temporal/re/svc/latent
  occ_*.h / dyn_occ_*.h / integrated_occ_*.h — single/dynamic(HMC fwd)/integrated
  cell_coupling_occu_cover.h / cell_coupling_occu_multiscale_cover.{cpp,h} / occu_coupling_shared.h — coupling specs + shared helpers
  occu_cover_ragged.h       — `Arms` (#185): the one-row-per-valid-visit predictor view every ragged occu_cover DIAGNOSTIC kernel assembles from (occu_cover_ploglik.cpp, occu_cover_diag.cpp). `make_arms()` = occ + det arms (what the CDF-limits kernel needs), `attach_cover()` adds the pos arm (loglik + PPC). NOT the fit kernels -- those are the cell-coupling specs above
  ms_occu_cover_spatial_nuts.cpp / abun_nuts.cpp / ms_abun_nuts.cpp / occu_cover_nuts.cpp — NUTS (#67/#41/#14; occu_cover non-spatial)
  nuts_engine.h            — shared run_tulpa_nuts driver for the in-tree FullGradFn targets
  nuts_field_block.h / nuts_field_hyper.h — the two non-centered areal field blocks. `_block` PINS the hypers at a nested-Laplace estimate and marshals one loading (abun/removal/distance/fp_occu/dyn_abun). `_hyper` SAMPLES them (#204, occu_cover): fixed basis `B1` + rho-dependent per-column weights + bounded transforms, so no leapfrog step re-decomposes anything; the pinned case is the same block with every hyper's coordinate absent, byte-identical to `_block`'s loading. A family moving from pinned to sampled swaps the header, not its eval
  community_chol.h         — shared log-Cholesky helpers (#14 non-centered, #67 centered) + `CommunityCholPri` / `community_chol_pri_read()` (#181): the log-Cholesky hyperprior scalars + the `pri` list keys, ONE declaration for all seven community NUTS targets. `MsOccuCoverPri` / the spatial-factor `PriScalars` INHERIT it and add their own fields; do not restate the three
  community_grid_pack.h    — `community_pack_grid()` (#181): the per-outer-grid-point pack shared by the community areal drivers (ms_occu_spatial.cpp, nmix_community_spatial.cpp). The state arm is reached by pointer-to-member and named by the caller ("psi" vs "lambda"), and a family with no boundary diagnostic passes a null member pointer. Include AFTER RcppEigen.h
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
