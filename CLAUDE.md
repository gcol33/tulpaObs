# CLAUDE.md

Guidance for Claude Code in this repo. Caveman speak: terse, telegraphic. File
paths + function names exact. (Compacted; if a detail is missing, read the source
file named.)

**This file has a hard 95k-char budget.** Claude Code caps ALL loaded
instruction files at 150k combined, and the global `~/.claude/CLAUDE.md` (~48k)
shares that cap. Keep RULES + CONTRACTS
here; put the MEASUREMENTS they rest on (fixtures, seeds, wall times, per-seed
recovery numbers, cost ratios) in `NOTES_measurements.md`, committed alongside and
Rbuildignored, and the MECHANICS behind a rule (kernel internals, block layouts,
derivations) in `NOTES_families.md`, same convention. Adding a paragraph of
numbers here costs someone else's context.

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
`latent()`; `share("id")` (share one realization state+detection). Ctors
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
wraps, supplies occ/det site marginal as `make_site`. AGHQ removes most of the
small-cluster sigma bias on both arms (`NOTES_measurements.md`). Det-arm
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
   `Sys.setenv(TULPAOBS_FAST = "1"); devtools::test()`. `skip_if_fast()` gates
   every fitting block; no-op when env var unset. Baseline assertion/skip counts
   + the #148 fix history: `NOTES_measurements.md`. Wall time NOT measured (box
   shared) -- time it yourself, do NOT trust old figures.

   **The switch is `TULPAOBS_FAST`, NOT the engine's `TULPA_FAST`.** Setting the
   engine spelling leaves `skip_if_fast()` INACTIVE, so the heavy blocks fall to
   `skip_on_cran()` instead and the run reports a clean tail having fitted
   nothing. Both spellings produce a green suite; only one of them ran it.
   A run is verified by its COUNTS, not its tail: `tests/expected-counts.csv`
   records assertions + skips per (file, tier) and `.github/scripts/run-tests.R`
   fails a run that asserts fewer or skips more than recorded (#302). Files
   absent from it are unchecked -- it holds what has actually been measured.

   Regenerate the manifest in the SAME commit as a change that legitimately
   moves counts; the checker's own message says so. Its first catch was
   `test-cover-perarm.R` at 72 against a recorded 107 after #295 turned two
   fit-based blocks into compile assertions -- a true positive wanting a
   regeneration, not a fix.

   Regenerate it from a MEASURED run, never by hand:
   `TULPAOBS_TEST_OUT=out Rscript .github/scripts/run-tests.R`, then
   `Rscript .github/scripts/record-counts.R out/timings.csv --tier=smoke`
   (`--dry-run` to see the delta first). The full tier records the same way
   from the `tier3-timings.csv` that `aggregate-tier3.R` publishes, one
   completed shard at a time. The recorder REFUSES a run carrying any failure
   or error -- a recorded count is the floor later runs are held to, so taking
   it from a red run pins the breakage -- and drops rows for files not in
   `tests/testthat` or that neither asserted nor skipped (a 0/0 row cannot
   fail, so it reads as coverage while checking nothing).

   **It is ADD-ONLY: it never rewrites an existing row without `--update`.**
   The two numbers are BOUNDS in opposite directions -- `assertions` a FLOOR
   (fails when a run asserts fewer), `skipped` a CEILING (fails when a run
   skips more) -- so a row is sound only if it holds in EVERY supported
   environment, and one run measures one. A box with more Suggests installed
   asserts more and skips less, so recording from it moves BOTH bounds inward
   and every leaner box fails on a healthy suite. Lowering is worse: a block
   that stops executing asserts fewer, which is what this
   manifest exists to catch, so auto-taking the minimum would erase the finding
   and hand back a green suite. Neither direction is decidable from counts
   ("fewer" is a lean environment or a dead block, identically), so differences
   on existing rows are REPORTED for a person and left alone; `--update` is for
   when you know why a count moved and are recording it in the same commit as
   the change that moved it.

   **Which box owns a floor is PER TIER, and the two answers differ on
   purpose.** Smoke floors come from the LEANEST box: smoke is cheap and runs
   everywhere, so a floor that holds without the optional Suggests is what
   keeps it honest. Full floors come from CI, deliberately: the full tier is a
   ~25h on-request dispatch that in practice only ever runs there, so a floor
   measured anywhere else describes a run nobody performs. That asymmetry is a
   decision, NOT an inconsistency to tidy up.

   26 of 286 files carry an environment- or outcome-dependent skip; the full
   partition and the betareg worked example are in `NOTES_measurements.md`. The
   rule that has to be here: **4 of those gate on a FIT OUTCOME**
   (`fit$converged` and friends), which a Suggests list does NOT predict and
   two sweeps on one box CANNOT bound -- fixed seeds make them deterministic
   per (platform, BLAS, engine), the #153 axis, so a floor recorded where the
   fit converged fails on a box with an IDENTICAL package set. **For those four
   the guard-fails note is TRUE and recording a lower floor is the wrong
   response** -- a lean box seeing it has found something.

   A stale `expect_error()` string cannot be found by grepping `R/`: nearly
   every message is `sprintf`/`paste`-composed, so no contiguous literal exists
   to match. A static audit over the suite returned 165 "unfindable" literals,
   almost all false. The instrument for that class is RUNNING the block, which
   is why the counts manifest and a finishable tier are the fix and a scanner is
   not.
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
`pkgload::load_all()` the same tree. Harmless pure-R; here `src/tulpaObs.dll` ~165MB,
so N workers recompile into one `src/` and race the build artifacts -> corrupted DLL
registration (`must specify DLL via a "DLLInfo" object`), which testthat wraps in
`cli_abort(..., parent = msg$error)` -- the real error is in that `parent`, NOT the
printed message (`rlang::cnd_message(e, inherit = TRUE)` / `e$parent$message`). A
silent 10-min hang w/ zero output = the same race landing differently: a partial DLL
load blocks in the loader w/o emitting its startup handshake -> `queue$poll(Inf)`
waits forever.

Recipe = install ONCE, then call testthat directly (never via
`devtools::test()`), every worker loading the built package instead of
recompiling:

```r
devtools::install(quick = TRUE)
testthat::test_dir("tests/testthat", package = "tulpaObs",
                   load_package = "installed")
```

`load_package = "installed"` is what closes #151: each worker's
`library(tulpaObs)` = one read of one built DLL, warm or cold `src/`, so no
worker recompiles into the shared tree.

**TWO different parallel mechanisms; neither is evidence about the other.**
`test_dir()` reads `Config/testthat/parallel` and dispatches files to testthat's
own **callr worker pool** -- and `test_check()` IS `test_dir(load_package =
"installed")`, so `R CMD check` takes that pool too. `.github/scripts/run-tests.R`
takes NEITHER: it holds a `parallel::makeCluster` pool and hands each worker ONE
file via `test_file()`, whose `parallel` formal defaults FALSE -> serial inside
the worker. So `TESTTHAT_PARALLEL: true` in `smoke.yaml` / `full-recovery.yaml`
reaches nothing (only `find_parallel()` reads it, and nothing on their path calls
it), and `R-CMD-check.yaml` sets it `false` on push/PR deliberately.

**The callr pool is where #289 lands, and `safe at any worker count` is NOT
established for it.** A worker exits `0xC0000005`; testthat reports `R session
crashed with exit code -1073741819` against whatever file that worker held,
aborts the run, and leaves every unreached file unrun. Intermittent, different
file each time. #151 established only that a worker no longer recompiles the
DLL.

**It is an R bug, NOT ours -- do not debug the package for it.** The fault is
`wcsstr` reading past a non-NUL-terminated `FILE_NAME_INFO.FileName` in R's
`R_is_redirection_tty()` (`src/main/sysutils.c`); no tulpaObs/tulpa frame is on
the stack. A testthat worker is `Rterm` with fd 0/1 wired to processx pipes ->
workers die at the very START of a file, the file varies, serial never shows it.
Fixed upstream 2026-08-11 (r-source `a2066dd40`, PR19104); no released R has it
yet. **Until the next R release: never run the parallel tier while anything is
BUILDING the package** -- that takes it from rare to common -- or run serially.
Ruled out by measurement, do not re-chase: binaries, resource exhaustion,
out-of-range vector/Eigen indexing, file ordering, callr (`NOTES_measurements.md`
for the numbers). The SERIAL route (`TESTTHAT_PARALLEL=false`) is clean; the
`makeCluster` route is neither implicated nor cleared.

Slow test -> pair `skip_if_fast()` + `skip_on_cran()` atop any multi-seed fit /
NUTS block (`tests/testthat/helper-speed.R`). C++ recompiles ccache-backed; only
killed/partial build needs `pkgbuild::clean_dll()`.

### CI (#149)

`.github/workflows/`:

- `R-CMD-check.yaml` -- push/PR + `workflow_dispatch`, NO cron. Checks a **built
  tarball**, not `load_all()`, so a missing NAMESPACE export surfaces (#147
  shipped `abun()` unexported precisely because `load_all()` resolves internals
  regardless). `--no-manual` (dev non-ASCII in Rd). ubuntu on push;
  ubuntu+windows+macOS on `workflow_dispatch`. `TESTTHAT_PARALLEL` is `false` on
  push/PR so the gate is deterministic and `true` on dispatch, which is the ONLY
  place CI reaches testthat's callr worker pool (#289); `test_check()` hardcodes
  `load_package = "installed"` -> not subject to #151 either way.
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
that = tier 3, blew the job cap. Green baseline counts: `NOTES_measurements.md`.

**Linux runner != local box; two test classes feel it.** Both surfaced on the
first CI run, neither reproduces on Windows. Magnitudes + the #153 case history:
`NOTES_measurements.md`.

- *Exact float equality across two BLAS call shapes.* Batched ploglik
  (`R/diagnostics.R`) builds eta w/ one GEMM; the R oracle builds it per draw as
  GEMV. Different routine -> different accumulation order -> eta differs in the
  last ULP, and a marginal summing many terms through `exp`/`lgamma` carries it.
  Reference BLAS collapses both to the same naive ordering -> equality held on
  Windows ALONE. Assert a tolerance AND report the magnitude: bare
  `expect_true(a == b)` cannot separate 1 ULP from a real bug, and the two want
  opposite responses. Thread-count invariance IS exact (same GEMM, same kernel)
  -> stays `expect_identical`.
- *Single-seed recovery under an aggregate tolerance.* `expect_equal(coefs,
  truth, tolerance = t)` scores the WHOLE vector -> a one-sided bias in one
  coefficient sits under `t` on one platform, over it on another. That was #153
  (a community slope biased high on nearly every seed, on BOTH platforms,
  tolerance hiding it). Test flips -> measure the bias over many seeds BEFORE
  touching the tolerance; widening it deletes the signal.

  Two things make that tolerance looser than it reads -- estimand scored
  against the wrong constant, and `tolerance`'s relative/absolute switch --
  take out BOTH before a community recovery number means anything; mechanism +
  the fix for each: `NOTES_families.md`.

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
`.nested_psi_quantiles()`). Calibrated on a 10x10 icar/bym2 fixture
(`NOTES_measurements.md`). Old tulpa w/o `fitted_eta_var` -> NA intervals.

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

Per-family mechanics (removal / distance / fp_occu / dyn_abun: kernels, arms, RE, areal + detection-arm fields, NUTS scope): `NOTES_families.md` "Observation families: per-family mechanics".

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
- **Latent-N truncation is PER SITE, and guarded.** The sum runs on `[max(y_i), K]`:
  the lower end is always the site's own maximum, only the ceiling is shared -> a
  shared ceiling makes every site pay for the largest count ANYWHERE in the data, and
  cost is linear in the state count, so one heavy-tailed species-site drags every
  other site's evaluation up by that ratio. `K_max = NULL` (default) caps each site at
  `max(y_i) + 100` (`.nmix_truncation`, `R/nmix_site_marginal.R`) = the SAME headroom
  the site holding the global maximum already had; an explicit `K_max` keeps its
  documented meaning (hard global truncation) and is NEVER capped, since a caller
  raising it is compensating for the one regime headroom cannot see. Threaded as a
  `headroom` argument (`< 0` = no cap) through `nmix_precompute_site` -- the single
  place `K_hi` is set, so every caller inherits it. Do NOT give `compute_nmix_site` a
  headroom argument: its arity IS the `CountKernelFn` function-pointer contract the
  shared count-NUTS / count-Laplace drivers take, and a defaulted 7th parameter
  silently breaks that conversion (`abun_nuts.cpp`); capping callers build the cache
  themselves instead.
  **The cap is verified, NOT assumed.** Exact where the posterior over N decays inside
  its window, wrong where it does not; no fixed headroom survives `p -> 0`
  (never-detected count ~`Poisson(lambda (1-p)^J)`). The guard is NOT the boundary
  mass `nmix_laplace()` warns on -- that bounds the error where the fit STOPPED, and
  the optimiser's PATH can run through the truncated region while boundary mass reads
  clean. What IS tested: the score at the fitted coefficients must agree under the
  capped AND the shared ceiling to `.NMIX_SCORE_TOL` (1e-4); failing that the window
  widens 4x and the fit is redone, escalating to uncapped -> a guarded fit is NEVER
  worse than the shared-ceiling fit, at worst it costs the extra fits. Wired into
  `nmix_laplace()`, `nmix_laplace_re()` + the `ms_abun() + latent()` path (scores at
  the predictor the fit ran on, field + factor offsets included -- a latent surface can
  push a site above what its own counts suggest). NB / zero-inflated keep the shared
  ceiling: the check runs a Poisson marginal and would understate a heavier tail. The
  dedicated C++ areal community path (#12) + the removal / distance kernels stay
  uncapped -- a capped oracle beside an uncapped field solve in ONE fit is worse than
  neither.
  **A capped fit matching its uncapped twin is NOT evidence the guard ran.**
  `.nmix_community_score_gap()` sits behind a `tryCatch(error = NA)` and
  `is.finite(NA)` is FALSE -> an ERRORING guard and a PASSING guard both leave the fit
  untouched and look identical from outside; assert on the gap itself. What it buys is
  confined to heavy-tailed seeds, so this is a tier-3 fix, NOT a way to shrink
  `test-ms-abun-factor.R` -- its lever is `factor.starts`. Numbers in
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
  by `summary` like `log_r`, with no arm in `tidy()`); `fit$zi_omega`
  = omega, `fit$zero_inflated`. Box-bounds only the pathological corners
  (`logit_omega`, `log_r`) so a degenerate seed can't diverge on the ZINB
  zero-source ridge; betas unbounded (no interior bias). ZIP recovers cleanly
  (betas/omega, coverage nominal); ZINB weakly identified (structural omega vs
  NB overdispersion both absorb zeros), recovers only in a higher-count regime.
  Gated at `.dispatch_abun`
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
- **Scalar nuisance AGHQ blocks floored at 3 nodes** (`.TOBS_MIN_SCALAR_NQUAD`,
  `.nmix_scalar_nquad()`, #234): `log_r` (NB) + `logit_omega` (ZI). A 2-node
  Gauss-Hermite rule places 2 nodes and has no freedom left for curvature, so on a
  non-Gaussian 1-D posterior the marginal can come out arbitrarily sharp -- a
  converged fit, an ordinary point estimate, no warning, and every downstream Wald CI
  / SBC rank inheriting a collapsed SE. Converged at 3, so a floor NOT a clamp -- and
  the floor beats a smaller `n_quad` (a coefficient block may run at plain Laplace, a
  block whose marginal SE is REPORTED may not). Also a bias fix, not just a seed
  rescue.
- **Scalar community variance components GATED against their boundary (#250).** A
  collapsed component is invisible on its own fit (converged, no warning) and the
  community mean it scales inherits a shrinking interval. Gate = the component's OWN
  uncertainty, NEVER a cut on `sigma_hat`: `.tobs_aghq_variance_boundary()`
  (`R/re_aghq.R`) = boundary-aware Wald `1/SE(log sigma)` off tulpa `re_par_se`,
  declining on a non-scalar block; `.tobs_warn_variance_boundary()` = ONE warning per
  fit, wired `nmix_laplace_re.R` (`sigma_log_r`, `sigma_omega`). Opt-in PC prior
  (`control$logr.sigma.prior` / `omega.sigma.prior`); default pure ML. NOT the whole
  guard: a failed solve leaves sigma at its INITIAL value, caught by the solve status
  (#281).
- **`mu_log_r` Wald calibration (#235/#280/#285).** "~24% too narrow" does NOT survive:
  put the seed block's realized draw at its expectation AND rebuild the SE at the
  simulated sigma -> calibrated at EVERY group count; the `S = 18` spike is the DRAW,
  not the estimator. Real: `sigma_hat` attenuates (monotone in S), SE built on it ->
  `?ms_abun`. A coverage floor here CANNOT be a tripwire (the SE at TRUE sigma still
  misses -- draw excess untouched), so `test-ms-abun-nb-rs-coverage.R` asserts
  unbiasedness, normality, a gross floor + a WIDE draw-corrected scale band. Evidence
  (+ what was refuted en route) + numbers: `NOTES_families.md` /
  `NOTES_measurements.md`.
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
field recovers); icar/bym2 also sample + centre via the #71 sum-to-zero
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

### Stats accessors: ONE layout on every fit (#332)

`coef()` = flat named vector, names `<arm>_<term>` exactly as `vcov()` /
`confint()` / `summary()` rownames; `coef(fit, arm = "psi")` = one arm,
unprefixed. `tidy()` = `arm, term, estimate, std.error, conf.low, conf.high`
(arm `NA` for a coordinate in no arm: `log_r`, `log_sigma_pos`, AR1 hypers).
Arms + their terms come from `.tobs_fit_arms()` (process_info `coef_names`,
detection arm also owns `visit_<cov>`; multi-arm fits `presence`/`positive`/
`class`), split by `.tobs_split_terms()` on DECLARED terms, never prefix alone
(`mu_log_r` beside arm `mu`; `f_sp1_sp2` beside `f_sp1`). `glance()` = `nobs, df,
logLik, n_fixed, n_samples, n_divergent, mean_accept, converged` (+ outer grid,
pareto) via `.tobs_glance_layout()`; multi-arm fits too. `summary()` = one data
frame on every method (no `summary.cover_fit`). `simulate()` = ALWAYS a list of
`nsim` (`sim_1`..), `"seed"` attr, stats:::simulate.lm RNG protocol; a family
handler `.tobs_simulate_<key>(object, nsim)` returns the list, never unwraps at
nsim = 1. `residuals()` keeps its LEVEL contract `list(occ, det)` -- `det` is
scored against `z * p`, not a `fitted()` quantity, so it is not renamed.

### Diagnostic doors = S3 methods, NEVER a second name

ONE verb per diagnostic, owned by whoever owns the CONCEPT; tulpaObs registers a
`tobs_fit` method + re-exports the generic (`R/reexports.R`), so `library(tulpaObs)`
alone reaches every door and no session masks anything. Roster:
`waic()`/`loo()` = **loo**'s generics (so `loo_compare()` reads a `tobs_fit`;
`loo.tobs_fit` returns a real `psis_loo` via `.tobs_loo_one()`, shared w/
`tobs_stack()`; `waic.tobs_fit` returns `loo::waic()`'s object, NOT a
`tulpa_criteria` -- read `$estimates`, loo deprecates `$waic`, #333); `sbc()` `pit_residuals()` `test_uniformity()` `test_dispersion()`
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

`n-L` = nested_laplace. Support matrix only: per-row detail (mechanics, gates, test files, what recovers) is in `NOTES_families.md` "Feature roster: per-row detail". `.tobs_family_methods` (`R/tobs_helpers.R`) stays the authority on which method a family accepts.

| Feature | Laplace | NUTS |
|---|---|---|
| Single-season occupancy | Yes | Yes |
| Dynamic (HMM) | Yes | Yes |
| Multi-season AR1 year effect (`t_occu`) | — | — (`pg_gibbs` only) |
| Community single-season (`ms_occu`) | Yes | Yes |
| Community dynamic (`ms_dyn_occu`) | Yes | Yes |
| Community integrated (`ms_int_occu`) | Yes | Yes |
| Integrated multi-source | Yes | Yes |
| Multi-season integrated occupancy | Yes | — |
| JSDM (`jsdm`) | Yes | Yes |
| Cover hurdle (joint) | Yes | Yes |
| Cover hurdle spatial coef fields (`\|\|` / `\|`) | n-L | — |
| Cover hurdle arm-specific fields (single-arm `to`) | n-L | — |
| Joint occu + cover | Yes | Yes |
| occu + cover + areal field + per-group RE | n-L | — |
| occu + cover independent cover-arm field (single-arm `to="positive"`) | n-L | — |
| occu + cover shared field + arm deviation (`share(residual=)`) | n-L | — |
| occu + cover independent detection-arm field (single-arm `to="detection"`) | n-L | — |
| occu + cover obs-arm RE (detection / pos) | n-L | Yes |
| occu + cover + detection / cover-arm RE | n-L | — |
| Community joint occu + cover | Yes | Yes |
| Spatial-factor community occu+cover (JSDM) | Yes | Yes |
| Multiscale occu + cover | Yes | Yes |
| Count / relative-abundance GLMM | Yes | — |
| Community spatial-factor (`ms_count`) | n-L | — |
| Community latent factor (`ms_count`) | Yes | — |
| Community SVC (`ms_count`) | n-L | — |
| Community count + areal (`ms_count`) | n-L | — |
| Community count (`ms_count`) | Yes | Yes |
| Count + areal / SPDE / GP | n-L | — |
| N-mixture (`abun`) | Yes | Yes |
| N-mixture + areal | n-L | Yes |
| Community N-mixture | Yes | Yes |
| Community N-mixture + areal | n-L | Yes |
| Community N-mixture latent factor (`ms_abun`) | Yes | — |
| Community N-mixture spatial-factor (`ms_abun`) | n-L | — |
| N-mixture + grouped RE | Yes | — |
| Community distance (`ms_distance`) | Yes | — |
| Community distance latent factor (`ms_distance`) | Yes | — |
| Community distance spatial-factor (`ms_distance`) | n-L | — |
| Removal (Pois/NB) | Yes | Yes |
| Distance (Pois/NB) | Yes | Yes |
| False-positive occupancy | Yes | Yes |
| Presence + nominal class | Yes | — |
| Royle-Nichols occupancy | Yes | — |
| Time-to-detection occupancy | Yes | — |
| Double-observer abundance | Yes | — |
| Joint distance + removal | Yes | — |
| Open-population distance sampling | Yes | — |
| Multi-species co-occurrence occupancy | Yes | — |
| Open N-mixture (Dail-Madsen) | Yes | Yes |
| Spatial ICAR/BYM2/NNGP | — | Yes |
| Spatial + dynamic | — | Yes |
| Nested-Laplace (areal) | n-L | — |
| NA-response prediction | n-L | — |
| Formula RE (intercept) | Yes | Yes |
| Formula RE (uncorr slope) | Yes | Yes |
| Formula RE (corr slope) | Yes | Yes |
| Formula RE on detection | Yes | Yes |
| All S3 methods | Yes | Yes |
| Diagnostics | Yes | Yes |
| Simulation | Yes | Yes |
| Spatial prediction | — | Yes |
| Components | Yes | Yes |

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

**Areal-field hypers are SAMPLED, not pinned (#204).** A fit conditioned on the outer
grid's own point estimate cannot serve as an independent reference for that grid, so
under `nuts` + an areal term on psi the field SD, the bym2/car mixing rho and the copy
amplitude alpha are sampled alongside the field. Every areal kind's loading factors as
a FIXED basis w/ hyper-dependent column weights, so a leapfrog step costs a rescale,
never a re-decomposition (`src/nuts_field_hyper.h`, R mirror `.ochf_*` in
`occu_cover_nuts.R`) -- which is why car_proper rho is NOT the O(n^3) per-step
Cholesky the issue scoped it as. Each hyper rides a bounded transform (no wall,
Jacobian in the target) and its **prior is FLAT over the WARM FIT'S OWN outer-grid
span** = the measure the nested-Laplace grid integrates against, so
`control$sigma.grid`/`alpha.grid`/`rho.car.grid` move BOTH backends; alpha's grid `0`
atom is not HMC-representable (positive nodes only), icar pins `rho=1`, and an axis the
grid pinned to one node stays pinned. `control$fixed.hyper=TRUE` restores #74
conditioning, byte-identical to the old loading (pinned = the degenerate configuration
of the same block, not a second path). `fit$nuts$sampled_hyper` / `$fixed_hyper` =
CHARACTER vectors, empty when nothing pinned (deliberately a type change from the old
`fixed_hyper=TRUE`, so a stale `isTRUE()` read fails loudly) + `$fixed_hyper_values`;
`fit$hyper_draws` cols `sigma|rho|alpha|field_sd`. **`field_sd`** = the geo-mean
marginal SD the block implies at that draw = the ONLY field-scale summary comparable
across kinds (the three normalise their precisions differently), so it is what a
simulation truth is stated in. `inv_metric` MUST be sized `n_raw + n_hyper`: the engine
takes a short metric pointer w/o a length check and bym2's `n_raw = 2n-1`.
**Sampled terms reach the criteria (#211/#215).** A structured term the scorer cannot
see is scored at ZERO. `.tobs_occu_cover_components()` returns the sampled field (per
SITE `field_occ`/`field_pos`) and the sampled obs-arm RE (per VISIT
`off_det`/`off_pos`) as OFFSETS beside the coefficient draws, and all four diagnostics
fold them in -- ploglik (WAIC/LOO/CPO), PPC, PIT/LOO-PIT -- because the per-visit
offset enters the SHARED `Arms` view (`src/occu_cover_ragged.h`), so ONE change reaches
all three kernels and dense == compact by construction (#185). A 0-column matrix = "arm
carries none" -> null pointer -> the no-offset path byte-identical. Cell-aggregated
cover scores one cover row per detected UNIT, so a per-visit offset errors there w/ a
pointer. The grid-integrated (`nested_laplace`) route reaches the same criteria; an
occupancy-arm RE (#56) is per SITE, so its group codes sit on `model$re_psi` and the
builder adds the per-group draw to `field_occ`.
**Copy axis (#287):** `control$alpha.n[.trend]` = RESOLUTION of the engine's alpha
axis; `alpha.grid` STATES nodes. One per block. See `NOTES_families.md`.

**A SECOND (SVC / trend) field samples too (#214)**: the block is a LIST, each field
carrying its own basis, site->node map, per-site design WEIGHT (`field_weight`; absent
= the intercept field) and its own sampled `(sigma, rho, alpha)` -- two fields share no
hyper. The warm fit is the multi-block coupled path and FORCES `integration = "grid"`:
above 3 axes the engine switches to a mode-centred CCD star whose column range is a
design radius, not an integrated span, and the sampler reads each axis's span as its
flat prior's support. Reported: `fit$trend_field`/`trend_fields` (named by weight
column), `fit$trend_field_draws`, per-block suffixed hypers (`sigma_trend`,
`alpha_trend`, `field_sd_trend`, indexed when several) in `hyper_draws` /
`sampled_hyper` / `fixed_hyper`, and `fit$spatial$field_suffix`/`field_weights` so
`.tobs_occu_cover_sampled_field()` sums every block's loading into the criteria.
Loading algebra, per-kind bases, validation protocol + grid thinning:
`NOTES_families.md`. `test-occu-cover-nuts-ic.R` + `test-occu-cover-nuts-svc.R`.
Correlated `|` (one free-Sigma MCAR block), temporal + RE still gated -> n-L.
group_var maps sites>cells; predict() needs the joint object (non-spatial laplace AND
nuts both error w/ pointer); sampled-field (estimated-variance) route =
`ms_occu_cover()` factor (tulpa#67).
**Spatial default** (`nested_laplace`,
`R/occu_cover_joint.R`): `joint` engine via `tulpa_nested_laplace_joint()` w/
`occu_cover_{lognormal,beta}` cell-coupling spec (tulpa#32) — 3-arm joint
nested-Laplace, outer-grid over `(sigma, alpha)`, per-cell occupancy mixture
closed-form derivs drive inner Newton. Much faster than v3_nested, completes at
sizes v3_nested does not. Lognormal + beta recovery
(`test-occu-cover-joint-coupled.R`); status `"working"` (#96). Shared-field occ
SLOPE Wald CI mildly anti-conservative small-N (NUTS non-spatial calibrated). Outer Pareto-k diagnostic (`control$diagnose.k`) defaults OFF
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
dense grid flattened site-major, visit-ascending (the order the dense `rowSums`
accumulates). It also derives `n_valid` + `any_det`, so the per-site detection summary
has ONE definition. All three consumers -- pointwise loglik, PPC, PIT/LOO-PIT CDF
limits -- go through it, and the three kernels assemble per-draw predictors from one
shared `Arms` view (`src/occu_cover_ragged.h`) => dense == compact TO THE BIT. Reading
the dense grid instead is what made `cpo()`/`ppc()` error on every compact fit while
`waic()` worked. The aggregated (mean/median/latent) PPC keeps the padded grid +
`cpp_occu_cover_ppc_agg` -- aggregation is dense-only by gate. Cover density gates on
`detected AND finite` everywhere (a detected visit may carry NA cover,
missing-at-random); the PPC used to score that NA and returned `fit.y = NA` for every
draw. Kernel/consumer detail: `NOTES_families.md`.

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
random intercept(s) and slopes on `detection=`/`positive=` (`(1|g)`/`re(g)`, crossed
`(1|g)+(1|h)`, nested `(1|g/h)`, uncorr `(x||g)`/`(0+x|g)`, corr `(1+x|g)`), per-VISIT
grouping (one code per (site,visit)), composing w/ the required psi field on the
nested_laplace joint path. Intercept -> one scalar iid block; uncorr slope -> one
weighted iid block per coef (tulpa `svc_weight` = `Z[,c]`); corr -> one `miid`
free-Sigma block (tulpa#114, needs tulpa>=0.0.39 -- the DESCRIPTION floor enforces it,
no gate). Slope covariates STANDARDIZED to unit SD so the fixed Sigma grid is
scale-invariant; BLUP/sigma + predict draws back-transformed `/scale` (cor
scale-free). **Key fix (#102)**: the det arm's `field_coef`=1 when `model$re_det` is
present so the iid block scatters; the shared field stays off detection by its
`spatial_idx=0` sentinel -- NO tulpa engine change. Reports
`sigma_re_p`/`sigma_re_pos` (+`_<var>` when >1 term/arm, `_<coef>` per coef,
`cor_re_*_<ci>_<cj>`) + BLUPs; `fit$re` = flat list keyed by arm (lone term) or
`"<arm>:<var>"` (crossed). predict SUMS a term's offsets on its arm, slope coef
weighted by `newdata[[coef_name]]` (intercept=1), unseen level -> 0 = pop mean. Gates:
obs RE needs `nested_laplace`; pos RE needs `cover_aggregate="none"`; not composed w/
correlated MCAR / latent cover / batch. Corr slope / crossed grow the grid ->
`control$integration="ccd"`. sigma carries binary small-cluster attenuation (lower
bound); BLUPs recover (cor>0.5). Knobs `re.sigma.grid.p`/`.pos`,
`re.logchol.grid.p`/`.pos`. Parse / design / block-layout / postprocess mechanics:
`NOTES_families.md`. `test-occu-cover-obs-re.R`.

### `ms_occu_cover()` detail

Community `occu_cover()` (`R/ms_occu_cover.R`); per-species coef RE w/ Gaussian
community covariances across psi/p/cover arms, shared dispersion. Latent z integrates
closed-form per species-cell (reuses `.occu_cover_site_ll`); per-species deviations by
in-tree pure-R Laplace-EM — arrowhead joint Newton (RE Schur-folded, analytic grads
`.occu_cover_eta_grad`) + closed-form community-cov M-step. Community-mean SEs =
marginal observed info (Louis 1982). Beta, lognormal, AND identity-Gaussian
(`ms_occu_cover("gaussian")`, #127: delta-normal magnitude, `mu = eta`, shared residual
`sigma_pos`; recovers to lognormal parity). Non-spatial Laplace ONLY (structured term
on any arm errors + pointer). Recovery + coverage (`test-ms-occu-cover.R`);
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

**MAR cover, same rule as the twin** (#262): a detected visit (`y=1`) w/ NA cover keeps
its detection term + drops only `f_pos` -- SHARED
`.occu_cover_validate_pos_values()` builder, `isfinite(y_pos)` gate in the
cell-coupling spec, the NUTS target AND the R marginal. A detected plot factorises, so
under `laplace` the psi/theta/p estimates are the FULL-DATA ones (assert
`expect_equal`, not "close"); under `nested_laplace` the shared field couples them, so
close only. The dispersion pre-fit must filter `is.finite(pos_vals)` -- an unfiltered
`sd(log(pos_vals))` is where an NA reaches the starts.
`test-occu-multiscale-cover-mar-cover.R`.

**Detected-plot score written ONCE** (#270): `mscale_det_plot_block()`
(`occu_coupling_shared.h`), beside its `mscale_nodet_cell()` sibling; the coupling spec
and the NUTS target each supply eta/y accessors + `emit_p`/`emit_pos` sinks, so neither
allocates and each keeps its own row layout. Cover-arm dispatch via `PosPolicyAccess<P>`
(spec, compile-time) / `PosCodeAccess` (sampler, runtime code) -- do NOT collapse those
two into one runtime switch, it would put a branch in the Laplace inner loop.

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
`n_iter - n_warmup` draws). Under `pg_gibbs` it means the OPPOSITE
(TOTAL sweeps, warmup taken out of it) -- see `engine_defaults.R`.

**Convergence diagnostics on EVERY NUTS path** (#174): ONE writer
`.tobs_nuts_attach_convergence()` (`R/nuts_chains.R`) fills the record
`summary.tobs_fit` + `print.tobs_fit` read -- `convergence$parameter/rhat/
ess_bulk/ess_tail` -- plus scalars `fit$max_rhat`/`fit$min_ess`. tulpa owns the
estimator (`tulpa::diagnostics()`: rank-normalized split-Rhat, bulk ESS, 5%/95%
tail-indicator ESS, Vehtari 2021); tulpaObs only hands it per-chain matrices +
names. Wired at shared choke points, one per family group (list:
`NOTES_families.md`). `convergence$parameter` MUST carry the names `summary()`
puts on its rows, else the record is present but unreadable -> the writer takes
`par_names` (the fit's `fixed_names`) + `cols` (the sampler coordinates those
name). Most community builders report moment-matched pseudo-draws around the
posterior mean, NOT the chain -> their record comes from the sampler's own
`rc$chains` at `par_cols`; computing it off `fit$draws` there = diagnosing a chain
never run. `converged` on a sampled fit = every reported split-Rhat < 1.01
(print's warn threshold), NOT the warm-start optimiser flag.
`test-nuts-convergence-contract.R` fits every family advertising `nuts` and
asserts the record resolves an Rhat for every `summary()` row; its name-set check
fails when a new NUTS family lands w/o a case. `.tobs_nuts_rhat_ess()` reads the
same table. PG-Gibbs writes the same record through the same writer, from the
per-chain matrices `.tobs_pg_summarize()` carries (#319); `converged` is `NA`
when the chains are too short for the estimator, never the old Rhat < 1.1 rule.

**Every `tobs()` fit passes one tail** (`R/tobs.R`, after dispatch):
`.tobs_default_field_eta_offset()` derives the state-arm field offset from the
field a fit REPORTS when its fitter recorded none (declines on any existing
offset slot, no field, or a map it cannot line up, e.g. `group_var`), then
`.tobs_attach_sampled_loglik()` fills `log_lik` ONLY where `logLik()` has no
finite value (gate = `logLik()`, NOT the slot: NUTS reports `mean(log_prob)`
with `log_lik` NA). `.tobs_eta_draws()` adds the offset itself, so every
process-major ploglik kernel is field-aware; custom-layout kernels call
`.tobs_add_eta_offset()` by hand, and the two never overlap (#318/#320).

**`control` knobs are read by EXACT key when another accepted key extends
them** (`control[["n.threads"]]`, not `control$n.threads` beside
`n.threads.outer`). `test-control-exact-keys.R` walks every namespace function's
parse tree and fails on a new prefix-shadowed `$` read (#322).

**Areal field on the occu NUTS path** (#142): `occu() + icar()/bym2()` under
`method="nuts"` EXPOSES its field. `occu_fit.cpp` emits `spatial_layout` (engine
ParamLayout) + names the columns (`spatial_field[i]`/`spatial_theta[i]`/
`log_tau_spatial`/`log_sigma_spatial`/`logit_rho_spatial`, no more `param[k]`);
`.tobs_areal_field()` (occu_fit.R) sets `fit$spatial_field` = centred per-cell
surface for icar/car_proper (level confounded w/ intercept -> centred). bym2
field = Riebler rho-mix of both blocks x graph scale factor -> columns named,
reconstruction left to the draws. icar field recovers
(`test-occu-areal-nuts-recovery.R`). car_proper NOT a wired occu-NUTS term
(errors at fit, pre-existing).

**SVC = two flavors behind ONE verb (#118, unified #146).** Both written under
`spatial()`; `model =` picks which. `spatial(~ 1 + w || cell, graph)` = areal,
`spatial(lon, lat, model = "svc", coefficients = )` = continuous. `svc()` stays as a
direct ctor (like `icar()`). The continuous form selects columns BY NAME
(`coefficients = c("(Intercept)", "elev")`), matched against the arm design at fit
time by `.tobs_svc_columns()` (`R/occu_svc.R`, single source of truth -- BOTH the
Laplace field builder `.tobs_svc_field_blocks()` and the NUTS packer in
`R/occu_fit.R` call it, so the two backends cannot drift). `indices =` (column
positions) is the lower-level form, still supported; giving both errors. Fits report
`fit$svc_indices` + `fit$svc_coefficients` on both backends. A bar with
`model = "svc"` errors with a pointer to the continuous form.
`test-svc-spatial-umbrella.R`. The two flavors are still DIFFERENT ENGINES with
different coverage -- do not conflate what they fit, only how they are spelled:
- **Areal spatially-varying coefficient** (the `svcPGOcc` analogue): a WEIGHTED
  areal bar `spatial(~ 1 + w || cell, graph)` on `occu()` via `nested_laplace`,
  rerouted through the joint direct-grid engine (`.tobs_fit_occu_joint`, #81).
  RECOVERY-TESTED: `test-occu-spatial-svc-recovery.R` +
  `test-occu-svc-joint-recovery.R` recover the known intercept + trend surfaces
  (`fit$spatial_field`/`fit$trend_field`, `cor > 0.75`), the SD hyperparameters and
  the coefficients. This arm arrives as a `spatial` term, not `svc`.
- **Continuous NNGP `svc()` term** (`spatial(lon, lat, model="svc",
  coefficients=)`, or direct `svc()`): wired to single-season `occu()` on NUTS
  (`populate_svc`, `src/populate_helpers.h` -> `data.svc_data`) AND on the
  deterministic backends (#143). On NUTS the eta-assembly + NNGP prior + gradient
  live in the compiled UPSTREAM tulpa engine, NOT tulpaObs. Surface exposed as
  `fit$svc_field` (n_obs vector / n_obs x n_svc matrix of posterior means, per-draw
  on `attr(., "draws")`), sliced by position from `fit$svc_layout`; block named
  (`svc_w[i,j]`, `log_sigma2_svc[j]`, `log_phi_svc[j]`). NEEDS
  `prior_range = c(r0, alpha)` (PC prior, `P(range < r0) = alpha`) -- tulpa ships the
  anchors unset and refuses without them; neither package defaults them.
  **RECOVERY-VALIDATED** (#119, tulpa 0.0.82, `test-occu-svc-nngp-recovery.R`):
  divergences fell to 0 every seed and phi recovered onto truth, after two upstream
  fixes -- a Uniform-behind-a-wall range prior (gcol33/tulpa#144) and an improper
  half-Cauchy on the SVC marginal SD (nothing bounded sigma above); fixing the SD
  prior is what pulled phi on, the two being ends of the GP ridge. Surface cor did
  NOT move -- information-bounded at these settings, not sampler-bounded, so the
  calibration test asserts divergences + phi + sigma and deliberately NOT cor. tulpa
  CANNOT make this measurement: `svc()` is a tulpaObs term and no tulpa-side fit
  reaches the NNGP SVC path -- how #144 survived there. Numbers in
  `NOTES_measurements.md`.
  **Laplace backends (#143, `R/occu_svc.R`)**: `occu() + svc()` also fits under
  `method="laplace"` / `"nested_laplace"`. K surfaces = latent field blocks on the psi
  logit riding the SHARED areal-BFGS driver (`.tobs_areal_bfgs_fit`,
  `R/areal_bfgs.R`); the Vecchia precision is assembled in R from the term's OWN
  neighbour structure w/ the compiled kernel's kernels/jitter/variance floor, so both
  backends integrate the SAME density, asserted == tulpa `cpp_test_svc_nngp_twins` to
  1e-8 (mechanics: `NOTES_families.md`). Hypers (sigma, phi) grid-integrated on both
  routes (`laplace` == `nested_laplace` here) -> `fit$svc_hyper`; surface ->
  `fit$svc_field` (NUTS naming). Surface cor matches the NUTS path on the same truth
  (information-bounded, NOT backend-bounded). `fitted()` adds the surface in-sample via
  `model$occ_eta_offset`; `predict(newdata=)` does NOT krige to new locations (as on
  NUTS). Gated on occu(): detection-arm svc, a spatial/temporal/re term alongside svc,
  `pg_gibbs` -- all error w/ pointer. `test-occu-svc-laplace-recovery.R` +
  `test-svc-guard.R`.
  **Observation families (#144, `laplace`/`nested_laplace`)**: `removal()`,
  `distance()`, `fp_occu()`, `dyn_abun()` carry `svc()` too, with NO family-specific
  code -- those four already ride `.tobs_areal_bfgs_fit`, and an svc surface IS just
  another latent block on the arm their `eval(theta, offset)` already exposes, so the
  wiring is `.tobs_svc_field_blocks()` + `.tobs_build_field_spec(svc=, X_svc=)` +
  `.tobs_attach_field_results(svc=, has_spatial=)` (mechanics: `NOTES_families.md`).
  Composes with an areal and/or temporal field on the same arm. Surfaces load on the
  STATE arm only (log lambda / psi); a detection-arm areal field alongside svc errors
  (`.tobs_check_svc_arm()` -- the driver exposes ONE `grad_eta`, so the surfaces would
  otherwise be fit against the detection arm). The N-mixture families
  (`abun`/`ms_abun`) do NOT get it: their areal path is the C++ count-spatial driver,
  not this one. NUTS still errors everywhere but single-season occu(). Surface cor in
  `NOTES_measurements.md`. `test-svc-families-recovery.R`. The driver also returns
  `res$eta_offset`, the marginalised per-observation offset the blocks jointly load;
  the family wrappers read it instead of re-deriving each block's site map, which also
  fixed the temporal-only fp_occu path (it indexed a length-`n_t` field by site).

## Performance

Penalized EM (`laplace`, default) trades speed to break the psi-p identifiability
ridge at small J; Gibbs/MI add a Rubin-pooled correction on top. Wall-time comparison
vs `laplace_gibbs`/`nuts`/inlaocc/spOccupancy: `NOTES_measurements.md`.

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
  `tulpa::tulpa_iter_progress()` (R/progress_iter.R).
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
    **Every sampler knob is a `NULL` formal filled by `.tobs_fill_sampler(environment(), engine, ...)` as the fitter's first statement (#188)** -- a fitter restating a literal has a DEAD formal wherever a caller forwards explicit values (`.tobs_fit_model()` does), so its own default never runs. `test-engine-defaults.R` asserts structurally that no fitter carries a literal + each calls the filler -> a new sampler family cannot reopen it.
    Divergences RECORDED (`.TOBS_SINGLE_SPECIES_NUTS` / `_LAPLACE`, via `.tobs_single_species_defaults(engine)`): the `.tobs_fit_model()` entry keeps its own `sigma.beta` / `adapt.delta` / `seed` -- each the value the single-species recovery/coverage tests were calibrated against, so the override belongs to the ENTRY, not the nine families passing through (values: `NOTES_families.md`).
    **`n.iter` means the OPPOSITE thing under `pg_gibbs`**: TOTAL sweeps, warmup taken out of it (kept = `length(seq.int(n.warmup+1, n.iter, by=n.thin))` = 1500); under `nuts` it is the KEPT count and the run = `n.iter + n.warmup`. Stated in the table's `pg_gibbs` block, pinned by a test recomputing the kept count from the profile.
    **`sd.load` (1.0) + `re.lkj` (1.5) = laplace rows (#189)**; `n.quad` NOT one number, deliberately -- `.TOBS_NQUAD_ROUTES` / `.tobs_n_quad(route)` enumerates seven marginals, `?tobs` lists them per route, and there is no single default most routes do not use.
  occu_categorical.R        — presence + nominal K-class hurdle (#106); Bernoulli presence arm + baseline-category multinomial logit class arm over tulpa `cpp_multinomial_logit_terms`
  occu_cover_dispatch.R     — formula-native cross-arm share()/spatial-field DAG coupling dispatch for occu_cover()
  occu_joint.R              — standalone occu() SVC-spatial-bar nested-Laplace path (#81): occu_cover()'s joint direct-grid engine with the cover arm removed
  cover_hurdle_joint.R      — the joint nested-Laplace cover fit (lognormal / beta), incl. occurrence-arm suff-stat aggregation `.cover_aggregate_occ`. One of the largest files in the package
  cover_nuts.R / src/cover_nuts.cpp — standalone cover() NUTS
  joint_substrate.R         — shared joint-fit draw / Pareto-k substrate (.tobs_joint_fit, .tobs_joint_draws*) across the occu_cover / cover nested-Laplace joint routes.
    **Outer-grid placement promoted (#187)**: `.tobs_promote_outer_grid(jf)` lifts `outer_grid_placement` ("fixed"/"auto_recentered"), `_recenter_attempts`, `_prior_added` and `_recenter_declined` (reason a "fixed" placement stayed fixed, tulpa#293) to `tobs_fit` top level, everywhere `.tobs_promote_pareto_k` is spliced, and `.tobs_glance_outer_grid()` adds two columns to BOTH glance methods. NOT gated on the grid having MOVED -- a declined recenter is exactly the case worth seeing.
  pg_gibbs_shared.R         — shared Polya-Gamma machinery (.tobs_pg_draw_beta, .tobs_pg_community_update, .tobs_pg_finalize_fit) behind EVERY pg_gibbs fitter
  stacking.R                — `tobs_stack()`, LOO-weighted predictive stacking across fitted tobs_fit objects of any family
  sbc.R                     — `sbc.tobs_fit()` (#207): posterior SBC, the `tobs_fit` method of `tulpa::sbc()`. Registry `.TOBS_SBC_REGISTRY` per family (simulate/refit/draws/+statistic); pooling, group_ids, arms, controls shared. Replicate on FRESH cells (block-diag graph) -> premise VERIFIED. occu_cover registered+verified; field SD not gated (#212). Numbers in NOTES_measurements.md
  aggregation_scan.R        — cell-size / NN-spacing / changepoint scan utilities
  field_offset.R            — the fitted latent field as a per-ARM eta offset (#254). `model$field_eta_offset` = one entry per process, in THAT process's own design row layout (removal's per-PASS detection arm included), written at the ONE dispatch tail `.tobs_finalize_family_fit()` (`.tobs_attach_model_eta_offset()`) -- after `fit$model <- model` drops whatever the fitter set on its autoscaled copy. EVERY door that rebuilds eta from `X_processes` reads it through `.tobs_eta_offset()` / `.tobs_add_eta_offset()` / `.tobs_sim_arm_block()` -- fitted / residuals / predict / simulate + the ploglik behind WAIC / LOO / DIC / CPO -- so a field cannot be scored at 0 again, and a field-free fit stays byte-identical. `.tobs_nuts_field_loglik()` (same tail) re-evaluates `log_lik`/`log_prob` on a SAMPLED-field fit through that ploglik kernel: the field is not a draw column, so the family's own posterior-mean marginal ran it at offset 0 and `logLik()`/AIC/BIC described a model WAIC did not score (grid-integrated paths already carry the field in their marginal, keep their grid-weighted value). `count()` + the community families keep their own named slots (`count_field_offset`, `occu_field_offset`, ...); the shared field -> per-site contribution reader is `.tobs_spatial_field_offset()` (was `.count_spatial_field_offset`). Writer list + simulate trick: `NOTES_families.md`
  community_em.R            — shared community Laplace-EM .tobs_community_em() (ms_occu/ms_dyn_occu/ms_int_occu/ms_count/jsdm/ms_abun-latent/ms_distance). Three OPTIONAL args, each defaulting to the previous behaviour byte-identically: `sp_info` (analytic per-species observed info; the FD fallback costs 2(P+G) marginal sweeps per species per Newton step -- supply it when the kernel exposes the Louis block), `init_b`/`init_Sigma` (warm start; a block-coordinate caller re-enters the EM once per outer pass, so a cold restart each time dominates)
  ms_{occu,dyn_occu,int_occu}.R      — community single/dynamic/integrated
  abun.R / abun_nuts.R               — nmix family + non-spatial NUTS (#41)
  count_spatial.R / count_methods.R  — areal count fitter .tobs_fit_count_spatial (#117); fitted/predict/residuals
  ms_count.R                — community count / relative-abundance GLMM (msAbund, #117); .tobs_fit_ms_count over shared community_em.R
  ms_count_nuts.R / src/ms_count_nuts.cpp — community count NUTS (msAbund NUTS, #117); in-tree C++ FullGradFn (reduced ms_abun_nuts: no detection/latent-N), R oracle .tobs_ms_count_nuts_logpost, warm-start from Laplace-EM
  ms_count_spatial.R        — community count + shared areal field (sfMsAbund) + SVC bar (svcMsAbund, #117/#118); block coordinate ascent (community EM offset <-> multi-field Poisson-ICAR), pure R
  community_latent.R        — SHARED latent-structure engine for EVERY community family (#119/#120/#121): one block-coordinate ascent (community EM w/ the latent as an offset <-> field / factor updates) + the areal Newton, the factor update, bym2/car_proper/spde hyper grids. A family supplies ONE callback `working(eta) -> list(score, curv)` (per-(site,species) score+curvature wrt an additive offset on the structured arm): Poisson `(y-mu, mu)`, occupancy two-state, Bernoulli `(y-psi, psi(1-psi))`. Field solve is `t(A) diag(w) A + tau Q`, so the site->node map slot takes an areal group_var incidence OR an spde barycentric projector unchanged. Adding a family to every latent route = one callback, not a new fitter. **The measured evidence behind every number below -- fixtures, seeds, wall times, per-family screens -- is in `NOTES_measurements.md`; the rules here are what it concluded.**

**Backtracking + guards.** Factor Newton (`.tobs_latent_factor_update`) backtracks and ridge-bumps singular curvature; non-finite guards are inline, not a named helper, and a non-finite `working()` score/curv `break`s the pass. `.tobs_latent_field_solve` has its OWN local `safe_solve()` (ridge retry for a singular Hessian only, Newton update unconditional) -- do not confuse the two, and there is no `safe_step()` anywhere in the repo. Mechanics: `NOTES_families.md`.

**Loadings by MARGINAL likelihood, NOT the joint mode (#153 -> #156).** Holding zeta at its joint mode makes `(zeta, lambda)` a joint-likelihood estimate with `Ns*Q` incidental params growing with the sample = Neyman-Scott, inconsistent: the site factors' estimation error lands in the fitted co-occurrence, lambda absorbs it, over-fit growing with Q/S. Fix = `.tobs_latent_factor_mmle()`, an EM on the SAME joint site marginal over all S*Q loadings (steps: `NOTES_families.md`; numbers: `NOTES_measurements.md`).

**ONE estimator, ONE state (#156).** `.tobs_latent_factor_update()` + `.tobs_latent_factor_scale()` run ONCE, outer pass 1, purely to INITIALIZE (the marginal's lambda-gradient vanishes at lambda=0, and the 1-D bracket is a global magnitude search the local EM cannot do). Running the joint-mode update every pass alongside the MMLE diverges both ways, so it never repeats.

**Offset by SCORE-MATCHING, NOT `zeta t(lambda)` (#156).** `.tobs_latent_factor_offset()` matches the plug-in score to the integrated one per cell, so it holds for ANY family with no link-specific derivation. `fit$model[[offset_slot]]` reads THIS, not `zeta t(lambda)`.

**Block-coordinate callers MUST warm-start `init_b`/`init_Sigma`.** Omitting cold-restarts every per-species deviation and both community covariances on every outer pass.

**`max.outer`: factor path 150, field path 25 by default.** Field reaches `tol` and breaks early; factor does NOT (alternates with the coefficient block along a slow mode) -- 25 leaves real community-mean bias. `factor.outer` is per-family, set from that family's own measurement: ms_count/jsdm/ms_occu 150, everything else 25 until measured. Do NOT globalize -- cost is not transferable.

**`factor.starts` (multi-start width) dominates a latent-N fit, NOT `max.outer`.** #157's basin escape runs K candidate starting directions on the first factor pass, each a full loading-EM to convergence. Per family, measured from that family's OWN recovery suite (a random seed screen alone is not enough -- the committed regression test #157 was built for caught a value ms_count's own screen missed): ms_abun 1, ms_occu 1, ms_count/jsdm + ms_distance at the driver default 8. NEVER copy a family's value -- cost and benefit both scale with how expensive one oracle eval is.

**NEITHER summary screens a fit alone.** `residual_cor` is row-normalized -> blind to a magnitude regression; `mag_ratio` is rotation-invariant -> blind to a direction regression. Screen on both.

**`n.quad` NOT threaded from any caller** -- driver default 5 is what every community latent fit actually uses; `control$n.quad` silently does nothing here.

TRAP surviving all of the above: magnitude MUST come from the JOINT marginal, never a per-species one (identifiability ridge between sigma and coef scale). Assert on `sqrt(tr(Sigma_res))` (rotation-invariant), as `test-ms-count-factor.R` / `test-ms-occu-factor.R` do
  ms_occu_field.R           — community occupancy SVC (svcMsPGOcc, #118); block coordinate ascent (community occ EM psi offset <-> two-state-marginal occupancy field solve), intercept + SVC field(s), pure R; plain intercept -> C++ ms_occu_spatial.R
  ms_abun.R / ms_abun_nuts.R         — community nmix + NUTS (#14)
  ms_abun_latent.R          — community nmix + latent() factors (lfMsNMix) / + shared field (spatial-factor); the ONLY new piece is the working oracle over the Royle marginal: score=grad_eta_lambda, curv=info_eta_lambda-var_N*score_wt_lambda^2 (the Louis (1982) (1,1) block = abundance curvature w/ the detection arm profiled out), which nmix_site_marginal() already exposes. Supplies `sp_info` (the design-sandwiched per-site Louis block) so the community EM skips its FD Hessian: 387s -> see table. Plain field w/o factors KEEPS the C++ #12 path
  ms_distance.R             — community binned distance sampling (msDS, #117) + latent() factors (lfMsDS) + shared field (sfMsDS). NO new C++: `cpp_distance_site_sweep` already returns log_lik/grad_lam/info_lam/var_N/swl, so the community EM reads its per-species score off it and the driver oracle is the SAME Louis formula as ms_abun. `.tobs_ms_distance_info_block()` assembles the per-species observed information from that sweep (#161) so the EM does not finite-difference it; the sign inside `v` differs from the N-mixture's and is asserted, not inherited -- see the spatial-factor row. Hazard-key log-shape = a community `global`, and that key keeps the FD fallback (its per-site cross terms are not exported). `simulate_ms_distance()` draws through `cpp_simulate_distance`, the kernel the likelihood integrates against -- a separate R-side quadrature simulates from a pi the model is not fit against and biases recovery. Poisson only; NUTS not wired
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
  community_spatial_em.h   — (#239) the Laplace-EM driver + field-structure layer (geometry/offset/log-prior/prior-gradient/sum-to-zero centering, icar/bym2/car_proper/spde) shared by the community areal drivers (nmix_community_spatial.cpp, ms_occu_spatial.cpp), templated on a per-family `SiteBlockFn` (one (species, site) cell -> log-lik/grad/curvature; nmix's reads nmix_kernel.h, occu's reads ms_occu_kernel.h). `run_community_spatial_grid`/`run_community_spatial_grid_spde` are the shared outer-grid drivers (areal vs continuous mesh); `community_pack_grid()` above still does the result packing. Only the per-site marginal cell function stays per-file
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
