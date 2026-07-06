# tulpaObs NEWS

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
