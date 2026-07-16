# dyn_occu + areal / SVC field under nested_laplace: field over-fit, no SEs

Status: FIXED. Two root causes, both resolved, both recovery-tested. dyn_occu +
icar and dyn_occu + an SVC bar (svcTPGOcc) now recover at nominal coverage, and
so do occu + icar and int_occu + icar.
Date: 2026-07-15 (found), 2026-07-16 (fixed).

Scope grew past the title: the same two root causes reached `int_occu` (which
this note first recorded as unaffected), and the single-season path -- the
"healthy" control throughout -- turned out to be untested rather than tested,
its only areal tests being shape smoke tests. Both now have recovery files.

## Root cause 1: the pseudo-binomial encoding (M)

`build_single_callbacks` / `build_dynamic_callbacks` encoded the state block as
`y = round(M*w)`, `n_trials = M` with **M = 1000** whenever the latent block was
AREAL, dropping to M = 4 only for a continuous SPDE block. A state row carries
**one** binary occupancy observation, so `n_trials = 1000` overstates its
information ~1000x: the field prior is swamped, pure between-cell binomial noise is
read as a real field, the field inflates, and the state slope inflates with it
through the logistic conditional-vs-marginal factor `sqrt(1 + 0.346 sigma^2)`.

The SPDE branch already documented this exact failure ("the field over-fits, and
the occupancy slope inflates ... keeping the per-site effective sample size O(1)").
The adjacent claim that areal blocks "tolerate the sharp M = 1000 encoding" was an
untested assertion and is false.

**M = 1 now applies to any nested latent block** -- one pseudo-trial per site, the
site's real information content.

`int_occu` was asserted here to inherit it "(same `.tobs_laplace_nested` encoding)".
It did NOT -- that claim was never measured and is false; see "Follow-ups" below.
The M rule was copy-pasted across the three `build_*_callbacks` and the integrated
copy had no latent branch at all. It is now one `.tobs_state_M()` /
`.tobs_encode_state_block()` that all three call.

M was swept on the real fits, 12-20 seeds (`_run_m_final.R`, `_run_m_single.R`),
truth slope 0.5 / f0 sd 1.0 / f1 sd 0.8:

| arm           | M    | slope  | coverage | f0 sd  | f1 sd  |
|---------------|------|--------|----------|--------|--------|
| dyn + icar    | 4    | 0.6232 | 0.83     | 1.3750 | -      |
| dyn + icar    | 1    | 0.5253 | 1.00     | 0.8618 | -      |
| dyn + SVC     | 4    | 0.7413 | 0.58     | 1.8212 | 1.8877 |
| dyn + SVC     | 1    | 0.4872 | 0.92     | 0.9231 | 0.8601 |
| single        | 4    | 0.5669 | 0.95     | 0.9295 | -      |
| single        | 1    | 0.5207 | 0.95     | 0.9294 | -      |

Monotone in M in every arm; M = 1 is uniformly best and NO arm regresses.

**A prediction that failed, recorded because it is presumably why the large M was
chosen:** M = 1 was expected to wreck the slope, since `round(w * 1)` is 0/1 and the
fractional resolution the inflation exists for is lost. It does collapse -- and the
slope IMPROVES anyway. The information-content error dominates the rounding loss.

**M = 4 was an intermediate wrong answer.** It was landed first, `#11` was marked
complete, and the leftover inflation (field 0.67 vs truth 0) was written up as a
permanent "bounded weakness" and nearly attributed to identifiability. It was not:
it was the same over-stated likelihood, undiluted. The SVC grid dump
(`_run_svc_grid.R`) is what broke that -- both tau axes pinned on the grid's
smallest-tau boundary (weights 0.756 / 0.446 on level 1 of 9) while the truths
(tau 1.0, 1.5625) sat comfortably INSIDE the grid, i.e. the marginal wanted fields
even larger and only the grid edge stopped it.

## Root cause 2: the SE gate

`use_louis` (R/laplace_helpers.R) was gated on `model_type == "single"`, so the
dynamic state block fell through to `.se_from_laplace_fit()`, which finds no
`H_beta` on a nested-Laplace fit and returns NA -- every dynamic nested fit reported
NA state SEs while p / gamma / epsilon were finite and correctly scaled (they are
ordinary tulpa_laplace blocks; only the state block goes through the nested engine).

The Louis identity is not single-season-specific: the state arm's complete-data
score is `x_i (z_i - psi_i)` in both families, so
`I_obs = X' diag(psi(1-psi) - w(1-w)) X` with `w = E[z_i | y]`. For a dynamic fit the
state arm is psi1 and its latent is z_1, whose smoothed posterior mean is the
season-1 weight column. (`em_result$weights` is [n_sites x n_seasons] for dynamic and
the helper length-checks against `nrow(X_occ)`, so it must be handed `w[, 1]`.)

## Final state (20 seeds, `_run_coverage_gate.R`)

| arm                  | coverage | mean_se/sd_est | field sd | truth |
|----------------------|----------|----------------|----------|-------|
| single, no field     | 0.95     | 1.13           | 0.1900   | 0     |
| single, field        | 0.95     | 1.19           | 0.9294   | 1.0   |
| dynamic, no field    | 0.95     | 1.08           | 0.1640   | 0     |
| dynamic, field       | 0.95     | 0.98           | 0.9174   | 1.0   |

Dynamic now matches single-season. svcTPGOcc, 12 seeds
(`_run_dyn_svc_seeds2.R`): slope 0.5095 (bias 0.0095, z 0.19), coverage 0.92, field
cor f0 0.942 (min 0.915) / f1 0.916 (min 0.814).

Progression on svcTPGOcc across the session: slope 1.2401 -> 0.9735 -> 0.5095;
z 4.84 -> 3.04 -> 0.19; field cor 0.518/0.461 -> 0.747/0.679 -> 0.942/0.916.

## Tests

- `test-dyn-occu-areal-recovery.R` (8) -- field shrinks / recovers, slope, SEs finite,
  coverage.
- `test-dyn-occu-svc-recovery.R` (7) -- both field surfaces (median AND min over
  seeds), slope, calibrated intervals.
- `test-dyn-occu-svc.R` (12) -- structural: binder slot, two distinct fields.
- `test-occu-areal-recovery.R` (13) -- single-season: interior-field surface, the
  grid mode off both boundaries, null-field shrinkage, slope + SEs.
- `test-int-occu-areal-recovery.R` (9) -- integrated: field shrinks / recovers,
  shared slope, SEs finite.
- Smoke 2898/0/0 before and after.

These replace reliance on `test-nested-laplace-families.R:41-64`, whose only
assertions were class, a type string, and `length(fit$spatial_field)` -- on a fixture
containing NO field. It passed for the entire time the path was broken.

## Summary

`dyn_occu()` + an areal field (`icar()` / SVC bar) under `method = "nested_laplace"`
does not estimate the field precision. It invents a field where none exists, roughly
doubles a field that does, inflates the season-1 slope as a mechanical consequence, and
returns no standard errors at all.

This is a **shipped, documented combination** ("Nested-Laplace (areal) | icar/bym2/car
(+temporal/iid) on occu/int_occu/dyn_occu"), not a new feature. It is not a regression
from any change in the current working tree.

## Evidence

12 seeds per arm, one thing varied at a time. Truth: slope = 0.5; "field sd" is the
realised sd of the fitted surface, whose truth is 0 in the no-field arms.
Runners: `_run_dyn_svc_ladder.R`, `_run_svc_nullfield.R`, `_run_dyn_icar_check.R`,
`_run_single_icar_scope.R`.

| arm | model  | data     | fit             | slope  | z     | SEs   | field sd (truth) |
|-----|--------|----------|-----------------|--------|-------|-------|------------------|
| A   | dyn    | no field | none, `laplace` | 0.5659 | +1.49 | 12/12 | --               |
| B   | dyn    | no field | SVC bar         | 0.8483 | +4.75 |  0/12 | --               |
| C   | dyn    | field    | SVC bar         | 0.9735 | +4.84 |  0/12 | --               |
| D   | single | no field | SVC bar         | 0.5692 | +2.06 | 12/12 | 0.17  (0)        |
| E   | single | field    | SVC bar         | 0.4896 | -0.19 | 12/12 | 0.90/0.85 (1.0/0.8) |
| F   | dyn    | no field | plain `icar()`  | 0.6409 | +2.96 |  0/12 | **1.35 (0)**     |
| G   | dyn    | field    | plain `icar()`  | 0.6545 | +2.60 |  0/12 | **2.46 (1.0)**   |
| H   | single | no field | plain `icar()`  | 0.5914 | +1.80 | 12/12 | 0.17  (0)        |
| I   | single | field    | plain `icar()`  | 0.6177 | +2.12 | 12/12 | 0.91  (1.0)      |

Broken: B, C, F, G -- every `dynamic` arm carrying a field.
Healthy: A (dynamic, no field), D/E/H/I (all single-season).

The `+0.07 .. +0.12` slope bias in the healthy arms is the fixture's baseline
small-sample Laplace bias (arm A, no field at all, carries +0.066). The defect is the
5x excess on top of it, not that number.

## Mechanism

The fitted field variance is grossly over-stated. The logistic conditional-slope factor
`sqrt(1 + 0.346 * sigma^2)` then predicts the observed slope inflation almost exactly:

- arm F: sigma 1.35 -> 0.5 * 1.28 = **0.64**; observed **0.6409**
- arm B: implied sigma ~2.3 -> 0.5 * 1.68 = **0.84**; observed **0.8483**

So the slope bias is a symptom. The disease is field-precision estimation on the
dynamic path. The prediction was registered before arms F/G were run, so it is not a
post-hoc fit to the numbers.

Both symptoms (field variance not shrinking; `se_finite = 0/12`) point the same way:
`.se_from_info()` (`R/laplace_helpers.R:348`) returns `NA` of length p when the
observed-information matrix is NULL or non-invertible, so the dynamic path is not
assembling the curvature that the hyperparameter grid and the SEs both depend on.

## Why single-season is unaffected

`.tobs_occu_reroute_to_joint()` (`R/em_nested_laplace.R:803`) returns FALSE unless
`model_type == "single"`. A single-season SVC bar is rerouted to the joint direct-grid
engine (the recovery-tested #81 path) and never touches the generic EM nested-Laplace
path. Arms D/E therefore did not test the generic path at all -- H/I were added for
exactly that reason, and single-season plain `icar()` is healthy regardless of route.

## Why it shipped

The only test of this combination is a shape smoke test
(`tests/testthat/test-nested-laplace-families.R:41-64`):

```r
expect_s3_class(fit, "tobs_fit")
expect_identical(fit$nested_laplace$multi_prior[[1]]$type, "icar")
expect_equal(length(fit$spatial_field), n_sites)
```

Class, a type string, and a length. Its simulation contains **no field**
(`zmat[, 1] <- rbinom(n_sites, 1, plogis(0.2 + 0.4 * elev))`), so it fits an ICAR to
field-free data -- precisely arm F, where the field comes back at sd 1.35 against a
truth of 0 -- and never inspects a value. `max.iter = 6L` caps it at smoke depth.

Plumbing tested; method never tested.

## Hypotheses tested and REFUTED (do not re-run these)

Six attempts, each eliminated by measurement. Recorded so the next session starts
past them.

1. **E-step field-blindness.** `.tobs_laplace_nested()` threads `latent_prior` into
   `build_single_callbacks` but not into the dynamic / integrated ones, so their
   E-step computes `psi1 = plogis(X beta)` with no field. Threading it in (the
   obvious fix, and what the in-code comment predicts) made the null-data field
   WORSE: sd 1.35 -> 1.92, slope z +2.96 -> +3.00. REVERTED.
2. **M = 1000 pseudo-binomial inflation.** Identical in both paths -- read both:
   `occ_block = list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ)` with
   M = 1000 for an areal (non-SPDE) latent prior. Not a difference.
3. **Response saturation (`y_occ = 0`).** Single-season at J = 8/12/20 encodes
   `y_occ = 0` exactly for every undetected site (w ~ 1e-4 .. 1e-9) and still
   returns sigma 0.2199, field_sd 0.1667, psi_x 0.5042 -- identical across J.
   Saturation is real and harmless. `_run_visit_saturation.R`.
4. **group_var / field length.** Both paths build 40 nodes on the 40-cell graph
   with identical `spatial_idx` (length 240, max 40). The field is the same object.
   `_run_field_length.R`.
5. **E-step weights / encoding.** Driven at the TRUE parameters on the SAME z1, the
   cell-mean pseudo-proportion SD -- the only thing the ICAR sees -- is 0.1812
   (single) vs 0.1766 (dynamic), ratio 0.97, both matching the truth's 0.1841.
   The engine is handed the same signal. `_run_yocc_dump.R`.
6. **Field-awareness as the stabiliser.** Turning single-season's field-aware
   E-step OFF changes nothing: sigma 0.2215 -> 0.2214, field_sd 0.1425 both ways.
   So field-awareness is NOT what keeps single healthy, and the comment in
   `build_single_callbacks` ("without this the EM converges to the fixed-effect-only
   fixed point and the field cannot track the data") does not hold for the areal
   case. `_run_single_blind.R`.

7. **EM convergence / the 25-iteration cap.** `occu_fit.R:406` caps the nested EM at
   `min(max.iter, 25L)` while `.tobs_fit_model` asks tol = 1e-4, damping = 0.7, so
   the fit reports `converged = FALSE`. Calling the driver directly at tol = 1e-4 /
   damping = 0.7 (dyn + icar, null field):

   | max_iter | converged | sigma  | field_sd | psi1_int | psi1_x |
   |----------|-----------|--------|----------|----------|--------|
   | 25       | FALSE     | 1.8257 | 1.5455   | -0.9723  | 0.5195 |
   | 100      | TRUE      | 1.8257 | 1.5473   | -0.9731  | 0.5196 |
   | 400      | TRUE      | 1.8257 | 1.5473   | -0.9731  | 0.5196 |

   Converged by 100 and sigma identical to 4dp at every budget. Not a convergence
   failure; the cap is cosmetic here. `_run_sigma_iters.R`.
8. **A degenerate / mis-built hyperparameter grid.** Both prior blocks are identical
   (type icar, spatial_idx len 240, n_spatial_units 40, same adjacency); the grid is
   built inside tulpa, not carried in the block. `_run_grid_dump.R`.

   RETRACTED over-read: sigma's reported SD is **0.0006**, not 0. An earlier dump
   printed `0.0000` under `%.4f` on a different fit and I read that formatting
   artefact as "sigma is pinned / not estimated". It is estimated -- the posterior
   is just ~100x sharper than single's (0.0006 vs 0.0616) and concentrated on the
   wrong value. sigma being identical to 4dp across iteration budgets is a
   CONSEQUENCE of that sharpness, not evidence of a frozen parameter.

9. **A mis-built / too-narrow hyperparameter grid.** Read the grid off the real fits
   (`_run_theta_grid.R`; it lives on `fit$nested_laplace$occ_fit$theta_grid`, NOT on
   `fit$nested_laplace` -- an earlier dump looked one level too shallow and wrongly
   concluded the grid was invisible from R). Both paths get the SAME grid:
   `theta_names = "b1.tau"`, 9 cells, tau in [0.3, 30], i.e. sigma in [0.1826, 1.8257].
   Only the WEIGHTS differ, and they land on opposite ends:

   | path    | peak cell | tau  | weight   | theta_mean | sigma  |
   |---------|-----------|------|----------|------------|--------|
   | single  | 9         | 30.0 | 0.603    | 23.75      | 0.2206 |
   | dynamic | 1         | 0.3  | 0.999998 | 0.30       | 1.8257 |

   Grid construction is fine; sigma is reachable. The marginal simply selects the
   wrong end.
10. **Bistability / damping picking the basin.** Sweeping damping 0.3 / 0.5 / 0.7
    through the real driver: single stays on cell 9 and dynamic on cell 1 at every
    value, identical to 4dp (`_run_damping_basin.R`). Not a basin problem.

    NOTE: a hand-rolled EM trace (`_run_em_trace.R`) showed single walking to
    tau = 0.3 under damping 0.7. That trace does NOT reproduce the real single path
    (the real one lands on cell 9 at every damping) and must not be built on. It is
    a broken replication of `.tobs_laplace_nested`, cause unidentified; it is kept
    only as a warning not to trust it.

## Where the evidence now points

Everything upstream of the engine is identical -- same block, same prior, same
field, same signal. The divergence is in the EM FIXED POINT.

Dynamic settles at `psi1_(Intercept) = -0.9887` (truth -0.619) with `sigma = 1.8257`.
Those two are not independent errors: they are a trade along the
(intercept, field-variance) ridge. Marginally,
`E[plogis(b0 + f)] ~ plogis(b0 / sqrt(1 + 0.346 sigma^2))`, and
`plogis(-0.9887 / sqrt(1 + 0.346 * 1.8257^2)) = plogis(-0.673) = 0.338` -- against a
true marginal occupancy of 0.35. So the dynamic fit reproduces the marginal
occupancy correctly and splits it WRONGLY between intercept and field: it buys a
large field by lowering the intercept, at no cost in marginal fit.

Single-season resolves that ridge (sigma 0.22, sd 0.0616). Dynamic does not, and
reports `sigma` sd 0.0000 -- the outer grid has collapsed onto one cell, i.e. it is
confident in the wrong sigma. Since the block handed to the engine is identical,
the difference must be in how the marginal likelihood over the hyperparameter grid
is evaluated for the dynamic state block.

That is the grid log-marginal -- exactly where the checklist said to start, and the
one thing not yet derived. Eight code-level guesses failed because the defect is not
visible from the R call sites; every input to the engine agrees between the two paths.

The sharpest remaining clue is the WIDTH of the sigma posterior, not its location:
single 0.2206 +- 0.0616, dynamic 1.8257 +- 0.0006 -- a posterior ~100x sharper on
the arm that is wrong. Both encode the state block at M = 1000, so both inflate the
apparent information by the same factor; only the dynamic one collapses. A marginal
that is sharply confident in the wrong sigma is the signature of a log-marginal
missing a term that varies with sigma (a normalising constant / log-determinant),
not of a mis-built grid: a wrong grid would misplace sigma without making the
posterior tight.

NEXT STEP (do this BEFORE any more code): derive the grid log-marginal for the
HMM-forward state marginal and compare it, term by term, against the single-season
one that works. Ten hypotheses have now failed, each costing a full fit-and-measure
cycle. Do NOT guess at an eleventh.

## Sharpest statement of the defect

Both fits are SELF-CONSISTENT FIXED POINTS of the same EM on the same grid:

- single : tau = 30  (grid ceiling), beta = (-0.93, 0.49), field sd 0.12
- dynamic: tau = 0.3 (grid floor),   beta = (-0.99, 0.52), field sd 1.47

Dynamic is not failing to converge -- it converges, stably, to a different and wrong
solution. Its weight on the wrong cell is 0.999998. The marginal over tau is
selecting the opposite end of an identical grid for a state block that, evaluated at
the truth, is near-identical to single's (cell-mean SD 0.1766 vs 0.1812).

`log pi(y|theta) ~ log pi(y|xhat,theta) + log pi(xhat|theta) - 0.5 log|H(theta)|`.
Only a term that varies with tau can flip which end wins. That is what has to be
derived and compared against the single-season path, which is the one thing this
session never did.

## A separate finding worth acting on regardless

**Single-season's sigma is a grid-ceiling artifact.** tau = 30 is the grid's MAXIMUM
and sigma = 0.1826 the smallest field those 9 cells can express; single reports
0.2206 with 60% of the weight on that last cell. Against a truth of sigma = 0 it
looks right because the grid floor sits near the truth, not because the marginal
identified anything. Every "single-season is healthy" statement in this session
rests on that. The scope bound still holds (single recovers slopes and SEs, and is
unaffected by this bug), but a null-field fixture cannot distinguish "correctly
shrunk" from "pinned at the grid ceiling" -- so a recovery test for the areal path
must use a field the grid can actually represent in its INTERIOR, and should assert
the peak grid cell is not on a boundary.

## Options

1. **Fix.** Find why the dynamic path does not deliver the observed information the grid
   and SEs need. Deep: needs the grid log-marginal for the HMM-forward state marginal
   worked out before touching code.
2. **Gate.** Reject `dyn_occu` + a structured field under `nested_laplace` in
   `.tobs_family_methods` with a pointer, matching the package rule that unsupported
   combinations error rather than silently downgrade. A silently biased fit with no SEs
   is the worst available outcome, so this is the floor.

Not decided here -- gating removes a documented feature and that is the maintainer's call.

## Follow-ups: both now run (2026-07-16)

### `int_occu()` + areal: did NOT inherit the fix -- three defects, all fixed

Recorded above as "inherits it (same `.tobs_laplace_nested` encoding)". That was
wrong, and was never measured. `int_occu` shared the disease and none of the cure:

1. `build_integrated_callbacks`'s `m_step_encode` had NO `latent_prior` branch --
   `if (is.null(spatial_occ)) M <- 1000L else M <- 4L`. A nested areal block arrives
   with `spatial_occ = NULL` (the prior is attached to `occ$prior` upstream), so the
   nested case took the **M = 1000** branch. The single-season fix keyed off
   `nested_block <- !is.null(latent_prior)`; integrated never got that branch.
2. `.tobs_laplace_nested()` called `build_integrated_callbacks(model, spatial = NULL)`
   -- no `latent_prior`. The field-aware E-step added to it was DEAD CODE.
3. `use_louis` was `model_type %in% c("single", "dynamic")`, so integrated fell through
   to `.se_from_laplace_fit()` -> NA state SEs, exactly the dynamic symptom.

Integrated shares one psi across sources, so `em_result$weights` is already the
per-site `w = E[z_i | y]` the Louis identity needs -- no column slice (unlike dynamic).

Measured after (8 seeds, `probe_int_occu_areal.R`, truth slope 0.5, sources J = 4/3):

| arm        | slope  | SEs | field sd (truth) | coverage |
|------------|--------|-----|------------------|----------|
| null field | 0.4888 | 8/8 | 0.1783 (0)       | 0.88     |
| field      | 0.5830 | 8/8 | 0.9758 (1.0)     | 0.88     |

Test: `test-int-occu-areal-recovery.R`.

### The grid-ceiling artifact: acted on, no code change

CONFIRMED by dumping the grid (`probe_grid_dump.R`). One axis, `b1.tau`, 9 cells:

```
tau  : 0.3000 0.5335 0.9487 1.6870 3.0000 5.3348 9.4868 16.8702 30.0000
sigma: 1.8257 1.3691 1.0267 0.7699 0.5774 0.4330 0.3247  0.2435  0.1826
```

- null field  -> peak cell **9** (tau 30, the ceiling), weight **0.603**, sigma 0.2206.
- interior field (marginal sd 1.0) -> peak cell **6** (tau 5.33), weight **0.309**,
  weight spread over cells 4-9, field cor 0.943.

So single-season DOES identify the precision when the truth is on the grid; the null-
field number was the grid floor. No code fix: for a truth of 0 the marginal piling on
the smallest expressible field is correct behaviour given the grid, and the slope/SE
scope bound holds. The fix is the TEST design, per this note's own prescription --
`test-occu-areal-recovery.R` uses an interior field and asserts the peak cell is off
both boundaries. It deliberately does NOT assert `peak == 9` on the null arm: that
would lock the grid floor in as expected behaviour and fail if the tau axis is ever
widened.

CORRECTION to this note's earlier arithmetic: `tau` is the ICAR CONDITIONAL precision
on increments, NOT `1/sd^2` of the realised surface -- a marginal field sd of 1.0 fits
at sigma ~= 0.40 (tau ~= 5.3) on this chain graph, which is why the interior peak is
cell 6 and not cell 3. Recovery is asserted on the field surface (sd, cor), never on
the reported `sigma`.

### `svcTPGOcc`

Now recovery-tested (`test-dyn-occu-svc-recovery.R`), not structural-only:
slope 0.5095, coverage 0.92, field cor 0.942 / 0.916 over 12 seeds.
