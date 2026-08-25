# Measurement notes

Measured evidence behind decisions recorded tersely in `CLAUDE.md`. Numbers here are
what was actually run: fixture, seeds, wall time, and what moved. `CLAUDE.md` keeps the
resulting rule; this file keeps the measurement it rests on.

## Community latent-structure engine (`R/community_latent.R`)

Issues #119 / #120 / #121 (shared engine), #153 -> #156 (loadings by marginal
likelihood), #157 (multi-start basin escape), #164 / #166 (per-family `factor.starts`).

the factor Newton (`.tobs_latent_factor_update`) DOES backtrack: local `ascend()` halves the step until the penalized objective improves and holds the previous iterate if it never does, and `nstep()` ridge-bumps a singular curvature. Non-finite guards are inline, NOT a named helper — `if (all(is.finite(Dz)))` / `(Dl)` skip a bad step, and a non-finite `working()` score/curv `break`s the pass. (The field solve `.tobs_latent_field_solve` has its OWN local `safe_solve()`, a ridge retry for a SINGULAR Hessian only; its Newton update is unconditional. Do not confuse the two — there is no `safe_step()` anywhere in the repo.) A latent-count marginal (nmix/distance) can return non-finite curvature far from the mode; the NaN used to surface only later as a non-finite `sd()` in the rescale. Found by a 6-SEED loop, never by a single fit. SHARED latent-structure engine for EVERY community family (#119/#120/#121): one block-coordinate ascent (community EM w/ the latent as an offset <-> field / factor updates) + the areal Newton, the factor update, bym2/car_proper/spde hyper grids. A family supplies ONE callback `working(eta) -> list(score, curv)` (per-(site,species) score+curvature wrt an additive offset on the structured arm): Poisson `(y-mu, mu)`, occupancy two-state, Bernoulli `(y-psi, psi(1-psi))`. Field solve is `t(A) diag(w) A + tau Q`, so the site->node map slot takes an areal group_var incidence OR an spde barycentric projector unchanged. Adding a family to every latent route = one callback, not a new fitter. **Loadings by MARGINAL likelihood, not the joint mode (#153 -> #156).** The factor update holds zeta at its joint mode, so `(zeta, lambda)` is a joint-likelihood estimate w/ `Ns*Q` incidental params growing w/ the sample -- Neyman-Scott, and the joint mode is INCONSISTENT. The site factors' estimation error lands in the fitted co-occurrence and lambda absorbs it, which is why #153's scalar could not close it: rescaling a direction finds the magnitude right FOR THE FITTED DIRECTION, hence too high for the true one. Over-fit grows w/ Q/S (6 seeds/cell, ms_count N=160 Poisson): S=7 Q=2 **1.435** (worst seed 2.229) / S=14 Q=3 1.076 / S=14 Q=2 1.064 / S=28 Q=2 1.057 / S=14 Q=1 1.014 -- monotone in Q/S but NOT proportional (flat ~1.06 across the middle of the range; most of the S=7 excess is 2 seeds of 6). Fix = `.tobs_latent_factor_mmle()`: EM on the SAME joint site marginal, ascended over all S*Q loadings. E-step = posterior `p(z_i|y_i)` on the grid the marginal already builds (`.tobs_latent_joint_grid()`, split out of `.tobs_latent_joint_marginal()` so both readers share ONE grid); M-step = per-species Qk-dim weighted Newton w/ the nodes as design rows, backtracked on the expected complete-data ll (Dempster/Laird/Rubin 1977). Gradient FD-validated (rel 6e-4, cor 0.9999998). 30 paired fits across the Q/S cells: pooled |mag-1| **0.163 -> 0.065**, worst-cell max 2.23 -> 1.10, residual cor UP in every cell. On the issue's own 16-seed fixture (seeds 201-216): ms_count_factor mag median 1.060 -> **0.991** / mean 1.092 -> 1.029, community intercept z **-3.2 -> +0.82** (bias gone), slope z 1.2 -> 1.56, rescor median 0.991; ms_occu_factor slope z **2.6 -> 0.71** (the #156 symptom), intercept -0.014 -> -0.017 (z -1.2 -> -2.23, same size, tighter sd), mag median 1.021 max 1.212. ONE regression survives: ms_count seed 215 goes 1.532 -> 1.654, a bad DIRECTION basin (marginal 31 nats below reachable, `residual_cor` still reads 0.90) -- **#157**, and rescaling the start does NOT escape it (measured: moves it 0.001 at 3x cost).

**ONE estimator, ONE state (#156).** `.tobs_latent_factor_update()` + `.tobs_latent_factor_scale()` now run ONCE, on outer pass 1, purely to INITIALIZE -- the marginal's lambda-gradient vanishes at lambda=0 so the EM cannot start from the zero init, and the 1-D bracket is a GLOBAL magnitude search the local EM cannot do. Running the joint-mode update every pass alongside the MMLE diverges BOTH ways: write the refined lambda back over the update's state and its Newton regrows the magnitude while the refinement shrinks it (5.3e3 x truth, residual cor 0.01); keep the two separate and the EM conditions on an offset whose spread is attenuated, the update grows lambda to cover the shortfall, and each pass compounds (measured S=7: |lambda| 1.98 -> 2.94 -> 4.29 -> ... -> 925 -- and because the community EM absorbs the inflated offset into the coefficients the runaway pair is LOCALLY SELF-CONSISTENT, so nothing downstream rejects it).

**Offset by SCORE-MATCHING, not `zeta t(lambda)` (#156).** The driver hands the coefficient update ONE point offset, and no plug-in reproduces the integrated objective through a nonlinear link: posterior means carry too little latent variance (log-link Jensen -> community intercept +0.165), scores rescaled to unit variance carry too much (magnitude back to 1.70x w/ a 4.92x tail). `.tobs_latent_factor_offset()` solves `score(eta+off) = E_z[score(eta+lambda_s'z)]` per cell (scalar Newton, `d score/d eta = -curv`), which makes the plug-in + integrated STATIONARY CONDITIONS identical for ANY family w/ no knowledge of the link; reduces to `lambda'zhat + v/2` on a Poisson log link. Validated w/ the community EM removed (true loadings, per-species GLM, 10 seeds): int -0.0098 / slope -0.0006, vs +0.3143 for a zero offset. `fit$model[[offset_slot]]` reads THIS, not `zeta t(lambda)`, so fitted()/WAIC see the predictor the coefficients were fit against.

**Block-coordinate callers MUST warm-start `init_b`/`init_Sigma` (#156).** `.tobs_community_em()` has accepted them since #119 but only `ms_abun_latent.R` passed them -- every other latent caller cold-restarted all per-species deviations AND the community covariances on every outer pass. Wiring them into ms_count / ms_occu / ms_dyn_occu cut the factor fit **2.5x** (96s -> 39s per 6 seeds) with the answer unchanged to 4 decimals.

**Factor path needs `max.outer` 150; field path 25 (#156).** A field block reaches `tol` and breaks early. The factor block does NOT: it alternates w/ the coefficient block along a SLOW mode (the per-species intercepts and the offset's per-species level absorb the same latent level), the per-cell offset change decaying ~2%/pass, so `tol` 1e-4 on the STEP would need ~300 passes -- and the step is ~1/50 of the REMAINING error, so the criterion reads converged long before it is. Stopping at 25 leaves a real community-mean bias: int **+0.0613** (25) / +0.0237 (60) / **+0.0001** (150) / -0.0038 (400). Driver resolves `max.outer = NULL` -> `factor.outer` when `has_factor` else 25; an explicit `control$max.outer` still wins. Callers pass the raw `control[["max.outer"]]`, NOT `%||% 25L`. **`factor.outer` is per-family and each family sets it from its OWN measurement** -- ms_count/jsdm 150 (the curve above), ms_occu 150 (its 16-seed recovery was measured AT that budget), everything else stays 25 until measured. Do NOT globalize it: the cost is not transferable either. ms_count nets 13.8s -> 22.5s (1.6x) because the 2.5x warm start pays back part of the 3.5x longer loop, but `ms_abun_latent` already warm-started, so 150 costs it the full 6x and pushed `test-ms-abun-factor.R` past **85 min on one file** before it was reverted to 25. `max.outer = 60` is the measured middle for ms_count (11.2s/fit, FASTER than shipped, intercept +0.024 = the size the shipped fit already carried) if suite time ever matters more than the last of the bias. The MMLE's inner `em.iter` (default 10) is NOT a speed lever: on ms_abun cutting it 10 -> 2 saved ~20% of a ~1000s contended fit and pushed the loading magnitude 1.93 -> 2.03 (under-shrinks) with `res_cor` dropping 0.994 -> 0.993 -- the per-outer-pass community EM refit over latent N is the cost, not the loading EM, so leave em.iter at 10. **That last conclusion predates the #157 multi-start and no longer holds** -- see `factor.starts` below.

**`factor.starts` (multi-start width) dominates a latent-N fit, NOT `max.outer`.** #157's basin escape runs K candidate starting directions on the FIRST factor pass -- cosine + principal-factor init + `.tobs_latent_factor_random_starts(k = 6L)` = 8 -- and each runs a FULL loading-EM to convergence, because the raw scale-search values were measured not to rank the same as the converged ones. Its cost was measured on `ms_count` ("+30s against an ~90s fit"), whose oracle evaluates a closed-form density. `ms_abun`'s oracle sums over the latent N per species-site, so the same 8 candidates cost far more there. Measured on `test-ms-abun-factor.R`'s own fixture (N=80, S=8, Q=2, seed 4, tulpaObs 0.0.179 / tulpa 0.0.101, idle box): `max.outer` 1 -> **19.7 min**, 3 -> 21.5, 25 (default) -> 30.1, with `resid_cor` 0.983/0.992/0.994. So a later outer pass costs **0.4-0.9 min** while the FIRST pass -- the only one that runs the multi-start -- costs **19.7 min**. Pass 1 also runs the warm start (`nmix_laplace_re(max_iter = 100L)`, a full community EM), so the two were separated directly on that same fixture (K_max=130): **warm start 2.5s**, full fit at `factor.starts` 1 -> **442.2s**, at 8 -> **892.3s**. The warm start is not the cost; the 7 extra candidates are **450s, half the fit**. What they buy on that seed is nothing measurable -- `mag_ratio` 0.9539 (1 start) vs 0.9535 (8), `res_cor` 0.9942 both. Successive halving was considered and rejected on these numbers: it makes a component cheaper that should instead be smaller, and it prunes on truncated candidate values that #157 measured as mis-ranking against their converged ones. NB the 30.1-min figure above and the clean 14.87-min fit are the SAME fixture and settings; the first sweep shared the box with another session's R job, so treat the `max.outer` row costs as contended and the ratios (not the absolutes) as the usable part. The share is also fixture-specific -- at N=40/S=6/seed 300 the same 7 candidates were only 109s (16%) -- because candidate cost scales with `K_max = max(y)+100` (the Royle marginal is O(K) per site). The nested_laplace spatial-factor fit is 34.3 min at `max.outer = 1` alone. At defaults the file is ~4.5-5h (block3 ~30 min + block5 6x30 + block4 ~70), which is the 6.6h abort once contended. Driver takes `factor.starts` (default 8L, byte-identical to the pre-knob path: cosine, eigen attempted iff >=2, then `factor.starts - 2` random), threaded family -> `control$factor.starts` via the `block_coordinate` group. Set it per family from that family's OWN measurement, exactly as `factor.outer` is. **`ms_abun` is now measured and set to 1** (`R/ms_abun_latent.R`, family default; `control$factor.starts` still overrides): 16 seeds at N=80/S=8/Q=2 at ONE start produced no magnitude failure at all (widest 1.25x, none reaching the 1.3x flag), and six seeds re-fit at 8 -- the three worst by `res_cor` plus three healthy controls -- did not move (largest `res_cor` change **0.0027**, on a CONTROL, moving DOWN; largest `mag_ratio` change 0.0145) while costing a consistent **2.0-2.3x**. So on this family the 7 extra candidates buy nothing measurable. Do NOT copy the value: both cost and benefit scale with how expensive one oracle evaluation is. **NEITHER summary screens a fit alone.** `residual_cor` is row-normalized and blind to a magnitude regression (the ms_count seed at 1.53x truth still read 0.93), so score `sqrt(tr(Sigma_res))` against truth -- but `mag_ratio` is rotation-invariant and equally blind to a DIRECTION regression, which is how the 16-seed screen's `mag_ratio > 1.3` criterion reported "0 flagged" while seed 314 sat at `res_cor` **0.0769** with a healthy 0.956 magnitude. Screen on both. Seed 314 is NOT a basin the multi-start escapes: it reads 0.0769 at one start and 0.0772 at eight, so it is a hard fixture, and a starts-only remedy was never going to reach it.

**`ms_occu` is now measured and set to 1** (`R/ms_occu_field.R`; gcol33/tulpaObs#164). A 10-seed screen (smaller than `ms_abun`'s 16+6 -- a session budget, treat as indicative) at N=250/S=16/Q=2/J=5 (median 16.1s/fit) had 1/10 flagged by `mag_ratio`; reseeding the worst 4 + 4 healthy controls at 8 starts cost 3.5-4.3x and on 4 of the 8 reseeded seeds landed on the EXACT SAME loadings (`d_mag = d_res = 0.0000`) -- the eight starting directions converging to one basin rather than escaping it -- and the flagged seed reproduced identically at 8 starts too (a hard fixture, not a basin escape, same character as `ms_abun`'s seed 314). Confirmed against the full `test-ms-occu-factor.R` suite (NOT_CRAN, every existing recovery/coverage assertion) at the new default: 0 failures.

**`ms_count`/`jsdm` are NOT set from a random screen -- #157's own adversarial seed caught it.** A 10-seed screen at N=100/S=10/Q=2 (median 6.1s/fit) found no seed where 8 starts reliably beat 1 (2/10 flagged by `mag_ratio`, reseeding cost 1.35-3.25x and moved `res_cor_err` with no consistent direction), so the family default was set to 1 the same way as `ms_occu`. Running the full `test-ms-count-factor.R` suite at that default caught what the random screen missed: `test-ms-count-factor.R:163` ("escapes the #157 direction basin on seed 215") is a COMMITTED regression test for the exact adversarial case #157's multi-start was built for, and it failed at `factor.starts = 1` -- `mag_ratio` regressed to 1.65, the documented PRE-#157-fix value, against the test's own `< 1.40` gate. A 10-random-seed screen has no way to reproduce a specific known-hard seed unless it happens to draw it; `ms_count` genuinely needs (at least some of) the multi-start on seeds like 215, unlike `ms_abun`'s 16-seed screen which found no such case at all. **Reverted to the driver default (8, unset in `R/ms_count_spatial.R`)** rather than shipping a plausible-looking but falsified number. Lesson for the next family measured this way: run the family's OWN existing recovery suite (not just a random screen) at the candidate default BEFORE setting it -- a documented adversarial seed is exactly the case a random screen is least likely to catch. `ms_count`/`jsdm` folded into gcol33/tulpaObs#166 for a redo that checks known seeds first.

**`ms_distance` measured (#166, post-#165) -- kept at the driver default of 8, same outcome as `ms_count`/`jsdm`.** Its single-fit cost (15-35 min even on a small fixture) made the original in-session screen intractable; landing gcol33/tulpaObs#165 (per-fit quadrature caching) first did NOT fix that -- the quad rebuild #165 removed was measured at <0.3% of a single kernel call, so it could never have explained the multi-minute fit times (see #165's closing comment / #167 for the real O(n_sites * K_max) cost). A lower-abundance fixture (`mu_lambda = log(8)` vs the default `log(30)`) was tried to shrink the screen and did NOT reliably help either -- one seed came back at 1781s, slower than the 887s default-abundance reference, because `K_max = 3*R_max+100` is driven by order statistics across all site x species draws, not the mean. Ran the full protocol anyway as a detached multi-hour Scheduled Task (checkpointed per-fit, survived past session boundaries) on the same N=40/S=6/Q=2 fixture: 10-seed screen at `factor.starts=1` found 2/10 flagged by `mag_ratio` outside [1/1.3, 1.3] (seed 5 at 1.3213, seed 8 at 1.3931) -- the SAME 2/10 flag rate `ms_count`/`jsdm`'s screen found. Reseeding the worst 4 by `res_cor_err` (seeds 4,1,5,2) plus 4 healthy controls (9,7,10,6) at `factor.starts=8`: seed 5's `mag_ratio` moved from 1.3213 to 1.1610 (crossing back under the flag threshold -- a genuine rescue, the largest movement in the reseed set) while every healthy control stayed flat. Seed 8 -- also flagged, in fact the MOST deviant seed in the screen -- was never reseeded, because the selection criterion sorted by `res_cor_err` alone and seed 8's `res_cor_err` (0.1193) wasn't in the bottom 4; a real gap in this run (should have also selected worst-by-`mag_ratio`-deviation, not just worst-by-`res_cor_err`, per the "screen on both" rule above), left open rather than closed by assumption. Cost ratio (`8starts`/`1start`) came back a median 1.19x post-#165, far short of the 2.35x measured pre-#165 on this same fixture -- and not because of the caching fix (see #165), just per-seed `K_max` variance dominating. Given a real rescue case in-hand and an unverified equally-flagged seed, **kept at the driver default of 8, unset in `R/ms_distance.R`** (already its as-shipped behavior -- `factor.starts` there was never overridden), rather than shipping 1 on the strength of a screen with a known gap. Full log + per-fit RDS checkpoints (seed x starts) in `dev_notes/factor_starts_166/` (gitignored, local only).

**`n.quad` is NOT threaded from any caller** -- the driver's default 5 is what every community latent fit actually uses, so `control$n.quad` silently does nothing on this path. Not a defect for the magnitude (ARGMAX stable to <0.4% vs n=21, `test-community-latent-quad.R`) NOR for the offset (5 nodes == 21 to 8e-4), but do not read a passed `n.quad` as having taken effect.

TRAP that survives all of the above: the magnitude must come from the JOINT marginal. Species s' OWN marginal reduces exactly to 1-D in `||lambda_s||` but does NOT identify it -- one Bernoulli/site is no replication, and a normal-mixed logit ~ a rescaled logit `plogis(eta/sqrt(1+0.346 sigma^2))`, so sigma and the coef scale are confounded along a ridge (fitted that way the scale slides to the grid floor, 0.22x truth). Cross-species co-occurrence, visible only in the joint integral, is what pins it. `residual_cor` is row-normalized and so CANNOT see a magnitude regression -- the seed carrying 1.53x truth still reported 0.93; assert on `sqrt(tr(Sigma_res))` (rotation-invariant), which is what `test-ms-count-factor.R` / `test-ms-occu-factor.R` now do

## Latent-N truncation guard (`abun()` / `ms_abun()`, `R/nmix_site_marginal.R`)

Contract + rules are in `CLAUDE.md` under "N-mixture abundance". Measurements:

- **Why the shared ceiling hurts.** One species-site drawing `y = 2248` (mean count
  8.15) took `K_max` to 2348 and made the other 639 sites evaluate 2344 states each
  instead of ~100. Cost is linear in the state count, so the fit slowed by that ratio
  and a single fixture seed ran 4x its siblings.
- **Boundary mass does not detect it.** A fixture passed the boundary-mass warning at
  4.6e-06 while sitting 0.032 nats below the uncapped optimum and 0.57 away in the
  coefficients. The score-gap guard separates the same two states by eight orders of
  magnitude: 1e-10 at the uncapped optimum, 7.4e-02 at the capped one.
- **Guard is live, not just wired.** Measured on a fitted `ms_abun() + latent()`
  (`dev_notes/_probe_guard_live.R`): gap exactly 0 at the shipped headroom, 93.5 with
  the window collapsed to the single state `max(y_i)`, and 137.0 collapsed WITH the
  fit's own factor offset. Finite on every call, six orders above `.NMIX_SCORE_TOL`.
- **End to end** (`dev_notes/_probe_kmax_fit.R`, `ms_abun() + latent(2)`): the
  pathological seed (`max(y) = 2248`, N=80 S=8) went from running past 108 min
  UNFINISHED to completing in 52.2 min, `res_cor` 0.955. A paired small fixture
  (N=40 S=6, ceiling 259) ran 1496.1s uncapped vs 788.6s capped, **1.90x**, with
  `mag_ratio` and `res_cor` agreeing to all six printed digits.
- **Buys nothing on light tails.** `test-ms-abun-factor.R`'s own seed 4 has
  `max(y) = 30`, so the shared ceiling was already nearly per-site (9832 -> 8080
  states, 1.2x by construction) and the guard spends part of that back.
## Continuous NNGP `svc()` recovery (#119, tulpa 0.0.82)

`test-occu-svc-nngp-recovery.R`, seeds 1/2/3/11 at N=150, J=6, p=0.6, truth
phi=0.25 / sigma=1.3:

- divergences 72-83% -> **0 on every seed**
- phi ~4 -> 0.14-0.23; sigma 1.06-2.31
- surface cor 0.66 -> 0.73 mean (did NOT move materially; information-bounded at
  these settings, not sampler-bounded)

Laplace backends (#143) on the same truth at N=150, J=6, p=0.6: surface cor
0.78 / 0.60 / 0.83 on seeds 1/2/3, matching the NUTS path's 0.76 / 0.60 / 0.81.
Measured surface cor for the observation families (#144): `removal()` ~0.97,
`fp_occu()` ~0.78 at N=120.
## Adaptive copy-axis boundary fixture (#194, tulpa 0.0.163)

`test-cover-hurdle-adaptive-grid.R`, `simulate_d3_like()` at N=300, n_s=25, BYM2 chain
graph, beta positive arm, `alpha_true = 1.5` pinned at the top node of
`control$alpha.grid = c(0, 0.5, 1.0, 1.5)`, `sigma.grid = c(0.3, 0.6, 0.9)`,
`rho.grid = c(0.5, 0.7, 0.9)`, seeds 3401-3420. Engine: tulpa 0.0.163 (commit
0820e93) in a private library; tulpaObs 0.0.190.

Four arms, 20 seeds each, 80 fits in 497 s total (5.6-6.6 s per fit):

| arm | alpha coverage | mean ci_hi | sd ci_hi | max alpha node | triggered axes |
|---|---|---|---|---|---|
| adaptive, default 7-node phi | 1.00 | 2.2227 | 0.3034 | 3.375 (20/20) | `alpha,phi_pos` (20/20) |
| fixed, default 7-node phi    | 1.00 | 1.7328 | 0.0120 | 1.5 (20/20)   | none (20/20) |
| adaptive, 13-node phi        | 1.00 | 2.2698 | 0.3734 | 3.375 (20/20) | `alpha` (20/20) |
| fixed, 13-node phi           | 1.00 | 1.7333 | 0.0107 | 1.5 (20/20)   | none (20/20) |

- **The coverage gap does not open.** All four arms cover 20/20, so the fixed-vs-adaptive
  coverage comparison the file used to assert is not evidence about the adaptive path at
  this placement. With the truth exactly ON the top node the upper side of the fixed
  arm's quantile interval cannot miss: since gcol33/tulpa#353 the outer cell contributes
  its own half-width to the reported support, and before that change the quantile clamped
  AT 1.5 so `truth <= ci_hi` held as an equality. The 0.53 -> 0.83 gain recorded when the
  adaptive path landed (commit d054bf3) was measured on a Wald summary, which the file
  replaced with the quantile CI.
- **What separates the arms is the upper edge.** Fixed: full range 1.6832-1.7375 over 20
  seeds, i.e. the axis's own outer-cell geometry. Adaptive: 1.8250-3.2263, following the
  data. Paired, the adaptive upper edge is larger on 20/20 seeds, minimum ratio 1.0506.
- **Truncation certificate.** The fixed arm puts mean 0.876 (13-node phi; min 0.206) /
  0.888 (7-node phi; min 0.029) of its posterior weight on the top node 1.5.
- The 13-node phi pin makes the adaptive pass fire on the copy axis ALONE; on the 7-node
  default it also densifies `phi_pos`.

At `alpha_true = 0` (seed 3101, same pins) the boundary trigger does not fire on the copy
axis (`triggered_axes` = `phi_pos`), CI `[0, 0.449]`, median 0. The var-of-means
consistency pass -- a separate mechanism, reported separately -- does place copy-axis
cells there (down to 2e-10 and up to 8e+04), which is why the no-op assertion reads the
triggered-axes list rather than the node set.

Prune ladder on the same fixture (seed 101, N=400, `alpha_true = 1.0`, 5x4x4x7 grid, 560
cells): `prune.tol` 1e-4 / 1e-3 / 1e-2 / 3e-2 / 1e-1 / 2e-1 prunes 539 / 543 / 548 / 553 /
556 / 558 cells with no safety fallback, and moves the alpha posterior mean by 8e-05 /
3e-04 / 3e-05 / 1e-02 / 5e-02 / 1.5e-01. At 0.5 the engine's cheap-screen argmax check
fires and the fit falls back to the full grid, returning a bit-identical result -- so
`prune.tol = 0.5` is not a usable perturbation and `0.2` is.

## Cover-hurdle fixtures re-anchored off the retired copy axis (#196-#200)

Engine: tulpa 0.0.163 (commit 0820e93) built into a private library from a clean
`git archive`; tulpaObs 0.0.191 installed against it. Probes in
`dev_notes/issue196/` (gitignored). Every number below is a run, not a carry-over.

### What the #196 fit actually integrated

`test-sla-cover-joint.R`'s vanishing-sigma test, as shipped
(`sigma.grid = c(0.01, 0.02, 0.03)`, `rho.grid = 0.5`, no copy pin, adaptive grid
at its default TRUE): **138 cells**, axes `sigma, rho, alpha, phi_pos`, copy axis
`0 0.1 0.234 0.376 0.548 1.282 1.961 3 3.622 7.021 13.22 16.43` (12 nodes after
refinement) and `phi_pos` 13 nodes over [2, 300]. So the comment's "3-point
(sigma, sigma_pos) grid at scale ~0.02" named neither the node count nor the
amplitude: the largest cover-arm amplitude `alpha * sigma` any cell reached was
`16.43 * 0.03 = 0.49`.

### #196, vanishing-amplitude collapse (10 seeds, 105-114)

Final pins `sigma.grid = c(0.02, 0.08, 0.15)`, `rho.grid = 0.5`,
`alpha.grid = c(0, 0.05, 0.20)`, `adaptive.grid = FALSE` -> exactly 63 cells,
max amplitude 0.03. Truth `sigma = 0.05`, `alpha = 0.01`, both off-node.

| statistic | vanishing (sigma .05, alpha .01) | control (sigma .5, alpha 1) |
|---|---|---|
| max abs skew_pos, joint | 0.0757 | 0.0409 |
| max abs skew_pos, separate | 0.0633 | 0.0957 |
| max abs skew_occ, joint | 8.63e-04 | -- |
| max abs(beta_pos joint - sep) / se_sep | **0.2024** | min **0.6594** |
| max abs(beta_occ joint - sep) / se_sep | 0.0822 | min 0.0101 |
| max abs(se_pos ratio - 1) | 0.1878 | min 0.1116 |

Reads: the skewness DIFFERENCE does not separate the two regimes. On the
neighbouring pin `sigma.grid = c(0.02, 0.05, 0.1)` (same alpha pin, same seeds,
adaptive off) the vanishing arm's worst per-coefficient difference is 0.1281 and
the non-vanishing control's is 0.1273 -- so the retired `tolerance = 1e-1` on it
was not evidence of the collapse, and on that grid it fails at the very seed the
test runs (0.1281). The cover-arm posterior gap DOES separate, on all 20
paired fits, and the occurrence-arm gap does not (control minimum 0.0101). Hence
the shipped assertion is the cover arm alone at 0.30 se.

Large-N band, `alpha.grid = c(0.5, 1.0, 1.5)`, `adaptive.grid = FALSE`, seeds
104-108: at N = 1000 max abs skew_occ 0.00002-0.00200, max abs skew_pos
0.00199-0.01292 (seed 104: 0.00010 / 0.00444). At N = 120: skew_occ
0.00114-0.00760, skew_pos 0.01756-0.12827 (seed 104: 0.12827).

### #197 / #198, on-node vs off-node pinned axes

Paired arms, same seeds, same fixtures, only the pins moved.

| fixture | centred pins | off-node pins |
|---|---|---|
| lognormal `sigma_pos`, 10 seeds: mean rel err / max rel | 0.0384 / 0.1400 | 0.0346 / 0.1392 |
| lognormal slopes, 15 seeds: abs bias occ / pos / pooled 95% coverage | 0.0062 / 0.0158 / 0.933 | 0.0081 / 0.0162 / 0.933 |
| beta `phi_pos`, 10 seeds: mean rel err / max rel | 0.0298 / 0.2095 | 0.0246 / 0.2043 |
| `beta_pos_0` coverage, 20 seeds | 0.95 (19/20), mean 0.4096 | 0.95 (19/20), mean 0.4098; with the copy axis also pinned off-node, 0.95, mean 0.4037 |

So centring the pinned axes on the truth was worth less than a percentage point
on every one of these, and nothing here rests on it. Recorded because a reader
could not previously tell, not because it moved a number.

Break perturbations (each rewritten band, shown failing):

- `sigma_pos` bands 0.08 / 0.20: simulate the same 10 seeds at 0.55, score against
  0.4 -> 0.2942 / 0.3594. (Against its own truth 0.55: 0.0588 / 0.1060.)
- slope bands 0.05 / 0.05: simulate at slopes (0.9, 0.45), score against
  (0.7, 0.3) -> 0.1900 / 0.1144. Note the fixture's coverage arm reads
  `sim$truth`, so moving the simulated slopes moves the interval's target with
  them and does NOT exercise the floor -- scored against the unshifted (0.7,
  0.3) the same 15 seeds cover 0.367.
- the 0.85 coverage floor: narrow the interval the seeds are scored in, 1.96 se
  -> 0.4 se.
- `phi_pos` bands 0.06 / 0.25: simulate at phi = 12, score against 30 -> 0.6165 /
  0.6537. (Against its own truth 12: 0.0412 / 0.1888.)
- `phi_pos_sd < 8`: 2.1975 at n_pos 316 (the fixture), 2.6580 at 151, 5.0595 at
  83, 11.6546 at 42.
- `beta_pos_0` coverage floor 0.80: 0.95 as shipped; 0.55 with the demean stripped
  out of the simulator (0.60 under the old centred pins). The header's 0.43-0.47
  is INLAabun's d3 sweep, a different fixture -- the number measured here is 0.55.

### #199, multi-block bands that could not fail

Grid: `sigma c(0.2, 0.45, 0.9)`, `rho c(0.5, 0.85)`, `alpha c(0.4, 1.0, 2.5)`,
`tau c(1, 4, 16, 64)`, `rho_temporal 0.6`, `sigma_re c(0.06, 0.2, 0.7)`,
`phi c(6, 15, 38, 95)`, `adaptive.grid = FALSE` -> 864 cells, ~3 s per fit.
Truth A is the fixture's own; truth B raises alpha and both non-spatial SDs,
lowers the field SD, and makes the beta arm far more disperse.

Seed 7001: A `sigma .8933 alpha 1.0013 tau 51.50 re .2395 phi 38.00`;
B `sigma .4480 alpha 1.5717 tau 8.771 re .7000 phi 6.009`. All five orderings
also hold on seeds 7011-7014.

Single-component reversals (B with one parameter put back to A's value), each
failing exactly its own rule:

| reversal | sigma | alpha | tau | re | phi |
|---|---|---|---|---|---|
| B full | pass | pass | pass | pass | pass |
| sigma_year -> 0.3 | **FAIL** | pass | **FAIL** | pass | pass |
| sigma_obs -> 0.25 | pass | pass | pass | **FAIL** | pass |
| phi_b -> 30 | pass | pass | pass | pass | **FAIL** |
| sigma -> 0.6 | **FAIL** | pass | pass | pass | pass |
| alpha -> 1.2 | pass | **FAIL** | pass | pass | pass |

The `sigma_year` reversal takes down the spatial rule as well: the year effect and
the areal field are partly confounded at N = 400 with 6 seasons, so raising the
year SD absorbs field amplitude.

Per-seed variability of the block moments at this fixture size is large -- over
seeds 7001/7011-7014 truth A gives `sigma` 0.386-0.894, `alpha` 1.001-2.399,
`tau` 8.18-57.68 -- which is why the shipped assertions are orderings against a
paired truth rather than single-seed recovery bands.

## Sampled field hypers on the occu_cover NUTS + areal path (#204)

`dev_notes/probe_204_calib.R`, 12 seeds (7001-7012) at side 8 (64 cells), J = 5,
lognormal cover, truth `sigma = 0.7` / `alpha = 1.0`, 2 chains x 1000 kept draws
after 800 warmup, run at the then-default `adapt.delta = 0.9`.

Coverage of the 95% credible interval, 12/12 seeds on every line:

| field | alpha cov | mean alpha | field_sd cov | mean field_sd | div total (max) | hyper Rhat max | min ESS | field cor | s / fit |
|---|---|---|---|---|---|---|---|---|---|
| icar | 1.00 | 1.053 | 1.00 | 0.943 | 0 (0) | 1.023 | 153 | 0.797 | 24.4 |
| car_proper | 1.00 | 1.030 | 1.00 | 1.012 | 68 (13) | 1.015 | 178 | 0.786 | 30.3 |
| bym2 | 1.00 | 1.036 | 1.00 | 0.871 | 3 (2) | 1.030 | 155 | 0.764 | 28.1 |

`alpha` is essentially unbiased (0.03-0.05 high). `field_sd` runs 0.17-0.31 high:
its posterior is right-skewed at 64 binary occupancy sites, so the posterior MEAN
sits above the bulk -- `fit$nuts$hyper_median` is the summary to quote against a
truth. The warm nested-Laplace point estimates over the same seeds
(`sigma` 0.918 / 1.010 / 0.928, `alpha` 1.020 / 1.002 / 0.978) sit inside the
sampled intervals, which is the like-for-like check; the shipped gate is coverage
of the TRUTH, since agreeing with the deterministic backend is the circularity
#204 exists to remove.

### Re-measured over 40 seeds after the Cholesky field draw (#279)

`dev_notes/probe_ocsn_hyper.R`, the same fixture and controls, seeds 7001-7040 on
both field types, at the shipped `adapt.delta = 0.99` for this path. The
simulator's ICAR draw changed in 8c67903, so these are different realisations of
the same truth at the same seeds.

| field | n | mean alpha | \|err\| | mean field_sd | \|err\| | median field_sd | alpha cov | field_sd cov | div | s / fit |
|---|---|---|---|---|---|---|---|---|---|---|
| icar | 12 | 0.974 | 0.026 | 1.067 | 0.367 | 1.009 | 1.00 | 1.00 | 0 | 47 |
| icar | 20 | 1.062 | 0.062 | 0.940 | 0.240 | 0.870 | 1.00 | 1.00 | 0 | 46 |
| icar | 40 | 1.095 | 0.095 | 0.933 | 0.233 | 0.853 | 0.97 | 0.97 | 0 | 43 |
| car_proper | 12 | 0.964 | 0.036 | 1.166 | **0.466** | 1.108 | 1.00 | 1.00 | 0 | 63 |
| car_proper | 20 | 1.053 | 0.053 | 1.042 | 0.342 | 0.967 | 1.00 | 1.00 | 0 | 62 |
| car_proper | 40 | 1.076 | 0.076 | 1.011 | 0.311 | 0.935 | 0.97 | 0.97 | 0 | 60 |

The bold cell is the CI failure: `expect_lt(abs(mean(fsd_mean[ok]) - sig_true),
0.45)` at 12 seeds reads 0.466 here and 0.470 on the Linux runner. It reproduces
locally and is deterministic, so it is a property of a code state, not a flake.

**The band was not wrong; 12 seeds was too few for it.** `field_sd` is weakly
identified at 64 binary occupancy sites: a single fit's posterior mean ranges
0.45 to 1.87 with a per-seed SD of 0.344 (icar) / 0.362 (car_proper), so the mean
over 12 seeds carries a standard error of 0.099 / 0.105 and over 20 seeds 0.077 /
0.081. The old gate sat about 1.3 SE above the car_proper bias, which the fixed
seed set can cross on its own: seeds 7008-7013 are a run of high draws (1.51,
1.30, 1.44, 1.50, 1.58, 1.12 against a 40-seed mean of 1.011). Bootstrapping the
40 measured values, a random 12-seed set clears the 0.45 band 90% of the time and
a 20-seed set 95%; on the actual seed prefix the value is 0.342 at 20 and 0.298
at 30. The fix raises `n_seeds` to 20 and leaves both bands where the earlier
measurement put them.

Robustifying the summary instead does NOT work, and was measured rather than
assumed: an across-seed median of the per-fit medians is *less* stable at this
size (bootstrap 0.80 at n = 12 for car_proper against 0.90 for the mean form),
because the committed seed set is extreme on every summary of it.

Two things changed between the table above and this one -- the field realisations
(8c67903) and `adapt.delta` (0.9 then, 0.99 now, the step-size fix in the section
below) -- so the two are not like-for-like and the agreement between the old
12-seed means (0.943 / 1.012) and the new 40-seed means (0.933 / 1.011) is
recorded as descriptive, not as evidence that the estimator is unchanged. The
divergences that motivated that section are gone at the shipped setting: 0 across
all 80 fits, against 68 (max 13) for car_proper at 0.9.

The interval coverage, which is the calibration claim the block exists to make,
is unaffected throughout: 1.00 on both hypers and both field types at 12 and 20
seeds, 0.97 at 40.

### The smoke block in the same file is fragile the same way

`dev_notes/probe_ocsn_smoke.R`, seeds 5001-5040, beta cover arm at phi = 20, the
block's own controls (1200 kept after 800 warmup, 1 chain), 47 s per fit.

| n | worst of beta 1-6 | max bias 1-6 (gate 0.35) | bias log_phi (gate 0.45) |
|---|---|---|---|
| 4 | psi_occ_cov1 | 0.161 | 0.080 |
| 8 | psi_(Intercept) | 0.096 | 0.050 |
| 12 | psi_occ_cov1 | 0.068 | 0.072 |
| 20 | psi_occ_cov1 | 0.119 | 0.040 |
| 40 | psi_occ_cov1 | 0.104 | 0.012 |

Field cor 0.775 mean / 0.510 min over 40 seeds, 2 divergences in total.

This block averages 4 seeds, and the occupancy coefficients carry a per-seed SD
of 0.35 (`psi_(Intercept)`) and 0.35 (`psi_occ_cov1`) against detection at 0.16
and cover at 0.06-0.10. A 4-seed mean therefore carries a standard error of 0.18
against a 0.35 band -- one standard error of headroom, and it is the assertion
that failed on the pre-8c67903 draws (run 32672430330, `all(bias[1:6] < 0.35)`)
while the field_sd band above passed, the two swapping sides when the
realisations moved. Raised to 12 seeds, where the standard error is 0.10 and the
measured worst coefficient is 0.068. Bands unchanged.

Together the two changes add about 20 minutes to this file (14 for the hyper
block, 6 for the smoke block) against a shard that has run 2.5-3.4 h under a
350-minute cap.

### car_proper divergences are a step-size artifact, and adapt.delta clears them

`Q(rho) = D - rho W` approaches the intrinsic (rank-deficient) limit as rho -> 1,
and an ICAR-simulated field pushes the sampled rho there (posterior means
0.78-0.97 against bounds [0.5, 0.99]), so the field's near-null direction
stretches against the psi intercept. `dev_notes/probe_204_adapt.R`, seeds
7002 / 7009 / 7010 / 7001:

| adapt.delta | divergences | rho mean | field_sd mean | s / fit |
|---|---|---|---|---|
| 0.80 | 4 / 1 / 16 / 13 | 0.887 / 0.957 / 0.931 / 0.974 | 0.978 / 1.086 / 1.004 / 1.403 | 37-70 |
| 0.95 | 6 / 0 / 2 / 0 | 0.888 / 0.959 / 0.929 / 0.977 | 1.041 / 1.120 / 1.025 / 1.443 | 51-101 |
| 0.99 | 0 / 0 / 0 / 0 | 0.886 / 0.954 / 0.929 / 0.976 | 1.013 / 1.090 / 1.010 / 1.477 | 86-147 |

The posterior does not move (rho identical to three decimals across the three
settings), so the divergences were the step size, not a region the chain was
missing. Hence the `occu_cover_spatial` row in `.TOBS_FAMILY_DEFAULTS`
(`adapt.delta = 0.99`, roughly 2x wall time).

### Cost of the parameterisation car_proper rho did NOT need

Reading `Q(rho) = D - rho W` literally makes a sampled rho a per-leapfrog dense
Cholesky. Measured on square lattices:

| n | dense `chol` | one-off eigen | chol at 1e5 leapfrog steps |
|---|---|---|---|
| 400 | 0.010 s | 0.09 s | 0.3 h |
| 1024 | 0.124 s | 1.27 s | 3.4 h |
| 2025 | 0.974 s | 9.69 s | 27.1 h |

The eigenbasis of the symmetrically normalised adjacency removes it entirely:
`Q(rho) = D^{1/2}(I - rho Lambda)D^{1/2}` gives a rho-independent `B1 = D^{-1/2}U`
with per-column weights `(1 - rho lambda_j)^{-1/2}`, so the one-off eigen above is
the whole cost and every step is the same matvec the pinned block already did.

### Shipped verification run

`NOT_CRAN=true`, `TULPAOBS_FAST` unset, against the shipped defaults (so
`adapt.delta = 0.99` on this path): `test-occu-cover-spatial-nuts.R` 62.4 min and
`test-occu-cover-nuts.R` 1.7 min, both 0 fail / 0 error / 0 skip. That run
includes the 24-fit coverage loop above and the pre-existing recovery blocks, so
`max(divergences) <= 5` per fit holds at the shipped adaptation target across all
24 sampled-hyper fits.

### Field-SD scale, and why `field_sd` is the reported quantity

`.tobs_field_load(adj, "icar", tau = 1, ...)` gives the field geo-mean marginal SD
`sqrt(scale_q)` at `sigma = 1`: measured 0.7426 (36 cells) / 0.7777 (64 cells),
matching `.occu_cover_icar_scale` = 0.5515 / 0.6048. The simulator's `f` carries
geo-mean marginal variance 1, so the fitter's `sigma` and the simulator's are NOT
the same number, and the bym2 structured block normalises differently again
(`compute_bym2_scale` = 2.6949 / 2.7732). `field_sd` -- the geo-mean marginal SD
the block's own covariance implies at a draw's hypers -- is the one summary all
three kinds share, so it is what the truth is stated in.

## A second (SVC) field on the occu_cover NUTS path (#214)

The sampler carries one block per coupled field. Everything below was measured on
Windows, R 4.6.0, against the shipped defaults unless a control is named.

### The target, on two blocks

Hand-built pairs of blocks (intercept + one weighted by a standardised per-site
column), 16 cells, J = 4, every hyper sampled, one theta drawn off the Laplace
mode. Same battery the single block was held to in #204:

| field | C++ vs R oracle, lp | grad | analytic vs central FD, max | hyper coords only | worst hyper | prior flat in t |
|---|---|---|---|---|---|---|
| icar | 0.0e+00 | 1.8e-15 | 7.6e-10 | 2.4e-10 | `sigma` (block 1) | 6.9e-16 |
| bym2 | 0.0e+00 | 8.9e-16 | 7.2e-10 | 1.6e-10 | `sigma` (block 1) | 1.4e-15 |
| car_proper | 3.6e-15 | 7.1e-15 | 7.8e-10 | 2.1e-10 | `sigma` (block 2) | 1.4e-15 |

The FD check is not vacuous: the hyper coordinates carry scores of 0.59-0.80 at
that theta. Prior flatness is `sd/mean` of the implied density in the axis
coordinate at `raw = 0`, over `u` in [-2, 1.5], per hyper of BOTH blocks.

Byte identity, one field: the top-level spec spelling and the one-element
`field_blocks` list give `identical()` lp and gradient, in C++ and in the R
oracle. End to end, a whole one-field fit (16 cells, J = 4, icar + `copy()`,
2 chains x 200 draws) is `identical()` before and after on `means`, `sds`,
`spatial_field`, `hyper_mean`, `log_lik` and the raw draw / hyper-draw matrices.

### Recovery, two icar fields

6 seeds (4001-4006), 36 cells, J = 5, truth `sigma = 0.8` / `alpha = 1.0` and
`sigma_trend = 0.7` / `alpha_trend = 0.9`, 2 chains x 800 kept draws after 600
warmup, ~26 s per fit at the shipped defaults:

| quantity | mean over seeds | truth |
|---|---|---|
| cor(intercept field, truth) | 0.637 | - |
| cor(trend field, truth) | 0.461 | - |
| `field_sd` | 0.817 | 0.8 |
| `field_sd_trend` | 0.872 | 0.7 |
| `alpha` | 0.996 | 1.0 |
| `alpha_trend` | 0.838 | 0.9 |

0 divergences on every seed, max Rhat 1.013, cover-dispersion bias 0.090. The two
amplitudes separate cleanly, which is the property a shared hyper would break.
Surface correlation is the weak part and stays weak by construction: one binary
occupancy observation per node, and the trend surface is seen only through its
weight column -- the per-seed spread is 0.42-0.82 (intercept) and 0.06-0.83
(trend). At 49 cells / J = 6 over 4 seeds the intercept surface improves (0.659)
and the trend does not (0.446), so the gates are set from the 36-cell fixture the
test runs. The same loop before the defaulted amplitude axis was thinned (324
warm cells rather than 144, 44 s per fit) gave 0.627 / 0.468 and field SDs
0.820 / 0.888: thinning an axis over its own span moves nothing but the cost.

### Why the warm fit forces a tensor grid

Each field adds its own `(sigma, alpha[, rho])` axes, so two icar fields take the
warm joint to 4 axes -- past the engine's CCD threshold, where it designs a
mode-centred star instead of a tensor. The sampler reads each axis's realised
span as the support of that hyper's flat prior, and a star's column range is a
design radius: on a 30-cell chain fixture with requested axes `sigma` [0.2, 1.5]
and `alpha` [0.2, 2.0], the CCD warm produced `b1.alpha` nodes from -1.23 to 4.70
(25 cells), the derived bounds became [1.49, 4.70], and the fit sampled alpha to
a posterior mean of 2.79 -- outside the axis the caller asked for. With
`integration = "grid"` the same fixture gives 81 cells and bounds exactly
[0.2, 1.5] / [0.2, 2.0].

The tensor is a product over blocks, so a DEFAULTED axis is thinned to three
nodes spanning the same range when a second field is present: the sampled prior
is defined by the span alone, so it is unchanged, while the warm fit drops from
900 cells (5 sigma x 6 alpha, squared) to 144. Provenance is the `auto_grid()`
mark, not whether the argument arrived: `copy(spatial())` with no amplitude hands
this path the DEFAULT alpha axis through `control$alpha.grid`, and reading that
as a user pin put the two-field warm grid over the engine's 2048-cell cap
(car_proper 2916 cells, bym2 7776 -- the latter because the trend block's own
`alpha.grid.trend` and bym2's engine-defaulted mixing axis were both unthinned).
Thinned on the mark, every areal kind fits at defaults on a 16-cell fixture:
icar 1.4 s, bym2 1.6 s, car_proper 0.8 s. An axis the caller chose
(`alpha = grid(c(...))`, `control$sigma.grid`) keeps its nodes.

## Posterior SBC on the coupled `occu_cover` (`R/sbc.R`, #207)

`sbc()` runs the posterior experiment of `tulpa::sbc()` (gcol33/tulpa#380,
Sailynoja et al. 2026 Algorithm 2) on a fitted `tobs_fit`. Everything below was
measured on tulpa 0.0.195 / tulpaObs 0.0.193, Windows, R 4.6.0.

### Fixture and cost

Chain adjacency, `N = 50` cells, `J = 6` visits, lognormal cover, shared `icar()`
field on the occurrence arm copied onto the cover arm, `phi_pos` on the outer
grid. `n.sim = 100`, `n.draws = 1000`, `n.ref = 200`, `seed = 0`. Observed fit
0.3 s; one augmented refit on the pooled 100 cells 0.75 s; the whole 100-sim run
**204.6 s** on the pinned grid (9 x 11 x 7 = 693 outer cells) and **~105 s** on
the smaller defaulted grid. Both premises report `verified` -- the pooling guard
and the fresh-groups guard, the latter only because `group_ids()` is supplied.

### What the observed fit recovers

Truth `beta_occ = (0.2, 0.6)`, `beta_p = (0.4, -0.5)`, `beta_pos = (-1.386,
0.3)`, `sigma = 0.8`, `alpha = 1`, dispersion 0.4. Fit: 0.312 / 0.545, 0.381 /
-0.606, -1.359 / 0.276, sigma 0.511, alpha 0.615, phi_pos 0.408.

### The calibration read, pinned grid, n = 100

Uniformity p-values on the `posterior` (reported) arm against the `narrow`
control (the same draws with their SD divided by 1.25, same fits, same
simulations):

| quantity | posterior ks | posterior p | narrow p |
|---|---|---|---|
| `psi_(Intercept)` | 0.063 | 0.890 | 0.060 |
| `psi_occ_cov1` | 0.094 | 0.579 | 0.052 |
| `p_(Intercept)` | 0.086 | 0.795 | 0.195 |
| `p_det_cov1` | 0.123 | 0.176 | 3.5e-04 |
| `pos_pos_cov1` | 0.130 | 0.133 | 1.8e-03 |
| `alpha` | 0.051 | 0.521 | 3.3e-16 |
| `pos_(Intercept)` | 0.141 | 2.8e-03 | 4.9e-03 |
| `sigma` | 0.718 | 0 | 0 |
| `sigma_pos_field` | 0.790 | 0 | 0 |
| `disp` | 0.085 | 1.2e-06 | 2.0e-09 |
| `log_lik` | 0.804 | 0 | -- |

The gated set is the first six: `min(p)` = **0.133** on the reported posterior
against **3.3e-16** on the mis-scaled control, fourteen orders of magnitude
apart, which is why the test gates on `min(p_unif)` either side of 1e-3 rather
than on the in-band indicator. Each band holds at 0.95 *within one quantity*, so
requiring six of them at once fails ~26% of the time on a perfectly calibrated
algorithm; the p-value gate is the multiplicity-aware form of the same read.

`alpha` carries the demonstration on its own: p = 0.52 as reported, 3.3e-16 once
its width is distorted. It is the DERIVED quantity the deliverable reports, and
it is read per draw off the outer grid (cells sampled by their own normalized
weight), so what is scored is its grid-marginalized posterior.

### The departures, and the mechanism behind them

`sigma` puts **81 of 100** PIT values in the top decile and `sigma_pos_field`
**87 of 100**: the truth is almost always ABOVE the augmented posterior, i.e.
pooling shrinks the shared field SD. `log_lik` is 90 of 100 in the top decile.
`alpha` -- the RATIO of the two field SDs -- is unaffected, and so are the
covariate slopes.

That pattern is what an improper direction in the pooled field prior produces,
and the improper direction is there. Posterior SBC needs the replicate on FRESH
cells, so the pooled areal graph is block-diagonal, i.e. DISCONNECTED. Measured
on three pooled fits, the two blocks' field means come back exactly equal and
opposite (-0.0126 / +0.0126, -0.0140 / +0.0140, +0.0546 / -0.0546): the fit
applies ONE global sum-to-zero, not one per connected component. `solve()` on
the pooled precision with a single global constraint is singular (reciprocal
condition 5.6e-18), so the block-level contrast has no prior at all. The
generator draws each block centred, the model leaves the contrast free, and the
shared SD absorbs the difference. `.occu_cover_icar_scale()` is NOT the culprit:
it returns 4.53051 on a 30-node chain and 4.53051 on two of them, identical, as
it must for a block-diagonal precision.

The occurrence-arm field SD, the cover intercept (which confounds with the field
level) and the joint-statistic rank are therefore NOT gated. Follow-up:
gcol33/tulpaObs#212, one sum-to-zero per connected component. Repros in
`dev_notes/repro_212_pooled_icar_constraint.R` and
`dev_notes/repro_212_constrained_precision.R`.

### Why the hyperparameter grid is pinned in the fixture

With the DEFAULT (auto-recentred) grids, the observed fit and the augmented
refits integrate different supports: measured over three refits, `alpha` differed
on 1 of 3, `phi_pos` on 3 of 3 (extra nodes interpolated around the mode), and
one refit recentred `sigma` from [0.100, 3.000] onto [0.031, 0.503] entirely.
The truth is then an atom of one support scored against a predictive on another.
Pinning `control$sigma.grid` marks it a USER grid, which is not recentred. It
moves `alpha` from outside the band (ks 0.359, p = 0) to inside (ks 0.051,
p = 0.52) and leaves the slopes where they were, so the unpinned `alpha` read was
measuring grid placement.

### Why the cover dispersion has to be on the outer grid

Left off it, the joint engine fixes the dispersion per data set, so the observed
fit and every refit hold it at a different value and the replicate is generated
at one while scored under another. `sbc()` warns and drops it from the
scored set (its posterior is a point mass, and `disp` took ONE unique value
across 100 draws in that configuration). With `control$phi.grid.pos` it is
estimated and scored.

## EM detection block carries the E-step weight (`.tobs_encode_det_block`)

Issue #208 (tulpaObs), enabled by gcol33/tulpa#385 (tulpa >= 0.0.184; measured
against installed 0.0.195, the DESCRIPTION floor).

**The engine really does consume `weights` on the spatial route.** Checked
directly against `tulpa::tulpa_laplace(spatial = <spde spec>)`, not taken from
the issue: an integer weight vector reproduces literal row replication (the same
design re-rowed through `.tobs_spde_broadcast_spec()`) to max |diff| **2.8e-17**
on both fixed effects; `weights = NULL` is bit-identical to `weights = 1`; a
fractional weight vector moves the mode by 0.0088, so it is not being silently
dropped. Pinned by `test-em-det-block-weights.R`, which fails if the engine ever
reverts.

**What `round(w * nd)` / `round(w * nv)` was destroying.** The rounding is not
symmetric between successes and trials, because the rows it erases are exactly
the low-weight ones and a low-weight row is (almost always) a row with no
detection. On the integrated fixture below (N = 300, 2 sources, E-step weights
median 0.062): **161 of 300 rows** collapsed to a `(0, 0)` row, effective trials
`sum(round(w nv))` = 265.0 against the exact `sum(w nv)` = 285.3 (**-7.1%**),
while effective successes were unchanged (132 vs 132). Dropping 7% of the trials
and none of the successes biases p_hat UP, and the occupancy arm compensates by
biasing psi DOWN. Both show up in the recovery table.

**Recovery, old encoding vs new, 8 seeds each.** Truth: p intercept -0.2,
p slope 0.4, psi slope 0.7; smooth Matern-like field `0.9 cos(3 lon) sin(3 lat)`
on logit(p).

| | | old (count-scaled) | new (weights) |
|---|---|---|---|
| single-season, N = 500, J = 3 | p intercept bias / MAE | +0.151 / 0.156 | **-0.047 / 0.069** |
| | p slope bias / MAE | +0.042 / 0.061 | **+0.014 / 0.050** |
| | psi slope bias / MAE | -0.057 / 0.079 | **-0.019 / 0.061** |
| | field cor | 0.814 | 0.815 |
| integrated, N = 500, 2 sources | p intercept bias / MAE | +0.202 / 0.202 | **-0.004 / 0.135** |
| | p slope bias / MAE | -0.030 / 0.112 | **+0.032 / 0.095** |
| | psi slope bias / MAE | -0.070 / 0.080 | **-0.030 / 0.056** |
| | field cor | 0.740 | 0.792 |

Every mean bias and every MAE improves; the field correlation is a wash on the
single-season arm and better on the integrated one. The direction is the one the
mechanism predicts (the p intercept, which the erased trials hit hardest, moves
from +0.15 / +0.20 to roughly zero). The residual -0.047 on the single-season p
intercept is the ordinary Laplace-EM attenuation, a third of what the rounding
was adding on top of it.

**The non-SPDE branch does not move.** It already passed `weights = w[dk]`, so
it is the control: block encodings compared against HEAD's
`build_single_callbacks` / `build_integrated_callbacks` (sourced from
`git show HEAD:R/laplace_callbacks.R`) over four weight vectors each are
`identical()`, and full fits agree bit-for-bit on `means`, `sds` and `vcov`
across five fixtures including the `laplace_gibbs` correction path.

**The integrated detection-arm field is not reachable from `tobs()`.**
`.tobs_build_integrated()` (`R/occu.R:247`) rejects any structured term on a
per-source detection formula, and process membership comes only from which
formula a term is written in, so `spatial_det` inside
`build_integrated_callbacks()` is always NULL in a public fit. The integrated
numbers above were produced by handing the callbacks a mesh spec with
`shared = c(FALSE, TRUE)` and driving `tulpa::tulpa_em_laplace()` directly. The
SPDE coverage matrix at the top of `R/laplace_helpers.R` calls that cell "wired
(per source)", which describes the fitter, not the parser.

## Component count of an areal graph (`.tobs_check_graph()`, #212)

Measured on tulpa 0.0.195 / tulpaObs 0.0.194, Windows, R 4.6.0, against the
engine's own `log_prior_icar` through `tulpa:::cpp_test_log_prior_icar` (which
derives the component partition from the adjacency it is passed, so it reads the
same partition a fit does).

### The constraint was already per component

Two disjoint 50-node chains, `tau = 1.7`, each component's values centred:

| read | measured | one global constraint would give |
|---|---|---|
| pooled log-prior minus sum of the two independent ones | `0.000e+00` | non-zero |
| `+c` on component 1, `-c` on component 2, `c = 0.10` | `-0.850000` | `0` |
| same at `c = 0.25` | `-5.312500` | `0` |
| same at `c = 0.50` | `-21.250000` | `0` |

The charged values are exactly `-0.5 * tau * (n1 + n2) * c^2`. The shift keeps
the GLOBAL sum at zero, so a single global sum-to-zero cannot see it at all.

Two reads the equal-size fixture cannot make, both exact: a 70-node plus a
30-node chain is still additive to `-5.7e-14`, and permuting the node ordering
of that graph leaves the log-prior unchanged to `5.7e-14`. The second is what
pins component MEMBERSHIP rather than only the count -- a contiguous
equal-split assumption reproduces neither.

An isolated node (a size-1 component) is charged `-0.5 * tau * v^2`, i.e. the
augmentation leaves it a proper `N(0, 1/tau)` effect and does NOT pin it to
zero, which a hard sum-to-zero on one node would do. That is the engine's
behaviour in the abstract; tulpaObs's own areal cover paths reject an isolated
node before fitting (`.occu_cover_icar_Q()`, `.occu_cover_adj_to_csr()`,
`R/occu_cover_nuts.R:903`, `R/ms_occu_cover_spatial.R:345`), and that decision
is unchanged -- the component report does not restate it.

### What the "one constraint" reading actually measured

`fit$spatial_field` is demeaned globally per field block
(`.occu_cover_demean_fields()`), so on two EQUAL-size components the reported
block means are exactly equal and opposite whatever the prior did. Three seeds
of the pooled `occu_cover` fixture, 50 + 50 cells:

| seed | RAW block means | reported block means | raw sum |
|---|---|---|---|
| 707 | `+0.026417 / -0.026866` | `+0.026642 / -0.026642` | `-0.000449` |
| 708 | `-0.013865 / +0.013130` | `-0.013497 / +0.013497` | `-0.000735` |
| 709 | `+0.060253 / -0.060560` | `+0.060407 / -0.060407` | `-0.000307` |

`reported == raw - mean(raw)` holds on every seed. The raw component means are
not forced equal and opposite; the reported ones are, by the demean. The
global demean shifts every cell by the same constant, so it preserves the
between-component contrast exactly (seed 707: raw `0.053283`, reported
`0.053284`) and only removes the grand mean.

### Single-component bit-identity

A connected 50-node chain `occu_cover` fit, run with the component check in
place and then with `.tobs_check_graph()` restored to its pre-change body:
`identical()` is TRUE on `means`, `sds`, `vcov`, `spatial_field`,
`log_marginal`, grid `modes`, `weights`, `theta_grid`, `log_det_Q` and
`converged`, and the connected graph emits zero messages.

## SBC registry beyond `occu_cover` (`R/sbc.R`, #207)

Eight single-response families registered and run end to end through
`sbc()` at `n.sim = 100`, `n.draws = 1000`, `n.ref = 200`, `seed = 0`,
`controls = "narrow"` (`bad.factor` 1.25). Both premises report `verified` on
every one. `min p_unif` is over the arm coefficients; the joint `log_lik` arm is
listed separately because on one family it carries the whole read.

| family | wall | min p_unif (coefs) | log_lik arm | narrow control min |
|---|---|---|---|---|
| `occu` | 10 s | 0.437 | 0.855 | 5.4e-4 |
| `count` (poisson) | 1.5 s | 0.049 (`mu_x`) | 0.464 | 3.6e-3 |
| `abun` | 63 s | 0.096 | 0.266 | 4.9e-4 |
| `royle_nichols` | 60 s | 0.235 | 0.012, OUTSIDE | 7.9e-5 |
| `occu_ttd` | 3.5 s | 0.299 | 0.787 | 8.3e-4 |
| `fp_occu` | 4 s | 0.077 | 0.481 | 2.4e-3 |
| `removal` | 145 s | 0.134 | 0.240 | 4.7e-4 |
| `distance` | 75 s | 0.291 | 0.291 | 1.0e-3 |

Two reads chased across seeds rather than reported once:

- `count`'s `mu_x` = SEED NOISE. At `n.sim = 300`, seeds 0/7/21 give 0.027 /
  0.511 / 0.126. One low read in nine at three seeds is what a calibrated
  algorithm looks like; do not re-tune on it.
- `royle_nichols`'s `log_lik` = REPRODUCIBLE departure: 0.0116 / 0.0413 / 0.0288
  at seeds 0/7/21, outside the band every time, while all four coefficient
  marginals sit inside. The joint-statistic arm doing what it exists for -- the
  Gaussian pseudo-posterior gets each marginal about right and the lambda/r
  dependence wrong. Tracked separately; NOT an argument for loosening a
  coefficient assertion.

Refusals built in rather than approximated: a fit carrying a structured term
(its field is a latent shared across sites that theta does not hold -- that is
the coupled `occu_cover` route), and a fit carrying a visit-level observation
design (the `single` generator assembles detection from the site-level design
only, so a visit-level column would be scored by the refit and absent from the
data it sees).

## SBC field-SD pile-up: generator scale + outer-grid atomicity (#213)

Two independent causes behind 81/100 and 87/100 top-decile ranks on the two
field SDs while every coefficient and `alpha` stayed uniform.

**Cause 1, the generator's units (the pile-up).** The joint nested-Laplace
engine hands its ICAR block the RAW `Q = D - W` at `tau = 1` and carries the
amplitude in the arm scale; no Sorbye-Rue scaling exists on that path
(`scale_factor` is bym2-only). So the fit's `sigma` multiplies a field of
geo-mean marginal SD `sqrt(scale_q)`, while `.occu_cover_draw_icar_field()`
returns a field normalised to 1. One common factor on both arm SDs; `alpha`,
their ratio, unaffected.

Measured on chain(60), 12 seeds, data generated at geo-mean SD 0.8:
`sigma_fit / 0.8 = 0.250`, 95% CI **0.179-0.350**. Predicted `1/sqrt(scale_q)`
= 0.332 sits inside; `c = 1` is excluded by ~4x. The eta-scale field itself was
always fine (sd ratio 0.88, cor 0.97) -- the LABEL was wrong, not the surface.
`sqrt(scale_q)` = 1.74 / 2.13 / 2.75 / 3.01 for chains of 20 / 30 / 50 / 60. A
pooled 2-component graph carries the same constant as one component.

**Hypothesis 2 (grid truncation) EXCLUDED**, not argued away: over 100
simulations the truth sits at the top node 0/100, bottom node 3/100, and 80/100
at node 5 of 9, interior. `outer_grid_placement = "fixed"`,
`recenter_declined = "axis_pinned"` (correct: the fixture pins the axis).

**Demean lead measured, not the answer**: on the pooled 2-component fixture the
two component levels are +-1.40e-02 against a field SD of 2.95 = 0.48% of the
field's own scale.

Before / after at n.sim = 100, seed 0, on the fixture's ORIGINAL 9/7-node grids:

| quantity | top decile | mean PIT | p_unif |
|---|---|---|---|
| sigma | 81 -> 18 | 0.865 -> 0.605 | 0 -> 4.9e-04 |
| sigma_pos_field | 87 -> 14 | 0.923 -> 0.579 | 0 -> 0.050 |
| log_lik | 90 -> 18 | 0.940 -> 0.558 | 0 -> 0.083 |
| pos_(Intercept) | 16 -> 11 | 0.582 -> 0.525 | 2.8e-03 -> 0.813 |

**Cause 2, outer-grid ATOMICITY (the residual).** `sigma`'s posterior IS its
axis -- the outer grid is its entire support -- so the axis's resolution is the
resolution of the predictive the truth is ranked against. At 9 nodes 1000 draws
carry 6 distinct values and the rank ECDF is a 6-step function, which a
continuous-uniform read scores as a departure whatever the fit does. At 21 sigma
/ 17 phi nodes every quantity is inside: sigma 0.714, sigma_pos_field 0.697,
alpha 0.511, disp 0.819, log_lik 0.491, min over all 0.056. Cost: acceptance run
3.6 -> ~13 min.

**Isolating control** (fine grid, generator constant forced back to 1): sigma
42/100 top decile p = 0, sigma_pos_field 52/100 p = 0. The GENERATOR moved them,
not the resolution -- the two causes are separable and both are real.

Consequence for anyone reading a spatial `occu_cover` fit: on the joint
nested-Laplace path `sigma` is a RAW amplitude, while `simulate_occu_cover()`
and the #204 NUTS `field_sd` state geo-mean marginal SD. Compare a fit against a
simulator on the FIELD (sd / cor), never by reading `sigma` off both sides. Also
`fit$spatial_field` is the LATENT field: the eta-scale surface is
`sigma * spatial_field`.

## Outer Pareto-k diagnostic cost (occu_cover joint, tulpa#118, issue #101)

`joint` engine (`R/occu_cover_joint.R`) vs `v3_nested`: ~150-300x faster at N=100,
completes at N=200+ where v3_nested does not. Outer Pareto-k diagnostic
(`control$diagnose.k`, default OFF): 84-98% of joint-fit wall time, since it
re-solves the inner Laplace at `k.samples`=200x on the full field vs the grid's
~30-70 points. Per-phase profiling: the binding per-solve cost is the
per-Newton-iteration Hessian/grad SCATTER (beta-arm digamma/trigamma fill,
73-83% of a solve), NOT the factorize (flat ~0.5ms, 8-12%, not super-linear up
to 1100 cells) -- corrected an earlier assumption that factorization dominated.
tulpa#118 sped the diagnostic 2-4x (Shamanskii reuse via `.K_DIAG_REFRESH` ->
grad-only scatter; loosened inner tolerance `.K_DIAG_TOL`=1e-4; near-neighbour
batch order); k-hat stayed byte-stable, externally validated against
`loo::psis`/`posterior::pareto_khat` on real EVA ratios. Knobs
`tulpa.kdiag.{refresh,tol,reorder,capture}`.

## royle_nichols() SBC joint log_lik arm (gcol33/tulpaObs#219)

Investigated per the issue's own checklist; no code defect found, closing as
documented approximation behaviour.

`royle_nichols()` fits by `optim(..., method = "BFGS", hessian = TRUE)` over the
closed-form Poisson-sum marginal (`.tobs_fit_royle_nichols` -> the shared
`.tobs_bfgs_marginal_fit`, `R/laplace_helpers.R:260`), reporting `vcov = solve(
observed information)` -- a single Gaussian (Laplace/Fisher) approximation to the
posterior at its mode. `tobs_sbc()`'s coefficient marginals pass on this fit
(seeds 0/7/21, min p_unif 0.235) while the joint `log_lik` arm rejects hard
(0.0116 / 0.0413 / 0.0288). Checked the issue's three items:

1. **Where the vcov comes from, and whether lambda/r is the discordant block.**
   Confirmed it is `solve(hessian)` at the BFGS mode, no numerical shortcuts.
   Profile-likelihood check on `r_x` (fix at a grid, re-optimize the other 3
   coordinates, compare `2*(profile_nll - mle_nll)` against `qchisq(.95, 1)` and
   against the Wald quadratic `((v - mle)/sd)^2`) at N=100, seeds 0/7/21: the
   profile 95% CI is 11-23% NARROWER than the Wald CI (seed 0: 0.979 vs 1.279;
   seed 7: 1.024 vs 1.147; seed 21: 0.693 vs 0.776) -- the 1-D marginal is close
   to quadratic, mildly conservative, not badly broken. A crude importance-
   resampling check (draw `N(mode, 4*V)`, reweight by the exact marginal
   log-lik) is far more telling on the JOINT: effective sample size collapsed to
   1-25 out of 20000 draws across the three seeds, meaning almost all mass of a
   4x-variance-inflated Gaussian proposal lands far from where the true
   likelihood surface actually sits. A well-approximated unimodal posterior does
   not do that to a 4x-inflated proposal; a posterior concentrated along a
   lower-dimensional CURVED manifold does. The reweighted correlations moved
   from the Gaussian's |rho| ~ 0.3-0.9 toward |rho| ~ 0.9-1.0 in every seed
   (numerically unreliable at that ESS, but directionally consistent across all
   three) -- lambda and r trade off along a curved ridge the linear-correlation
   Gaussian cannot represent, exactly the mechanism the issue named.
2. **Whether it survives at higher counts.** Same profile-vs-Wald check at
   N=1000 (seed 0): gap narrows from 23% to 13% (profile 0.258 vs Wald 0.289).
   Shrinking with N is what a finite-sample curvature artifact should do, not
   what a permanently mis-specified vcov would do.
3. **Whether `abun` shares the signature.** It does not (issue's own measurement:
   `log_lik` p = 0.266, inside the band). `abun`'s detection arm is a plain
   logit-linear p; `royle_nichols`'s is `p_site = 1 - (1 - r)^N`, N-dependent and
   nonlinear in the pair (lambda, r). That link is the one structural difference
   between a family whose joint SBC arm passes and one whose does not, and it is
   exactly the kind of link that turns a linear correlation into a curved one.

Conclusion: individual coefficient marginals are close to quadratic and mildly
conservative (matches the SBC coefficient pass); the joint dependency between
`lambda` and `r` follows a curved ridge that ONE Gaussian cannot capture, and the
joint `log_lik` SBC statistic is doing exactly what it exists to catch. Not a
bug in the vcov computation, not a truncation artifact (`K_max` unchanged
between checks), and not something a tolerance change should paper over --
`royle_nichols()` has no `nuts` method to sample the exact joint instead, so this
is a documented limitation of the `laplace`-only backend, not a defect to close
by code change.

## BYM2 scale factor: Riebler marginal-variance constant (#232)

The constant is `s = geomean(diag(Q^+))` for the intrinsic ICAR precision
`Q = D - W`, NOT `geomean(non-zero eigenvalues of Q)` -- the two are different
numbers and the gap grows with the graph. `s` is what `INLA::inla.scale.model()`
applies. Measured against it (`inla.scale.model(Q, constr = list(A = matrix(1,
1, n), e = 0))`, factor read off `Qs[1,1] / Q[1,1]`), and against tulpa's
`compute_bym2_scale()` (unexported, reach it with `getFromNamespace()`), which
returns the ENGINE-side reciprocal `1 / sqrt(s)`:

| lattice | `.bym2_scale()` | `inla.scale.model` | old eigenvalue mean | `.bym2_engine_scale(s)` vs tulpa |
|---|---|---|---|---|
| 5x5   | 0.516386 | 0.516386 | 2.646529 | 1.391595 vs 1.391595 |
| 10x10 | 0.644879 | 0.644879 | 2.831882 | 1.245263 vs 1.245263 |
| 20x20 | 0.765027 | 0.765026 | 2.984944 | 1.143304 vs 1.143304 |

The old constant is not the reference under a reciprocal or a square root
either: `1/sqrt(eigmean)` FALLS with graph size (0.6147 / 0.5942 / 0.5788)
where `s` RISES (0.5164 / 0.6449 / 0.7650), so the two cross near a 10x10
lattice. An absolute gap at ONE graph size is therefore not a separating
property -- `test-bym2-scale.R` asserts the size TREND plus a ratio range,
after an absolute-gap assertion failed at 10x10 (0.051) while holding at
5x5 (0.098) and 20x20 (0.186).

**Recovery** (`dev_notes/probe_232_bym2_constant.R`, nmix areal BYM2, 35 cells,
8 seeds, truth `sigma = 1`, `rho = 0.7`):

| constant | sigma | rho | field cor |
|---|---|---|---|
| Riebler `s` | 0.921 (bias -0.079) | 0.595 (bias -0.105) | 0.897 |
| old eigenvalue mean | 1.429 (bias +0.429) | 0.726 (bias +0.026) | 0.896 |

`sigma` is the decisive axis: the old constant inflated the field SD by 43% and
did it on 8 seeds out of 8 (per-seed range 1.13-1.91, never once below truth),
while the corrected one lands within 8% and straddles truth (7 of 8 below, mild
attenuation). Mechanism: the old constant shrinks the structured block by
`sqrt(2.773 / 0.605) ~ 2.1`, and the fit buys that back through BOTH `sigma`
and `rho`.

`rho` is mildly WORSE in the mean under the correct constant (-0.105 vs +0.026)
and this is not a regression to chase: per-seed `rho` scatters 0.28-0.80 (new)
and 0.40-0.89 (old) at 35 cells, so the SE over 8 seeds is ~0.065 -- the new
bias is ~1.6 SE, the old ~0.4 SE, neither decisive at this sample. Correcting
the scale removes the `sigma` compensation and leaves `rho` carrying its own
weak identification at this graph size. A `rho` verdict needs more cells, not a
different constant.

**Field cor is 0.897 vs 0.896** -- identical. Every BYM2 recovery fixture in the
package asserts field cor, so a green suite is not evidence about this constant
in either direction. That is why the issue asked for `sigma`/`rho`.

Verified at close: `test-bym2-scale.R` (46/46), `test-nmix-spatial-bym2.R`
(pass; its one warning is a `K_max` truncation notice on the fixture, unrelated).
NOT swept: the other ~57 test files mentioning bym2, most of which only exercise
term parsing or gates. `test-occu-cover-field-sd-units.R` is the one whose name
says it asserts on the axis this moves; the smoke and full-recovery tiers cover
the rest.

## ms_abun NUTS per-(species, site) latent-N ceiling (#233)

CONTRACT (CLAUDE.md points here rather than restating it; that file is at its
150k cap). This is a DIFFERENT mechanism from the Laplace headroom guard in the
`abun()` truncation section -- do not conflate them.

- The ceiling is set PER (species, site) by `.tobs_ms_nmix_nuts_kmax()` from an
  exact tail quantile (`qpois` / `qnbinom`, `.MS_NMIX_NUTS_TAIL_TOL` 1e-12) at
  `.MS_NMIX_NUTS_INFLATE` x that cell's WARM-MODE lambda. NOT a `max(y_i)`
  -relative headroom: `max(y_i)` understates N by exactly the detection rate, so
  a fixed headroom is tightest where detection is lowest = the wrong direction.
- The ceiling is FIXED at the warm-start mode and held for the whole chain, so
  the inflation factor is an EXCURSION margin, not a precision knob. A cell whose
  lambda wanders past its ceiling loses real posterior mass -- a wrong answer,
  not a slow one. `.tobs_ms_nmix_kmax_check()` re-reads the realised excursion
  off the draws on every fit (`fit$nuts$K_site_check`), so the margin is measured
  rather than assumed to have held.
- Threaded to the kernel as `K_site` (`nmix_precompute_site(K_site=)`).
  `compute_nmix_site`'s arity is UNTOUCHED -- it IS the `CountKernelFn` function
  -pointer contract the shared count-NUTS / count-Laplace drivers take, and a
  defaulted extra parameter silently breaks that conversion.
- C++ AND the R oracle REJECT a ceiling below a cell's own `max(y)` rather than
  silently raising it; raising it would hide a caller bug.
- The R oracle slices the SAME per-cell ceiling (`.tobs_ms_abun_nuts_marginals`
  takes `K_site`), which is what keeps it a real check of the capped target
  rather than of a different one.

Fixture: `simulate_ms_abun(n_species = 8, N = 40, J = 4, ...)`, seeds 201-220,
`n.iter = 300`, `n.warmup = 300`, Poisson (the simulator's default -- so every
number in this section is the Poisson path). Wall time on this box is +-2.6x
noisy; the state counts are byte-identical across runs, so the conclusions rest
on those.

Summed latent states per fit, one shared ceiling vs one ceiling per cell:

| | shared `K_max` | per-cell | saving |
|---|---|---|---|
| mean over 20 seeds (margin 8) | -- | -- | 2.04x |
| mean over 20 seeds (margin 4) | -- | -- | 3.09x |
| seed 211 (margin 8) | 128761 | 27934 | 4.61x |
| seed 203 (margin 8) | 96669 | 33080 | 2.92x |
| seed 212 (margin 8, mildest) | 36815 | 25411 | 1.45x |

The point is the flatness, not the ratio: `K_site_mean` stays in 77.8-108.3 over
all 20 seeds while the shared `K_max` swings 118-407. Under one shared ceiling a
seed's per-step cost is set by its single heaviest cell; per-cell resolution
decouples it, so the ratio is largest exactly where the shared ceiling was worst.

Excursion margin (`.MS_NMIX_NUTS_INFLATE`). The ceiling is fixed at the warm-start
mode and held for the whole chain, so it has to cover where the sampler GOES; a
cell wandering past it loses real posterior mass, which is a wrong answer rather
than a slow one. Max sampled/warm-mode lambda ratio over the 20 seeds: **3.93**
(seed 208; 3.76 seed 212, 3.68 seed 207). That figure is IDENTICAL at margin 4 and
margin 8, as it must be -- the ceiling caps the latent-N sum and never constrains
lambda -- so it measures the posterior, not the setting. Margin 8 is therefore a
2x safety factor over the worst measured excursion. Cost of the headroom, from the
warm-start-only probe across all 20 seeds:

```
  inflate  4 : 3.09x fewer states  (100% of the inflate-4 cost)
  inflate  6 : 2.42x fewer states  (128%)
  inflate  8 : 2.04x fewer states  (151%)
  inflate 12 : 1.64x fewer states  (186%)
  inflate 16 : 1.44x fewer states  (210%)
```

The probe's 2.04x predicted the end-to-end 20-fit result exactly. Coverage 0.975
at both margins (gate 0.85 pooled), divergences 0-6 per seed.

Negbin saves much less, and the reason is structural rather than a fixture
accident: the ceiling is an exact tail quantile at `.MS_NMIX_NUTS_TAIL_TOL` 1e-12,
and the NB tail is far heavier than the Poisson one, so each cell's own ceiling
lands near the shared one. On the `test-ms-abun-nuts.R` fixture at margin 8:

```
  poisson  K_max 136  K_site 36..136 (mean  95.4)  69% of cells < K_max  1.45x
  negbin   K_max 174  K_site 59..174 (mean 153.5)  31% of cells < K_max  1.14x
```

That is why `test-ms-abun-nuts.R` asserts that the ceiling BINDS (>25% of cells
strictly below `K_max`, total states strictly down) and not a states ratio: the
ratio tracks the excursion margin and the mixture's tail weight, both tuning, so a
pinned ratio would need re-deriving whenever either moved. Correctness is the
separate assertion -- capped == uncapped to 1e-8, and C++ == R oracle to 1e-9 WITH
the ceiling present, both of which the binding check keeps non-vacuous.

## Scalar nuisance AGHQ order floored at 3 (`ms_abun(mixture = "negbin")`, #234)

Fixture verbatim from `test-ms-abun-nb-rs.R:78-93`: `simulate_ms_abun(n_species = 18,
N = 100, J = 5, mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3), sd_lambda = 0.4,
sd_p = 0.35, mixture = "negbin", size = 5, sigma_logr = 0.5)`, seeds 501-520, fitted
`method = "laplace"`, `control = list(n.quad = 3L)`, tulpaObs 0.0.235 / 343f425.
Truth `mu_log_r = log(5) = 1.60944`. One fit is 12-59 min at 2 threads (mean 27.6),
measured 10-wide.

**The 2-node rule is not a coarse setting, it is below the order that represents the
integrand.** Holding `n.quad = 3` and moving ONLY `n.quad.scalar` -- the axis the
collapse is on -- seed 515:

```
n.quad.scalar   mu_log_r        se_log_r        z         covered
2               1.8494206734    0.0110225163    +21.772   no
3               1.7231670271    0.1886464461     +0.603   yes
5               1.7231670271    0.1886464461     +0.603   yes
9               1.7231670271    0.1886464461     +0.603   yes
```

3 / 5 / 9 are identical to all ten recorded decimals, not merely close: a cliff
between 2 and 3, then a flat plateau. So this is NOT an accuracy/cost tradeoff along
that axis and there is no reason to expose 2. Across the other 19 seeds the SE ranges
0.090-0.186 (mean 0.139), so 0.011 is off-scale, not a tail.

**Mechanism.** A 2-node Gauss-Hermite rule places two nodes and has no freedom left to
represent curvature; where the 1-D `log_r` posterior is not near-Gaussian the marginal
it returns can come out arbitrarily sharp. Seed 515 is exactly that regime -- one
species with `max(y) = 122` against seventeen ordinary ones (`K_max` 222; the healthy
control seed 516 tops out at 56), `sigma_log_r` 0.730 against 516's 0.373. The fit
does not complain because from the optimizer's side nothing went wrong.

`vcov` at seed 515 vs 516 (coordinate 5 is `log_r`), `sqrt(diag)`:

```
                    515 (2 nodes)   516 (2 nodes)
lambda_(Intercept)      0.158789        0.122560
lambda_abund_cov1       0.124808        0.090861
p_(Intercept)           0.175050        0.077682
p_det_cov1              0.070857        0.088662
log_r                   0.010901        0.122535
eigen min               1.065e-04       5.916e-03
eigen max               3.072e-02       1.546e-02
condition               2.884e+02       2.613e+00
```

The four coefficient coordinates are ordinary in both -- 515's are if anything WIDER.
Only `log_r` moves, its variance 126x smaller (1.1883e-04 vs 1.5015e-02), and it lands
essentially ON the matrix's minimum eigenvalue while nearly decoupled from everything
else (off-diagonals ~1e-4). So the marginal is not going singular; inverted, it says
the likelihood is extraordinarily peaked in `log_r`.

**Two candidate mechanisms were measured and REFUTED, do not revisit them.**
`sigma_log_r` does NOT collapse to the boundary (the `sigma_omega` failure the PC prior
at `R/nmix_laplace_re.R` addresses): it reads 0.730 against a target of 0.5, and the
per-species spread is recovered almost exactly (`sd(log r_s)` 0.6734 against that
seed's realised truth 0.6715 -- BETTER than the healthy control, 0.2899 against 0.4224).
Only the SE on the community mean is wrong, so a PC prior on that variance is the wrong
lever. Second, `n.quad` alone is NOT an escape hatch, because the scalar order was
clamped to `min(n_quad, n_quad_scalar)`: the family's own documented default
`n.quad = 5` still ran the `log_r` block at 2. That clamp is why the fix is a FLOOR
that overrides `n_quad`, not a raised default.

**Cost.** The tensor grid is `prod_b n_quad_b^dim_b` = `n_quad^4 * n_quad_scalar` for
this model (`p_lambda` = 2, `p_p` = 2, `log_r` scalar):

```
                          n.quad = 3   n.quad = 5
n.quad.scalar = 2              162         1250
floor 3                        243         1875   (1.5x)
scalar order = n_quad          243         3125   (2.5x)
```

So the floor gives back only part of the per-block-quadrature saving -- 1.5x on one
axis at `n.quad = 5`, where abandoning per-block quadrature entirely would be 2.5x.

**It is a bias fix as well as a seed rescue.** All 20 seeds, `n.quad = 3`:

```
                      n    bias      sd(mu)   mean(se)   ratio   coverage
nqs = 2               20   +0.0578   0.1740   0.1328     1.310   15/20 = 0.750
nqs = 2, excl 515     19   +0.0482   0.1733   0.1392     1.245   15/19 = 0.789
nqs = 3               20   +0.0160   0.1687   0.1348     1.251   16/20 = 0.800
nqs = 3, excl 515     19   +0.0109   0.1717   0.1320     1.301   15/19 = 0.789
```

Bias drops 77% (+0.0482 -> +0.0109) on the nineteen seeds that never had the collapse,
and EVERY seed's SE moves under the floor -- median -4.6%, range -10.6% to -0.5%
excluding 515, which moves +1611%. The 2-node rule was distorting all twenty fits
mildly and one of them catastrophically. Residual bias afterwards is +0.016 against a
sampling sd of 0.169, about 0.09 SE.

**This does NOT turn `test-ms-abun-nb-rs.R:94` green and the threshold was not
touched.** Coverage goes 15/20 -> 16/20 = 0.800 and the assertion is
`expect_gt(mean(covered), 0.8)`; 0.800 is not greater than 0.8. What remains is a pure
multiplicative miss on the SE with no bias left (`mean(z)` +0.104, t p = 0.73), carried
in #235. Max `|z|` over all 20 seeds is 2.421 and scaling every SE by 1.2512 gives
20/20, so no seed misses for any other reason.

**Threshold power, for whoever sizes #235.** The assertion needs >= 17 of 20:

```
true coverage 0.95 -> P(pass) 0.984
              0.93 -> 0.953
              0.90 -> 0.867
              0.85 -> 0.648
              0.80 -> 0.411
```

At genuinely nominal coverage the test still fails about 1 run in 60, and at 20 seeds
it cannot separate 0.90 from 0.95. Seed count is the only lever, against 27.6 min per
fit.

Probes: `dev_notes/nbrs_bisect/` (`dissect_s515_nq3.txt`, `dissect_s516_nq3.txt`, the
`nq5x_*` / `def5_*` arms).

## mu_log_r Wald calibration is conditional on sigma_log_r (`ms_abun(negbin)`, #235)

The residual after #234's quadrature floor. Filed as a ~24% understatement of the
`mu_log_r` SE; it is not that, and the correction matters because a uniform
inflation would over-widen the ~90% of fits that are already right.

Three arms, all `n_quad = 3` / `n_quad_scalar = 3`, fixture
`simulate_ms_abun(N = 100, J = 5, size = 5, sigma_logr = 0.5, ...)`, seeds
501-520, truth `mu_log_r = log(5)`:

```
 S    n   mean(z)   t p     k_hat   chi2 p   coverage        sd/mean(se)
 8   19   +0.376    0.172   1.182   0.116    17/19 = 0.895   0.970
18   20   +0.104    0.733   1.311   0.024    16/20 = 0.800   1.251
36   20   +0.376    0.114   1.060   0.315    18/20 = 0.900   0.895

pooled 59 fits:  k_hat 1.189  95% CI [1.008, 1.450]  chi2 p 0.020
                 mean(z) +0.284 (t p 0.066)   coverage 51/59 = 0.864
```

**NOT a finite-sample-in-S effect.** `k_hat` peaks at the shipped fixture instead
of decaying, and the raw bias does not decay either (+0.0681 / +0.0160 / +0.0431
at S = 8 / 18 / 36). The positive test that hypothesis predicted fails.

**NOT a t-df effect either.** Re-reading the same z values against `t(S-1)`
rather than 1.96 moves pooled coverage only 0.864 -> 0.881 and still rejects
nominal (binomial p = 0.027 against 0.95). A Hartung-Knapp-style df correction
does not account for it.

**NOT a uniform scale miss.** The whole excess sits in four fits:

```
                        k_hat    share of sum(z^2) removed
all 59                  1.189
drop 1 largest |z|      1.142    9.3%
drop 2                  1.095    18.0%
drop 3                  1.046    26.6%
drop 4                  1.003    33.6%
```

and the robust scale of the body is BELOW 1 (IQR-based 0.885, MAD-based 0.931).
z passes normality (Shapiro p = 0.48, KS vs N(0,1) p = 0.18) but those tests have
little power at n = 59, so they do not separate "mild uniform inflation" from
"calibrated body plus a few large errors" on their own. The conditioning below
does.

**It is `sigma_log_r` collapsing toward its boundary.** `sigma_logr` is recorded
per fit in the S = 8 and S = 36 arms (n = 39; the S = 18 runner predates that
column). Truth 0.5.

```
spearman(|z|,             sigma_log_r) = -0.470   p = 0.0028
spearman(|mu_hat - truth|, sigma_log_r) = -0.353   p = 0.0279
spearman(SE,              sigma_log_r) = +0.373   p = 0.0198

                     n    k_hat   mean |err|   mean SE   |err|/SE   coverage
sigma_log_r < 0.30   5    2.117     0.2287      0.1144     2.00      2/5  = 0.400
sigma_log_r >= 0.30  34   0.884     0.1055      0.1590     0.66     33/34 = 0.971
```

Split is stable for any threshold in 0.20-0.35 (the same five fits). Conditional
on the variance component being recovered the Wald interval is calibrated and
slightly conservative. The two failure modes fire together -- the point estimate
is 2.2x worse AND the interval 28% narrower -- which is why the pooled statistics
read as a clean multiplicative miss. It is not an SE artefact of a smaller sigma:
the raw error carries its own significant correlation.

This is the `sigma_omega` boundary collapse already documented at
`R/nmix_laplace_re.R`, on the sibling block that comment deliberately leaves at
pure ML. The S = 8 arm's hard failure (seed 502, `singular marginal Hessian or
non-finite optimum`) is the end state of the same collapse; that arm's
`sigma_log_r` minimum is 0.012. `theta_cov` marginalising over the variance
components does not protect against it -- marginalising accounts for uncertainty
in `sigma` GIVEN the data, not for a marginal Hessian evaluated where `sigma_hat`
sits against the boundary.

**Why no gate ships yet.** Detecting the collapse from the fit needs the variance
component's own uncertainty, `log(sigma)` against 0. `tulpa_re_aghq()` computes
the joint inverse `V` over `c(theta, log-Cholesky Sigma)` and returns only its
top-left block, so `SE(log sigma_log_r)` is discarded and neither `V` nor
`opt$hessian` is on the return to rebuild it from. Filed as gcol33/tulpa#418.
Without it the only available gate is an absolute cut on `sigma_hat`, which is a
number with nothing behind it and does not transfer off this fixture. The
measured conditioning is documented on `?ms_abun` instead.

**Threshold power, unchanged from #234's note.** `test-ms-abun-nb-rs.R:94` asserts
`mean(covered) > 0.8` on 20 seeds; at genuinely nominal coverage that still fails
about 1 run in 60, and at this seed count it cannot separate 0.90 from 0.95. The
threshold has no measurement behind it (written in 79e5eb7 with none recorded)
and is deliberately left untouched here: what it should assert depends on whether
the collapse is fixed or only reported, which is a decision about what the
package claims.

Probes: `dev_notes/nbrs_bisect/` (`Ssweep_S{8,36}_*.csv`, `fix3b*_nq3_*.csv`),
runners `dev_notes/_nbrs_S_probe.R`, analysis `dev_notes/_nbrs_nuts_analyse.R`.
A NUTS reference posterior on the S = 8 fixture was started and abandoned at 2 of
20 seeds (86 min per fit, ~5x per-seed spread, projecting 4-6 h);
`NUTSref_S8_*.csv` holds the two, width ratios 1.023 and 1.197.

### PC prior on the log_r block: measured, and NOT made the default (#235)

`logr_sigma_prior = c(1, 0.05)` (the omega block's own default), 20 seeds at
S = 8 / N = 100 / J = 5, `n_quad = 3` / `n_quad_scalar = 3`, paired seed for seed
against the pure-ML arm. 19 pairs -- seed 502 has no ML row, see below.

```
                       n    k_hat   bias      mean sigma   coverage
pure ML               19   1.182   +0.0681     0.432       17/19 = 0.895
PC prior              19   1.126   +0.0649     0.427       17/19 = 0.895

collapsed under ML     4   2.124   +0.1767     0.117 -> 0.208   2/4 both arms
healthy  under ML     15   0.752 -> 0.780      0.516 -> 0.485  15/15 both arms
```

**Coverage does not move: 17/19 either way, and the same two seeds miss.** The
prior is graded purely by proximity to zero, which is what a PC prior does and
also why it is mistargeted here:

```
seed  sigma_ml -> sigma_pc     z_ml -> z_pc        covered
502   FAILED   -> 0.385        (none) -> +0.784    rescued: ML hit a singular
                                                   marginal Hessian and errored
518   0.0116   -> 0.180        -1.321 -> -0.924    covered both
515   0.0612   -> 0.175        +2.684 -> +2.354    MISS -> MISS
513   0.1981   -> 0.238        +1.154 -> +1.141    covered both
507   0.1982   -> 0.238        +2.788 -> +2.695    MISS -> MISS
```

It bites hard at sigma ~ 0.01-0.06 (an order of magnitude lift) and barely at
sigma ~ 0.2 -- and sigma ~ 0.2 is where the misses live. Both misses improve and
neither crosses back.

**What it does buy.** Seed 502 under pure ML does not return a bad interval, it
`stop()`s with `singular marginal Hessian or non-finite optimum`; under the prior
it converges and covers. That is the failure mode the curvature at the boundary
is for, and it is worth having on its own -- but it is a robustness gain, not a
calibration one.

**What it costs.** The healthy 15 are not free: `k_hat` 0.752 -> 0.780, SEs a
median 4.2% narrower, sigma a median 0.031 lower, and a paired shift in
`mu_log_r` of -0.0061, 95% CI [-0.0092, -0.0029], p = 0.001. Small, but
systematic and in the wrong direction on the subgroup that was already
calibrated.

**Therefore the default stays `NULL` (pure ML).** Coverage unchanged, a
significant if small bias added to 15 of 19 fits, and the benefit confined to the
boundary cases. Turning it on by default would trade a documented failure on a
few fits for an undocumented shift on all of them. It ships opt-in
(`control$logr.sigma.prior`) and is worth recommending where a negbin community
fit at few species reports a `sigma_log_r` near zero or fails outright.

Deliberately NOT done: tuning `U` / `alpha` until coverage improves. The only
data available to tune against is the same 20 seeds used to judge the result, so
a strength chosen that way is fitted to the test rather than measured.

Probe: `dev_notes/nbrs_bisect/PCprior_S8_s*.csv`, runner `_nbrs_S_probe.R` (10th
argument, "U;alpha"), comparison `_nbrs_pcprior_compare.R`.

## Scalar nuisance AGHQ node floor (`ms_abun(mixture="negbin")` / ZI, #234)

Rule is in `CLAUDE.md` under "Community / multispecies N-mixture". Measurements:

- **The collapse.** A 2-node rule on `mu_log_r` gave SE 0.0110 against 0.1886 at
  `n_quad >= 3` (seed 515) -- a **17x** collapse, on a fit that converged with an
  ordinary point estimate and no warning.
- **Converged at 3, not clamped.** `n_quad` 3 / 5 / 9 agree to 10 decimal places.
- **Cost.** 1.5x on that axis at `n_quad = 5`.
- **Bias, not only the one seed.** Mean bias on `mu_log_r` +0.0482 -> +0.0109 over
  the 19 seeds that never collapsed.
- **Not sufficient for the gate.** `test-ms-abun-nb-rs.R:94` still reads 16/20 =
  0.800 against `> 0.8`; see the #235 note below.

## `mu_log_r` Wald calibration conditions on `sigma_log_r` (#235)

Rule is in `CLAUDE.md` under "Community / multispecies N-mixture". Measurements
(n = 59 fits pooled; conditioning subset n = 39, S = 8 + 36, truth `sigma` 0.5):

- **Not a uniform scale miss.** 4 of 59 fits carry 34% of `sum(z^2)`. Dropping them
  takes `k_hat` 1.189 -> **1.003**, and the body's robust scale is below 1
  (IQR 0.885, MAD 0.931). The filed "~24% too small" does not hold.
- **Conditioned on the recovered variance component.** `sigma_log_r >= 0.30` ->
  **33/34 = 0.971** covered, k 0.88. `sigma_log_r < 0.30` -> **2/5**, k 2.12, mean
  |err| 2.2x worse and SE 28% narrower -- both fire together, which is what makes
  the pooled statistics read as a clean scale miss.
- **Threshold-free.** spearman(|z|, sigma) -0.470, p = 0.003; spearman(|err|, sigma)
  -0.353, p = 0.028 -- so the effect is in the ERROR too, not an SE artefact alone.
- **Refuted en route.** Finite-sample-in-S: k peaks at S = 18 and the bias does not
  decay. t(S-1) df correction: coverage 0.864 -> 0.881, still p = 0.027 vs nominal.
- **Why no gate ships.** A detector needs `SE(log sigma)`, and `tulpa_re_aghq()`
  discards the joint inverse's variance-component block (`gcol33/tulpa#418`).

## Grouped-RE AGHQ debias, occupancy arms (`R/re_aghq.R`)

Rule is in `CLAUDE.md` under "RE both engines". Measured at n = 8 groups:

- Occupancy arm: sigma bias ~18% (raw EM variance component) -> ~4% (AGHQ).
- Detection arm: ~70% attenuation -> ~1%, with 88-96% interval coverage.

## Community variance under `pg_gibbs` vs the Laplace-EM

Rule is in `CLAUDE.md`'s "What works" rows: the Laplace-EM community SD is a
documented lower bound, the PG-Gibbs posterior is calibrated. SD reported as the
posterior MEDIAN (robust to the variance skew). Measured:

| family | arm SDs, truth | Gibbs | Laplace-EM |
|---|---|---|---|
| `ms_occu` | 0.6 / 0.4 (`sd_psi` / `sd_p`) | 0.667 / 0.417 | 0.595 / 0.342 |
| `jsdm` / `ms_count("bernoulli")` | 0.7 / 0.5 (`sd_mu`) | 0.674 / 0.439 | 0.635 / 0.401 |
| `ms_int_occu` | 0.5 (`sd_psi`) | ~0.53 | attenuates |

## AR1 year effect `rho` at few seasons (`t_occu()`, #124)

Rule is in `CLAUDE.md`'s `t_occu` row: `rho` is weakly identified at few seasons
and is reported, not asserted tightly. Even GIVEN the true year effects, `rho_hat`
for truth 0.6 climbs ~0.08 (T=8) -> 0.43 (T=20) -> 0.59 (T=200) -- a property of
the short AR1 series, not of the sampler.

## Per-species Hessian for `ms_distance` (#161)

Rule is in `CLAUDE.md`'s "Community distance spatial-factor" row. Before
`.tobs_ms_distance_info_block()`, the community EM finite-differenced the
per-species Hessian at `2U` sweeps per species per Newton step, each sweep summing
over the latent N: `test-ms-distance.R` ran **5.12 h** at tier 3 = **88%** of
`full-recovery.yaml`'s 350-min cap on its own. Assembling the block from what
`cpp_distance_site_sweep` already returns made the community EM **2.3x** faster
with the fit UNCHANGED to ~1e-13.

## Structured terms missing from the criteria (#211 / #215)

Rule is in `CLAUDE.md` under "`occu_cover()` detail": a term the scorer cannot see
is scored at ZERO. How far the scores moved once each term reached the kernels:

- Sampled (NUTS) route: elpd moved **31 nats** on a detection-arm RE and **549
  nats** on a sampled field.
- Grid-integrated (`nested_laplace`) route: `elpd_waic` moved **11.5 nats** on a
  detection-arm RE at `sigma` 1.1 and **61.5 nats** on a cover-arm RE.
- Occupancy-arm RE (#56, per site, via `model$re_psi`): `elpd_waic` **9.34**,
  `elpd_loo` **9.41**, at `sigma_re` 1.57. A fit carrying none stays byte-identical.

## Sampler-default divergences removed at #188

Rule is in `CLAUDE.md` under `engine_defaults.R`. Two shipped defaults were
proven artifacts rather than calibrated values, and were deleted:

- `n.iter` 2000 -> 1000 on the occu / cover / occu_cover NUTS paths. Commit
  8975470 changed `n.iter` from TOTAL to kept-post-warmup on exactly those paths
  and left the literal, so a default that had always kept 2000 - 1000 = 1000
  silently began keeping 2000. That commit's message records the behaviour change.
- occu `pg_gibbs` 2000/1000 -> the shared 3000/1500. It kept 1000 where every
  `ms_*` sibling kept 1500, with no recorded reason;
  `test-occu-pg-gibbs{,-spatial}.R` 29/29 green after.

Symptom that located the dead formals: `abun(method = "nuts")` ran
`.tobs_fit_model()`'s `n.iter = 2000` rather than its own `1000L` -- measured 2050
iterations, 2000 kept.

## #279 -- areal field simulators drew from a LAPACK eigenbasis

**Degeneracy of the fixture graphs.** Repeated eigenvalues of the ICAR precision
(gap < 1e-9 relative): 8x8 grid **31 of 64**, 9x9 grid **39 of 81**.

**The realisation was basis-dependent.** Rotating inside the repeated blocks
gives an exact eigendecomposition of the same `Q` (residual 5.1e-15, against
5.3e-15 for the original). Same seed, two bases:

| | eigen draw | Cholesky draw |
|---|---|---|
| `cor(f_basisA, f_basisB)` | **-0.0137** | **+1.0000** |
| via `MASS::mvrnorm` (nmix route) | +0.1363 | -- |
| `sd` each side | 0.9352 / 0.9352 | -- |

**The Cholesky replacement is the same distribution** (9x9 graph, 200k draws):
empirical covariance vs the centred `Q^+` off by 0.0045 (eigen) / 0.0063 (chol),
Monte Carlo error; geometric-mean marginal SD target 0.791238, eigen 0.791519,
chol 0.791398.

**It was not the fit.** `ms_occu_cover` BYM2 fixture, Windows: joint Hessian min
eigenvalue +1.20 against max 2.6e4 over all 15 EM iterations with `chol()`
succeeding at eps = 0 every time; the phi profile unimodal at every M-step with
`optimize()` matching a 960-point grid argmax to 4 dp; the EOF warm start's
leading singular value clearing the second by 13% and unmoved (`|cor|` 1.000000,
min over 200) by perturbations at 1e-15; the inner `optim(BFGS)` returning
convergence 0 at 127 gradient evaluations against its 400 cap. Perturbing the
objective and gradient by a relative 1e-15 -- the size of a BLAS reordering --
moves `phi_w` from 0.756105 to 0.758952 / 0.758586, a 0.003 move. The nmix BYM2
fit does not move at all: `scale_factor` 0.604775456220684, and perturbing it by
1e-15, 1e-12, 1e-9 and 1e-6 leaves the slope error at 0.00605, the winning cell
at (sigma 0.866, rho 0.3) and its weight at 0.2379.

### Seed counts behind the two multi-seed assertions

Both thresholds were tuned to one realisation of an arbitrary eigenbasis, and
the basis pinning is what hid how seed-fragile they were.

**`test-occu-cover-pos-field.R`**, cover-arm field SD, truth 0.6, 40 seeds:
median / mean 0.5480, sd 0.0932, min 0.3315, max 1.0324, median bias **-8.7%**;
quantiles (0, .05, .25, .5, .75, .95, 1) = 0.332, 0.441, 0.546, 0.548, 0.550,
0.598, 1.032. **2 of 40** outside the asserted band (0.4, 0.85). Quartiles span
0.004, so the median over 6 seeds is stable; `min`/`max` over 6 held with
probability 0.95^6 ~ 0.735, i.e. a ~26% failure rate on any platform. Band kept,
scored on the median; per-seed guard (0.1, 2.0) sits clear of the measured
extremes.

**`test-ms-occu-cover-spatial.R`** K = 2, `cor(F_hat, F_true)`, 12 seeds:

```
2024 0.6481 | 77 0.5365 | 1 0.5010 | 2 0.5532 | 3 0.6244 | 11 0.6348
  42 0.7110 | 101 0.7312 | 202 0.7890 | 303 0.8107 | 404 0.8432 | 505 0.5996
mean 0.6652  median 0.6414  sd 0.1120  min 0.5010  max 0.8432
```

**4 of 12 below the asserted 0.6** -- the threshold sits inside the
distribution, so asserting it on each of two seeds held with probability near
0.67^2 ~ 0.45. Scored on the MEAN, not the median: the spread is symmetric with
no outliers (0.501 to 0.843), which makes the mean's standard error the tighter
statistic -- 0.112/sqrt(12) = 0.032, clearing 0.6 by ~2 SE, where the median's
1.253 sd/sqrt(n) = 0.040 clears it by only ~1. (The pos-field case is the
opposite shape -- a tight core with rare escapees -- and takes the median.)
12 seeds where the block ran 2, so the block costs ~6x its old wall time.

### #279 follow-on: nested detection-RE BLUP recovery (`test-occu-cover-obs-re.R:388`)

Surfaced by the Cholesky field draw, NOT caused by it. The two draws have the
same distribution (verified directly, `dev_notes/probe_279_cov_equiv.R`): at
40,000 draws on the 9x9 fixture both empirical covariances converge to
`Q^+ / scale`, max relative error 0.0135 (Cholesky) vs 0.0187 (eigen), and the
two differ from each other (max abs 0.0346) by no more than each differs from
the target -- Monte-Carlo noise, not a distributional shift. Geometric-mean
marginal SD 1.0004 vs 1.0007.

`cor_s` = cor(BLUP, truth) for the nested `region:site` level, 20 seeds each
(`dev_notes/probe_279_obsre_seeds.R`, and the same with the pre-fix eigen draw
restored via `assignInNamespace`):

| draw            | seeds 1:5 mean | 20-seed mean | median | range       |
|-----------------|----------------|--------------|--------|-------------|
| pre-fix (eigen) | 0.6195         | 0.5424       | 0.558  | 0.170-0.840 |
| post-fix (chol) | 0.5176         | 0.5063       | 0.510  | 0.217-0.772 |

20-seed means differ by 0.036 against a difference-SE of ~0.054 -> inside noise,
as the covariance equality requires. The asserted floor of 0.6 was above the
estimator's mean under BOTH draws; pre-fix it cleared on a 5-seed mean by 0.0195,
about a quarter of that mean's own standard error (~0.076). The eigenbasis pinned
which draw each seed produced, so the luck was reproducible and read as stable.

Threshold moved 0.6 -> 0.40 on a 20-seed mean (SE ~0.038, so ~3 SE of headroom),
plus a per-seed `max(cor_s) > 0.5` gross-regression guard. The parent level is
untouched at 0.55 (20-seed mean 0.786 post-fix, 0.841 pre-fix).

The parent/nested gap (~0.79 vs ~0.52) is information content, not a fitting
defect: at `K = 3` sites per region a nested group carries ~4.5 sites of binary
detections and its BLUP is shrunk hard. Estimator left alone.

## `occu_multiscale_cover()` NUTS coverage, and why the assertion was pooled

Behind `test-occu-multiscale-cover-nuts.R` (3). Fixture as committed: `n_cells =
70`, `plots_per_cell = 4`, `visits_per_plot = 4`, `positive = "beta"`, `phi = 12`,
**`sigma = 0`** (no areal field in the truth), `seed = 600 + s`; NUTS `n.iter =
600`, `n.warmup = 600`, `n.chains = 2`, `seed = 100 + s`. 40 seeds, 0 errors,
0 divergences, 4.7 s per seed (one NUTS fit plus one Laplace fit), tulpaObs
47b728f / tulpa 0.1.20.

| coef | truth | bias | mean se | sd(est) | sd/se | coverage |
|---|---|---|---|---|---|---|
| `psi_x_cell` | 0.600 | -0.044 | 0.297 | 0.331 | 1.11 | **0.850** |
| `psi_(Intercept)` | 0.200 | +0.085 | 0.277 | 0.267 | 0.97 | 0.925 |
| `pos_(Intercept)` | -0.847 | +0.003 | 0.042 | 0.040 | 0.96 | 0.950 |
| `pos_x_cov` | -0.300 | +0.003 | 0.042 | 0.042 | 1.01 | 0.950 |
| `theta_(Intercept)` | 0.500 | +0.017 | 0.211 | 0.236 | 1.12 | 0.950 |
| `p_(Intercept)` | 0.300 | -0.033 | 0.119 | 0.104 | 0.88 | 0.975 |
| `p_x_pdet` | -0.400 | +0.009 | 0.121 | 0.097 | 0.80 | 0.975 |
| `theta_x_plot` | 0.400 | +0.029 | 0.199 | 0.204 | 1.03 | 0.975 |

Pooled 302/320 = **0.944**, exact 95% CI [0.913, 0.966], which brackets nominal.
Every `sd(est)/mean(se)` lies in [0.80, 1.12], so no arm carries an SE miss.

**The old `min(cover) >= 0.80` at 8 seeds was a coin flip.** A per-coefficient
rate over 8 seeds is a small binomial and the minimum over eight of them
compounds; at the coverages above

```
P(pass)  n= 8 seeds  0.464     <- the committed form
         n=20        0.812
         n=40        0.862
         n=60        0.893
```

so raising seeds does not fix that form. It failed when 8c67903 moved the
simulator's RNG stream by one normal (`rnorm(sum(keep))`, N-1 = 69, to
`rnorm(N)` = 70): the field is discarded at `sigma = 0` but the stream position
is not, so the fixture landed on a different draw and the coin came up the other
way. Not a regression; the "42 files re-draw their data" cost of #279.

Replaced by pooled coverage `>= 0.85` (false-failure rate ~0 at every size tried)
plus a per-coefficient floor of 0.60 as a gross-regression guard (false-failure
rate 3e-5 at 40 seeds), and `n_seeds` 8 -> 40. Margins on the other assertions in
the block at 40 seeds: `max|bias|` 0.085 against 0.30, `max|NUTS - Laplace|`
0.042 against 0.15, divergences 0 against 5.

`psi_x_cell` at 0.850 (34/40, CI [0.702, 0.943]) is genuinely under nominal --
0.95 is excluded, 0.80 is not. It is the cell-level occupancy slope competing
with the one-node-per-cell field the model fits whether or not the truth carries
one, the same confounding spPGOcc documents. Reported, not tuned away.

## Divergence rate on the spatial-factor K = 2 NUTS path (`ms_occu_cover`)

`test-ms-occu-cover-spatial.R`, "NUTS recovers the community means under the
constrained K = 2 path", asserts one fit's divergence fraction below 0.25. The
first full-recovery run that ever reached the block returned 0.271 on the Linux
runner. The same seeds on Windows return 0.016, so the assertion was not
reproducible locally and the number could not be read as a regression without
measuring what the quantity does.

Fixture: 6x6 grid ICAR, S = 8, J = 4, K = 2, `sd_occ` 0.5, `sd_load` 1.1,
`sigma_pos` 0.4; sampler 2 chains x 500 kept after 300 warmup. 8 seeds per arm,
Windows, tulpaObs 0.0.237 / tulpa 0.1.20, about 70 s a fit.

**The rate is a property of the data draw, not the sampler seed.**

| varied | n | mean | sd | median | range | over 0.25 |
|---|---|---|---|---|---|---|
| sampler seed (data seed 4) | 8 | 0.027 | 0.018 | 0.019 | 0.009-0.049 | 0 |
| simulation seed (sampler seed 7) | 8 | 0.112 | 0.167 | 0.058 | 0.018-0.515 | 1 |

Holding the data and moving the sampler moves the fraction by a factor of 5;
holding the sampler and moving the data moves it by a factor of 29, with a right
tail that crosses 0.25 on its own. A single fit against a fixed band is therefore
scoring the draw, and the platform shifts it by about as much as a change of
data seed does.

**The step size is the cause, and raising it is a correctness fix rather than a
threshold move.** Same fixture, `adapt.delta` swept on the test's own data seed
and the two worst draws of the sweep:

| data seed | 0.95 | 0.99 | 0.999 |
|---|---|---|---|
| 4 (the test's) | 0.016 | 0.018 | 0.012 |
| 104 | **0.515** | 0.043 | 0.007 |
| 105 | 0.134 | 0.032 | 0.003 |

The community-mean correlation rises with it on the draws that were diverging --
seed 104 goes 0.915 -> 0.970 and seed 105 0.962 -> 0.968 -- so at 0.95 the
divergences were biasing the posterior, not just counting against a band. Cost
on the test's own seed is 75 s -> 83 s at 0.99 and 137 s at 0.999, which is what
puts the block at 0.99 rather than 0.999: 0.99 already clears the band by a
factor of six on the worst draw measured.

The 0.25 band is untouched. Whether 0.99 also clears the Linux 0.271 is a
prediction from the mechanism and from the 12x reduction it produces on the
Windows draws in the same regime, not a measurement -- the next full-recovery
run is what tests it.

## Two occu_cover fixtures after the Cholesky field draw (#279)

`8c67903` draws every areal field through a Cholesky instead of a LAPACK
eigenbasis, so a fixed seed no longer produces the realisation it used to. The
estimators are untouched; two fixtures that assert on one realisation landed on
the wrong side of their bands with the new draws. Neither band moved.

### `test-occu-cover-nl-ic.R`: the fixture's true RE, 0.9 -> 2.0

Every criteria claim in the file rests on the fit having found a substantial
random effect. At `sigma_re` truth 0.9 the grid returns 0.30 to 1.57 over 10
seeds (arm p: mean 0.895, median 0.900; arm pos: mean 0.986, median 1.022), so
the premise fails on its own on roughly a fifth of draws, and seed 1 came back
0.301 against the 0.4 the block asks for.

A fit whose RE came out near zero is one WAIC is right to penalise. On seeds 1
and 6 the offsets roughly triple the effective parameter count for almost no
gain in lppd -- `p_waic` 25.7 against 9.9, and 28.9 against 10.9 -- so elpd
falls. On seeds where the RE is substantial `p_waic` DROPS when the offsets are
carried (15.4 against 19.1) and elpd rises 12 to 30 nats. The direction of the
criteria claim is conditional on the fitted variance, not on the draw count:
raising the posterior draws 400 -> 1500 -> 4000 leaves the negative gains
negative to within a nat.

At truth 2.0, over 6 seeds on each arm, the worst margin is 2.8x the band:

| arm | min sigma_re (band 0.4) | min d_lppd (5) | min d_elpd (5) | min WAIC gain (3) | min LOO gain (3) |
|---|---|---|---|---|---|
| p | 1.12 | 13.35 | 14.39 | 14.57 | 13.67 |
| pos | 1.23 | 12.83 | 15.56 | 15.78 | 15.35 |

1.5 also passes everywhere but leaves 4.04 against a band of 3 on one draw.

The draw count in the first block goes 400 -> 1500 for a separate reason: it
compares a per-group mean OVER DRAWS against the BLUPs at a tolerance of 0.05,
and that comparison's Monte Carlo error is 0.068 at 400 draws and 0.026 at 1500.
That one IS draw-limited, and it is the only assertion in the file that is.

`zero offset == no offset`, the structural invariant the file exists for, held
on 20 of 20 fits throughout, before and after.

### `test-occu-cover-predict.R`: 40 cells and 6 visits -> 120 and 10

`cor(delta_p, delta_true)` tracks the trend surface it is read off almost
exactly, and at 40 chain cells that surface is not identified: over 16 seeds its
correlation with truth runs 0.122 to 0.856 and the predicted change clears its
0.5 band on 11 draws of 16. Averaging does not rescue it -- bootstrapped
P(mean > 0.5) is 0.82 at 12 seeds and 0.85 at 16.

Widening the outer grids was tried first and rejected. The fixture integrates
over `sigma.grid = c(0.5, 1.0)` and `alpha.grid = c(0, 0.5)` while its truth is
sigma 0.8 and alpha 0.6, so the copy amplitude cannot reach the value the data
were drawn at; bracketing both moved the correlation from 0.556 to 0.584 and
cost the interval coverage in the same block, 16 of 16 draws over 0.8 falling to
11 of 16 with a minimum of 0.175.

More data is what identifies it:

| fixture | cor mean | cor min | over 0.5 | trend-field cor min | coverage over 0.8 |
|---|---|---|---|---|---|
| N=40 J=6 | 0.556 | 0.182 | 11/16 | 0.122 | 16/16 |
| N=80 J=6 | 0.661 | -0.027 | 12/16 | 0.114 | 16/16 |
| N=80 J=10 | 0.635 | -0.171 | 13/16 | -0.136 | 16/16 |
| N=120 J=10 | 0.851 | **0.532** | **16/16** | 0.527 | 16/16 |

The file runs in 6 s at the larger fixture.
