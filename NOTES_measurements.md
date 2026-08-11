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
