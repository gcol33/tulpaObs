# Plan: phi on the outer joint hyperparameter grid (tulpaObs#7 candidate)

## Why

The joint cover-hurdle beta path treats `phi` as a profiled scalar: it is
estimated once on the positive subset (no spatial) by `.prefit_beta_phi()`
and held fixed while `tulpa_nested_laplace_joint()` integrates the spatial
hyperparameters. After the joint fit, `.refit_beta_phi_postfield()` does
a 1D Brent refit of `phi` against the posterior-weighted linear predictor,
correcting the field-leakage that motivated issue #5.

This is the right fix at moderate `n_pos` (D7 cell B median 320: bias
~−2%). At thin `n_pos` (D7 cell B sparse cell, median 46) the same path
leaves a residual −40% bias on `phi_pos`. The cause is **not** small-
sample MLE bias of the 1D Brent step — controlled probes prove otherwise
(see *Probe evidence* below). It is upstream posterior shrinkage of
`(sigma, rho, alpha, w_s)` toward the prior at small `n_pos`. The post-
hoc refit cannot self-correct because it is conditional on a `mu` whose
variance has already been compressed.

The principled fix is to put `phi` on the outer joint hyperparameter grid
alongside `(sigma, rho, alpha)` (resp. `(tau, alpha)` for ICAR,
`(tau, rho_car, alpha)` for CAR_proper) and let the marginal likelihood
weight it. With `phi` part of the outer integration, mean-variance
balance is enforced by the data, not by a plug-in.

## Probe evidence (already in tree)

* `tulpaObs/dev_notes/probe_beta_phi_small_sample.R` — known-mu 1D MLE
  bias at `phi=30`: `+9%` at `n=25`, `+4%` at `n=50`, `+1%` at `n=100`.
  Slightly upward, not downward. Bootstrap correction kills it.
* `INLAabun/example/validation/.probe_d7_phi_origin.R` — at D7 cell B
  config (median `n_pos = 46`):
  * posterior-mu refit: median `phi = 17.7` (−42%)
  * truth-mu refit:     median `phi = 30.1` (+7.8%)
  * bootstrap-corrected posterior-mu refit: median `phi = 17.1` (−44%)
  * posterior linpred variance: 0.30 (truth 0.64)
  * `alpha_hat` median 0.82 (truth 1.0)

## Design

### R-side (tulpaObs, ~1 day)

In `family_cover_hurdle.R::fit_cover_hurdle_joint_nested`:

1. Drop `.prefit_beta_phi` as a hard plug-in. Replace with a sensible
   default grid in log-space, e.g.
   `phi_grid = exp(seq(log(2), log(300), length.out = 5))`. Honor
   `control$phi_grid` user override.
2. Outer-loop over `phi_grid`: for each `phi_k`, call
   `tulpa_nested_laplace_joint()` with that value passed to the arm spec.
3. Collect per-`phi_k` log marginal likelihoods (the kernel already
   returns one per outer grid point); combine across `phi_k` × spatial
   grid into a single weighted-mean / weighted-mode output.
4. Compute posterior summaries of `phi`: weighted mean, weighted sd,
   marginal posterior over the grid for plot/sanity.
5. Drop `.refit_beta_phi_postfield()`. Keep it as a helper for users who
   want a post-hoc field-corrected refit, but it is no longer called
   from the joint path.

### tulpa-side (kernels, ~½ day)

No C++ changes required for option **A** (R outer-loop calls the
existing kernel `n_phi` times, ~5 calls). The kernel already accepts
`phi` per arm in `arms_list[k]$phi` and integrates the rest. Cost:
`n_phi × kernel_time`. At `n_phi = 5` and current kernel ~3 s on D7
cell B, that's ~15 s per fit — acceptable.

Option **B** (push `phi` inside the kernel as another axis) saves
maybe 30% by amortizing the per-grid setup, but adds an axis to the
C++ inner loop. tulpaObs is a production package, not a POC, so the
C++ work is justified *if* profiling at realistic problem size shows
the R outer loop is the bottleneck. Decision rule: ship option A
first, profile on the recovery suite + a panel-scale benchmark
(`n_pos > 10^3`), and promote to option B if the outer loop accounts
for >30% of fit time at production scale.

### Verification (~½ day)

* Recovery test at `n_pos ~ 46` (D7 cell B sparse): require
  `|phi_hat - 30| / 30 < 0.20` across 10 seeds.
* Recovery test at `n_pos ~ 320` (D7 cell B dense): no regression.
* Coverage test on a `phi` SE / 95% interval derived from the marginal
  weights: ≥0.85 across 30 seeds at `n_pos ~ 100`.
* D7 sparse re-run on the INLAabun side: confirm bias collapses to
  within FCN expectation (+/- 10% at `n_pos ~ 46`).

## Scope of change

* `tulpaObs/R/family_cover_hurdle.R` — 30–60 lines edited / added.
* `tulpaObs/tests/testthat/test-cover-hurdle-nested-joint-recovery.R`
  — un-bracket the small-`n_pos` recovery, add coverage test.
* `tulpaObs/R/predict.cover.R` (or equivalent) — propagate `phi`
  posterior into predictive draws (currently uses the point estimate).
* `INLAabun/example/validation/d7_sparse_positive.R` — no change; the
  rerun just consumes the new output.
* No tulpa-package change for option A.

Total: ~1–2 days end-to-end including the probe rebuild and the D7
re-run.

## Open questions

* Default `phi_grid` placement. The default
  `exp(seq(log(2), log(300), 5))` covers
  `phi ∈ {2, 6.2, 19, 58, 300}` — wide enough to bracket realistic cover
  dispersions. If the posterior weight piles on the boundary, the
  control message should ask the user to widen.
* Prior on `phi`. Right now there is no explicit prior (it's a profiled
  point estimate). Putting it on a grid implicitly imposes a uniform-
  prior weighting on the grid points. A Gamma(2, 0.1) prior on `phi`
  (mean 20, sd ~14) would be more honest and is what most beta-
  regression Bayesian packages use. Worth checking what `INLA::inla()`
  defaults to for the equivalent setting.
* Whether ICAR / CAR_proper backends need the same treatment. They use
  `tau` not `(sigma, rho)`, but the structural argument is identical:
  the spatial precision is integrated, the dispersion `phi` is plugged.
  Same fix applies, same scope of change.

## Not in scope

* Re-deriving the closed-form bias correction for the 1D Brent step
  (Cox-Snell `(E[U''] - 2 E[U^3]) / (6 I^2)`). Already shown to be
  empirically irrelevant at the regime of interest.
* Parametric bootstrap bias correction on the post-hoc refit. Same
  reason.
* Pushing `phi` into the C++ kernel as an internal axis (option B
  above) in *this* plan. Deferred to a follow-up once option A is
  shipped and profiled — see decision rule under *tulpa-side* design.
