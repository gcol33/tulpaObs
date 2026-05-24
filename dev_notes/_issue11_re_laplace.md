# Deterministic random effects under Laplace (gcol33/tulpaObs#11)

**Status: RESOLVED.** Formula random effects are fit by the default
`engine = "laplace"` (iid intercept + uncorrelated slopes) and by NUTS (all
forms), with named parameters, structured `ranef()`, and `coef()` carrying the
visit-level detection coefficients. Forms the deterministic path cannot fit
error toward NUTS instead of being silently dropped.

## What was wrong

`.tobs_fit_model()` computed `structs$re` but only the NUTS branch consumed it;
`method == "laplace"` called `.tobs_laplace()` without `re`, so the default
engine silently ignored random effects. `nested_laplace` fit intercept RE but
rejected slopes with a message that pointed back at the (also-dropping) Laplace
path.

## Engine constraints (verified in the tulpa source / by probe)

* `tulpa_laplace()` finds the latent mode at a **fixed** sigma and uses a
  **diagonal** RE precision (`fit_laplace.R`): no variance-component estimation,
  no correlated (Cholesky) covariance. Correlated slopes live only in the NUTS
  sampler.
* `tulpa_laplace()` **rejects a fractional binomial response** (`y in (0,1)` with
  `n = 1` gives garbage betas — see `probe_laplace_weights.R`). The direct
  two-row weighted encoding (`(y=1, w)`, `(y=0, 1-w)`) is exact in theory but the
  RE-aware Newton solve does **not converge** (`probe_two_row_debug` showed beta
  creeping 0.007 -> 0.03 over 500 iters).
* Observation `weights` ARE applied correctly to the fixed-effect block
  (`probe_laplace_weights.R` case A matches `glm(weights=)` to 7e-05).

## Method (R/em_laplace_re.R)

A variance-component EM fused with the occupancy missing-data EM:

1. **E-step** — `psi_i = logit^{-1}(X_i beta + Z_i b)` (RE mode fed back in);
   `w_i = P(z_i = 1 | y_i)`.
2. **M-step (occupancy)** — pseudo-binomial M-trick (`y = round(w*M)`, `n = M`)
   with the RE prior **rescaled to `sigma/sqrt(M)`** so the penalty scales with
   the M-inflated data term and the penalized MAP of `(beta, b)` is unchanged
   (argmax is scale-invariant). Mode-find via `tulpa_laplace(re_list = ...)`.
3. **VC update** — `sigma^2_c <- mean_g(b_gc^2 + Var(b_gc | y))`, with
   `Var(b | y)` the diagonal of the RE block of the joint Hessian inverse,
   recomputed at **natural** (n=1) scale via Schur. PQL-style (complete-data
   working weights), so sigma carries the usual small-cluster downward bias.
4. **M-step (detection)** — weighted binomial (site level).

**SE caveat that bit us:** the M-step `H_beta` is M-inflated, so it is unusable
for the occupancy fixed-effect SE (~sqrt(M)=~31x too small). The SE is computed
separately at natural scale as the **Louis observed info**
(`psi(1-psi) - w(1-w)`) marginalised over the RE block via Schur
(`.tobs_re_occ_fixed_se`). A recovery-test assertion guards against the inflated
SE regressing back in.

## Naming + BLUPs (R/re_effects.R)

NUTS lays RE params out TYPE-BLOCKED after fixed+visit:
`[all log_sigma][all chol][all z]`, each section iterating terms; `z` is
group-major within a term. BLUPs are reconstructed per draw as
`b_{g,c} = sigma_c * (L z_g)_c` (`L = I` uncorrelated) and summarised. The
deterministic path builds the same `re_effects` shape directly.

## Validation

* `probe_vc_em.R` — end-to-end VC-EM on simulated occupancy (iid intercept +
  uncorrelated slopes): beta/p recovered, sigma shows the expected mild PQL
  bias, BLUPs correlate with truth.
* `tests/testthat/test-re-laplace-recovery.R` — deterministic recovery + a
  NUTS-ballpark comparison (sigma within ~60%, BLUP cor > 0.85) + safety errors.
* `tests/testthat/test-re-bar-recovery.R` (NUTS) still passes; the type-blocked
  BLUP reconstruction matches its manual recipe.

## Not done (deferred feature, not a bug)

Correlated random slopes under Laplace would need a Cholesky-factored covariance
in tulpa's Laplace engine (currently diagonal-only). Scoped to NUTS; the
deterministic path errors with a clear pointer rather than approximating.
