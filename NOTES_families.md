# NOTES_families.md

Deep implementation notes moved out of `CLAUDE.md` to keep it inside its
character budget. `CLAUDE.md` keeps the rules, contracts, gates and pointers;
this file keeps the mechanics behind them. Same convention as
`NOTES_measurements.md`: committed alongside, Rbuildignored.

## occu_cover() obs-arm RE: parse, design, blocks, postprocess (#102/#103)

**Per-group RE on detection / cover arm (#102 intercept, #103 crossed/nested)**:
random intercept(s) on `detection=`/`positive=` (`(1|g)`/`re(g)`), per-VISIT
grouping (one code per (site,visit)), composing w/ the required psi field on the
nested_laplace joint path. `.occu_cover_obs_re_parse` (occu_cover.R) strips ALL re
terms off the obs formula BEFORE copy-extraction + design build (rejects other
structured terms; copy/re allowed), returns `$terms` (LIST of specs,
crossed/nested) + `$has_slope`. `.occu_cover_obs_re_design` resolves each term's
per-(site,visit) codes site-major from `data` (site-level, broadcast) or `visits`
(visit-level), levels from observed visits ONLY, builds slope `Z` (intercept +
covariate cols), attaches per-term LIST `model$re_det`/`re_pos`. **Slopes (#103,
tulpa>=0.0.39)**: NO gate -- the DESCRIPTION floor enforces the engine. Per term:
intercept -> one scalar iid block; uncorr slope -> one weighted iid block per coef
(tulpa `svc_weight` = `Z[,c]`, intercept col all-ones = scalar iid); corr
(`(1+x|g)`) -> one `miid` block (tulpa#114: `Q=I` mcar, `n_fields`=n_coefs,
`field_weight`=Z cols, free Sigma log-Cholesky). **Slope covariate STANDARDIZED**
to unit SD in `.occu_cover_obs_re_design` (`coef_scales`, intercept scale 1) so the
fixed Sigma grid is scale-invariant; BLUP/sigma + the predict draws
back-transformed `/scale` to natural units (cor scale-free). miid grid =
`.occu_cover_miid_logchol_grid` (p=2 principled compact: SYMMETRIC rho incl 0 +
strong +/-, log-spaced SD) to stay under the engine's 2048 outer-grid cap, knob
`control$re.logchol.grid.p`/`.pos`. `re_descs` = ONE desc per TERM
(block_start/n_blocks span). **Key fix (#102)**: the det arm's `field_coef`=1 when
`model$re_det` present so the iid block scatters; field skipped by
`spatial_idx=0`. Postprocess: per term gather its blocks' latent cols (uncorr =
n_coefs iid blocks; corr = one miid, latent coef-major `(c-1)*ng+g`) ->
`[n_groups x n_coefs]` BLUP (centred per coef); sigma per coef from `b<P>.sigma`
(uncorr) or marginalized from `b<P>.L<ij>` log-Cholesky (corr) -> per-coef `sigma`
+ `cor` matrix on `re_terms`. Names `.occu_cover_re_sigma_names` (base
`sigma_re`/`_p`/`_pos`, `_<var>` when >1 term/arm) + `_<coef>` per coef +
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

## occu_cover() NUTS: sampled field hypers, criteria offsets, second field (#204/#211/#215/#214)

**Hypers SAMPLED, not pinned (#204).** A fit conditioned on the outer grid's own
point estimate cannot serve as an independent reference for that grid. Every areal
kind's loading factors as a FIXED basis w/ hyper-dependent column weights, so a
leapfrog step costs a rescale, never a re-decomposition (`src/nuts_field_hyper.h`, R
mirror `.ochf_*` in `occu_cover_nuts.R`):
`z = sigma*(B1 %*% (s1(rho)*raw1) + s2(rho)*raw2)`. icar `B1` = sum-to-zero
eigen-loading of intrinsic `Q` (#71), `s1=1`, no `raw2`; bym2 same `B1`,
`s1=sqrt(rho/sf)`, `raw2` iid w/ `s2=sqrt(1-rho)` (Riebler); car_proper
`Q(rho)=D-rho W` in the eigenbasis of the symmetrically normalised adjacency ->
`B1 = D^{-1/2}U` FIXED, `s1_j=(1-rho lambda_j)^{-1/2}` -- which is why car_proper rho
is NOT the O(n^3) per-step Cholesky the issue scoped it as. Each sampled hyper rides
`t = t_lo + (t_hi-t_lo)*plogis(u)`, `value = inv_link(t)` (log for sigma/alpha, logit
for rho) -- bounded, so no wall, Jacobian in the target.
**Prior = flat in `t` over the WARM FIT'S OWN outer-grid span** (`fit$theta_grid`
column range) = the measure the nested-Laplace grid integrates against, so
`control$sigma.grid`/`alpha.grid`/`rho.car.grid` move BOTH backends; flat prior +
change of variables = normalised `log(e)+log(1-e)`, `e=plogis(u)`. alpha's grid `0`
atom is not HMC-representable -> bounds take the positive nodes only. icar pins
`rho=1` (intrinsic precision has no mixing param); an axis the grid pinned to one node
stays pinned. `fit$nuts$sampled_hyper` / `$fixed_hyper` = CHARACTER vectors, empty
when nothing pinned (deliberately a type change from the old `fixed_hyper=TRUE`, so a
stale `isTRUE()` read fails loudly) + `$fixed_hyper_values`; `fit$hyper_draws` cols
`sigma|rho|alpha|field_sd`. **`field_sd`** = geo-mean marginal SD the block implies at
that draw = the ONLY field-scale summary comparable across kinds (the three normalise
their precisions differently), so it is what a simulation truth is stated in (the
simulator's `f` carries geo-mean marginal variance 1, Sorbye-Rue, truth = `sigma`).
`control$fixed.hyper=TRUE` restores #74 conditioning, byte-identical to the old
loading (pinned = the degenerate configuration of the same block, not a second path).
Validated C++ == R oracle, analytic == central FD on every hyper coordinate, and the
prior verified flat in `t` by holding `raw=0` (field vanishes -> the target IS the
prior; catches a missing Jacobian, which the gradient check cannot). `inv_metric` MUST
be sized `n_raw + n_hyper`: the engine takes a short metric pointer w/o a length check
and bym2's `n_raw = 2n-1`.
**Sampled terms reach the criteria (#211/#215).** A structured term the scorer cannot
see is scored at ZERO. `.tobs_occu_cover_components()` returns them as OFFSETS beside
the coefficient draws: per SITE `field_occ`/`field_pos` (#204's sampled field off
`fit$field_draws` + the per-draw `alpha`; the v3 route still reads `field_table`), per
VISIT `off_det`/`off_pos` (#205's sampled obs-arm RE, `re_draws` -> `sigma_re*z` mapped
through the view's `flat_idx`). All four diagnostics fold them in -- ploglik
(WAIC/LOO/CPO), PPC, PIT/LOO-PIT -- b/c the per-visit offset enters the SHARED `Arms`
view (`src/occu_cover_ragged.h`), so ONE change reaches all three kernels and dense ==
compact by construction (#185). A 0-column matrix = "arm carries none" -> null pointer
-> the no-offset path byte-identical. Cell-aggregated cover scores one cover row per
detected UNIT, so a per-visit offset errors there w/ a pointer. The grid-integrated
(`nested_laplace`) route reaches the same criteria: `.tobs_joint_draws()` returns its
RE latents on `bundle$re` in the layout the offset builder already reads. An
occupancy-arm RE (#56) is per SITE, so the dispatcher stores its group codes on
`model$re_psi` (the counterpart of `model$re_det`/`re_pos`) and the builder adds the
per-group draw to `field_occ`; `.occu_cover_spatial_fields()` also carries the term's
`var` + factor `levels` through, so `fit$re$psi`/`ranef()` label the grouping like the
obs arms and `predict(newdata=)` matches it. `test-occu-cover-nuts-ic.R`.
**A SECOND (SVC / trend) field samples too (#214).** The block is a LIST
(`hyper_field_build_list()`, `src/nuts_field_hyper.h`): each field carries its own
basis, site->node map, per-site design WEIGHT (`field_weight`; absent = the intercept
field) and its own sampled `(sigma, rho, alpha)` -- two fields share no hyper. Site i
loads `sum_b w_b(i) z_b[cell(i)]` on psi and `sum_b alpha_b w_b(i) z_b[cell(i)]` on
cover; the three places that loading is written are
`hyper_field_site_{value,offsets,score}()`, so the eta assembly and the score cannot
express it differently. Spec spelling `field_blocks` (list); a one-block fit is
`expect_identical` on lp+grad AND on a whole fit's means/sds/field/draws. Layout per
block `[raw, sampled hypers]`, blocks back to back, RE blocks after -- so `n.iter`-for-
`n.iter` the one-field vector is unchanged. **The warm fit is the multi-block coupled
path** (`multi = TRUE` arms + one `icar/bym2/car_proper` block per field w/
`svc_weight` + a copy spec per block, `alpha.grid` / `alpha.grid.trend`), and it FORCES
`integration = "grid"`: above 3 axes the engine switches to a mode-centred CCD star
whose column range is a design radius, not an integrated span, and the sampler reads
each axis's span as its flat prior's support. A DEFAULTED axis is thinned to 3 nodes
over the SAME span when a second field is present (the prior is defined by the span
alone -> unchanged; the tensor is a product over blocks). Reported:
`fit$trend_field`/`trend_fields` (named by weight column), `fit$trend_field_draws`,
per-block suffixed hypers (`sigma_trend`, `alpha_trend`, `field_sd_trend`, indexed when
several) in `hyper_draws` / `sampled_hyper` / `fixed_hyper`, and
`fit$spatial$field_suffix`/`field_weights` so `.tobs_occu_cover_sampled_field()` sums
every block's loading into the criteria. Both surfaces recover w/ 0 divergences;
numbers in `NOTES_measurements.md`. `test-occu-cover-nuts-svc.R`. Correlated `|` (one
free-Sigma MCAR block), temporal + RE still gated -> n-L. group_var maps sites>cells;
predict() needs the joint object (non-spatial laplace AND nuts both error w/ pointer);
sampled-field (estimated-variance) route = `ms_occu_cover()` factor (tulpa#67).
## joint_substrate.R: outer-grid placement + defaulted-grid marking (#187/#186)

    **Outer-grid placement promoted (#187)**: `.tobs_promote_outer_grid(jf)` lifts `outer_grid_placement` ("fixed"/"auto_recentered"), `_recenter_attempts`, `_prior_added`, `_recenter_declined` (reason a "fixed" placement stayed fixed, tulpa#293) to `tobs_fit` top level, spliced everywhere `.tobs_promote_pareto_k` is (occu_cover postprocess, occu_joint, occu_multiscale_cover_joint) + the cover decode. NOT gated on the grid having MOVED -- a declined recenter is exactly the case worth seeing; an inert one invisible across a whole batch is what filed it. `.tobs_glance_outer_grid(g, x)` adds the two columns in BOTH `glance.tobs_fit()` + `glance.tobs_multiarm_fit()`; latter terminal for `cover_fit` (class order cover_fit/tobs_multiarm_fit/tobs_fit) -> never reaches the former.
    **Defaulted grids declare themselves (#186, needs tulpa >= 0.0.132)**: engine auto-recenter decides axis PROVENANCE, not field presence -- moves an axis that is absent, `auto_grid()`-marked, or exactly equal to its own default; anything else = user pin. tulpaObs writes a grid on EVERY joint fit -> unmarked default reads as a pin, rescue goes inert. `.tobs_default_{sigma,alpha,bym2_rho}_grid()` now RETURN `tulpa::auto_grid(...)` (they are the layer that chose the values; a user grid never passes through them); `.tobs_mark_auto(x, auto)` re-applies the mark wherever a site reshapes a defaulted grid, since `sort()`/`[`/`c()`/`as.numeric()`/`expand.grid()` all drop the attribute -- cover arm-specific tau translation, occu_cover pos-arm tau, copy `alpha_grid`s, RE `sigma_grid`s, EM-path bym2/ar1 pairings. Verified end to end: defaulted axis reports `declined = "grid_not_collapsed"`, a `control$sigma.grid` one `"axis_pinned"`

## svc(): observation-family wiring over the shared areal-BFGS driver (#144)

  **Observation families (#144, `laplace`/`nested_laplace`)**: `removal()`,
  `distance()`, `fp_occu()`, `dyn_abun()` carry `svc()` too, with NO family-specific
  code -- those four already ride `.tobs_areal_bfgs_fit`, and an svc surface IS just
  another latent block on the arm their `eval(theta, offset)` already exposes, so the
  whole wiring is `.tobs_svc_field_blocks()` (single source of truth for the term's
  validation + hyper grid) + `.tobs_build_field_spec(svc=, X_svc=)` appending one
  NNGP block per `indices` entry AFTER the areal / temporal blocks +
  `.tobs_attach_field_results(svc=, has_spatial=)` slicing the trailing blocks into
  `fit$svc_field`/`svc_hyper`/`svc_field_arm`. Composes with an areal and/or temporal
  field on the same arm. Surfaces load on the STATE arm only (log lambda / psi); a
  detection-arm areal field alongside svc errors (`.tobs_check_svc_arm()` -- the
  driver exposes ONE `grad_eta`, so the surfaces would otherwise be fit against the
  detection arm). The N-mixture families (`abun`/`ms_abun`) do NOT get it: their areal
  path is the C++ count-spatial driver, not this one. NUTS still errors everywhere but
  single-season occu(). Surface cor in `NOTES_measurements.md`.
  `test-svc-families-recovery.R`. The driver also returns `res$eta_offset`, the
  marginalised per-observation offset the blocks jointly load; the family wrappers
  read it instead of re-deriving each block's site map, which also fixed the
  temporal-only fp_occu path (it indexed a length-`n_t` field by site via
  `res$field_mean[map]`).

## svc(): deterministic-backend field blocks and the Vecchia precision (#143)

  **Laplace backends (#143, `R/occu_svc.R`)**: `occu() + svc()` also fits under
  `method="laplace"` / `"nested_laplace"`. K surfaces = latent field blocks on the
  psi logit -> rides the SHARED areal-BFGS driver (`.tobs_areal_bfgs_fit`,
  `R/areal_bfgs.R`); two new pieces only: `.tobs_svc_nngp_field()` (continuous NNGP
  block w/ optional per-site design weight, continuous sibling of
  `.areal_field_car(weight=)`) + `.tobs_occu_svc_marginal()` (exact two-state
  occupancy marginal, Fisher-identity gradients `w-psi` / `w(y-p)`, FD-validated).
  Vecchia precision `Q=(I-A)'D^-1(I-A)` assembled in R (`.tobs_nngp_precision`) from
  the term's OWN neighbour structure w/ the compiled kernel's kernels/jitter/variance
  floor -> both backends integrate the SAME density, asserted == tulpa
  `cpp_test_svc_nngp_twins` to 1e-8. Hypers (sigma, phi) grid-integrated on both
  routes (`laplace` == `nested_laplace` here) -> `fit$svc_hyper`; surface ->
  `fit$svc_field` (NUTS naming). Surface cor matches the NUTS path on the same truth
  (information-bounded, NOT backend-bounded). `fitted()` adds the surface in-sample
  via `model$occ_eta_offset`; `predict(newdata=)` does NOT krige to new locations (as
  on NUTS). Gated on occu(): detection-arm svc, a spatial/temporal/re term alongside
  svc, `pg_gibbs` -- all error w/ pointer. `test-occu-svc-laplace-recovery.R` +
  `test-svc-guard.R`.

## occu_cover(): the per-visit view behind every diagnostic (#185)

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
## community_latent.R: why loadings come from the marginal, not the joint mode (#153/#156)

**Loadings by MARGINAL likelihood, NOT the joint mode (#153 -> #156).** Factor update holds zeta at its joint mode -> `(zeta, lambda)` is a joint-likelihood estimate with `Ns*Q` incidental params growing with the sample = Neyman-Scott, inconsistent. Site factors' estimation error lands in the fitted co-occurrence and lambda absorbs it, over-fit growing with Q/S. Fix = `.tobs_latent_factor_mmle()`: EM on the SAME joint site marginal over all S*Q loadings (E-step = posterior `p(z_i|y_i)` off `.tobs_latent_joint_grid()`; M-step = per-species Qk-dim weighted Newton, backtracked on expected complete-data ll). Numbers in `NOTES_measurements.md`.

## ms_abun(negbin): the mu_log_r interval, and the variance-component gate (#235/#250/#280/#285)

**The interval is calibrated; the filed miss was the seed block's DRAW (#280 -> #235 ->
#285).** `mu_log_r` is a POPULATION mean. Each seed draws `S` log-dispersions around it
and the fit sees only those, so the error splits `est - mu = (est - mu_real) +
(mu_real - mu)` with `mu_real` the seed's own realized species mean. The second term is
`N(0, sigma_logr^2 / S)`, belongs inside the interval, and the SE does carry it -- but
measured over ~20 seeds it carries most of the spread, so an ordinary fluctuation in the
DRAW reads as a miscalibrated SE. Put the draw at its expectation AND rebuild the SE at
the simulated sigma and the interval scale is calibrated at EVERY group count (`S = 8 /
18 / 36`); the `S = 18` spike #280/#285 filed is that seed block drawing its species
means ~18.5% wide, not the estimator. Numbers in `NOTES_measurements.md`.

Two readings this replaces, both from the same data. #280's "~24% too narrow, uniformly"
is wrong: a handful of fits carried a third of `sum(z^2)`. #235's successor reading --
that the miss is CONDITIONAL on `sigma_log_r`, healthy-sigma fits covering at nominal and
near-boundary fits not -- is the right shape for the boundary phenomenon the gate below
addresses, but it is not what drives the pooled `k_hat` at `S = 18` either. REFUTED en
route, each against the same fits: finite-sample-in-`S`, a `t(S-1)` df correction, and
heavy tails at `S = 18` (top-3 share of `sum(z^2)` came in BELOW what a normal expects).

**What IS real: `sigma_hat` attenuates and the SE is built on it.** Monotone in `S`, so it
shrinks with group count rather than indicating a defect at any one. Documented on
`?ms_abun`; not gated.

**Why the coverage test asserts what it does.** A coverage floor cannot be a defect
tripwire here: rebuilding the SE at the TRUE sigma -- strictly better than any estimator
delivers -- still misses on those seeds, because their draw excess is untouched by it. So
`test-ms-abun-nb-rs-coverage.R` asserts unbiasedness, normal `z`, a gross floor, and the
draw-corrected interval scale in a deliberately WIDE band (its own null is `[0.81, 1.19]`
at 19 seeds). Tightening it needs MORE SEEDS, not a smaller number: separating `1.28` from
noise takes `n >= 40`, and any future arm needs its realized draw reported beside its
`k_hat`.

**The boundary gate (#250).** A collapsed variance component is invisible on the fit it
comes on -- the optimizer converges, the point estimate is ordinary, nothing warns -- while
the community mean it scales carries an interval that shrinks with it. So the gate has to
be the component's OWN uncertainty; an absolute cut on `sigma_hat` is a number with nothing
behind it and does not transfer between fixtures.

A 1x1 covariance block's integration coordinate IS `log(sigma)` -- the log-Cholesky diagonal
of a 1x1 factor -- so `tulpa_re_aghq()`'s `re_par_se` is `SE(log sigma)` with no transform,
and the delta method gives

    W = sigma_hat / SE(sigma_hat) = sigma_hat / (sigma_hat * SE(log sigma_hat))
      = 1 / SE(log sigma_hat)

with `sigma_hat`'s own scale cancelling out. `sigma = 0` sits on the boundary of the
parameter space, where the one-sided null is the 50:50 mixture `0.5 chi^2_0 + 0.5 chi^2_1`
(Self & Liang 1987; Stram & Lee 1994) rather than `chi^2_1`; for `W >= 0` that mixture gives
`P(W > c) = 1 - Phi(c)`, so the boundary-aware critical value and the ordinary one-sided
normal quantile coincide at `qnorm(1 - alpha)` (`.TOBS_VC_BOUNDARY_ALPHA` = 0.05). Nothing
here is fitted to a fixture. `re_par_se` is a block of the SAME joint inverse Hessian
`theta_cov` is the top-left block of, so this reads the fit's own curvature and costs no
solve. A block carrying several coordinates has correlations in it and no single SD to
test, so `.tobs_aghq_variance_boundary()` declines (`available = FALSE` + a `reason`)
rather than reporting one diagonal as though it stood alone -- the record comes back in
ONE shape whatever happened, so a caller reads fields rather than branching on NULL.
`.tobs_warn_variance_boundary()` raises ONE warning naming every failing component, so a
fit with two collapsed blocks does not raise two warnings a reader has to correlate.

Upstream: `re_par_se` / `re_par_layout` / `joint_cov` arrived in tulpa `12b641d`
(gcol33/tulpa#418), first released v0.1.18; the DESCRIPTION engine floor is already above
it, so there is no availability gate to write.

**The gate is NOT the whole guard.** A failed per-species posterior solve returns the
community dispersion block at the values it STARTED from -- `sigma_log_r` exactly the
initial 0.5, `mu_log_r` a hair below `log(r_init)`, `converged = TRUE`, and the tightest
`log_r` SE of the sweep. A small-sigma detector reads that as healthy, precisely because
the initial value is not small. What catches it is the engine refusing an optimum carrying
its failure sentinel plus the per-species solve status the fit now carries (#281) -- both
read in the coverage loop, and both counted, so the denominator is visible rather than
assumed.

**The PC prior is opt-in.** `control$logr.sigma.prior` / `control$omega.sigma.prior` add
curvature at the boundary; the DEFAULT stays pure ML, so no fit changes unless asked.
## occu_multiscale_cover(): MAR cover + the shared detected-plot score (#262/#270)

**MAR cover, same rule as the twin** (#262): a detected visit (`y=1`) w/ NA cover keeps
its detection term + drops only `f_pos`. Builder = the SHARED
`.occu_cover_validate_pos_values()` (NA sentinel + zero-fill), gate = `isfinite(y_pos)`
in the cell-coupling spec, the NUTS target AND the R marginal
(`.occu_mscale_cover_nonspatial_ll`); the ploglik kernel already had it. A detected
plot factorises -> the cover factor is ADDITIVE, so under `laplace` the psi/theta/p
estimates are the full-data ones (assert `expect_equal` there, not just "close"); under
`nested_laplace` the shared field couples them, so close, not equal. Dispersion pre-fit
(`occu_multiscale_cover{,_joint}.R`) filters `is.finite(pos_vals)` -- an unfiltered
`sd(log(pos_vals))` is where an NA reaches the starts. `test-occu-multiscale-cover-mar-cover.R`.

**Detected-plot score written ONCE** (#270): `mscale_det_plot_block()`
(`occu_coupling_shared.h`), beside its `mscale_nodet_cell()` sibling; the coupling spec
and the NUTS target each supply eta/y accessors + `emit_p`/`emit_pos` sinks, so neither
allocates and each keeps its own row layout. Cover-arm dispatch via `PosPolicyAccess<P>`
(spec, compile-time) / `PosCodeAccess` (sampler, runtime code) -- do NOT collapse those
two into one runtime switch, it would put a branch in the Laplace inner loop.
`want_logdisp=false` keeps beta's digamma out of the grid-integrated paths.

## community_latent.R: backtracking, the initialise-once estimator, the offset solve

**Backtracking + guards.** Factor Newton (`.tobs_latent_factor_update`) backtracks: local `ascend()` halves the step until the penalized objective improves, holds previous iterate if never; `nstep()` ridge-bumps singular curvature. Non-finite guards are inline (`if (all(is.finite(Dz)))` / `(Dl)`), not a named helper; non-finite `working()` score/curv `break`s the pass. Field solve `.tobs_latent_field_solve` has its own local `safe_solve()` (ridge retry for a singular Hessian only, Newton update unconditional) -- do not confuse the two, there is no `safe_step()` anywhere in the repo.
**ONE estimator, ONE state (#156).** `.tobs_latent_factor_update()` + `.tobs_latent_factor_scale()` run ONCE, outer pass 1, purely to INITIALIZE -- the marginal's lambda-gradient vanishes at lambda=0, and the 1-D bracket is a global magnitude search the local EM cannot do. Running the joint-mode update every pass alongside the MMLE diverges both ways, so it never repeats.
**Offset by SCORE-MATCHING, NOT `zeta t(lambda)` (#156).** `.tobs_latent_factor_offset()` solves `score(eta+off) = E_z[score(eta+lambda_s'z)]` per cell (scalar Newton) -- plug-in and integrated stationary conditions match for ANY family, no link-specific derivation; reduces to `lambda'zhat + v/2` on a Poisson log link. `fit$model[[offset_slot]]` reads THIS, not `zeta t(lambda)`.

## field_offset.R: writers and readers of the per-arm eta offset (#254)

field_offset.R            — the fitted latent field as a per-ARM eta offset (#254). `model$field_eta_offset` = one entry per process, in THAT process's own design row layout (removal's per-PASS detection arm included). Written by `.tobs_set_field_eta_offset()` (areal-BFGS `.tobs_attach_field_results()`, `.tobs_nuts_field_attach()`, the two C++ areal branches abun/removal) + `.tobs_attach_model_eta_offset()` at the ONE dispatch tail `.tobs_finalize_family_fit()` -- after `fit$model <- model` drops whatever the fitter set on its autoscaled copy. Read by `.tobs_eta_offset()` (vector), `.tobs_add_eta_offset()` ([n_draws x n_rows] ploglik), `.tobs_sim_arm_block()` (simulate: the offset rides in as a design COLUMN w/ its coefficient pinned at 1, so the C++ simulators are untouched and a field-free fit is byte-identical). EVERY door that rebuilds eta from `X_processes` goes through these -- fitted / residuals / predict / simulate + the ploglik behind WAIC / LOO / DIC / CPO -- so a field cannot be scored at 0 again. `.tobs_nuts_field_loglik()` (same tail) re-evaluates `log_lik`/`log_prob` on a SAMPLED-field fit through that ploglik kernel: the field is not a draw column, so the family's own posterior-mean marginal ran it at offset 0 and `logLik()`/AIC/BIC described a model WAIC did not score (grid-integrated paths already carry the field in their marginal, keep their grid-weighted value). `count()` + the community families keep their own named slots (`count_field_offset`, `occu_field_offset`, ...) written by `.tobs_latent_attach_field()`; the shared field -> per-site contribution reader is `.tobs_spatial_field_offset()` (was `.count_spatial_field_offset`, count_spatial.R)

## engine_defaults.R: single-species divergences + the n.quad route table

    Divergences RECORDED (`.TOBS_SINGLE_SPECIES_NUTS` / `_LAPLACE`, via `.tobs_single_species_defaults(engine)`): `.tobs_fit_model()` entry keeps `sigma.beta = 10` (vs community 5, on laplace/nested_laplace too), `adapt.delta = 0.8` (vs 0.9), `seed = 42` (vs 1). Each = value the single-species recovery/coverage tests were calibrated against; override belongs to the ENTRY, not the nine families passing through -> one constant, not nine family rows.
    **`sd.load` (1.0) + `re.lkj` (1.5) = laplace rows (#189)**; `n.quad` NOT one number, deliberately -- `.TOBS_NQUAD_ROUTES` / `.tobs_n_quad(route)` enumerates seven marginals (`re_aghq` 9, `ms_nmix` 1, `ms_nmix_scalar` 2, `ms_occu_cover` 5, `cover_latent_beta` 15 vs `cover_latent_lognormal` 1 (closed form, needs none), `community_latent` 5). `?tobs` lists per route, no single default most routes do not use

## nuts_chains.R: where the convergence record is written from (#174)

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
`summary()` row; its name-set check fails when a new NUTS family lands w/o a
case. `.tobs_nuts_rhat_ess()` (`list(rhat, ess)`
accessor for the `fit$nuts` full-coordinate block + PG-Gibbs summariser) reads
the same table. PG-Gibbs still keeps rhat/ess at `fit$rhat`/`fit$ess` ONLY,
where summary/print do not look.

## Copy-axis resolution (`control$alpha.n`, #287, tulpa >= 0.2.6)

The alpha axis carries prior structure (atom at 0 + log slab over [0.1, 3]), so
`control$alpha.grid` / `share(alpha = grid(...))` STATE its nodes and restate that
structure with them; `control$alpha.n[.trend]` states a RESOLUTION -- the engine
re-reads its OWN axis w/ n SLAB nodes (axis length n+1), atom + bounds unchanged.
Needed because the alpha axis does NOT densify when the donor `sigma.grid` does, so
outer-grid quadrature ESS saturates on it (`NOTES_measurements.md`). ONE resolver
`.tobs_alpha_axis()` (`R/joint_substrate.R`) -> the two engine-facing fields
(`alpha_grid`, `alpha_n`), shaped by `.tobs_alpha_copy_spec()` (multi-block driver) /
`.tobs_alpha_field_coef()` (single-block pos arm); EVERY joint route reads it --
cover, occu_cover, occu_multiscale_cover, MCAR + coupled-trend blocks included. Both
spellings on ONE block = error (`.tobs_check_alpha_control()`, raised in the
DISPATCHER so the message names the knob typed); a `share()` that STATES nodes beside
`alpha.n` = error (`.tobs_check_alpha_copy()`), a bare `share(spatial())` composes (it
asks for the default axis, and on occu_cover now records NULL rather than resolving
the default nodes, so the resolution reaches it). A block with no `share()` is pinned
`alpha = 0` and has no axis to resolve -- stated nodes win there. NUTS resolves the
axis to NODES instead (`.tobs_alpha_nodes()`): the sampled alpha's flat prior takes
the realised node set's span as its support.
