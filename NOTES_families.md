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

## Community recovery tolerance (#153/#155) -- why the naive tolerance is looser than it reads

CLAUDE.md's CI section states the rule (score a multi-seed community recovery
against `colMeans(bs)` and assert absolutely, not against the nominal constant
with a relative `tolerance`). Two things make the naive form loose:

- *Estimand = wrong constant.* The community simulator draws per-species coefs
  around a POPULATION mean -> the seed's realized mean sits `beta_sd/sqrt(S)`
  off that constant, so scoring vs the constant spends most of the budget on
  draw noise. Score vs `colMeans(bs)` (seed's own realized mean) -> pure
  estimator budget. `.jsdmc_sim` returns `beta_real`; most other community
  simulators still return only the nominal constant (#155).
- *`tolerance` silently switches scale.* `all.equal.numeric` = relative while
  target exceeds tolerance, absolute below -> the same nominal `tolerance`
  means a much wider absolute band on a big coefficient than a small one.
  Assert absolutely (`expect_lt(abs(est - truth), tol)`) when the budget is
  meant to be uniform.

## Observation families: per-family mechanics

Moved verbatim out of `CLAUDE.md` Architecture.

- **`removal`** (#39): sequential depletion, pass k sees `N - sum_{l<k} y_l` trials
  = depleting-binomial = multinomial removal; latent N to `K_max`. Shares the
  count-marginal Laplace (`src/marginal_count_laplace.h`) + NUTS
  (`src/marginal_count_nuts.h`) w/ abun; per-site `src/removal_kernel.h` over shared
  `accumulate_count_moments`/`fill_nb_dispersion` (`src/nmix_kernel.h`). Areal
  icar/car_proper/bym2 on the abundance arm via `nested_laplace` (#51):
  `compute_removal_site` returns the same `NMixSiteResult` moments as the Royle
  kernel -> reuses the templated count-spatial driver
  (`src/nmix_count_spatial_driver.h`, byte-identical to nmix) via
  `src/removal_spatial.cpp` + `R/removal_spatial.R`; one unit/site, K_max=per-site
  total. **Detection-arm areal field (#114)**: `detection=~icar()/car_proper()/
  bym2()` loads a spatially-varying capture logit on eta_p via the shared areal-BFGS
  driver (`.tobs_fit_removal_spatial_bfgs(det_arm=TRUE)`); removal's detection
  design is per-PASS, so the per-pass grad_eta_p is `rowsum`-summed to a per-site
  field gradient + the per-site offset expanded via site_idx. `fit$spatial_field_arm`
  labels it. NUTS carries the field on the abundance arm only (detection-arm ->
  nested_laplace). spde + NUTS+detection-arm gated.
- **`distance`** (#38): latent N in covered region, per-bin counts multinomial over
  `(bin 1..B, undetected)`, `pi_b = int_bin g(x;sigma) f(x) dx` (half-normal/hazard
  key, line/point transect density `f`); bin integrals + 1st/2nd eta-derivs by
  Gauss-Legendre quad. Reuses the count-marginal core (`src/nmix_kernel.h`); det arm
  = site-level `log sigma` + optional scalar hazard shape (`src/distance_quad.h`/
  `src/distance_kernel.h`), own Laplace driver (`src/distance_laplace.cpp`) + NUTS
  (`src/distance_nuts.cpp`). K_max default `3*max(rowSums)+100`. Laplace grouped RE
  on the abundance arm (one grouping factor, dim<=3) via the shared count AGHQ
  (`DistanceGroupedOracle` over `CountGroupedOracle`; theta = count layout
  `[beta_lambda|beta_sigma|log_r?]`). **hazard-key grouped RE (#114)**: the scalar
  log-shape eta_b is not a per-site design column, so rather than a second global
  theta slot it is PROFILED over the AGHQ log-marginal (`.tobs_distance_re_aghq`
  optimises eta_b in an outer loop, each candidate a full AGHQ fit at a FIXED shape
  the oracle carries via `key_code`/`eta_b`); the profile log-shape row/col is
  inserted into the AGHQ vcov (`.tobs_distance_insert_shape_vcov`, profile SE,
  off-diag 0). Pois + NB; det-arm RE gated (couples a site's bins through the latent
  N). Areal icar/car_proper/bym2 on the abundance arm via `nested_laplace` (#51,
  `R/distance_spatial.R` over shared areal-BFGS driver `R/areal_bfgs.R`): the
  distance marginal = bin-multinomial (NOT N-mixture, so not the count-spatial
  driver) but exposes an analytic per-site gradient (`cpp_distance_site_sweep` over
  `compute_distance_site`), so BFGS + FD-Hessian recovers observed info.
  **DETECTION-arm areal field (#114)**: `detection=~icar()` loads on the per-site
  detection scale eta_sigma (`det_arm` routes offset->eta_sigma,
  grad_eta=sw$grad_sig); works under BOTH the half-normal AND the hazard key (field
  on eta_sigma, log-shape eta_b threaded as a global in the same eval). Pois + NB;
  NUTS+spatial gated.
- **`fp_occu`** (#40): Miller 2011 multistate false-positive occupancy, `y in
  {0,1,2}` (none/ambiguous/certain); certain detections only at occupied sites
  identify it. Latent z summed (2-state); 4 logit arms psi/p11/p10/b
  (`fp_formula`/`b_formula`, default `~1`). Laplace = analytic-grad BFGS, vcov = inv
  of -FD-Jacobian of the analytic grad (`src/fp_occu_kernel.h`, gradient only). NUTS
  `src/fp_occu_nuts.cpp`. Laplace grouped RE on the psi OR p11 arm (one factor,
  dim<=3) via pure-R `make_site` AGHQ (`.tobs_fp_occu_re_aghq`, no native oracle,
  branches on arm): holding the other arms fixed makes the 2-state marginal
  `psi*A+(1-psi)*B` a fn of one scalar offset/site (psi shifts the mixture weight;
  p11 shifts the occupied emission `A`, uniform over a site's visits), both
  closed-form d1/d2. RE both arms at once rejected; p10/b never carry structured
  terms; NUTS samples psi-arm intercept RE only. Det-arm RE params named
  `sigma_p<t>_*` (p11 process). Areal icar/car_proper/bym2 on the psi arm via
  `nested_laplace` (#51, `R/fp_occu_spatial.R` over `R/areal_bfgs.R`
  `.tobs_areal_bfgs_fit`, shared w/ dyn_abun): BFGS over the two-state marginal
  (`cpp_fp_occu_total_log_lik` analytic grad) + CAR prior, FD-Hessian observed info;
  one unit/site. Occupancy fields more weakly identified than count (one binary
  site/node). **Detection-arm areal field (#114)**: `detection=~icar()` loads a
  spatially-varying true-positive logit on the per-site eta_p11 (per-site detection
  design, so no aggregation); p10/b never structured. NUTS carries the field on the
  psi arm only.
- **`dyn_abun`** (#37): Dail-Madsen open N-mixture, `N_1~Pois/NB(lambda)`,
  `N_t=Binom(N_{t-1},omega_{t-1})+Pois(gamma_{t-1})`, `Binom(N_t,p)` obs. Latent N
  sequence summed by exact HMM forward over 0..K_max; analytic grad by forward-mode
  diff (`src/dyn_abun_kernel.h`). 4 arms lambda/p/omega/gamma
  (`omega_formula`/`gamma_formula`, default `~1`). NUTS `src/dyn_abun_nuts.cpp`.
  K_max default `max(count)+40` (forward ~cubic in K). Pois OR negbin init
  (`mixture="negbin"`, trailing `log_r`). **Season-varying omega/gamma (#80)**: the
  transition from season t-1 to t uses interval-(t-1) vital rates, so a covariate
  carried as an `[n_sites x (T-1)]` matrix column of `data` drives the dynamics. The
  kernel takes interval-indexed `eta_omega`/`eta_gamma` (length T-1) with
  per-interval forward-mode gradients (direction iv born at its own transition,
  propagated through later ones); the scalar overload broadcasts -> byte-identical
  constant-rate path. Binder `.tobs_dyn_abun_arm_design` unrolls to long-form
  `[(site x interval) x p]` only when a covariate is such a matrix column; both
  laplace + NUTS scatter the per-interval score through it.
  `simulate_dyn_abun(beta_omega=, beta_gamma=)` season-varying truth. NUTS+temporal
  still gated. Laplace grouped RE on the initial-abundance arm (one factor, dim<=3):
  the RE shifts only eta_lambda -> the per-site marginal factorises
  `L(eta_lambda)=sum_n1 pi_n1(eta_lambda) c(n1)`, conditional `c(n1)=P(all
  data|N_1=n1)` (O(K^2 T) HMM BACKWARD pass `compute_dyn_abun_init_weights`)
  precomputed ONCE per make_site; each AGHQ node = O(K) dot
  (`cpp_dyn_abun_init_loglik`). Pure-R `make_site` AGHQ (`.tobs_dyn_abun_re_aghq`),
  Pois+NB; omega/gamma never structured. **Detection-arm RE (#82)**: a p-arm RE
  shifts eta_p, which enters EVERY season's obs pmf, so `c(n1)` cannot be
  precomputed -- each AGHQ node re-evaluates the full O(K^2 T) forward marginal via a
  closed-form SECOND-ORDER eta_p forward-mode pass (`compute_dyn_abun_p_curv` ->
  `cpp_dyn_abun_p_loglik`: transition/initial are p-free, so `d obs`/`d2 obs` are
  the only source terms, propagated + renormalised per season; logL/d1/d2
  FD-validated). NUTS p-arm RE (`re_arm` 0=lambda|1=p): non-centered offset to
  eta_p, grad via `grad_eta_p`; both-arm RE rejected (AGHQ integrates one arm).
  Det-arm RE params named `sigma_p<t>_*` / `log_sigma_p_*`. Shared pmf helpers
  (`da_obs_season_pmf`/`da_recruit_pmf`/`da_binom_pmf_row`) back the forward gradient
  kernel, the backward `c` pass, AND the eta_p second-order pass. Areal
  icar/car_proper/bym2 on the initial-abundance arm via `nested_laplace` (#51,
  `R/dyn_abun_spatial.R` over `R/areal_bfgs.R`): BFGS over the exact forward-HMM
  marginal (`cpp_dyn_abun_total_log_lik` analytic grad) + CAR prior, FD-Hessian
  observed-info Laplace integrated over tau[,rho]; one unit/site, Pois+NB.
  **Detection-arm areal field (#114)**: `detection=~icar()` loads a spatially-varying
  detection logit on the per-site eta_p (applied across every season, so no
  aggregation); omega/gamma never structured. NUTS carries the field on the
  initial-abundance arm only (detection-arm -> nested_laplace).

## Feature roster: per-row detail

Moved verbatim out of the `CLAUDE.md` "What works (tested)" table, which keeps the support matrix only. `n-L` = nested_laplace.

| Feature | Laplace | NUTS | Notes |
|---|---|---|---|
| Single-season occupancy | Yes | Yes | parity w/ inlaocc. **`method="pg_gibbs"`** (#126, spOccupancy `PGOcc`): a REAL Polya-Gamma Gibbs chain over the exact posterior (Polson-Scott-Windle, tulpa `cpp_rpg`), NOT the stochastic-EM `laplace_gibbs`; z|. Bernoulli, then omega~PG conjugate-Gaussian beta updates on both arms (`R/occu_pg_gibbs.R`, `.tobs_fit_occu_pg_gibbs`). Posterior matches the Laplace observed-Fisher fit; Rhat/ESS calibrated, nominal coverage. Wired via the method route table + `.tobs_control_allow(pg_gibbs="sampler")`. **spPGOcc** (`occu() + icar()` + pg_gibbs, `.tobs_fit_occu_pg_gibbs_spatial`): intrinsic ICAR field on the psi logit, jointly updated w/ beta as a GMRF (dense precision `t(W)Omega W + blkdiag(B0inv, tau Q)`, `W=[X|I]`); conjugate Gamma tau; field sum-to-zero'd each sweep, level moved to the intercept (else eta leaks). NB: a site-level occupancy covariate spatially-confounds w/ the saturated 1-node-per-site field (known property, spOccupancy too) -> recovery target = field+intercept+detection. icar only (bym2/car follow-ups). `test-occu-pg-gibbs.R` + `test-occu-pg-gibbs-spatial.R` |
| Dynamic (HMM) | Yes | Yes | colonization/extinction. **Season-varying gamma/epsilon (#124)**: a `[n_sites x (T-1)]` matrix column of `data` drives interval-indexed rates (`colonization = ~ cov`), the dyn_abun #80 recipe on the colext forward; `.tobs_interval_arm_design` (shared w/ dyn_abun) unrolls long-form `[(site x interval) x p]`. E-step forward-backward uses per-interval transition matrices (constant rates recycle -> byte-identical to pre-#124); M-step encodes the transition arm as a WEIGHTED LOGISTIC (response `xi/n_mat`, weight = origin-state prob `n_mat`) NOT the per-site M=1000 pseudo-binomial -- per interval that makes each row near-separable, so the inner Newton returns ~0 slope. Exact-marginal refine via the pure-R season-varying forward `.tobs_dyn_occu_marginal_nlp_sv` (`R/dyn_occu_marginal.R`; the cpp `cpp_occu_dynamic_ploglik` reads one rate/site) -- escapes the EM local optimum AND calibrates SEs. **Season-varying DETECTION (#124)**: a `[n_sites x T]` matrix column on `detection` drives per-SEASON `logit p_it` through the period-agnostic unroller (`.tobs_period_arm_design` + season wrapper `.tobs_season_arm_design`, single source of truth): emission reads `p_mat[i,t]`, M-step (+ `hard_encode`) one detection row per `(site,season)`, refine + gate carry `det_season_varying`; `.tobs_dynamic_smoothed_z` (`fitted()$z`) indexes detection per season + transitions per interval. `simulate_dyn_occu(beta_gamma=, beta_epsilon=, beta_det_season=)`. Laplace/gibbs/mi only; NUTS+SV gated (C++ forward per-site). 20-seed recovery + coverage (`test-dyn-occ.R`) |
| Multi-season AR1 year effect (`t_occu`) | — | — | `t_occu()` (#124, spOccupancy `tPGOcc`): NOT colext -- per-`(site,season)` `z_{i,t}~Bern(psi_{i,t})`, `logit psi = X beta + eta_t` w/ a SHARED AR1 year effect `eta_t = rho eta_{t-1} + w_t`, NO colonization/extinction. Given the year effects the seasons factorise -> EXACT Polya-Gamma Gibbs (`method="pg_gibbs"`, the engine spOccupancy uses; `R/t_occu.R`, `.tobs_fit_t_occu_pg_gibbs`): draw z, then joint `(beta_occ, eta)` as ONE GMRF update (dense precision `blkdiag(t(Xf)Omega Xf + B0inv, Z'Omega Z + Q/sigma^2)` w/ stationary AR1 precision `Q` = `.t_occu_ar1_Q` as the year-effect prior, as spPGOcc's ICAR precision is for a spatial field), then `beta_p`, then AR1 hypers (`sigma^2` conjugate IG, `rho` on a grid over the AR1 log-density). Year effect sum-to-zero'd each sweep, level moved to the intercept; AR1 hypers update on RAW `eta_raw` BEFORE centering (else centering distorts the increments -> rho collapses). `fit$temporal_field` = posterior-mean year effect; recovers it plus occ/det coefs + `sigma`. AR1 `rho` WEAKLY IDENTIFIED at few seasons -- `rho_hat` climbs toward truth only as T grows, a property of the short series NOT the sampler -> reported, not asserted tightly. `y` = 3D `[n_sites x n_seasons x max_visits]` or list of per-season matrices; `simulate_t_occu()`, `test-t-occu.R`. v1 = site-level occ/det covariates, `pg_gibbs` only (season-varying covariates, areal psi field, NUTS = follow-ups) |
| Community single-season (`ms_occu`) | Yes | Yes | per-arm community RE, shared community Laplace-EM (`R/community_em.R`, `R/ms_occu.R`); msPGOcc. **`pg_gibbs`** (#115/#126, `R/ms_occu_pg_gibbs.R`): hierarchical PG Gibbs -- per-species PG-augmented conjugate coef updates + conjugate community mean + Inverse-Gamma (near-Jeffreys) community variance draws (diagonal per-arm covariance). Gives a CALIBRATED community-VARIANCE posterior: Laplace-EM `sd_psi`/`sd_p` attenuate (documented lower bound), Gibbs recovers them (sd = posterior MEDIAN, robust to the variance skew). `test-ms-occu-pg-gibbs.R`. NUTS non-spatial samples community means/deviations/covariances jointly (#69, `R/ms_occu_nuts.R`); shared areal field icar/bym2/car_proper/spde on the occ arm via n-L (#75, `R/ms_occu_spatial.R`, sfMsNMix analogue; icar/bym2/car_proper and count/N-mixture's spde arm share ONE Laplace-EM driver, `src/community_spatial_em.h`, #239). NUTS+field -> n-L. **SVC (`svcMsPGOcc`, #118)**: a varying-coefficient bar `spatial(~ 1 + w || cell, graph)` on the occ arm -> intercept + SVC field(s) by BLOCK COORDINATE ascent (`R/ms_occu_field.R`, `.ms_occu_field_solve` = two-state-marginal occupancy field Newton + ICAR prior). Both fields recover; icar only, and a plain intercept field stays on the C++ path (no regression). `test-ms-occu-field.R` |
| Community dynamic (`ms_dyn_occu`) | Yes | Yes | per-species psi1/p RE + shared gamma/eps; HMM-forward; `R/ms_dyn_occu.R`. **NUTS** (#115, `R/ms_dyn_occu_nuts.R`, `src/ms_dyn_occu_nuts.cpp`): non-spatial community sampler over the exact per-(species,site) HMM-forward marginal via an in-tree C++ FullGradFn -- community means, per-species psi1/p deviations, both per-arm covariances AND the shared gamma/eps globals jointly, non-centered `b_{s,arm}=C_arm z_{s,arm}` (the two shared transition arms carry no RE; per-arm eta scores from `.ms_dyn_occu_fb_vec`, shared verbatim w/ the stMsPGOcc field fitter). Byte-exact vs the R oracle, warm-started at the Laplace-EM mode, 0 divergences, de-attenuates the community variance the EM under-reports. Non-spatial only. `test-ms-dyn-occu-nuts.R`. **`pg_gibbs`** (#115/#126, `tMsPGOcc`, `R/ms_dyn_occu_pg_gibbs.R`): msPGOcc community PG machinery + a 2-state HMM forward-filter backward-sample latent step (SHARED gamma/eps from the aggregated 0->/1-> transitions across all species). Calibrated community-variance posterior; shared gamma/eps recover tightly, sd_psi1 recovers. Constant transitions, site-level detection. `test-ms-dyn-occu-pg-gibbs.R`. **Shared areal field on psi1 (stMsPGOcc, #123)**: `~ 1 + icar(graph=adj)` under `nested_laplace` -> `R/ms_dyn_occu_spatial.R`. KEY MATH: psi1 sets ONLY the initial HMM mixing weight, so the per-(species,site) marginal is LINEAR in psi1 -- `L=(1-psi1)A+psi1 B`, A/B = forward likelihood cond on season-1 state 0/1 -- the IDENTICAL mixture to the single-season `ms_occu` field oracle, so its Louis score/curv (`score=r-psi1`, `curv=psi1(1-psi1)-r(1-r)`, `r=psi1 B/L`) carry over verbatim (only A/B differ). Routed through the shared block-coordinate driver (`R/community_latent.R`, `.tobs_community_latent_ascent`). Kernels VECTORIZED over sites + ANALYTIC Fisher-identity gradient (`grad_psi1=X'(w1-psi1)`, `grad_p=X'(sum_t w_t(ndet-nvalid p))`, `grad_gamma=X'(col_y-gamma col_n)`) so the community EM skips its O(U^2) FD Hessian (a per-site-loop version was ~10x too slow). icar only; field + community means recover. **svcTMsPGOcc (#123)**: a weighted bar `spatial(~ 1 + w || cell, graph)` adds a shared covariate-weighted field beside the intercept field -- NO new code, the psi1 oracle already returns per-site/species score+curv and the K-field weighted-ICAR solve is the SAME `community_latent.R` machinery as svcMsAbund. `simulate_ms_dyn_occu(field=, trend=)`. NUTS/bym2/car follow-ups |
| Community integrated (`ms_int_occu`) | Yes | Yes | per-species psi + per-source det RE; multi-source two-state marginal; `R/ms_int_occu.R`. **NUTS** (#115, `R/ms_int_occu_nuts.R`, `src/ms_int_occu_nuts.cpp`): multi-source generalisation of the `ms_occu`/`ms_dyn_occu` samplers -- community means, per-species occupancy / per-source detection deviations, and the D+1 independent per-arm covariances jointly over the exact multi-source two-state per-(species,site) marginal via in-tree C++ FullGradFn (per-arm non-centered `b=Cz`, NO shared globals), warm-started at the Laplace-EM mode; reuses the shared `.ms_ocs_*` epilogue (#128). Byte-exact vs R oracle, 0 divergences, de-attenuates the community variance the Laplace-EM under-reports. Non-spatial only. `test-ms-int-occu-nuts.R`. **`pg_gibbs`** (#115/#126, `R/ms_int_occu_pg_gibbs.R`): msPGOcc generalized to D detection arms -- per species draw the single latent z (occupied if any source detects, else Bernoulli on the pooled occupied-undetected mass), PG-conjugate `beta_psi_s` + D per-source `beta_pd_s` (each at that species' occupied covered sites), conjugate community mean + IG variance per arm. Calibrated community-VARIANCE posterior, community means recover; `ms_community` layout matches the Laplace fit (`Sigma_/sd_/coef_/blup_<arm>`). `test-ms-int-occu-pg-gibbs.R` |
| Integrated multi-source | Yes | Yes | shared psi |
| Multi-season integrated occupancy | Yes | — | `dyn_int_occu()` (#122, spOccupancy `tIntPGOcc`): dynamic occupancy (multi-season HMM psi1/gamma/eps) x integrated occupancy (per-season emission pooling S detection sources). Pure-R two-state HMM forward w/ multi-source pooled emission (`R/dyn_int_occu.R`, the `dyn_occu_marginal.R` forward generalised), optim BFGS + observed-Fisher vcov, no new C++. `y`=list of S `[sites x visits x seasons]` arrays; `colonization=~`/`extinction=~` required (as dyn_occu); shared per-source detection design (`p_<src>` arms); `fit$means` names `psi1_*`/`gamma_*`/`eps_*`/`p_<src>_*`. **Partial season overlap**: a source absent at a (site,season) -> NA -> contributes nothing to that season's emission (nvalid=0); a (site,season) unobserved by EVERY source is marginalised (e0=e1=1) by the forward -> staggered surveys NA-padded to the common `[n_sites x max_visits x T]` grid (`simulate_dyn_int_occu(source_seasons=list(1:4, 3:6))`). Two boundary ANCHORS pin the reduction: one source (others all-NA) reproduces `dyn_occu()`; one season of data (`T=2`, season 2 all-NA) reproduces `int_occu()`. **stIntPGOcc**: `~ 1 + icar(graph)` loads a shared areal field on psi1 via the areal-BFGS driver (`.tobs_fit_dyn_int_occu_spatial`); field gradient = Fisher-identity psi1 score `w1 - psi1` (`.dio_fb`, FD-validated; also sped the non-spatial fit). **svcTIntPGOcc (#122)**: a `spatial(~ 1 + w || cell, graph)` bar adds an SVC field beside the intercept field -- the driver already takes a LIST of field blocks and scatters `w1 - psi1` to each, so the weighted block = ICAR field w/ a `w`-weighted loading (`.areal_field_car(weight=)`, byte-identical unweighted). `fit$spatial_field`=intercept, `fit$trend_field(s)`=SVC; both recover. `simulate_dyn_int_occu(field=, trend=)`. icar only. v1 = constant transitions, site-level detection; season-varying rates (the #124 recipe) + bym2/car_proper + NUTS = follow-ups. Full S3 + WAIC. `test-dyn-int-occu.R` + `test-dyn-int-occu-areal-recovery.R` |
| JSDM (`jsdm`) | Yes | Yes | `jsdm()` = COMMUNITY GLMM on observed presence/absence (#121): per-species coefs + Gaussian community covariance, NO detection/latent state = `ms_count()` w/ logit link -> SHARES that binder (`.tobs_build_ms_count(response="bernoulli")`, model_type `ms_count`), community EM, latent driver, NUTS target + S3. `latent(n)` = lfJSDM; + shared field = sfJSDM; icar/car_proper/bym2/spde field via n-L. NUTS = exact joint community posterior over the Bernoulli response (`MSC_BERN` in `src/ms_count_nuts.cpp`, byte-exact vs R oracle). **`pg_gibbs`** (#126): community logistic GLMM has NO latent state -> pure per-species conjugate coef update + community mean + near-Jeffreys IG community variance (`R/ms_count_pg_gibbs.R`, `.tobs_fit_ms_count_pg_gibbs`, tulpa `cpp_rpg`), shared w/ `ms_count("binomial")` (n=trials, jsdm n=1); recovers the community variance the Laplace-EM attenuates. Non-spatial, bernoulli/binomial only (poisson/negbin/gaussian rejected w/ pointer). `test-ms-count-pg-gibbs.R`. laplace_sla/gibbs/mi were single-block routes, dropped |
| Cover hurdle (joint) | Yes | Yes | `family_cover_hurdle.R`, `sla_cover_*`, NUTS `R/cover_nuts.R` + `src/cover_nuts.cpp` (no `laplace_gibbs`/`laplace_mi`). positive = beta/lognormal/lognormal_trunc/ordinal/`beta_oi`. `beta_oi` (#108) = one-inflated Beta: ceiling (cover=1) plots = a point mass (constant pi = ceiling share, binomial SE), interior Beta on (0,1); encode splits `is_pos` to interior, `enc$oi` carries pi, decode reports `pi_one`, predict conditional cover = `pi + (1-pi)*mu` (`.tobs_cover_mu`). `control$aggregate.occ` (ON, #48) collapses occ arm to Binomial suff-stat; `control$aggregate.pos` (ON beta arm, #49) collapses beta pos arm to grouped suff-stat (tulpa `slog_y`/`slog_1my`), errors on non-beta. Both byte-identical to per-plot |
| Cover hurdle spatial coef fields (`\|\|` / `\|`) | n-L | — | `spatial(~ 1 + w \|\| node, graph, to=)` independent (#61, two coupled ICAR blocks, per-field alpha) OR `\| ` correlated (#64, one separable-MCAR block sharing free Sigma). `\|` both-arm `to=c("presence","positive")` = copied to pos arm w/ one alpha (#64); `\|` single-arm `to="presence"`/`"positive"` = free-Sigma field on that arm alone, NO copy (#109, 0-sentinel `spatial_idx` on the other arm via `mc$to`, `copy=NULL`, `alpha_mcar`=NA). `\|` -> `.cover_build_mcar_spec`/`.fit_cover_hurdle_joint_mcar` (tulpa `type="mcar"` block, copy only when both-arm); reports `sigma_mcar`/`rho_mcar`/`alpha_mcar`. SLA on `\|` no-op. icar only |
| Cover hurdle arm-specific fields (single-arm `to`) | n-L | — | `spatial(~ 1 + w \|\| cell, graph, to="positive")` (or `"presence"`): separate single-arm calls = independent per-arm fields, NO cross-arm copy (#65). NO engine change -- per-arm `spatial_idx=0` makes the other arm's rows skip the block (tulpa `l_b>0` scatter guard), own precision grid-integrated. `.tobs_armspecific_bar_fields` (formula_terms.R) -> `enc$armspec` -> `.fit_cover_hurdle_joint_armspecific`; `armspec_blocks` carries per-block arm/slot/type and `.tobs_joint_draws_cover_armspecific` scatters each block onto its arm only. icar/car/car_proper AND bym2 (#107): a bym2 block is the non-copied length-2 (phi ICAR + iid theta) block on a paired (sigma,rho) grid, the draw projection reconstructing `z = sqrt(rho)*sf*phi + sqrt(1-rho)*theta` so predict/WAIC see the full mix. `\|\|` only (`\|` arm-specific undefined: copy-only). No mix w/ shared field/trend/temporal/re; one field per arm. SLA no-op |
| Joint occu + cover | Yes | Yes | `occu_cover()` — see below. NUTS non-spatial (in-tree FullGradFn over exact two-state marginal, beta/lognormal; `R/occu_cover_nuts.R`, `src/occu_cover_nuts.cpp`) AND spatial coupled areal field (#74, car_proper recovers; icar/bym2 sample + centre via the #71 sum-to-zero coupled field, #113) with its HYPERS SAMPLED (#204: field SD, bym2/car mixing rho, copy amplitude alpha — fixed basis, no per-step decomposition; `fit$nuts$sampled_hyper`/`$fixed_hyper` report per fit). Community sampled-field (per-species loadings) route = `ms_occu_cover()` factor |
| occu + cover + areal field + per-group RE | n-L | — | `occu_cover()` + icar/bym2 + `re(g)`/`(1\|g)` on psi; one iid RE block (#56); `sigma_re` + BLUPs; intercept RE only |
| occu + cover independent cover-arm field (single-arm `to="positive"`) | n-L | — | `occu_cover()` + `spatial(~ 1 + w \|\| cell, graph, to="positive")` on the occurrence formula (#110): NON-copied ICAR block(s) on the cover arm ALONE, decoupled from the occupancy field's alpha copy, so `delta_cover_cond` varies instead of collapsing when `alpha->0` (composes w/ the shared occupancy field, which psi + `delta_cover_exp` keep). Parse via `.occu_cover_spatial_fields` -> `spatial_info$pos_armspec`; the fitter appends the non-copied blocks (`spatial_idx` psi=0/p=0/pos=cell, own `tau_grid`, svc_weight for trend) after the occ fields, before RE blocks, and `ctx$field_specs` labels each block shared-vs-pos. Draw substrate reads `field_specs` (pos block amp_occ=0, amp_pos=1/sqrt(tau)) -> predict/WAIC add it to `field_pos` automatically. Reports `sigma_pos_field`/`sigma_pos_field_<col>` (from `b<k>.tau`), `fit$pos_field`/`pos_field_table(s)`. Per-visit cover (`cover_aggregate="none"`); icar only (bym2/car->icar); NOT w/ MCAR `\|`, latent cover RE, or batch. Static intercept field weakly ID'd vs the alpha copy; time-weighted trend cleanly ID'd. `test-occu-cover-pos-field.R` |
| occu + cover shared field + arm deviation (`share(residual=)`) | n-L | — | `occu_cover()` + `share(spatial(), residual = "full" | r)` in the positive formula (`R/occu_cover_residual.R`): cover arm carries `alpha*w + delta`. `"full"` = one latent/node, compiles to the SAME blocks as the #110 placed-bar spelling (`expect_identical`); what it adds is the IDENTIFIED alpha = projection of the cover-arm surface onto the occurrence surface, taken on the ETA scale (latents are unit-scale, amplitudes ride the outer grid -- adding them unscaled compares two units). `r` = r weighted `iid` blocks on a basis orthogonalized against the shared field (`.tobs_residual_basis`), latent K+r not 2K, orthogonality exact. Deviation SD PINNED: a free SD is one outer axis PER BLOCK. Forces `integration="grid"` -- the CCD counts pinned axes and builds a 2^(d-1) factorial. Costs one warm fit (reference + pin). **Rank REGULARIZES**: the sweep converges on the full-rank fit and that is a LOSS when the arms are confounded -- 12-seed table + the pinned-SD control separating rank from pin in `NOTES_measurements.md`. Gated: one per fit, not beside a placed cover-arm field, not with `alpha=0`, not with `|` MCAR, n-L only. `test-occu-cover-residual.R` |
| occu + cover independent detection-arm field (single-arm `to="detection"`) | n-L | — | `occu_cover()` + `spatial(~ 0 + time \|\| cell, graph, to="detection")` on the `detection` formula (or lifted `to="detection"`) (tulpa#140): spatially-structured detection prob. Same non-copied arm-specific block machinery as the cover-arm field (#110) -- `arm_field_blocks(af, "p")` sets slot 2 (`spatial_idx` psi=0/p=cell/pos=0). Enters via detection-arm `field_coef=1` (set when `det_armspec` present) so the block scatters onto the p rows; the shared occ field stays off detection by its `spatial_idx=0` sentinel -- identical to the #102 detection-RE mechanism, so NO tulpa engine change. Reports `sigma_p_field`/`sigma_p_field_<col>`. icar only; per-visit cover. `test-occu-cover-pos-field.R` |
| occu + cover obs-arm RE (detection / pos) | n-L | Yes | `occu_cover()` + RE on `detection=`/`positive=`; per-visit grouping. Intercepts: `(1\|g)`, crossed `(1\|g)+(1\|h)`, nested `(1\|g/h)` = N iid blocks (#102 single, #103 crossed/nested). Slopes (#103, tulpa>=0.0.39): uncorr `(x\|\|g)`/`(0+x\|g)` = weighted iid per coef; corr `(1+x\|g)` = miid free-Sigma block. `sigma_re_p`/`sigma_re_pos` (+`_<coef>`, `cor_re_p_*`) + BLUPs. **NUTS = INTERCEPTS, non-spatial (#205)**: each grouping = one shared `ReBlock` (`src/nuts_re_block.h`, reused not re-derived) trailing the field block, non-centered `b_g=sigma_re*z_g`, group SD SAMPLED (prior `N(0,1.5^2)` on `log sigma_re`, the observation-family width -- written ON the sampled coordinate, so NO change-of-variables term; asserted by reducing the target to the prior at `z=0`). Row-indexed, NOT site-indexed: one `re_group` code per `(site,visit)` row, 0 = padded/unseen level (engine's scatter sentinel). Designs from the SAME `.occu_cover_obs_re_parse`/`_design` as n-L (one dispatch-level call site). `fit$re` keys + `ranef()` match n-L; adds `sigma_median`/`sigma_draws`, `fit$nuts$re_sigma_rhat`. `fit$draws`/`vcov` stay the coefficient block (WAIC reads it -> scores RE at 0, as the sampled field already does). Off-block byte-identity `expect_identical`. Gated -> n-L: slopes, RE+areal field. `test-occu-cover-nuts-re.R` |
| occu + cover + detection / cover-arm RE | n-L | — | `occu_cover()` + `(1\|g)`/`re(g)` on `detection=` or `positive=` (#102); per-VISIT grouping iid block on that arm, composes w/ psi field; `sigma_re_p`/`sigma_re_pos` + BLUPs; `fit$re` per-arm list; `predict(type="detection")` adds BLUP offset, unseen->pop mean. Detection arm `field_coef=1` (not 0) so the iid block scatters, field still skipped via `spatial_idx=0` sentinel. pos RE needs `cover_aggregate="none"`. Intercept only; slope/correlated/non-spatial/NUTS gated |
| Community joint occu + cover | Yes | Yes | `ms_occu_cover()` — see below. positive = beta/lognormal/gaussian (`gaussian` #127, delta-normal `mu=eta`). **NUTS (non-spatial, #115 B7, `R/ms_occu_cover_nuts.R`, `src/ms_occu_cover_nuts.cpp`)**: joint-cover analogue of the `ms_occu`/`ms_int_occu` samplers -- three non-centered per-species arms (occ/p/pos), each w/ a log-Cholesky community covariance + ONE shared community log-dispersion scalar, over the per-(species,cell) two-state occu_cover marginal via in-tree C++ FullGradFn (WRAPS existing `occu_coupling_shared.h` cover-density kernels + `community_chol.h`, reuses `.ms_ocs_*` #128), warm-started at the Laplace-EM mode. Byte-exact vs R oracle, 0 divergences, REMOVES the Laplace community-variance attenuation (`var_attenuation` caveat). **Per-species dispersion RE** via `control=list(dispersion.re=TRUE)`: the shared cover log-dispersion becomes a 4th 1-D community arm `log_disp_s=mu_ld+sigma_ld*z_ld_s` (ms_abun log_r_s analogue), byte-exact vs the RE oracle, shared-disp path byte-identical (extracted per-cell `msoc_cell_sweep`); `fit$ms_dispersion$sigma_log_disp`/`log_disp_species`. Non-spatial lognormal/beta/gaussian; negbin N/A (continuous cover arm). `test-ms-occu-cover-nuts.R`. The shared-field spatial-factor variant stays lognormal-only |
| Spatial-factor community occu+cover (JSDM) | Yes | Yes | `ms_occu_cover()` + icar/car_proper/bym2 shared field, per-species loadings (tulpa#67). Laplace-EM (`R/ms_occu_cover_spatial.R`) + NUTS (`src/ms_occu_cover_spatial_nuts.cpp`). Cover-arm factor, `tobs_associations()`, per-species `predict()` maps |
| Multiscale occu + cover | Yes | Yes | `occu_multiscale_cover()` — 3-level cell/plot/visit + cover. THREE engines: `nested_laplace` (shared + SVC-trend coupled field) and `laplace` + `nuts` (both non-spatial, iid cells, field fixed at 0; `R/occu_multiscale_cover_nuts.R`, `src/occu_multiscale_cover_nuts.cpp`). Spatial = joint only. positive = beta/lognormal/gaussian (`gaussian` #127, threaded through all three engines via the shared `pos_*` code dispatch + `OccuMultiscaleCoverGaussianCoupling`) |
| Count / relative-abundance GLMM | Yes | — | `count(response=)` (spAbundance `abund`): GLMM on the observed response directly, NO detection, NO latent state = abundance analogue of `jsdm()`. Pois/negbin (log link) + gaussian (identity) + **binomial (logit, `k`-of-`n`, `trials=`; spOccupancy `svcPGBinom`, #125)**. ONE tulpa GLMM block (`build_count_callbacks`, `R/laplace_callbacks.R`); negbin size / gaussian variance by an outer dispersion loop in `.dispatch_count` (tulpa_laplace takes fixed phi), reported `fit$count_dispersion`; binomial has NO dispersion (variance pinned by `n`) -> takes the Poisson no-loop path + carries per-site `model$n_trials`. `simulate_count()`, `.tobs_ploglik_count` (WAIC), `count_methods.R` (fitted = expected successes `n*p`, predict newdata = per-trial prob, residuals binomial deviance/pearson). Community (msAbund) binomial = `ms_count("binomial")`; NUTS pending. Non-spatial binomial matches `glm()` MLE to ~1e-3. `test-count.R` |
| Community spatial-factor (`ms_count`) | n-L | — | `ms_count()` + `icar()` + `latent(n)` (spatial-factor `sfMsAbund`; #117): shared areal field AND latent factors on ONE formula, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}+sum_q lambda_{s,q} eta_{q,i}`. ONE block-coordinate loop runs BOTH field + factor updates; when both present the factor loadings are CENTERED across species (`sum_s lambda_{s,q}=0`) -> field owns the shared spatial mean, factors own between-species residual co-occ; both recover. UNIFIED fitter `.tobs_fit_ms_count_latent` (field-only/factor-only/both = one source of truth; old fitters = thin wrappers). `test-ms-count-factor.R` |
| Community latent factor (`ms_count`) | Yes | — | `ms_count()` + `latent(n)` (spAbundance `lfMsAbund`; #117): residual species co-occurrence via Q per-site latent factors + per-species loadings, `log mu_{s,i}=X_i(mu+b_s)+sum_q lambda_{s,q} eta_{q,i}`. BLOCK COORDINATE via the shared driver (`R/community_latent.R`): community EM w/ factor offset <-> factor update (alternating Newton eta\|lambda, lambda\|eta + unit-variance anchor). Loadings/factors identified up to ROTATION -> recoverable target = residual cov `Sigma_res=lambda lambda'` (`fit$ms_factor$residual_cov`/`residual_cor`/`loadings`/`factors`); residual cor recovers. `fitted()`/WAIC factor-aware via `model$count_factor_offset`. Poisson; not composed w/ areal field yet. `test-ms-count-factor.R` |
| Community SVC (`ms_count`) | n-L | — | `ms_count()` + `spatial(~ 1 + w \|\| cell, graph)` (spAbundance `svcMsAbund`; #117/#118): intercept field + one shared spatially-varying-coefficient field per covariate, `log mu_{s,i}=X_i(mu+b_s)+f0+w_i*f1`. Same block-coordinate scheme; the field solve generalizes to K covariate-weighted ICAR fields solved jointly (K x K sparse block system, per-field tau). `fit$spatial_field`=intercept field, `fit$trend_field(s)`=SVC field(s). Both recover. icar only. `test-ms-count-spatial.R` |
| Community count + areal (`ms_count`) | n-L | — | `ms_count()` + `icar()` (spAbundance `sfMsAbund`; #117): ONE shared areal field across species, `log mu_{s,i}=X_i(mu+b_s)+f_{u(i)}`. Fit by BLOCK COORDINATE ASCENT (`R/ms_count_spatial.R`) -- community Laplace-EM w/ the field as a per-site OFFSET (captured in the sp_ll closure, NO community_em.R change), alternated w/ a self-contained Poisson-ICAR field update (analytic sparse Newton `.ms_count_field_solve` + closed-form tau M-step). PURE R, no C++ (sidesteps the community EM FD-Hessian, which does not scale to an O(n_sites) field). Field informed by all species/site -> recovers cleanly; `fitted()`/WAIC field-aware via `model$count_field_offset`. Poisson + icar ONLY (overdispersed community count NOT identified vs a per-site field; bym2/car/group_var/bar follow-ups). `test-ms-count-spatial.R` |
| Community count (`ms_count`) | Yes | Yes | `ms_count(response=)` (spAbundance `msAbund`; #117): per-species GLMM on the observed response w/ Gaussian community hyperpriors on the coefficients, NO detection/latent state = community analogue of `count()`, abundance analogue of `ms_occu()`. Pois/negbin (per-species dispersion RE = second arm)/gaussian (per-species resid var, outer loop)/**binomial (logit `k`-of-`n`, `trials=`; community `svcPGBinom`, #125 -- laplace non-spatial ONLY, NUTS + shared field/latent gated w/ pointer; the community-mean intercept carries a small O(1/n_species) first-order-Laplace bias, absent at `trials=1`, slope unbiased)**. Reuses the shared community Laplace-EM (`.tobs_community_em`, `R/community_em.R`, PURE R) w/ count `sp_ll`/`sp_grad` (`R/ms_count.R`). `y` = `[n_sites x n_species]` or named list; `coef`=community means, `ranef`=per-species deviations, `simulate_ms_count()`, `.tobs_ploglik_ms_count` (WAIC). Gaussian/Pois community means unbiased; negbin slope carries mild first-order-Laplace attenuation (~10%, documented). **NUTS (Pois/negbin/gaussian, #117)**: family-aware in-tree C++ FullGradFn (`src/ms_count_nuts.cpp`, `R/ms_count_nuts.R`) = reduced `ms_abun` NUTS (no detection/latent-N), non-centered `b_{s,arm}=C_arm z_{s,arm}`, byte-exact vs R oracle, warm-started at the Laplace-EM mode, 0 divergences, NUTS==Laplace. negbin adds a per-species dispersion RE `log_r_s~N(mu_log_r,sigma_log_r^2)` as a SECOND community arm (block-diag chol, mirrors `ms_abun`) -> removes the Laplace attenuation; gaussian adds S FREE per-species `log_phi_s` (no community prior, weakly-informative N(pooled log-var, 2^2), matching the Laplace outer loop); `fit$ms_dispersion` carries `r_s`/`variance`. Missing (`NA`) site x species entries are dropped from the per-(species,site) data sum (NA-in-`y` IS the mask; C++ + R oracle skip identically -> byte-exact w/ NA present), matching the Laplace path's per-species `valid` subsets. `test-ms-count.R` + `test-ms-count-nuts.R` |
| Count + areal / SPDE / GP | n-L | — | `count()` + `icar()`/`car_proper()`/`bym2()`/`spde()`/`gp()` (spAbund; #117, bym2+spde+gp #116): areal OR continuous-mesh/NNGP field on the abundance formula. Response observed (no latent state) -> fit = ONE `tulpa_nested_laplace()` call over the count block w/ the field as its GMRF prior, NOT the occupancy EM. Grid-integrated FE (law of total cov via per-cell `keep_grid_hessians`) + field/sigma via shared `.tobs_nested_attach_field_summary`; dedicated `.tobs_fit_count_spatial` (`R/count_spatial.R`), reuses `.tobs_to_multi_block_prior`; `fitted()` field-aware in-sample. **bym2 (#116)**: the generic field summary (`R/laplace.R`) reads the `b<n>.sigma`/`b<n>.rho` grid axes + the (phi\|theta) mode slices to reconstruct `z=sqrt(rho/scale)*phi+sqrt(1-rho)*theta` (Riebler 2016), grid-averaged. SAME summary drives occu() etc -> a bym2 areal field there is reconstructed too (was dropped/NULL); icar/car_proper byte-identical. **spde (#116)**: the Matern field lives on the mesh nodes (`fit$spatial_field`, length `n_mesh`) and the barycentric projector `fit$spatial$tulpa_spec$A` maps it to sites; `.tobs_spatial_field_offset` has an A-projection branch so `fitted()`/`predict()` add the projected per-site field (areal fits carry no projector -> branch skipped, byte-identical). Gated `skip_if_no_tulpamesh()`. **gp (#116)**: `gp(lon,lat,prior_range=c(r0,alpha))` routes to tulpa's single-block `nngp` kernel via `.tobs_fit_count_gp`: marginal variance + range integrated on the kernel's OWN outer grid, field Schur-folded out -> grid-integrated FE + `fit$gp_hyper`, but NO per-cell map (`fit$spatial_field=NULL`; use spde() for a reconstructed field). `multiscale_gp()` NOT hosted by nested-Laplace -> errors w/ pointer. **Poisson OR binomial** -- Poisson has no dispersion; binomial variance pinned by `n` -> identified against a per-node field (`svcPGBinom`, #125; works at `trials=1` too, weak but NOT confounded). negbin size / gaussian residual variance NOT jointly identified w/ the field under the fixed-phi loop (degenerate: size->Inf, resid var->0) -> gated w/ pointer; bym2 + improper `car()` gated to icar/car_proper. `test-count-spatial.R` |
| N-mixture (`abun`) | Yes | Yes | `abun()` = laplace + nuts (non-spatial, #41) + nested_laplace (areal). Royle 2004; see the detail section below |
| N-mixture + areal | n-L | Yes | `abun()`+icar/bym2/car_proper; Pois/NB (r grid-int); grid-int cov (constrained intercept). NUTS+areal = car_proper ONLY (#51, fixed-hyper non-centered field); icar/bym2 NUTS+spatial gated. Field on the ABUNDANCE arm only: a `detection=~icar()` term is REJECTED on BOTH methods (#255), where the four observation families route it (#114) |
| Community N-mixture | Yes | Yes | `ms_abun()` (msNMix); per-species coef RE, in-tree Laplace-EM (`nmix_laplace_re`) -> `NMixCommunityOracle` AGHQ, Schur SEs; Pois + negbin. NUTS (#14) |
| Community N-mixture + areal | n-L | Yes | `ms_abun()`+icar/bym2/car_proper (sfMsNMix; #12); shared field + per-species RE; `nmix_community_spatial.cpp`; Pois/NB. NUTS (#73, car_proper only) = #14 non-centered community sampler + a SHARED fixed-hyper non-centered proper-CAR field on abundance (tau Q(rho) fixed at the #12 nested-Laplace estimate, raw ~ N(0,I), f=Linv raw; optional field block in `src/ms_abun_nuts.cpp`, FD-validated, field-off byte-identical to #14). 0 divergences, field cor ~0.97; Pois only. icar/bym2 also sample + centre via the #71 sum-to-zero shared-field loading (#113). Field on the abundance arm only; a `detection=~icar()` term is REJECTED on both methods (#255) |
| Community N-mixture latent factor (`ms_abun`) | Yes | — | `ms_abun()` + `latent(n)` (spAbundance `lfMsNMix`; #117): residual species co-occurrence on the ABUNDANCE arm via Q per-site factors + per-species loadings. Latent N still marginalises closed-form per species-site -> the whole latent structure sits on `eta_lambda`, and the family reduces to ONE working oracle over the Royle marginal (`R/ms_abun_latent.R`) driving the shared block-coordinate engine: `score = grad_eta_lambda`, `curv = info_eta_lambda - var_N * score_wt_lambda^2` = the Louis (1982) (1,1) block (abundance curvature, detection arm profiled out), which `nmix_site_marginal()` ALREADY exposes -> no new kernel. Oracle FD-validated; residual cor recovers well (counts information-rich, cleaner than occupancy). Poisson only (a negbin size = a second per-site dispersion, not identified vs a per-site latent). `test-ms-abun-factor.R` |
| Community N-mixture spatial-factor (`ms_abun`) | n-L | — | `ms_abun()` + `icar()`/`car_proper()`/`bym2()`/`spde()` + `latent(n)` (#117): shared field AND factors on the abundance arm; the centred loadings separate them. Same block-coordinate driver. A plain field with NO factors KEEPS the dedicated C++ path (#12) -- faster + already recovery-tested, deliberately not replaced |
| N-mixture + grouped RE | Yes | — | `abun()`+`(1\|g)`/`(x\|g)` either arm (#13); non-species grouping; Pois/NB; `NMixGroupedOracle`. Gated: RE+spatial, RE+visit-det, RE both arms |
| Community distance (`ms_distance`) | Yes | — | `ms_distance(key=, transect=, cutpoints=)` (spAbundance `msDS`; #117): per-species binned distance sampling w/ Gaussian community hyperpriors on the abundance (`log lambda`) + detection-scale (`log sigma`) coefs. `y` = `[n_sites x n_bins x n_species]` or named list. NO new C++ -- latent N marginalises closed-form per species-site and `cpp_distance_site_sweep` already returns `log_lik`/`grad_lam`/`info_lam`/`var_N`/`swl`, so the shared community Laplace-EM (`R/ms_distance.R`) reads its per-species score straight off it. Hazard-key log-shape = a community `global` (shared across species) via the EM `global` slot. Community means unbiased, Wald coverage nominal; single-species `distance()` control agrees. NB: a single seed shows a ~0.18 paired lambda-down/sigma-up shift that LOOKS like bias and is NOT -- one draw on the lambda/sigma ridge. Poisson only (negbin size not yet a per-species RE); NUTS not wired. `simulate_ms_distance()` draws through `cpp_simulate_distance` (a separate R-side quadrature would simulate from a pi the model is not fit against). `test-ms-distance.R` |
| Community distance latent factor (`ms_distance`) | Yes | — | `ms_distance()` + `latent(n)` (`lfMsDS`, #117). Factors ALONE = plain block-coordinate Laplace-EM, so `method="laplace"`; `.dispatch_ms_distance` REJECTS `nested_laplace` without a field ("needs a shared field ... use method = \"laplace\""). Same shared driver + oracle as the spatial-factor row below |
| Community distance spatial-factor (`ms_distance`) | n-L | — | `ms_distance()` + `icar()`/`car_proper()`/`bym2()`/`spde()` (`sfMsDS`, #117), optionally + `latent(n)`. Shared field REQUIRES `method="nested_laplace"` (`.dispatch_ms_distance` errors under `laplace`) -- same gate as the sibling `ms_count`/`ms_abun` rows. Field on the ABUNDANCE arm only; a `spatial()`/`latent()` term on the detection formula errors. Working oracle = the SAME Louis formula as `ms_abun` -- `score = grad_lam`, `curv = info_lam - var_N * swl^2` (both count marginals w/ `B_i = diag(info) - var_N v v'`; only the kernel supplying the pieces differs), FD-validated. Poisson only; temporal/re/svc not wired. **Per-species Hessian IS assembled (#161)**: `.tobs_ms_distance_info_block()` builds the per-site Louis block `B_i = diag(info_lam, info_sig) - Var(N_i|y) v v'` from what `cpp_distance_site_sweep` already returns; distance has ONE unit per site -> three `crossprod`s against diagonal weights, no per-site loop (N-mixture needs one only b/c its visit count varies). Fit UNCHANGED, community EM markedly faster -- the FD it replaced cost `2U` latent-N sweeps per species per Newton step and made `test-ms-distance.R` the most expensive file in the recovery tier. **Sign inside `v` is NOT the sibling's**: `v = (-swl, -vN_sig)` -- the distance kernel stores `vN_d` already negated (`-p_k/(1-p)`) where `nmix_site_marginal()`'s `v` carries `+p`. Invisible on the diagonal, wrong by ~2x on the cross block -> asserted against the finite difference, NOT reasoned across (`test-ms-distance-info.R` also asserts the cross block >10% of the diagonal so the agreement cannot be vacuous, and that the flipped convention is REJECTED). Half-normal only: under the hazard key `grad_b`/`info_b` come back already summed over sites + per-site cross terms are not exported, so that key keeps the FD fallback |
| Removal (Pois/NB) | Yes | Yes | `removal()` (#39) — see Architecture. **continuous NNGP `svc()` varying coefficients on the abundance arm under laplace/nested_laplace, composing with the areal/temporal blocks (#144)**; NUTS single intercept RE; Laplace grouped RE one arm; areal icar/car_proper/bym2 abundance arm; **detection-arm areal field on the capture logit (`detection=~icar()`, #114)**; areal+temporal AND temporal-only AR1/RW1/RW2/iid on abundance arm via shared areal-BFGS (#78/#114); NUTS+areal icar/car_proper/bym2 field on abundance arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| Distance (Pois/NB) | Yes | Yes | `distance(key=, transect=, cutpoints=)` (#38); `formula`=log lambda, `detection`=log sigma, `y`=`n_sites x n_bins` — see Architecture. **continuous NNGP `svc()` varying coefficients on the abundance arm under laplace/nested_laplace (#144)**; NUTS single abundance intercept RE; Laplace grouped RE abundance arm (half-normal AND hazard key -- hazard log-shape profiled over the AGHQ log-marginal, #114); areal field (half-normal + hazard key); areal+temporal AND temporal-only on abundance arm (#78/#114); DETECTION-arm areal field on log sigma (`detection=~icar()`, spatially-varying detection scale, half-normal AND hazard key, #114); NUTS+areal icar/car_proper/bym2 field on abundance arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| False-positive occupancy | Yes | Yes | `fp_occu()` (#40) — see Architecture. **continuous NNGP `svc()` varying coefficients on the psi arm under laplace/nested_laplace (#144)**; NUTS single psi intercept RE; Laplace grouped RE psi OR p11; areal psi arm; **detection-arm areal field on the p11 logit (`detection=~icar()`, #114)**; areal+temporal AND temporal-only on psi arm (#78/#114); NUTS+areal icar/car_proper/bym2 field on psi arm (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113) |
| Presence + nominal class | Yes | — | `occu_categorical()` (#106, `R/occu_categorical.R`): categorical analogue of the `cover()` hurdle -- presence Bernoulli + an UNORDERED class (multinomial logit / softmax), not a continuous magnitude. Arms factorise exactly: `P(y=0)=1-psi`, `P(y=k)=psi*p_k`, `p=softmax(X beta_class)`. Multinomial math = FD-validated tulpa kernel `cpp_multinomial_logit_terms` (`src/multinomial_logit.h`); non-spatial Laplace fit = vectorised R Newton over the same closed forms. laplace ONLY (native multi-process LikelihoodSpec path for spatial fields / NUTS = documented follow-up). `simulate_occu_categorical`, S3 (`print`/`coef`/`predict`). `test-occu-categorical.R` |
| Royle-Nichols occupancy | Yes | — | `royle_nichols()` (#116, unmarked `occuRN`): abundance-induced detection heterogeneity, `p_site=1-(1-r)^N`, `N~Pois(lambda)`. Pure-R closed-form Poisson-sum marginal over per-site suff-stats (k,n), optim BFGS + observed-Fisher vcov (`R/royle_nichols.R`). Site-level detection; laplace only. `simulate_royle_nichols`, WAIC. `test-royle-nichols.R` |
| Time-to-detection occupancy | Yes | — | `occu_ttd()` (#116, unmarked `occuTTD`): survey records TIME to first detection; `z~Bern(psi)`, `t\|z=1~Exp(rate lambda)` censored at `surveyLength`. Pure-R two-state marginal (`.ttd_site_loglik` over suff-stats k/S_det/Tcens), optim BFGS + observed-Fisher vcov (`R/occu_ttd.R`); same recipe as royle_nichols w/ a continuous censored-exponential emission. `y`=`[N x J]` TTD matrix (val<Tmax=detection, >=Tmax=censored, NA=not surveyed). Site-level rate; exponential only (Weibull shape + visit-varying rate + areal follow-ups). `simulate_occu_ttd`, full S3 + WAIC. `test-occu-ttd.R` |
| Double-observer abundance | Yes | — | `double_observer(type=)` (#116, unmarked `multinomPois` double-observer pi): `N~Pois(lambda)`. **`type="independent"`** (default): two independent observers `p1`/`p2`; Poisson-multinomial thinning makes the 3 observable cells (obs1-only/obs2-only/both) INDEPENDENT Poissons `n_c~Pois(lambda*pi_c)` -> closed-form marginal, NO latent-N sum. `y`=`[N x 3]` cell counts. **`type="dependent"`** (removal-style): PRIMARY observer records what it detects, SECONDARY only the primary's misses -> two cells `n1~Pois(lambda*ppri)`, `n2~Pois(lambda*(1-ppri)*psec)` (`.dobs_dep_site_loglik`). A single fixed primary = 2 cells for 3 params = a ridge, NOT identified; observer ROLE-SWAPPING (per-site `primary` in {1,2}) -> 4 cell means, identifies all three. `y`=`[N x 2]` (primary-detected, secondary-only) + `primary=`; single-observer `primary` warns. optim BFGS + observed-Fisher (`R/double_observer.R`). `formula`=log lambda, `detection`=shared design w/ separate p1/p2 arms; `fit$means` names `lambda_*`/`p1_*`/`p2_*`. Both types recover lambda + both detection probs. `simulate_double_observer(type=)`, full S3 + WAIC. `test-double-observer.R` |
| Joint distance + removal | Yes | — | `gdistremoval()` (#116, unmarked `gdistremoval`, Amundson 2014): SINGLE-SEASON (NOT the open-population `distsampOpen` HMM). `N~Pois(lambda)`; detected birds cross-classified by a distance band AND a removal period. Total detected = binomial thinning of N, Pois closed under thinning -> `ysum~Pois(lambda*pdist*prem)` CLOSED-FORM (no latent-N sum); band + period counts = two conditional multinomials (the `double_observer` Poisson-multinomial pattern). Half-normal band integrals CLOSED-FORM line + point (`.gdr_dist_cp`, == integrate to machine precision); depleting-removal `pi_k=r(1-r)^(k-1)` (`.gdr_rem_cp`). 3 site-level arms: log lambda (`formula`), log sigma (`detection`), logit r (`removal=~`). optim BFGS + observed-Fisher (`R/gdistremoval.R`, `.gdr_site_loglik`). `y`=`[N x Jdist]` band counts, `y_rem`=`[N x Jrem]` period counts (rowSums must match). Recovery clean, bias small, coverage nominal; closed form validated == unmarked brute-force K-sum, diff 0. `simulate_gdistremoval`, full S3 + WAIC, `test-gdistremoval.R`. v1 = halfnorm key, line/point, Pois, constant r, availability phi FIXED 1 (a single period does not ID phi); hazard key / NB-ZIP / multi-period phi arm = follow-ups (all closed-form under the same thinning) |
| Open-population distance sampling | Yes | — | `distsamp_open()` (#116, unmarked `distsampOpen`, Dail-Madsen + distance): open-population counterpart of `gdistremoval`. `N_1~Pois(lambda)`, `N_t=Binom(N_{t-1},omega)+Pois(gamma)`; each primary period the detected birds are distance-sampled into bins. Band allocation conditional on the period total factors OUT of the abundance HMM (the gdistremoval trick) -> marginal REUSES the validated `cpp_dyn_abun_total_log_lik` kernel fed **eta_p=logit(pdist)** + y=per-period detected totals (J=1) + per-period band multinomials. NO new HMM kernel. 4 site-level arms: log lambda (`formula`), log sigma (`detection`), logit omega (`omega=~`), log gamma (`gamma=~`). optim BFGS + ANALYTIC gradient (`.dso_grad`: `dyn_abun` returns grad_eta_*, the sigma block chains through pdist via a distance-integral FD; ONE kernel call/eval; FD-validated) + observed-Fisher (`R/distsamp_open.R`). `y`=`[n_sites x n_bins x n_seasons]`. K_max default = DETECTION-CORRECTED (`max(ntot)/pdist + infl*sqrt + 10`, infl=4 Pois / 5 NB-ZI, NOT 3*max -- the forward is cubic in K). **`mixture=`**: `negbin` threads the kernel's `use_nb`/`eta_logr` (trailing `log_r` coord, analytic `grad_eta_logr` FD-validated); `zip`/`zinb` = pure-R additive layer `.tobs_fit_distsamp_open_zip` over the composed per-site marginal (HMM `log_lik_site` + band multinomials), `omega*1{all bands 0}+(1-omega)*L_open`, per-site eta/sigma grads weighted by the structural-zero posterior w_i + ZI-logit score `(1-om)-w_i`, log_r via central FD; ZI coord `zi_logit`/`zi_omega` (NOT survival `omega`). **`dynamics=`**: `constant` (default) / `notrend` / `trend` / `autoreg` / `ricker` / `gompertz` = unmarked tp1..tp5; `.dso_dyn_meta(dynamics)` sets the active-arm layout + kernel code. `constant`/`notrend` keep the analytic-gradient path (notrend ties `gamma=(1-omega)*lambda`); the density-dependent forms use the value-only forward kernel `compute_dyn_abun_site_dyn` (`src/dyn_abun_kernel.h`, `int dynamics` 2..5, brute-force validated) exported as `cpp_dyn_abun_dynamics_log_lik`, driven by numeric-gradient BFGS (`.tobs_fit_distsamp_open_dyn`). `ricker`/`gompertz` estimate a carrying capacity `K` (log link, `omega=~` slot) + growth `r` (`gamma=~` slot); `trend` drops survival; `autoreg` gamma is per-capita; the density-regulated forms grow toward K so bound `K_max`. lambda/sigma/zi recover on every mixture + dynamics; omega/gamma/NB-size/K/r sit on the short-series ridge. `simulate_distsamp_open(mixture=,size=,zi=)`, full S3 + WAIC, `fitted()`/`predict()`/`simulate()` dynamics-aware (`.dso_draw_dyn`), per-param coverage nominal (`test-distsamp-open.R`). v1 = halfnorm key, line/point, site-level arms; season-varying detection, NUTS = follow-ups |
| Multi-species co-occurrence occupancy | Yes | — | `occu_multi()` (#116, unmarked `occuMulti`, Rota 2016): joint state `z in {0,1}^S` from a log-linear model w/ first-order (per species) + second-order (per pair = the INTERACTION) natural params; conditional per-species detection. Pure-R marginal enumerating the `2^S` states (`.occu_multi_site_loglik`), optim BFGS + observed-Fisher vcov (`R/occu_multi.R`). `y`=list of S `[N x J]` matrices or 3D array; `species=`. Shared covariate design (separate coefs per natural param), site-level detection, laplace only (per-param formulas, higher order, visit-level det, areal = follow-ups). NB: log-linear natural-param SLOPES trade off (weakly ID'd individually); MARGINAL psi + the interaction recover cleanly -> recovery targets those. `fit$means` names `f_<sp>`/`f_<sp>_<sp>`/`p_<sp>`; `fitted()$psi`=marginal occupancy. `simulate_occu_multi`, full S3 + WAIC. `test-occu-multi.R` |
| Open N-mixture (Dail-Madsen) | Yes | Yes | `dyn_abun()` (#37); y 3D `[n_sites x J x T]` — see Architecture. Pois/NB init; **ZIP/ZINB (`mixture="zip"/"zinb"`, #116)**: a structural-zero site is never occupied in any season -> all its counts 0, per-site marginal `omega*1{all y=0}+(1-omega)*L_DailMadsen` (the abun ZIP additive layer over the forward-HMM marginal); PURE-R `.tobs_fit_dyn_abun_zip`, analytic-gradient BFGS, FD-Jacobian observed-info vcov; ZI coord named `zi_logit` (NOT `omega_*` = the SURVIVAL arm); `simulate_dyn_abun(zi=)`; non-spatial laplace intercept-only omega, field/RE/NUTS stay Pois/NB w/ pointer (`test-dyn-abun-zip.R`). Season-varying omega/gamma via a `[n_sites x (T-1)]` covariate, interval-indexed forward kernel, all backends (#80). **continuous NNGP `svc()` on the initial-abundance arm under laplace/nested_laplace (#144)**; NUTS single intercept RE on init-abundance OR detection arm; Laplace grouped RE on either (#82, p-arm = per-node full-forward second-order eta_p pass); areal init arm; **detection-arm areal field (`detection=~icar()`, #114)**; areal+temporal AND temporal-only on the init-abundance arm (#78/#114); NUTS+areal icar/car_proper/bym2 on init-abundance (#72; intrinsic icar/bym2 via the #71 sum-to-zero reparam, #113); NUTS+temporal-only fixed-hyper ar1/rw1/rw2/iid on init-abundance (#114, same field block as areal, `field_map`=period; areal+temporal under NUTS gated to n-L) |
| Spatial ICAR/BYM2/NNGP | — | Yes | |
| Spatial + dynamic | — | Yes | |
| Nested-Laplace (areal) | n-L | — | icar/bym2/car (+temporal/iid) on occu/int_occu/dyn_occu. RECOVERY-TESTED on icar for all three (`test-{occu,int-occu,dyn-occu}-areal-recovery.R`) + SVC bar (`test-dyn-occu-svc-recovery.R`). State arm encodes at **M = 1** for ANY nested latent block (`nested_block <- !is.null(latent_prior)` in all three `build_*_callbacks`): a state row = ONE binary occupancy obs, so the old M=1000 (areal) / M=4 (spde) overstated its info M-fold, swamped the field prior, read between-cell noise as field, inflated the slope via `sqrt(1+0.346 sigma^2)`. Monotone in M, no arm regresses at M=1. `use_louis` (laplace_helpers.R) covers single/dynamic/integrated -- the state score `x_i(z_i-psi_i)` is family-generic; without it the state block falls to `.se_from_laplace_fit()` -> **NA SEs** (dynamic + integrated both shipped that way). dynamic needs `weights[,1]` (season-1 col), integrated's is already per-site. **A null-field fixture CANNOT test this path**: the grid `b1.tau` is 9 log cells [0.3,30] (sigma [1.83,0.18]), so truth sigma=0 is OFF-grid and the marginal pins on the last cell -- "shrunk" and "pinned" read identical. Test with an INTERIOR field + assert the peak cell off both boundaries. `tau` = ICAR CONDITIONAL precision, NOT 1/sd^2 (marginal field sd 1.0 -> sigma~0.40) -> assert on field sd/cor, never on `sigma`. `dev_notes/finding_dyn_nested_laplace_field.md` |
| NA-response prediction | n-L | — | `predict(type="state")`: all-NA single-season sites field-interpolated, calibrated from an exact-marginal bernoulli pass |
| Formula RE (intercept) | Yes | Yes | `(1\|g)`; variance-component EM, occ OR det arm (#11) |
| Formula RE (uncorr slope) | Yes | Yes | `(x\|\|g)`, `(0+x\|g)`, `(1+x\|\|g)`; either arm |
| Formula RE (corr slope) | Yes | Yes | `(1+x\|g)`; EM M-step consumes tulpa `cov_blocks`; either arm (#11) |
| Formula RE on detection | Yes | Yes | `detection=~(1\|g)`; separate det-arm RE block, AGHQ branches on arm |
| All S3 methods | Yes | Yes | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, `$.tobs_fit` |
| Diagnostics | Yes | Yes | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatial_range |
| Simulation | Yes | Yes | `simulate_occu/ms_occu/dyn_occu/int_occu/dyn_ms_occu/int_ms_occu/cover/cover_joint` |
| Spatial prediction | — | Yes | `tobs_predict_spatial` (IDW on field) |
| Components | Yes | Yes | `tobs_re`, `tobs_temporal`, `tobs_svc`, `tobs_latent`, `tobs_community_re`, `tobs_areal` |
