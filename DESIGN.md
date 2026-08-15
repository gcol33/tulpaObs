# DESIGN.md

Where a new thing plugs in. Read this BEFORE adding a family, a prior, a backend,
or a diagnostic.

Three repo docs, three jobs. Do not mix them:

- `CLAUDE.md` / `AGENTS.md` (byte-identical) = the ROSTER. What each family
  supports, per backend. 1261 lines, at its stated 150k budget.
- `NOTES_measurements.md` = the NUMBERS the roster rests on. Fixtures, seeds,
  wall times, per-seed recovery.
- `DESIGN.md` (this file) = the EXTENSION POINTS. Which single object to edit,
  and which duplication is deliberate.

This file stays short so it gets read in full. It carries no family roster and no
measurements. If you are about to add either here, it belongs in one of the other
two.

## The rule

**Before writing a second copy of anything, find the registry.** Every axis this
package extends along already has one. The failure mode this file exists to stop
is copying the nearest sibling family instead of the factored path, which is how
every open issue in the repo got filed.

| Adding | Edit this, and nothing else | Contract lives at |
|---|---|---|
| structured term (spatial / re / temporal / svc / latent / copy) | `.tobs_terms` | `R/formula_terms.R:1017` |
| family x backend support | `.tobs_family_methods` | `R/tobs_helpers.R:146`, validated `:427` |
| sampler default (`n.iter`, `adapt.delta`, `sigma.beta`, ...) | `.TOBS_ENGINE_DEFAULTS` | `R/engine_defaults.R:25`, filled `:237` / `:251` |
| SBC coverage for a family | `.TOBS_SBC_REGISTRY` | `R/sbc.R:2788`, single-species shortcut `:625` |

Each of those already states its own contract at the definition. `.TOBS_SBC_REGISTRY`
(`R/sbc.R:2783-2785`) is the clearest: *"A family is one row. Everything not named
here ... is the shared driver's."* That sentence is the standard for all four.

## Shared drivers

A new family reaches these by calling them, never by growing a private copy.

R side:

- `.tobs_community_em()` `R/community_em.R:148` -- the community Laplace-EM behind
  every `ms_*` family.
- `.tobs_community_latent_ascent()` `R/community_latent.R:1131` -- block-coordinate
  latent structure (areal field / factors) for every community family. A family
  supplies ONE `working(eta) -> list(score, curv)` callback.
- `.tobs_areal_bfgs_fit()` `R/areal_bfgs.R:466` -- areal / temporal / NNGP field
  blocks over any family exposing `eval(theta, offset)`.
- `.tobs_pg_draw_beta()` `R/pg_gibbs_shared.R:25` -- Polya-Gamma machinery behind
  every `pg_gibbs` fitter.
- `.tobs_nuts_attach_convergence()` `R/nuts_chains.R:389` -- the single writer of
  the Rhat / ESS record `summary()` and `print()` read. A NUTS path that skips it
  reports no convergence at all.
- `compute_bym2_scale()` `R/spatial.R:62` -- the Riebler scale factor.
- `.tobs_svc_columns()` `R/occu_svc.R:46` -- resolves `svc()` columns for BOTH
  backends, so Laplace and NUTS cannot disagree on which coefficients vary.
- `.occu_cover_visit_view()` `R/occu_cover_diag.R:420` -- the one length-V visit
  view every per-visit `occu_cover` diagnostic reads. Reading `model$y` /
  `model$valid` directly is what broke `cpo()` / `ppc()` on compact fits.

C++ side:

- `run_tulpa_nuts` `src/nuts_engine.h` -- the sampler driver behind every in-tree
  FullGradFn target.
- `community_chol_pri_read` `src/community_chol.h` -- log-Cholesky hyperprior
  scalars, one declaration for all seven community NUTS targets.
- `community_pack_grid` `src/community_grid_pack.h` -- per-outer-grid-point pack
  shared by the community areal drivers.
- `make_arms` `src/occu_cover_ragged.h` -- the predictor view every ragged
  `occu_cover` diagnostic kernel assembles from.

## Duplication that is deliberate

Do NOT "fix" these. Each was checked and each fails loudly if it drifts.

- **`[[Rcpp::export]]` entry points** (`cpp_*_nuts`, `cpp_*_joint_logpost`, ~19
  across `src/`). Unpack spec, size-check, call the shared driver. Rcpp needs one
  concrete signature per export. There is nothing to extract.
- **R/C++ mirrored `*_nuts_layout` pairs.** The same offset arithmetic in
  `R/ms_*_nuts.R` and `src/ms_*_nuts.cpp` is deliberate: each family's
  byte-exactness oracle test compares the R logpost against its C++ twin, so a
  wrong offset fails the suite rather than biasing a fit. (The R copies AMONG
  THEMSELVES are real duplication, tracked as #231.)
- **Per-family `.dispatch_*` bodies** `R/tobs_dispatch.R`. They look alike and the
  gates genuinely differ per family.
- **`max.iter` / `tol` at their call sites.** These are per-ROUTE Laplace-EM
  values, NOT sampler knobs, so they deliberately do not live in
  `.TOBS_ENGINE_DEFAULTS`. See that file's scope note.

## The tell

Real duplication in this repo has one reliable signature: **a sibling already
solved it.** The target shape is in the tree, and one caller did not reach it.

Worked examples, all four currently open:

- #228 -- `compute_bym2_scale()` `R/spatial.R:62` exists; `R/ms_occu_spatial.R:217`
  and an inline block at `R/nmix_laplace_spatial.R:452` reimplement it, with a
  different eigenvalue filter.
- #229 -- `.removal_spatial_prep()` `R/removal_spatial.R:38` exists; its twin
  `nmix_laplace_spatial.R` inlines the same preamble in every fitter. The defaults
  have already drifted (`tau_grid` 9 vs 7, two `rho_grid` spellings).
- #230 -- `.tobs_sbc_simple_entry()` `R/sbc.R:625` reduces 9 single-species
  families to one row each; the 7 community families hand-roll it.
- #231 -- `.ms_ocs_arm()` / `.ms_ocs_chol_dim()`
  `R/ms_occu_cover_spatial_nuts.R:506` / `:49` exist, and the consumer at `:512`
  states the intent outright: *"adding an arm is a row in that declaration rather
  than a private copy of this loop."* Four of six layouts still spell the
  triangular count inline.

The inverse also holds: if no sibling solved it, the shape is probably load-bearing.
Read two members of the group and ask whether the shared part is already behind a
driver before concluding anything.

## What a registry buys

Recorded, so the next refactor does not have to re-argue it:

- #226 (community `mu` and `b_s` draw jointly) is a few lines inside one helper.
  Landing it meant editing all seven `draws_ms_*` copies. All seven do call it
  today. The structure offers no guarantee the seventh gets the next one.
- The `nmix` / `removal` areal defaults drifted precisely because each fitter owns
  its own. From the code it is not visible which value was chosen and which was
  inherited.

## Adding a family, in order

1. Constructor in `R/obs_families.R`, S3 class `tobs_*`.
2. Row in `.tobs_family_methods` (`R/tobs_helpers.R:146`). No row means dispatch
   rejects it, which is the intended failure.
3. `.dispatch_<family>()` in `R/tobs_dispatch.R`, resolving controls via
   `.tobs_control_defaults()`.
4. Reach the shared drivers above. A private copy of one is a review finding.
5. Sampler knobs as `NULL` formals plus `.tobs_fill_sampler(environment(), engine)`
   as the first statement. `test-engine-defaults.R` asserts structurally that no
   fitter carries a literal.
6. NUTS path: attach convergence via `.tobs_nuts_attach_convergence()`.
   `test-nuts-convergence-contract.R` fails on a family that advertises `nuts`
   without a case.
7. Row in `.TOBS_SBC_REGISTRY` (`R/sbc.R:2788`); single-species families take
   `.tobs_sbc_simple_entry()`.
