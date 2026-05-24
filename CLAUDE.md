# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## tulpaObs — Hierarchical Latent-State Observation Models via tulpa

Bayesian occupancy / abundance / distance / removal / cover models built on the
[`tulpa`](https://github.com/gcol33/tulpa) inference engine. R package, C++17
backend via Rcpp/RcppEigen, depends on a sibling checkout of `tulpa` at
`../tulpa` (via `LinkingTo: tulpa`).

> **Public API:** `tobs()` + family constructors (`occu()`, `dyn_occu()`,
> `ms_occu()`, `int_occu()`, `jsdm()`, `abun()`, `ms_abun()`, `dyn_abun()`,
> `distance()`, `removal()`, `fp_occu()`, `cover()`). All S3 classes are
> `tobs_*` (`tobs_fit`, `tobs_model`, `tobs_family`, `tobs_spatial`,
> `tobs_temporal`, `tobs_re`, `tobs_svc`, `tobs_latent`, `tobs_priors`).
>
> **Structured terms live inside the formula** (lme4 / mgcv / INLA style),
> not as `tobs()` arguments. The registry in `R/formula_terms.R` maps a term
> name to a constructor: spatial `icar()`, `bym2()`, `car()`, `car_proper()`,
> `gp()`, `multiscale_gp()`, `spde()`; `re()`; `temporal()`; `svc()`;
> `latent()`; and `copy("id")` to share one realization across the state and
> detection predictors. The term constructors (`.tobs_term_*`) and the AST
> parser (`.tobs_parse_formula`, `.tobs_bind_formulas`) are internal — there
> are no exported `tobs_icar()` / `tobs_re()` / `tobs_spde()` constructors and
> no `spatial = ` / `temporal = ` / `re = ` arguments on `tobs()`. A term's
> process membership (which linear predictor it enters) is determined by which
> process formula it appears in; the fitter derives the `shared = c(occ, det)`
> vector from that via `.tobs_structures_from_model()`.
>
> `lme4` bar syntax is supported as sugar over `re()`: `.tobs_desugar_bars()`
> (in `R/formula_parse.R`) rewrites bar terms into the equivalent `re()`
> calls on the formula AST *before* `terms()` runs, so there is one parser and
> one term type. Supported forms (gcol33/tulpaObs#10): `(1 | g)`, `(x | g)`,
> `(1 + x + z | g)` (multi-slope, stacked via `cbind()`), `(x || g)`
> (uncorrelated), `(0 + x | g)` (slope-only, no group intercept), and crossed
> / nested grouping `(1 | g:h)`, `(1 | g/h)` (one `re()` per implied grouping
> factor). Multi-slope and nested are pure R-side (`build_re_spec()` computes
> `n_coefs` from the slope-matrix width); slope-only needs the per-term
> `re_has_intercept` flag in tulpa's engine (ABI 22).
>
> Random effects fit under both engines (gcol33/tulpaObs#11). NUTS fits every
> form (incl. correlated slopes). The default `engine = "laplace"` fits iid
> intercept RE and *uncorrelated* slopes (`(x || g)`, `(0 + x | g)`,
> `(1 + x || g)`) on the occupancy predictor of a single-season model via a
> variance-component EM (`R/em_laplace_re.R`, `.tobs_em_laplace_re()`): it wraps
> tulpa's fixed-sigma `tulpa_laplace()` in the occupancy missing-data EM (the RE
> mode is fed back into psi) and an EM/REML variance-component update, with the
> RE posterior covariance recomputed at natural scale via Schur. tulpa's Laplace
> uses a *diagonal* RE precision, so correlated slopes (`(1 + x | g)`), RE on
> detection, RE + spatial, RE + visit-level detection, and RE on non-single
> families error from `.validate_re_laplace()` with a pointer to NUTS rather
> than being silently dropped. Deterministic Laplace variance estimates for
> binary occupancy carry the usual small-cluster (PQL) bias; NUTS is the
> calibrated route. RE parameter naming + per-group BLUP reconstruction for both
> engines live in `R/re_effects.R` (`ranef()` / `coef()` overrides in
> `R/methods.R`).
>
> Legacy data-binder / engine entry are internal: `.tobs_build_model()`,
> `.tobs_fit_model()`, `.tobs_laplace()`.
> See `PLAN_tulpaObs.md` for the family roster, `R/obs_families.R` and
> `R/tobs.R` for the public dispatcher.

## Common commands

```r
# Rebuild C++ exports + roxygen after editing src/ or R/ docs
Rcpp::compileAttributes()
devtools::document()

# Install upstream tulpa from sibling checkout (required after tulpa ABI bumps)
devtools::install("../tulpa", quick = TRUE)

# Iterative dev
devtools::load_all()

# Full check (skip manual since we may have non-ASCII issues in dev)
devtools::check(args = "--no-manual")

# Run all tests
devtools::test()

# Run a single test file
devtools::test(filter = "sla-cover-joint")          # matches test-sla-cover-joint.R
testthat::test_file("tests/testthat/test-occu.R")

# Run a single test_that block
testthat::test_file("tests/testthat/test-occu.R", desc = "single fit recovers truth")
```

Convention: ad-hoc probes / reproducers / planning notes live in `dev_notes/`
(R files prefixed `_` are runners, `probe_*.R` are diagnostics, `repro_*.R`
are minimal reproducers — useful when surfacing upstream `tulpa` bugs).

## Architecture

Compositional builder: one C++ entry point (`cpp_occu_fit`) for NUTS.
Laplace via tulpa's EM+Laplace engine (`tulpa::tulpa_em_laplace`); the default
`method = "laplace"` is a penalized EM with Gaussian marginals (no post-EM
correction), where the fixed-effect prior is attached to each M-step block as a
per-block `beta_prior` that tulpa applies. MI/Gibbs corrections are opt-in via
`method = "laplace_gibbs"` / `"laplace_mi"`; the same prior threads into the
correction refits (gcol33/tulpa#27), so they are penalised like `"laplace"`
unless `priors = FALSE`. A second nested-Laplace path (`em_nested_laplace.R`,
`simplified_laplace.R`, the `sla_*` files) handles families where the standard
EM closed-form does not apply (cover hurdle joint, dynamic/integrated occupancy
under `method = "nested_laplace"`).

### Boundary: What lives here vs tulpa

- **tulpaObs owns**: family-specific likelihoods (`src/*_likelihood.h`),
  E-step weights, M-step encoding, family-specific S3 / diagnostics.
- **tulpa owns**: EM engine, MI/Gibbs correction, Rubin's pooling,
  generic S3, generic diagnostics, NUTS/HMC sampler.

When NUTS crashes for a component that has correct `populate_*` code here,
the bug is in tulpa's `hmc_sampler.cpp`, not tulpaObs (see "Wired but
blocked" below). File the issue against `gcol33/tulpa`, not this repo.

### Key design rules

- **Never pass `Rcpp::Nullable<T>` to header helper functions** — MinGW
  crashes. Unwrap in the `.cpp` and pass concrete types into headers.
- **Composition over registry** — families combine with spatial /
  temporal / re / svc / latent specs orthogonally; no per-combination
  branches.
- **tulpaObs defines likelihoods, tulpa handles structure.**
- **Dotted argument names (spOccupancy / base-R style), never underscores.**
  Prefer a single word (`visits`, `priors`, `control`); when a compound is
  unavoidable, separate with `.` not `_` — e.g. `n.iter`, `n.warmup`,
  `n.chains`, `n.thin`, `n.threads`, `adapt.delta`, `max.treedepth`,
  `max.iter`, `sigma.beta`. This governs `tobs()` arguments and every
  `control = list(...)` key (the control names are the user-facing vocabulary;
  they are splatted into `.tobs_fit_model()`). Internal-only helper formals may
  stay underscore, but anything a user types follows the dotted convention.

## What works (tested)

| Feature                 | Laplace | NUTS  | Notes                                   |
|-------------------------|---------|-------|-----------------------------------------|
| Single-season occupancy | Yes     | Yes   | Full parity with inlaocc                |
| Dynamic (HMM)           | Yes     | Yes   | Colonization / extinction               |
| Community / multi-spp   | Yes     | Yes   | Species-level RE                        |
| Integrated multi-source | Yes     | Yes   | Shared psi                              |
| JSDM                    | Yes     | Yes   | No detection process                    |
| Cover hurdle (joint)    | Yes     | —     | `family_cover_hurdle.R`, `sla_cover_*`  |
| Spatial ICAR/BYM2/NNGP  | —       | Yes   |                                         |
| Spatial + dynamic       | —       | Yes   |                                         |
| Spatial + community     | —       | Yes   |                                         |
| Formula RE (intercept)  | Yes     | Yes   | `(1 \| g)`; Laplace via variance-component EM (gcol33/tulpaObs#11) |
| Formula RE (uncorr slope)| Yes    | Yes   | `(x \|\| g)`, `(0 + x \| g)`, `(1 + x \|\| g)` |
| Formula RE (corr slope) | —       | Yes   | `(1 + x \| g)`; Cholesky cov is NUTS-only (Laplace errors -> NUTS) |
| All S3 methods          | Yes     | Yes   | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, `$.tobs_fit` |
| Diagnostics             | Yes     | Yes   | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatialRange, temporalCorr |
| Simulation              | Yes     | Yes   | `simulate_occu`, `simulate_ms_occu`, `simulate_dyn_occu`, `simulate_int_occu`, `simulate_dyn_ms_occu`, `simulate_int_ms_occu`, `simulate_cover`, `simulate_cover_joint` |
| Spatial prediction      | —       | Yes   | `tobs_predict_spatial` (IDW on field)   |
| Components              | Yes     | Yes   | `tobs_re`, `tobs_temporal`, `tobs_svc`, `tobs_latent`, `tobs_community_re`, `tobs_areal` |

## NUTS coverage status

A prior CLAUDE.md flagged `temporal`, multi-term `re`, `svc`, and `latent` as
crashing under NUTS. Smoke-tested 2026-05-20 (`dev_notes/probe_blocked_nuts.R`,
N=40, 50 iter / 25 warmup, single-season occupancy + `ms_occu` for latent) —
all four return `tobs_fit` objects without crashing. The C++ population path is
not segfaulting.

What this **does not** verify: correctness of gradients, posterior calibration,
convergence at production iter counts, or behavior under stacking (spatial +
temporal + multi-RE in one fit). No tests in `tests/testthat/` exercise these
components under NUTS — only Laplace / nested-Laplace paths are covered. Treat
the smoke test as "not blocked", not "validated".

## Performance

N=200 sites, single-season (J=3), measured 2026-05-24:
- tulpaObs `method = "laplace"` (prior-aware penalized EM, default): **~0.2 s**
- tulpaObs `method = "laplace_gibbs"` (prior-aware EM + Gibbs, n.gibbs=10): ~1.7 s
- tulpaObs `method = "nuts"`: ~13 s (historical)
- inlaocc: 0.7 s; spOccupancy MCMC: 0.9 s (historical reference)

The earlier "0.01 s" figure predates the prior-aware default driver; the
penalized EM trades a little speed for breaking the psi-p identifiability
ridge at small J. Gibbs/MI add the Rubin-pooled correction phase on top.

## File organization

```
R/
  tobs.R                   — tobs() public dispatcher + print.tobs_fit
  obs_families.R           — family constructors (occu, dyn_occu, ms_occu, …, cover)
  occu.R                   — internal .tobs_build_model() (single/dynamic/community/integrated/jsdm)
  occu_fit.R               — internal .tobs_fit_model() (Laplace default, NUTS fallback)
  occu_priors.R            — occu_priors() + print + prior-spec -> per-block beta_prior plumbing (.expand_prior/.prior_for_submodel/.attach_priors_to_blocks)
  laplace.R                — internal .tobs_laplace() + EM callbacks per model type
  em_nested_laplace.R      — nested Laplace EM for non-conjugate hyperpriors
  simplified_laplace.R     — SLA wrapper used by sla_*.R families
  sla_dyn_occu.R           — SLA path for dynamic occupancy
  sla_int_occu.R           — SLA path for integrated occupancy
  sla_cover_hurdle.R       — SLA path for cover hurdle (separate-Laplace)
  sla_cover_hurdle_joint.R — SLA path for cover hurdle (joint-Laplace)
  family_cover_hurdle.R    — .dispatch_cover() (two-Laplace hurdle), large
  sim_cover_hurdle.R       — cover hurdle simulators (incl. joint)
  formula_terms.R          — structured-term registry + constructors (.tobs_term_icar/bym2/car/gp/spde/re/temporal/svc/latent/copy), tobs_* print methods, .tobs_term_to_tulpa_spatial
  formula_parse.R          — AST parser: .tobs_parse_formula / .tobs_parse_processes / .tobs_resolve_terms / .tobs_bind_formulas
  spatial.R                — internal precompute helpers (adjacency_to_csr, compute_bym2_scale, compute_nngp_neighbors)
  methods.R                — S3 methods on tobs_fit, $.tobs_fit, predict_spatial, checkIdentifiability, tobs_priors
  diagnostics.R            — tobs_waic, tobs_ppc, tobs_test_*, tobs_pit_residuals
  data.R                   — tobs_format, tobs_data, summary/plot.tobs_data, simulators
  within_between.R         — within_between() covariate decomposition
  RcppExports.R            — generated, do not edit
src/
  occu_fit.cpp                — Unified C++ entry point
  populate_helpers.h          — populate_spatial/temporal/re/svc/latent
  occ_data.h                  — OccResponseData
  occ_likelihood.h            — Single-season occupancy likelihood
  dyn_occ_data.h              — DynOccResponseData
  dyn_occ_likelihood.h        — HMC forward algorithm
  integrated_occ_data.h       — IntegratedOccResponseData
  integrated_occ_likelihood.h — Integrated multi-source likelihood
  jsdm_likelihood.h           — JSDM (Bernoulli, no detection)
  RcppExports.cpp             — generated, do not edit
  Makevars.win                — CXX_STD=CXX17, OpenMP, -Wa,-mbig-obj (large-object MinGW workaround)
tests/testthat/                — 19 test files
vignettes/                     — cover-hurdle-*.Rmd, occupancy-spatial-spde.Rmd, occupancy-vs-inla.Rmd
dev_notes/                     — ad-hoc probes, reproducers, planning notes (gitignored history, but kept in-tree)
```

## Roxygen / Rd

After editing exported function docs or adding/removing `@export`:

```r
Rcpp::compileAttributes()   # only if src/ changed
devtools::document()
```

If `devtools::document()` regenerates an Rd that contains a non-Latin-1
Unicode char (e.g. `n⁴`, `≤`, `→`), CRAN's PDF manual build will fail.
Stick to ASCII or Latin-1 Supplement in roxygen; see global CLAUDE.md
"ASCII-Only in Roxygen/Rd" rule for the safe/unsafe character map.
