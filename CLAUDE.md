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
> (in `R/formula_parse.R`) rewrites `(1 | g)`, `(x | g)`, `(x || g)` into the
> equivalent `re()` calls on the formula AST *before* `terms()` runs, so there
> is one parser and one term type. Forms the engine can't hold in one `re()`
> block (slope-only `0 + x`, multi-slope, nested `g/h`) error with a pointer to
> an explicit `re()`.
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
Laplace via tulpa's EM+Laplace engine with MI/Gibbs correction. A second
nested-Laplace path (`em_nested_laplace.R`, `simplified_laplace.R`, the
`sla_*` files) handles families where the standard EM closed-form does
not apply (cover hurdle joint, dynamic/integrated occupancy under
`engine = "nested_laplace"`).

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

N=200 sites, single-season, p=0.4:
- tulpaObs Laplace+Gibbs: **0.01 s** (90× faster than spOccupancy)
- inlaocc: 0.7 s (after Gibbs fix)
- spOccupancy MCMC: 0.9 s
- tulpaObs NUTS: 12.8 s

## File organization

```
R/
  tobs.R                   — tobs() public dispatcher + print.tobs_fit
  obs_families.R           — family constructors (occu, dyn_occu, ms_occu, …, cover)
  occu.R                   — internal .tobs_build_model() (single/dynamic/community/integrated/jsdm)
  occu_fit.R               — internal .tobs_fit_model() (Laplace default, NUTS fallback)
  occu_priors.R            — prior constructors + print.tobs_priors
  laplace.R                — internal .tobs_laplace() + EM callbacks per model type
  em_laplace_penalized.R   — penalized EM (ridge / IRLS-style)
  em_nested_laplace.R      — nested Laplace EM for non-conjugate hyperpriors
  simplified_laplace.R     — SLA wrapper used by sla_*.R families
  penalized_irls.R         — IRLS solver for the penalized M-step
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
