# tulpaObs — Hierarchical Latent-State Observation Models via tulpa

> **API clean-slate complete.** Public surface is now `tobs()` + family
> constructors (`occu()`, `dyn_occu()`, `ms_occu()`, `int_occu()`, `jsdm()`,
> `abun()`, `ms_abun()`, `dyn_abun()`, `distance()`, `removal()`, `fp_occu()`,
> `cover()`). The legacy `occu()` data-binder and `occu_fit()` engine entry
> are now internal (`.tobs_build_model()`, `.tobs_fit_model()`, `.tobs_laplace()`).
> All `tulpaObs_*` S3 classes renamed to `tobs_*` (`tobs_fit`, `tobs_model`,
> `tobs_family`, `tobs_spatial`, `tobs_temporal`, `tobs_re`, `tobs_svc`,
> `tobs_latent`, `tobs_priors`). Spatial/component helpers renamed to `tobs_*`
> (e.g. `occu_icar` → `tobs_icar`, `occu_re` → `tobs_re`).
> See `PLAN_tulpaObs.md` for the family roster and `R/obs_families.R` /
> `R/tobs.R` for the public API.

Full-featured Bayesian occupancy modeling. Feature parity with inlaocc.

## Architecture

Compositional builder: one C++ entry point (`cpp_occu_fit`) for NUTS.
Laplace via tulpa's EM+Laplace engine with MI/Gibbs correction.

### Boundary: What Lives Here vs tulpa

**tulpaObs owns**: family-specific likelihoods, E-step weights, M-step encoding, family-specific S3/diagnostics.
**tulpa owns**: EM engine, MI/Gibbs correction, Rubin's pooling, generic S3, generic diagnostics.

See memory file `project_tulpa_boundary.md` for full details.

### Key Design Rules
- **Never pass `Rcpp::Nullable<T>` to header helper functions** — MinGW crash
- **Composition over registry**
- **tulpaObs defines likelihoods, tulpa handles structure**

## What Works (Tested)

| Feature | Laplace | NUTS | Notes |
|---|---|---|---|
| Single-season | Yes | Yes | Full parity with inlaocc |
| Dynamic (HMM) | Yes | Yes | Colonization/extinction |
| Community | Yes | Yes | Species-level RE |
| Integrated | Yes | Yes | Multi-source, shared psi |
| JSDM | Yes | Yes | No detection process |
| Spatial ICAR | — | Yes | |
| Spatial BYM2 | — | Yes | |
| Spatial GP (NNGP) | — | Yes | |
| Spatial + dynamic | — | Yes | |
| Spatial + community | — | Yes | |
| All S3 methods | Yes | Yes | coef, confint, vcov, logLik, nobs, fitted, residuals, simulate, predict, tidy, glance, ranef, update, summary, $ accessor |
| All diagnostics | Yes | Yes | WAIC, PPC, PIT, dispersion, zero-inflation, outliers, Moran's I, DW, variogram, spatialRange, temporalCorr |
| All simulation | Yes | Yes | simulate_occu, simMsOcc, simTOcc, simIntOcc, simTMsOcc, simIntMsOcc |
| Model comparison | Yes | Yes | compare_models, modelAverage, postHocLM |
| Data exploration | Yes | Yes | summary.occu_data, plot.occu_data |
| Spatial prediction | — | Yes | predict_spatial (IDW interpolation of spatial field) |
| Components | Yes | Yes | occu_re, occu_temporal, occu_svc, occu_latent, occu_community_re, occu_areal |

## What's Wired But Blocked by tulpa Bugs

These components have correct populate_* code in tulpaObs but segfault in tulpa's
NUTS engine via the cross-package API. Need to fix in tulpa's hmc_sampler.cpp.

- **Temporal** (AR1, RW1, RW2, IID): `n_temporal_params` is set but NUTS crashes
- **Multi-term RE**: `populate_re` builds correct `re_group_multi_flat` but NUTS crashes
- **SVC**: `populate_svc` fills `SVCData` but NUTS crashes
- **Latent factors**: `populate_latent` sets `latent_n_factors` but NUTS crashes

Root cause: tulpa's multi-process NUTS path only tested with spatial + legacy single-term RE.
The temporal/RE/SVC/latent gradient code paths haven't been exercised via `LikelihoodSpec`.

## Performance

N=200, single-season, p=0.4:
- tulpaObs Laplace+Gibbs: **0.01s** (90x faster than spOccupancy)
- inlaocc: 0.7s (after Gibbs fix)
- spOccupancy MCMC: 0.9s
- tulpaObs NUTS: 12.8s

## Building

```r
devtools::install("../tulpa", quick = TRUE)
devtools::load_all()
devtools::check(args = "--no-manual")
```

## File Organization

```
R/
  tobs.R            — tobs() public dispatcher + print.tobs_fit
  obs_families.R    — family constructors (occu, dyn_occu, ms_occu, …, cover)
  occu.R            — internal .tobs_build_model() (single/dynamic/community/integrated/jsdm)
  occu_fit.R        — internal .tobs_fit_model() (Laplace default, NUTS fallback)
  laplace.R         — internal .tobs_laplace() + EM callbacks per model type
  family_cover_hurdle.R — .dispatch_cover() (two-Laplace hurdle)
  components.R      — tobs_re, tobs_temporal, tobs_svc, tobs_latent, tobs_community_re, tobs_areal
  spatial.R         — tobs_icar/bym2/gp/multiscale_gp/spde (returns tobs_spatial)
  methods.R         — S3 methods on tobs_fit, $.tobs_fit, predict_spatial, checkIdentifiability, tobs_priors
  diagnostics.R     — waicOccu, ppcOccu (generic diagnostics inherited from tulpa)
  data.R            — occu_format, occu_data, summary/plot.occu_data, all simulation functions
src/
  occu_fit.cpp               — Unified C++ entry point
  populate_helpers.h          — populate_spatial/temporal/re/svc/latent
  occ_data.h                  — OccResponseData
  occ_likelihood.h            — Single-season occupancy likelihood
  dyn_occ_data.h              — DynOccResponseData
  dyn_occ_likelihood.h        — HMC forward algorithm
  integrated_occ_data.h       — IntegratedOccResponseData
  integrated_occ_likelihood.h — Integrated multi-source likelihood
  jsdm_likelihood.h           — JSDM (Bernoulli)
tests/testthat/               — ~132 tests across 7 test files
```
