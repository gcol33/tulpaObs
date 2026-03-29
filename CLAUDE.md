# tulpaOcc — Bayesian Occupancy Models via tulpa

Full-featured Bayesian occupancy modeling. Feature parity with inlaocc.

## Architecture

Compositional builder: one C++ entry point (`cpp_occu_fit`) for NUTS.
Laplace via tulpa's EM+Laplace engine with MI/Gibbs correction.

### Boundary: What Lives Here vs tulpa

**tulpaOcc owns**: likelihoods, E-step weights, M-step encoding, occupancy-specific S3/diagnostics.
**tulpa owns**: EM engine, MI/Gibbs correction, Rubin's pooling, generic S3, generic diagnostics.

See memory file `project_tulpa_boundary.md` for full details.

### Key Design Rules
- **Never pass `Rcpp::Nullable<T>` to header helper functions** — MinGW crash
- **Composition over registry**
- **tulpaOcc defines likelihoods, tulpa handles structure**

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

These components have correct populate_* code in tulpaOcc but segfault in tulpa's
NUTS engine via the cross-package API. Need to fix in tulpa's hmc_sampler.cpp.

- **Temporal** (AR1, RW1, RW2, IID): `n_temporal_params` is set but NUTS crashes
- **Multi-term RE**: `populate_re` builds correct `re_group_multi_flat` but NUTS crashes
- **SVC**: `populate_svc` fills `SVCData` but NUTS crashes
- **Latent factors**: `populate_latent` sets `latent_n_factors` but NUTS crashes

Root cause: tulpa's multi-process NUTS path only tested with spatial + legacy single-term RE.
The temporal/RE/SVC/latent gradient code paths haven't been exercised via `LikelihoodSpec`.

## Performance

N=200, single-season, p=0.4:
- tulpaOcc Laplace+Gibbs: **0.01s** (90x faster than spOccupancy)
- inlaocc: 0.7s (after Gibbs fix)
- spOccupancy MCMC: 0.9s
- tulpaOcc NUTS: 12.8s

## Building

```r
devtools::install("../tulpa", quick = TRUE)
devtools::load_all()
devtools::check(args = "--no-manual")
```

## File Organization

```
R/
  occu.R            — occu() constructor (single/dynamic/community/integrated/jsdm)
  occu_fit.R        — occu_fit() dispatcher (Laplace default, NUTS fallback)
  occu_output.R     — print/summary
  laplace.R         — EM callbacks per model type, build_laplace_fit
  components.R      — occu_re, occu_temporal, occu_svc, occu_latent, occu_community_re, occu_areal
  spatial.R         — occ_icar/bym2/gp/multiscale_gp/spde
  methods.R         — S3 methods, $.tulpaOcc_fit, predict_spatial, checkIdentifiability
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
