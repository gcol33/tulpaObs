# tulpaOcc — Bayesian Occupancy Models via tulpa

Single-season and multi-season Bayesian occupancy models using the tulpa engine.

## Architecture

tulpaOcc plugs occupancy likelihoods into tulpa's generic multi-process interface:
- **Process 0**: Occupancy (psi) — site-level covariates
- **Process 1**: Detection (p) — site-level covariates
- **Extra params**: Visit-level detection covariates (optional)
- **Response data**: Detection histories (sites × visits) via `OccResponseData`

## Current Status

- **Single-season occupancy**: working (MacKenzie et al. 2002)
- **Multi-season dynamic occupancy**: working (MacKenzie et al. 2003). HMM forward algorithm with colonization/extinction.
- **Community occupancy**: working. Multi-species with species-level random effects.
- **Spatial occupancy**: working. ICAR, BYM2, GP (NNGP), multi-scale GP on psi and/or p.
- Sampler: full NUTS via tulpa cross-package API (R_GetCCallable)
- Gradients: arena AD (reverse-mode) + forward AD wired.
- Cross-DLL arena sync: required on Windows/MinGW (static thread_local duplicated per DLL).
  Fix: `if constexpr` sync of `current_arena()` at likelihood entry point.
- TODO: Spatial on dynamic occupancy / community models
- TODO: Season-varying detection/colonization/extinction covariates
- TODO: Correlated species random slopes for community models

## Building

```r
# Requires tulpa installed first
devtools::install("../tulpa", quick = TRUE)
devtools::load_all()
devtools::check(args = "--no-manual")
```

## File Organization

```
R/
  occ.R             — occ() single-season model constructor
  fit.R             — occ_fit() with spatial support, print/summary methods
  spatial.R         — occ_icar(), occ_bym2(), occ_gp(), occ_multiscale_gp()
  dyn_occ.R         — dynOcc() multi-season model constructor
  dyn_fit.R         — dynOcc_fit(), print/summary methods
  community_occ.R   — communityOcc() multi-species constructor
  community_fit.R   — communityOcc_fit(), print/summary methods
src/
  occ_data.h           — OccResponseData (single-season + community)
  occ_likelihood.h     — Templated occupancy log-likelihood + AD helpers
  occ_fit.cpp          — Rcpp entry: single-season with spatial support
  dyn_occ_data.h       — DynOccResponseData (multi-season)
  dyn_occ_likelihood.h — HMM forward algorithm likelihood
  dyn_occ_fit.cpp      — Rcpp entry: dynamic occupancy
  community_occ_fit.cpp — Rcpp entry: community occupancy with species RE
tests/testthat/        — Unit and integration tests
```
