# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## tulpaObs — Hierarchical Latent-State Observation Models via tulpa

Bayesian occupancy / abundance / distance / removal / cover models built on the
[`tulpa`](https://github.com/gcol33/tulpa) inference engine. R package, C++17
backend via Rcpp/RcppEigen, depends on a sibling checkout of `tulpa` at
`../tulpa` (via `LinkingTo: tulpa`).

> **Public API:** `tobs()` + family constructors (`occu()`, `dyn_occu()`,
> `ms_occu()`, `int_occu()`, `jsdm()`, `abun()`, `ms_abun()`, `dyn_abun()`,
> `distance()`, `removal()`, `fp_occu()`, `cover()`, `occu_cover()`). All S3 classes are
> `tobs_*` (`tobs_fit`, `tobs_model`, `tobs_family`, `tobs_spatial`,
> `tobs_temporal`, `tobs_re`, `tobs_svc`, `tobs_latent`, `tobs_priors_spec`).
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
> form. The default `engine = "laplace"` fits iid intercept RE, uncorrelated
> slopes (`(x || g)`, `(0 + x | g)`, `(1 + x || g)`), and *correlated* slopes
> (`(1 + x | g)`) on EITHER the occupancy or the detection predictor of a
> single-season model via a variance-component EM (`R/em_laplace_re.R`,
> `.tobs_em_laplace_re()`): it splits the RE terms into an occupancy and a
> detection arm by their `shared` membership (`.tobs_re_split_arms()`) and wraps
> tulpa's fixed-covariance `tulpa_laplace()` in the occupancy missing-data EM
> (each arm's RE mode is fed back into psi / p at the E-step) with an EM/REML
> update of the per-term RE covariance `Sigma_k <- mean_g [b_g b_g' + Cov(b_g |
> y)]`. The per-group posterior covariance `Cov(b_g | y)` comes from
> `tulpa_laplace(return_re_cov = TRUE)$cov_blocks`: the occupancy arm scales it to
> natural scale by `M` (the M-step's pseudo-binomial inflation), the detection
> arm is a genuine weighted binomial (`M = 1`, prior at natural scale); a
> correlated term keeps the full `Sigma`, an uncorrelated term is projected to its
> diagonal each M-step. A single RE shared across BOTH predictors, RE + spatial,
> RE + visit-level detection, and RE on non-single families error from
> `.validate_re_laplace()` with a pointer to NUTS rather than being silently
> dropped. The raw EM variance components (sigma,
> correlation) carry the Laplace small-cluster bias for binary data (the glmer
> nAGQ=1 regime, not Breslow-Clayton PQL -- attenuated at small per-group n);
> the fixed-effect SEs are read at natural scale and do not. By default
> (`re.aghq = TRUE`, `control = list(n.quad = )`) the variance components are
> debiased after the EM by an adaptive Gauss-Hermite quadrature refinement on the
> exact per-group marginal (the nAGQ > 1 fix; single grouping factor on one arm,
> RE dim <= 3, falls back to the EM otherwise -- including when the RE is split
> across both arms). The generic quadrature engine is `tulpa::tulpa_re_aghq()`
> (quadrature grid, per-group mode, log-Cholesky covariance, LKJ penalty, joint
> optimization, marginal Hessian); `R/re_aghq.R` (`.tobs_re_aghq()`) is a thin
> wrapper supplying only the occupancy / detection site marginal as a
> `make_site` callback. The per-group marginal
> branches on the arm: an occupancy-arm RE moves `psi`, a detection-arm RE moves
> `p` (binomial-in-`p` site derivatives, FD-verified). Measured (occupancy arm):
> per-group-n = 8 sigma bias ~18% (EM) -> ~4% (AGHQ), matching NUTS; correlated
> sigmas to ~1%. Detection arm: the EM attenuates sigma ~70% (only occupied sites
> inform `p`) -> AGHQ to ~1% with 88-96% fixed-effect CI coverage. Detection-arm
> RE parameters are named for the detection process (`sigma_p<t>_*`, `re_p<t>_*`). A default LKJ(`re.lkj = 1.5`) penalty on each correlated block's
> correlation matrix regularizes a weakly-identified RE correlation off the +-1
> boundary toward 0 (untouching the marginal SDs; `re.lkj = 1` disables) -- on
> the recovery sim it removes every boundary hit while keeping rho near-unbiased;
> NUTS is available for a full posterior treatment of the correlation. RE
> parameter
> naming + per-group BLUP reconstruction for both
> engines live in `R/re_effects.R` (`ranef()` / `coef()` overrides in
> `R/methods.R`); the correlated-block off-diagonal is reported as a
> `cor_<g>_<ci>_<cj>` correlation.
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
unless `priors = FALSE`. A second **nested-Laplace** engine (`em_nested_laplace.R`,
`method = "nested_laplace"`) handles non-conjugate hyperpriors via a multi-block
latent prior. `.tobs_em_nested_laplace()` is a thin wrapper over
`.tobs_laplace(latent_prior = )`: the multi-block prior built from the formula's
`spatial` / `temporal` / `re` terms is attached to the state ("occ") M-step
block (`.tobs_laplace_nested()`), so tulpa's EM routes that block through
`tulpa_nested_laplace()`. It is wired for single-season, integrated, community,
and dynamic occupancy (one set of `build_*_callbacks` shared with the
single-Laplace path; `.tobs_state_block_dims()` maps each state row to its
spatial unit, so a community model shares a site-level field across species),
and the cover hurdle joint path uses `tulpa_nested_laplace_joint()` separately.
The single-season path also does **INLA-style NA-response prediction** with
**calibrated credible intervals**: sites with an all-missing detection history
(`.tobs_heldout_sites()`) are interpolated by the latent field, and
`predict(type = "state")` returns the per-site psi posterior (`psi`,
`psi_lower`, `psi_upper`, 95%) marginalised over the hyperparameter grid. The
intervals are calibrated by refining the EM field with one exact-marginal pass
(`.tobs_occu_state_marginal_fit()`): the latent occupancy state z is integrated
out, so each site is a Bernoulli on D = 1{>=1 detection} with mean
q*sigma(eta), where q = 1 - (1-p)^J (a held-out site has no visits -> q = 0 ->
dropped from the likelihood, interpolated by the field). This fits tulpa's
generic `family = "bernoulli"` with a per-observation probability scale
`det_prob = q` (gcol33/tulpa: scaled-Bernoulli family), so the converged Hessian
is the marginal curvature -- `fitted_eta_var` (a_i' H^{-1} a_i off the live
factor) is the calibrated per-cell predictive variance and the grid weights are
on the natural scale (no M-inflation). The EM's M-inflated pseudo-binomial,
which weights the data ~M times the prior and whose unit-trial Hessian is the
complete-data information, is kept only for the mode/detection estimate; reading
variance or grid weights off it under-covers and collapses the grid (the
refinement is the field analogue of `.tobs_occu_marginal_refine()` for the
fixed effects). The per-row eta posterior is a Gaussian mixture over grid cells
(`.nested_psi_mean()` Gauss-Hermite mean, `.nested_psi_quantiles()` mixture-CDF
quantiles). Measured held-out coverage ~1.0 (conservative), cor ~0.88, MAE ~0.11
on a 10x10 icar/bym2 grid (`dev_notes/probe_nested_ci_coverage.R`). Other model
types keep the EM occ fit (their NA-response mapping is not yet wired). Older
tulpa without `fitted_eta_var` reports `psi` only with NA interval columns.
Orthogonal to the engine, the **simplified-
Laplace skewness correction** (`simplified_laplace.R`, the `sla_*` files) is a
post-fit marginal refinement selected by `approx = "simplified_laplace"` (the
`*_sla` method names); it is computed for single / dynamic / integrated
occupancy and the cover hurdle, and no-ops to the Gaussian marginal (recording
`sla_status`) for community / jsdm.

Backend coverage is uneven across families and is enforced centrally:
`.tobs_family_methods` in `R/tobs.R` is the single source of truth for which
`method` each working family supports, and `tobs()` errors with a pointer to the
supported set rather than silently downgrading the engine. `nested_laplace` is
occu / int_occu / ms_occu / dyn_occu + cover; the `*_sla` skew variants on the
nested path are occu + cover only; the cover hurdle has no NUTS likelihood or
EM-correction engine (no `nuts` / `laplace_gibbs` / `laplace_mi`). `abun`
(N-mixture) supports `laplace` (non-spatial) + `nested_laplace` (areal spatial);
`ms_abun` (community / multispecies N-mixture) supports `laplace`; see below.
The remaining roster (`dyn_abun`, `distance`, `removal`, `fp_occu`) is
`status = "planned"` and errors via `.stop_planned_family()`.

### N-mixture abundance (`abun()`)

Royle (2004) N-mixture: `N_i ~ Poisson(lambda_i)`, `y_ij | N_i ~ Binomial(N_i,
p_ij)`. Unlike occupancy, there is **no EM** — the latent `N` marginalises out
in closed form (sum to `K_max`), so tulpa fits the marginal directly with
analytical gradients and observed-Fisher curvature. tulpaObs owns only the
family wiring (`R/abun.R`): the data binder produces `model_type = "nmix"`
(site-level `X_lambda`, long-form `y` / `site_idx` / `X_p` dropping NA visits),
and `.tobs_fit_nmix()` dispatches to the in-tree N-mixture Laplace fitter
(`R/nmix_laplace*.R`, `src/nmix_*.{cpp,h}`), wrapped by `build_nmix_fit()`. Pseudo-draws are MVN from the **joint** (lambda, p)
coefficient covariance, so derived quantities propagate the cross-arm
covariance. The abundance arm uses a log link: `compute_intercepts()` reads a
per-process `pi$link` (default `"logit"`, `"log"` for lambda) so the lambda
intercept back-transforms with `exp()`. `simulate_abun()` +
`tests/testthat/test-abun.R` cover point recovery and 95% CI coverage.

- **Non-spatial** (`method = "laplace"`): `nmix_laplace()` (package-internal).
  `vcov` is the marginal observed-Fisher inverse (full joint lambda/p block).
- **Areal spatial** (`method = "nested_laplace"`): an `icar()` / `bym2()` /
  `car_proper()` term on the abundance formula ->
  `.tobs_fit_nmix_spatial()` -> `nmix_laplace_{icar,bym2,car_proper}`
  (one spatial unit per site). The coefficient covariance is **grid-integrated**
  (law of total covariance over the hyperparameter grid): the spatial kernels
  return per-grid `cov_blocks` (the beta-block of each grid mode's joint
  `H^{-1}`), and the R wrapper assembles
  `V = sum_k w_k [cov_k + (m_k - mbar)(m_k - mbar)']` (see `.nmix_grid_vcov()`).
  For the rank-deficient intrinsic fields (ICAR, BYM2's structured component)
  the within-grid block is computed under the sum-to-zero constraint (a large
  `(sum field)^2` penalty in `nmix_spatial_beta_cov()` / `nmix_beta_cov_bym2()`)
  so the intercept variance reflects the constraint the mode sits under rather
  than the flat (intercept, field-mean) confounding of the improper prior.
- **Negative binomial** (`abun(mixture = "negbin")`): the family `mixture`
  (`"poisson"`/`"negbin"`) maps to the kernel's `"P"`/`"NB"` code and threads
  into both paths (`nmix_laplace(mixture=)` non-spatial, the spatial grid
  fitters via `common$mixture`). The NB adds an overdispersion `r`
  (`Var(N) = lambda + lambda^2 / r`). Non-spatially `log_r` is estimated jointly
  with the betas and returned as the trailing vcov coordinate; `build_nmix_fit()`
  carries it (`means`/`vcov`/`draws` gain a `log_r` column, named in
  `nmix_dispersion` with its delta-method `r_sd`) -- the autoscale unscaler
  leaves the trailing non-process coordinate untouched, and `coef()`/`confint()`/
  `vcov()` (which walk `process_info`) report the two arms only. On the
  areal-spatial path `r` is integrated over the outer grid alongside tau/rho/
  sigma, so it carries no `log_r` coordinate and is summarized as `r_mean`/`r_sd`
  in `nmix_hyper$r`. `simulate()` draws `N ~ NegBin(mu, r)` under NB. Matches
  `unmarked::pcount(mixture = "NB")`; `test-abun.R` covers point recovery,
  dispersion recovery, 95% CI coverage, and the spatial NB path.
- **Pending**: no NUTS path (no N-mixture HMC likelihood wired yet).

### Community / multispecies N-mixture (`ms_abun()`)

The spAbundance `msNMix` model: a per-species N-mixture with Gaussian community
hyperpriors on the per-species abundance and detection coefficients,
`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`, `beta_p_s ~ N(mu_p, Sigma_p)`.
The latent `N_{s,i}` integrates out per species-site in closed form (the same
kernel as `abun()`); the per-species coefficient deviations `b_s =
(b_lambda_s, b_p_s)` are the random effects, with the Gaussian community priors
pinning `(mu_lambda, mu_p)` as fixed effects (no sum-to-zero constraint, unlike
the intrinsic spatial fields). `y` is a 3D array `[n_sites x max_visits x
n_species]` or a named list of count matrices; `species =` is required.

tulpaObs owns the family wiring (`R/ms_abun.R`): the binder produces
`model_type = "ms_nmix"`, `.tobs_ms_nmix_longform()` flattens the 3D array to
the stacked `(y, site_idx, species_idx, X_p)` long form, and `.tobs_fit_ms_nmix()`
calls the in-tree **C++ Laplace-EM** `nmix_laplace_re()` -- the fast path. The
community fitter builds a native `tulpaObs::NMixCommunityOracle` (a subclass of
`tulpa::REGroupOracle` declared in `<tulpa/aghq_oracle.h>`) via
`cpp_nmix_community_oracle()`, then drives it either through the in-tree
`cpp_nmix_community_em()` EM driver (default, `n_quad = 1`) or via tulpa's
generic AGHQ engine `tulpa::tulpa_re_aghq(oracle = , ...)` for the `n_quad > 1`
variance-component debias. The fit does per-species coefficient mode-finding
(complete-data Fisher, PSD), a closed-form EM covariance M-step
`Sigma_k = mean_s[b_s b_s' + Cov(b_s | y)]`, and fixed-effect SEs from the
marginal observed-information **Schur complement** of the b-block (with the
`Var[N|y]` rank-1 coupling between the two arms; Louis 1982). The N-mixture
per-site kernel (`src/nmix_kernel.h`) caches its eta-independent `lgamma`
combinatorial terms (`compute_nmix_site_cached`), so the EM's repeated per-site
evaluations skip the lgamma recompute -- the single-shot `compute_nmix_site()`
Poisson path delegates to the same cached helper, so the single-species /
spatial fits are byte-identical (nmix regression suite passes).

`coef()` returns the community means; `ranef()` the per-species coefficient
deviations (long form: species / arm / term / estimate); `vcov()` / `confint()`
the community-mean covariance; `fitted()` / `simulate()` the per-species
`lambda` / `p` / counts. `simulate_ms_abun()` + `tests/testthat/test-ms-abun.R`
cover community-mean recovery, 95% CI coverage over 20 seeds, per-species
coefficient recovery, and the S3 surface.

- **Poisson only** for now. The native `NMixCommunityOracle` already carries a
  global negative-binomial size (`log_r` as the trailing theta coordinate) and
  the per-site kernel branches on `is.finite(r)`; the community fitter just
  needs the `mixture` flag plumbed through. An areal-spatial community field is
  the other open follow-up.

### Boundary: What lives here vs tulpa

- **tulpaObs owns**: family-specific likelihoods (`src/*_likelihood.h`,
  `src/nmix_*.{cpp,h}`), E-step weights, M-step encoding, family-specific S3 /
  diagnostics. For the RE AGHQ debias it owns the occupancy / detection
  per-site marginal (the `make_site` callback in `R/re_aghq.R`) and the
  community N-mixture native oracle (`NMixCommunityOracle` in
  `src/nmix_community_oracle.{h,cpp}`, exported to the engine as an
  `XPtr<tulpa::REGroupOracle>`).
- **tulpa owns**: EM engine, MI/Gibbs correction, Rubin's pooling, the
  callback-driven AGHQ variance-component engine (`tulpa_re_aghq()`, with the
  quadrature grid / log-Cholesky covariance / LKJ penalty / marginal Hessian),
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
| Joint occu + cover      | Yes     | —     | `occu_cover()`. **Non-spatial** (`method = "laplace"`, `R/occu_cover.R`): cell-level psi + per-visit detection + per-visit cover (beta or lognormal) on the exact two-state marginal. **Spatial v3** (`method = "nested_laplace"`, `R/occu_cover_nested.R`): nested Laplace with field z profiled out per outer (alpha, sigma) candidate; ICAR/besag (rho fixed) field shared across psi and cover arms via `alpha`. Lognormal positive arm only on spatial path (beta arm v4). Recovery: 20-seed lognormal v1 + 10-seed beta v1 + 10-seed v3 spatial; status `"experimental"`. v2's joint Laplace (`control = list(engine = "v2_joint")`) kept as a debug escape hatch (had a (z, alpha, sigma) ridge that v3 breaks). |
| N-mixture (Poisson/NB)  | Yes     | —     | `abun(mixture=)`; closed-form marginal via in-tree `nmix_laplace`, joint-vcov draws, calibrated CIs (`test-abun.R`). NB adds jointly-estimated `log_r`. NUTS pending |
| N-mixture + areal spatial| n-L    | —     | `abun()` + `icar()`/`bym2()`/`car_proper()`, `method="nested_laplace"`; Poisson or NB (size `r` integrated over the grid); grid-integrated coefficient covariance (constrained intercept), calibrated slope CIs |
| Community N-mixture     | Yes     | —     | `ms_abun()` (spAbundance `msNMix`); per-species coef RE with Gaussian community covariances; in-tree C++ Laplace-EM (`nmix_laplace_re`) driving a native `NMixCommunityOracle` via tulpa's generic AGHQ engine, Schur-complement mean SEs; Poisson only; recovery + 20-seed coverage (`test-ms-abun.R`). NUTS / negbin / spatial pending |
| Spatial ICAR/BYM2/NNGP  | —       | Yes   |                                         |
| Spatial + dynamic       | —       | Yes   |                                         |
| Spatial + community     | —       | Yes   |                                         |
| Nested-Laplace (areal)  | n-L     | —     | `method="nested_laplace"`: icar/bym2/car (+ temporal/iid) on occu / int_occu / ms_occu / dyn_occu |
| NA-response prediction  | n-L     | —     | `predict(type="state")`: all-NA single-season sites interpolated by the field (any prior incl. bym2, via engine `fitted_eta`), with calibrated 95% `psi_lower`/`psi_upper` from the exact-marginal `bernoulli`-family pass (held-out coverage ~1.0) |
| Formula RE (intercept)  | Yes     | Yes   | `(1 \| g)`; Laplace via variance-component EM, occupancy OR detection arm (gcol33/tulpaObs#11) |
| Formula RE (uncorr slope)| Yes    | Yes   | `(x \|\| g)`, `(0 + x \| g)`, `(1 + x \|\| g)`; either arm |
| Formula RE (corr slope) | Yes     | Yes   | `(1 + x \| g)`; Laplace fits the full cov via the EM M-step consuming tulpa `cov_blocks`; either arm (gcol33/tulpaObs#11) |
| Formula RE on detection | Yes     | Yes   | RE on the detection predictor (`detection = ~ (1 \| g)`); Laplace fits a separate detection-arm RE block, AGHQ debias branches on the arm |
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
  occu.R                   — internal .tobs_build_model() (single/dynamic/community/integrated/jsdm/nmix)
  abun.R                   — N-mixture family: .tobs_build_abun(), .tobs_fit_nmix() (-> tulpa nmix Laplace), build_nmix_fit(), nmix S3 helpers, simulate_abun()
  ms_abun.R                — community N-mixture: .tobs_build_ms_abun(), .tobs_ms_nmix_longform(), .tobs_fit_ms_nmix() (-> nmix_laplace_re()), build_ms_nmix_fit(), ms_nmix S3 helpers, simulate_ms_abun()
  nmix_laplace.R           — in-tree non-spatial N-mixture (Royle 2004) Laplace fitter (Poisson + NB)
  nmix_laplace_re.R        — in-tree community / multispecies N-mixture (msNMix), .nmix_re_oracle() helper
  nmix_laplace_spatial.R   — nmix_laplace_icar() / _bym2() / _car_proper() areal spatial fitters
  nmix_site_marginal.R     — per-site marginal exposed as a composable AGHQ RE callback
  occu_fit.R               — internal .tobs_fit_model() (Laplace default, NUTS fallback; routes model_type=="nmix" to .tobs_fit_nmix)
  occu_priors.R            — occu_priors() + print + prior-spec -> per-block beta_prior plumbing (.expand_prior/.prior_for_submodel/.attach_priors_to_blocks)
  laplace.R                — internal .tobs_laplace() + EM callbacks per model type
  em_nested_laplace.R      — nested Laplace EM for non-conjugate hyperpriors
  simplified_laplace.R     — SLA wrapper used by sla_*.R families
  sla_dyn_occu.R           — SLA path for dynamic occupancy
  sla_int_occu.R           — SLA path for integrated occupancy
  sla_cover_hurdle.R       — SLA path for cover hurdle (separate-Laplace)
  sla_cover_hurdle_joint.R — SLA path for cover hurdle (joint-Laplace)
  family_cover_hurdle.R    — .dispatch_cover() (two-Laplace hurdle), large
  occu_cover.R             — joint occu-detection + cover hurdle: family wiring (builder, dispatcher, simulator, non-spatial Laplace fitter).
  occu_cover_spatial.R     — v2 joint-Laplace spatial path (debug escape hatch via control$engine = "v2_joint"); shares the ICAR Q + Sorbye-Rue scale helpers with v3.
  occu_cover_nested.R      — v3 nested-Laplace spatial path (default for method = "nested_laplace"): inner Newton on z profiled out per outer (alpha, sigma) candidate; outer BFGS over ~10 params; ICAR/besag field shared across psi and cover arms via alpha. Lognormal positive arm only.
  sim_cover_hurdle.R       — cover hurdle simulators (incl. joint)
  formula_terms.R          — structured-term registry + constructors (.tobs_term_icar/bym2/car/gp/spde/re/temporal/svc/latent/copy), tobs_* print methods, .tobs_term_to_tulpa_spatial
  formula_parse.R          — AST parser: .tobs_parse_formula / .tobs_parse_processes / .tobs_resolve_terms / .tobs_bind_formulas
  spatial.R                — internal precompute helpers (adjacency_to_csr, compute_bym2_scale, compute_nngp_neighbors)
  methods.R                — S3 methods on tobs_fit, $.tobs_fit, predict_spatial, checkIdentifiability
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
