# tulpaObs — unrot + slugg audit

- Date: 2026-06-01 · branch `main` · HEAD `dc88f5e`
- Tools: `unrot scan` (AST clone / dead-code / scaling) + `slugg` grep-battery placeholder audit.
- Scope for findings: `R/`, `src/`. `dev_notes/`, `experimental/`, `tests/` triaged out (see packaging note P1).
- Line numbers drift; re-locate by symbol name.

---

## Part A — unrot (structural rot)

### P1 — packaging: `.Rbuildignore` excludes almost nothing (HIGH)
`.Rbuildignore` contains a single line (`^src/.*\.gch$`). `R CMD build` therefore ships the entire
`dev_notes/` scratch tree (≈60 probe/repro scripts) **and** `experimental/` in the CRAN tarball, and
the unrot scan is ~90% dev_notes clone noise as a result. Fix: add the standard set —
`^dev_notes$`, `^experimental$`, `^.*\.Rproj$`, `^\.Rproj\.user$`, `^README\.Rmd$`, `^docs$`,
`^_pkgdown\.yml$`, `^\.github$`, `^CLAUDE\.md$`, `^\.claude$`, `^\.vscode$` (mirror tulpa's).
This single change clears the bulk of the scan.

### Real clones (R/ + src/, scratch excluded)

| # | Sim | Location | Symbols | Note |
|---|-----|----------|---------|------|
| C1 | IDENTICAL 225 tok | `R/abun.R:216` `.tobs_nmix_re_split_arms` ≡ `R/em_laplace_re.R:256` `.tobs_re_split_arms` | byte-identical arm-splitting logic under two names in two files | **top dedup target** — one shared `.tobs_split_arms()`; the duplicate names are also a naming-drift hit |
| C2 | IDENTICAL 103 + 74 tok (3 fn each) | `src/nmix_community_oracle.cpp`, `src/nmix_community_spatial_oracle.cpp`, `src/nmix_re_oracle.cpp` — `grad_hess`, `newton_hess`, `theta_score` | the three native AGHQ oracles each carry an identical copy of the inner gradient/Newton/score kernels | violates the no-copy-paste-across-specialized-functions rule; extract `static inline` helpers in a shared `nmix_oracle_kernels.h`, call from all three |
| C3 | RENAMED 80% (3 fn, 1.3-1.7k tok) | `src/nmix_spatial.cpp:363,552` + `src/nmix_spatial_bym2.cpp:336` `cpp_nested_laplace_nmix_{icar,car_proper,bym2}` | the **largest** clone: three ~1500-token spatial nested-Laplace kernels copy-paste-specialized by spatial type | highest-risk divergence point; template over the Q-builder / log-det term |
| C4 | RENAMED 85% (3 fn) | `src/nmix_community_spatial.cpp:713,768,829` `cpp_nmix_community_spatial_{icar,car_proper,bym2}` | community spatial triplet, same pattern as C3 | same fix |
| C5 | IDENTICAL 127-146 tok (3 fn) | `R/formula_terms.R:118,131,159` `.tobs_term_{icar,bym2,car_proper}` | formula-term constructors for the CAR family are near-identical | shared `.tobs_term_car_family(type)` |
| C6 | IDENTICAL 96% (3 fn) | `R/formula_terms.R:179,201,310` `.tobs_term_{gp,multiscale_gp,svc}` | continuous-field term constructors duplicate scaffolding | shared scaffold |
| C7 | IDENTICAL 64-68 tok (4 fn) | `R/obs_families.R:292,314,336,358` `abun`, `ms_abun`, `dyn_abun`, `distance` | family-constructor boilerplate identical across observation families | one `.make_obs_family(name, ...)` factory |
| C8 | IDENTICAL 74 tok | `R/occu_priors.R:95` `occu_priors` ≡ `R/occu_priors.R:389` `cover_priors` | two prior constructors identical in the same file | parameterize by family |

### Scaling issues (CRITICAL/HIGH from unrot, R/+src/ only)
- `cpp_nested_laplace_nmix_*` (C3) and `cpp_nmix_community_spatial_*` (C4): the
  copy-paste-specialization families — adding a spatial type is O(1500 lines). Template/Q-builder
  dispatch is the structural fix and the single highest-leverage refactor in the package.
- `tobs_term_*` (14-fn family, C5/C6), `cpp_nmix_*_oracle` (C2): medium.

### Dead code
unrot reported little dead code in R/src after excluding scratch; the remaining hits are the native
oracle methods reached through the `LinkingTo: tulpa` AGHQ driver (external-pointer dispatch the
static call graph misses). Do not delete on this basis.

### Naming drift
`chain_adj`/`adj_chain`, `grid_adj`/`.grid_adj` — both live only in `dev_notes/` scratch (resolved by P1).
The one in-package drift is C1's `.tobs_nmix_re_split_arms` vs `.tobs_re_split_arms`.

---

## Part B — slugg (placeholders / fabricated quantities)

**House posture: HONEST.** Confirmed on three representative sites: the package drops with `-Inf`
weight, returns `NA`, or fails loudly rather than fabricating. No fabricated reported quantity found.

### Findings

| Sev | Location | What it is | Verdict |
|-----|----------|-----------|---------|
| OK (cite) | `R/laplace.R:1295` `.se_from_laplace_fit` | returns `rep(NA_real_, p)` when `H_beta` is unavailable (spatial-mesh fits), explicitly "instead of carrying a placeholder" | this is the house posture — honest NA-on-unavailable. The headline contrast. |
| OK | `src/nmix_spatial.cpp:661` "Q(rho) not PD: skip but record placeholder" | sets `log_marginals[k]=R_NegInf`, `converged=false`, `grad_norm=+Inf` | honest drop: a `-Inf`-weight grid cell, contributes zero to the integration. Not a fabricated value. |
| OK | `R/occu_cover_joint_coupled.R:138-167` `arm_psi$y/n_trials = rep(0,...)`, `arm_p$spatial_idx = rep(0L,...)` | structural pseudo-data | overwritten downstream: `coupled = TRUE` skips the per-obs scatter and the cell-coupling spec writes every derivative from the cell-level occupancy mixture. Canonical benign case. |
| OK | `R/family_cover_hurdle.R:969` `arm_pos$phi` "placeholder overridden per grid" | init scalar | overwritten per hyperparameter grid cell. |
| OK | `R/nmix_laplace_re_spatial.R:29`, `R/occu_cover_spatial.R:336` `rough`/`crude` | init heuristics | OLS-like warm starts feeding an EM that runs to convergence per grid point; the comment is accurate, the start value does not enter the reported result. |
| GATE | `R/abun.R:{145,150,155,258}`, `R/em_nested_laplace.R:{499,507}`, `R/occu.R:34`, `R/tobs.R:918`, `R/sla_int_occu.R:169`, `R/simplified_laplace.R:279` | `stop("... not yet supported / wired / implemented")` and `status: validity_failed` returns | absent-not-faked. The unimplemented combinations (dynamic community, SLA on spatial Sigma, temporal-on-Nmix) gate loudly. |

Note: "simplified_laplace" throughout is the **named method** (Simplified Laplace Approximation,
Rue et al. 2009), not a self-simplified stand-in — the grep `simplif` match is a false positive for
slop purposes.

No test exercises a fabricated path because none was found; the recovery suite operates on the honest
drop/NA branches.

---

## Fixes applied (2026-06-01)
- P1: `.Rbuildignore` rewritten to the standard set (`dev_notes`, `experimental`, Rproj, docs,
  pkgdown, .github, CLAUDE.md, .claude, .vscode, .lintr, tarballs) — scratch no longer ships.
- C1: `.tobs_nmix_re_split_arms` (abun.R) and `.tobs_re_split_arms` (em_laplace_re.R) now both
  delegate to one parameterized `.tobs_split_re_arms()` in `em_laplace_re.R`. Behavior preserved
  (same return-list names, `$arm` tags, second-arm `p<t>` group label). Parse-verified.

## Not applied (needs a verified cycle)
- C2 (oracle `grad_hess`/`newton_hess`/`theta_score`): the identical hits are 3-8 line dispatchers;
  the real shared computation (`eval_group`) is legitimately specialized per oracle. Low value,
  cross-class entanglement. Left as-is.
- C3 (`cpp_nested_laplace_nmix_{icar,car_proper,bym2}`, ~1.3-1.7k tok each): highest-value clone but
  the variation is structural (ICAR uses a 2D `(tau,r)` grid; CAR_proper a 3D `(tau,rho,r)` grid;
  BYM2 adds a mixing axis), so it is a multi-axis kernel-templating refactor, not a parameterization.
  A sign/index slip in a spatial precision is silent statistical bias, so this needs a build +
  spatial recovery-test gate before it can land. Tracked: gcol33/tulpaObs#28 (verified pass).
