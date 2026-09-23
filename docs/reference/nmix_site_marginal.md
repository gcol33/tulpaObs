# Per-site N-mixture marginal as a composable random-effect callback

Exposes the Royle (2004) N-mixture per-site marginal – the latent
abundance \\N_i\\ summed out in closed form – as a reusable building
block for integrating grouped random effects over the abundance and/or
detection linear predictors. It is the bridge between the marginal
kernel (which already differentiates through both arms) and a
random-effect integrator such as
[`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html).

Unlike
[`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md),
which fixes the coefficients and returns a point fit, this helper
returns a *closure object* that evaluates the marginal and its eta-level
derivatives at arbitrary linear predictors. A random-effect integrator
perturbs the predictors by \\Z b\\ per group and calls
[`eval()`](https://rdrr.io/r/base/eval.html) at each quadrature point;
the per-site value, gradient and observed- information block are
everything the per-group integral needs.

The abundance arm is per site (\\\log\lambda_i = \eta^\lambda_i\\) and
the detection arm is per visit (\\\mathrm{logit}\\p\_{ij} =
\eta^p\_{ij}\\); the two are coupled through the shared latent \\N_i\\.
The per-site marginal observed information in the eta coordinates
\\(\eta^\lambda_i, \eta^p\_{i1}, \dots, \eta^p\_{iJ_i})\\ is \$\$B_i =
\mathrm{diag}(I^\lambda_i, I^p\_{ij}) - \mathrm{Var}(N_i\mid y_i)\\ v_i
v_i^\top, \qquad v_i = (-w_i, p\_{i1}, \dots, p\_{iJ_i}),\$\$ where
\\I^\lambda_i\\, \\I^p\_{ij}\\ are the complete-data Fisher diagonal,
\\w_i\\ is the abundance score weight (`1` Poisson, \\1-q_i\\ NB) and
\\p\_{ij} = \mathrm{plogis}(\eta^p\_{ij})\\. The off-diagonal
\\\mathrm{Var}(N_i)\\w_i\\p\_{ij}\\ is the abundance/detection coupling
an integrator placing random effects on both arms must carry (Louis
1982). This is the eta-level form of the curvature
[`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md)
sandwiches with the design matrices.

The single- vs multi-arm `make_site` adapter that wires this into a
specific integrator is deliberately not built here – this object is the
arm-agnostic foundation it sits on.

## Usage

``` r
nmix_site_marginal(
  y,
  site_idx,
  X_lambda,
  X_p,
  mixture = c("P", "NB"),
  K_max = NULL,
  headroom = NULL,
  K_site = NULL
)
```

## Arguments

- y:

  Integer vector of observed counts, one entry per visit (long form).

- site_idx:

  Integer vector, same length as `y`, 1-based site index indicating
  which site each visit belongs to.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_p:

  Numeric matrix `[n_obs x p_p]` of detection covariates (long form, row
  order matches `y`).

- mixture:

  Abundance mixing distribution: `"P"` (Poisson, default) or `"NB"`
  (negative binomial). Under `"NB"` an extra dispersion parameter
  `log_r` (log NB size) is estimated jointly and reported with its
  standard error.

- K_max:

  Marginal-sum truncation. Defaults to `max(y) + 100` (matches
  [`unmarked::pcount`](https://ecoverseR.github.io/unmarked/reference/pcount.html)),
  applied per site (see `headroom`). The returned `boundary_weight` per
  site flags any site whose posterior over N puts non-trivial mass on
  its truncation; raise `K_max` if any such weight exceeds ~1e-4.

- headroom:

  Latent-N states summed above each site's own `max(y_i)`. `NULL`
  (default) derives it from `K_max`: an unset `K_max` caps each site at
  `max(y_i) + 100`, an explicit one truncates globally and uncapped. A
  negative value disables the per-site cap. Callers that resolve the
  ceiling themselves (the community fitters share one `K_max` across
  species) pass both.

- K_site:

  Per-site latent-N ceiling, an integer vector of length `n_sites`
  (`NULL` = none). The same ceiling `headroom` derives from `max(y_i)`,
  stated directly instead: a caller holding a fitted \\\lambda_i\\ keys
  each site's ceiling to that site's abundance scale, which `max(y_i)`
  understates wherever detection is low. Takes precedence over
  `headroom`, and must be at least each site's own `max(y_i)`.

## Value

An object of class `nmix_marginal`: a list of the validated data plus
closures

- `eval(eta_lambda, eta_p, r = Inf)` – evaluate at the abundance log
  predictor `eta_lambda` (length `n_sites`) and detection logit
  predictor `eta_p` (length `n_obs`, visit order matching `y`). Returns
  a list with `log_lik` (total), `log_lik_site`, `grad_eta_lambda`,
  `grad_eta_p`, `grad_theta`, the complete-data Fisher `info_eta_lambda`
  / `info_eta_p`, the abundance score weight `score_wt_lambda`, `p` (=
  `plogis(eta_p)`), `mean_N`, `var_N`, `boundary_weight`, and (NB) the
  dispersion-coupling pieces `info_theta` / `info_lambda_theta` /
  `cov_N_stheta` / `var_stheta`.

- `eval_beta(beta_lambda, beta_p, r = Inf)` – the same, with the
  predictors formed from the stored design matrices.

- `obs_info_block(s, ev)` – the \\(1+J_s)\times(1+J_s)\\ per-site
  marginal observed-information matrix \\B_s\\ (above) for site `s`,
  from an [`eval()`](https://rdrr.io/r/base/eval.html) / `eval_beta()`
  result `ev`. Coordinates are
  `(eta_lambda_s, eta_p over site s's visits in input order)`. Plus
  `n_sites`, `n_obs`, `p_lambda`, `p_p`, `mixture`, `K_max`,
  `obs_by_site` (per-site visit row indices).

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Louis, T. A.
(1982). Finding the observed information matrix when using the EM
algorithm. *JRSS-B* 44, 226-233.

## See also

[`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md)
for the fixed-effects fit,
[`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)
for the grouped random-effect integrator this feeds.
