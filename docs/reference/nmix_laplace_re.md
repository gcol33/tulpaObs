# Community / multispecies N-mixture by Laplace

Fits the community (spAbundance `msNMix`) N-mixture model: a per-species
Royle (2004) N-mixture with Gaussian community hyperpriors on the
per-species abundance and detection coefficients, \$\$N\_{s,i} \sim
\mathrm{Poisson}(\lambda\_{s,i}), \quad y\_{s,i,j} \| N \sim
\mathrm{Binomial}(N\_{s,i}, p\_{s,i,j}),\$\$ \$\$\log \lambda\_{s,i} =
X\_\lambda^{(i)} (\mu\_\lambda + b^\lambda_s), \quad \mathrm{logit}\\
p\_{s,i,j} = X_p^{(ij)} (\mu_p + b^p_s),\$\$ \$\$b^\lambda_s \sim N(0,
\Sigma\_\lambda), \quad b^p_s \sim N(0, \Sigma_p).\$\$

The latent abundances integrate out per species-site in closed form (the
shared N-mixture kernel exposed by
[`nmix_site_marginal()`](https://gillescolling.com/tulpaObs/reference/nmix_site_marginal.md));
the per-species coefficient deviations \\b_s = (b^\lambda_s, b^p_s)\\
are the random effects. Both solvers assemble the per-species marginal
in compiled code via a native oracle (no per-group round trip into R).
At `n_quad = 1` (the default, the joint Laplace / glmer `nAGQ = 1`) the
fit is a Laplace-EM (block-coordinate Newton mode + closed-form
covariance M-step + Schur- complement SEs) – the fast production path. A
higher `n_quad` routes the same native oracle through the shared
compiled AGHQ engine
([`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)),
replacing each species' Laplace integral with adaptive Gauss-Hermite
quadrature to reduce the small-cluster (few species) downward bias of
the community covariances `Sigma_lambda` / `Sigma_p` – at a
`n_quad^(p_lambda + p_p)` per-species grid cost. Each species' marginal
– value, abundance/detection score, and the per-site observed-
information block carrying the \\\mathrm{Var}(N_i\mid y_i)\\ abundance/
detection coupling – is supplied as the per-group oracle, so there is
one marginal/quadrature/covariance implementation across the package.
The Gaussian community priors pin \\(\mu\_\lambda, \mu_p)\\ as fixed
effects, so no sum-to-zero constraint is needed; their standard errors
come from the marginal observed-information Hessian.

The abundance mixing distribution is Poisson (`mixture = "P"`) or
negative binomial (`mixture = "NB"`). Under NB the dispersion is itself
a per-species random effect, `log_r_s ~ N(mu_log_r, sigma_log_r)` (i.e.
`r_s ~ LogNormal`): the per-species RE vector widens to
`b_s = (b_lambda_s, b_p_s, b_logr_s)`, the community log-dispersion
`mu_log_r` joins `(mu_lambda, mu_p)` as a fixed effect, and
`sigma_log_r` joins the community covariances as a third (scalar) block,
all integrated jointly by the same AGHQ engine. The per-species size is
\\r_s = \exp(\mu\_{\log r} + b^{\log r}\_s)\\. Poisson is the \\r \to
\infty\\ limit (no dispersion coordinate).

## Usage

``` r
nmix_laplace_re(
  y,
  site_idx,
  species_idx,
  X_lambda,
  X_p,
  n_sites,
  n_species,
  mu_lambda_init = NULL,
  mu_p_init = NULL,
  Sigma_lambda_init = NULL,
  Sigma_p_init = NULL,
  K_max = NULL,
  headroom = NULL,
  max_iter = 200L,
  optimizer = c("em", "joint_fd", "joint_grad"),
  mixture = c("P", "NB", "ZIP", "ZINB"),
  r_init = 10,
  sigma_logr_init = 0.5,
  omega_init = 0.2,
  sigma_omega_init = 0.5,
  n_quad = 1L,
  n_quad_scalar = .TOBS_MIN_SCALAR_NQUAD,
  lkj_eta = 1,
  sigma_beta = 100,
  omega_sigma_prior = c(1, 0.05),
  logr_sigma_prior = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer vector of counts, one entry per observed visit (long form, all
  species stacked).

- site_idx:

  Integer vector (same length as `y`), 1-based site index.

- species_idx:

  Integer vector (same length as `y`), 1-based species index.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates (shared
  across species; one row per site).

- X_p:

  Numeric matrix `[n_obs x p_p]` of detection covariates (long form, row
  order matching `y`).

- n_sites, n_species:

  Integer counts.

- mu_lambda_init, mu_p_init:

  Optional warm starts for the community means. Default: the column
  means of independent per-species
  [`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md)
  fits.

- Sigma_lambda_init, Sigma_p_init:

  Optional warm starts for the community covariances. Default: the
  (ridge-regularized) sample covariance of the per-species coefficient
  estimates.

- K_max:

  Marginal-sum truncation (default `max(y) + 100`), applied per site
  (see `headroom`).

- headroom:

  Latent-N states summed above each site's own `max(y_i)`. `NULL`
  (default) derives it from `K_max`: an unset `K_max` caps each site at
  `max(y_i) + 100`, an explicit one truncates globally and uncapped. A
  caller that resolved the ceiling itself passes both.

- max_iter:

  Optimizer iteration cap (default 200).

- optimizer:

  Outer optimize driver over the shared native oracle: `"em"` (default)
  is the fast Laplace-EM (block-coordinate Newton mode + closed-form
  covariance M-step + Schur SE), the exact `n_quad = 1` solver.
  `"joint_grad"` and `"joint_fd"` are the joint
  `(theta, log-Cholesky Sigma)` optimizers
  ([`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html))
  and both do the `n_quad > 1` AGHQ debias of the community covariances.
  `"joint_grad"` (the fast debias path) supplies the analytic
  Fisher-identity gradient – one group sweep per step, no per-coordinate
  objective re-solve – and requires `n_quad > 1` (at `n_quad = 1` use
  the EM). `"joint_fd"` finite-differences the objective (slower; the FD
  sweep re-solves every per-species mode per coordinate) and is kept for
  correctness / architecture validation and as the `n_quad = 1` joint
  reference.

- mixture:

  Abundance mixing distribution: `"P"` (Poisson, default), `"NB"`
  (negative binomial with a per-species dispersion random effect
  `log_r_s ~ N(mu_log_r, sigma_log_r)`), or their zero-inflated
  counterparts `"ZIP"` / `"ZINB"` (a per-species structural-zero random
  effect `logit_omega_s ~ N(mu_omega, sigma_omega)`; a share `omega_s`
  of a species' sites is structurally empty). None of `"NB"` / `"ZIP"` /
  `"ZINB"` has a closed-form EM, so each defaults `optimizer` to
  `"joint_grad"` and errors on `optimizer = "em"`. The default `n_quad`
  when unsupplied is `3` for `"NB"` and the zero-inflated `"ZIP"` /
  `"ZINB"` (the tensor grid is `n_quad^(p_lambda + p_p + [NB] + [ZI])`,
  so a coarser grid keeps it tractable; \#355).

- r_init:

  Initial community-mean negative-binomial size for the joint optimizer
  (`mixture = "NB"` only; default `10`, a moderate overdispersion
  start). The optimizer carries `mu_log_r = log(r_init)` as the
  `(p_lambda + p_p + 1)`-th fixed effect.

- sigma_logr_init:

  Initial standard deviation of the per-species log-dispersion random
  effect `log_r_s` (`mixture = "NB"` / `"ZINB"` only; default `0.5`).
  Seeds the scalar log_r covariance block.

- omega_init:

  Initial community-mean structural-zero probability for the joint
  optimizer (`mixture = "ZIP"` / `"ZINB"` only; default `0.2`). The
  optimizer carries `mu_omega = qlogis(omega_init)` as a trailing fixed
  effect.

- sigma_omega_init:

  Initial standard deviation of the per-species structural-zero-logit
  random effect `logit_omega_s` (`mixture = "ZIP"` / `"ZINB"` only;
  default `0.5`). Seeds the scalar omega covariance block.

- n_quad:

  Quadrature points per random-effect dimension passed to
  [`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)
  (default 1, the joint Laplace). A higher `n_quad` debiases the
  community covariances at a `n_quad^(p_lambda + p_p)` per-species grid
  cost, so keep it modest when the coefficient dimension is large. This
  sets the order of the correlated coefficient blocks (`lambda`, `p`);
  the scalar nuisance blocks use `n_quad_scalar`.

- n_quad_scalar:

  Quadrature points for the scalar nuisance random-effect blocks – the
  negative-binomial dispersion `log_r` and the zero-inflation
  `logit_omega` (default 3). These integrate a 1-D posterior, so
  [`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)'s
  per-block quadrature runs them at this coarser order rather than the
  full `n_quad`, trimming the `n_quad^(NB + ZI)` factor the nuisance
  axes would otherwise multiply into the tensor grid. Floored at 3 (and
  the floor wins over a smaller `n_quad`): the rule is converged there,
  and below it the reported community-mean SE on that block can collapse
  by an order of magnitude with no other symptom.

- lkj_eta:

  LKJ shape regularizing each *correlated* community covariance block's
  correlation off the boundary (default 1, no penalty); passed through
  to
  [`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html).
  Does not touch the marginal SDs.

- sigma_beta:

  Weak Gaussian ridge SD on the community means (default 100, i.e.
  `tau = 1e-4`, matching the other Laplace paths); stabilizes a
  weakly-identified community mean without materially shifting it.

- omega_sigma_prior:

  Penalized-Complexity prior `c(U, alpha)`
  (`P(sigma_omega > U) = alpha`) on the per-species structural-zero
  random effect SD (`mixture = "ZIP"` / `"ZINB"` only; default
  `c(1, 0.05)`, interior mode ~ 0.33). `sigma_omega` is the softest AGHQ
  direction and at few species can collapse to the boundary, flattening
  the marginal Hessian and attenuating the recovered SD; the weak prior
  adds curvature there (passed to
  [`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)'s
  `sigma_prior`) without biasing an identified fit. `NULL` restores pure
  ML on the structural-zero variance. Ignored for Poisson / NB.

- logr_sigma_prior:

  Penalized-Complexity prior `c(U, alpha)`
  (`P(sigma_log_r > U) = alpha`) on the per-species log-dispersion
  random effect SD (`mixture = "NB"` / `"ZINB"` only; default `NULL`,
  pure ML). `sigma_log_r` is the same shape of parameter as
  `sigma_omega` – one scalar variance over species – and settles near
  its lower boundary the same way at few species. At 8 and 36 species
  with a simulated `sigma_logr = 0.5`, `n_quad = 3` and
  `n_quad_scalar = 3`, fits recovering `sigma_log_r >= 0.30` covered
  `mu_log_r` 33/34 while those below covered 2/5.

  That split is a property of those two group counts and does not carry
  to a third. Re-measured at 18 species on the same fixture and seeds
  (19 fits), one fit sits below 0.30, the point-estimate error is
  uncorrelated with the recovered SD (Spearman -0.08), and the reported
  SE is very nearly a deterministic multiple of it (Spearman +0.98, R^2
  = 0.99). So what a threshold split picks up at 18 species is the SE
  side alone: a low recovered SD buys a proportionally narrower interval
  around an error that did not shrink with it. Both halves moving
  together is what the 8-and-36 pool showed and what is absent there.

  The SE itself is \\\mathrm{se}^2 = (\hat\sigma\_{\log r}^2 + c) / S\\,
  with `c` the per-species inverse information for `log_r_s`. Measured
  over 97 fits at 8 / 18 / 36 species, `c` is 0.130 / 0.133 / 0.141 –
  one constant across a factor of 4.5 in `S` – so the shape is right and
  `sigma_hat` is the input that is off. Its attenuation (0.418 / 0.448 /
  0.487 against a simulated 0.5) is monotone in `S` and passes straight
  into the interval.

  A 1.28x-too-narrow `mu_log_r` interval at 18 species on this fixture
  is mostly the seed block, not the estimator. `mu_log_r` is a
  population mean and each seed draws `S` log-dispersions around it, so
  the across-seed error carries a `sigma_logr^2 / S` term that the SE
  includes and that is itself measured on ~19 seeds: it supplies about
  two thirds of the spread, and the 18-species blocks drew it 18-21%
  wide (chi-square p = 0.12 and 0.03). Putting that draw at its
  expectation and rebuilding the SE at the simulated sigma gives a scale
  of 0.990 / 1.035 / 0.963 at 8 / 18 / 36 species. Use
  [`simulate_ms_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_abun.md)'s
  `truth$mu_log_r_real` to score against the seed's own realized species
  mean and keep the two apart.

  A fit now says whether its dispersion variance is distinguishable from
  zero at all: `fit$ms_dispersion$sigma_log_r_boundary` carries the
  boundary test and a fit that fails it warns.

  What the penalty reaches on this block appears to be the boundary
  rather than the calibration. Measured on 20 seeds of that fixture at 8
  species with `c(U, alpha) = c(1, 0.05)`: a `sigma_log_r` of 0.01-0.06
  lifts by an order of magnitude and a seed that stopped at a singular
  marginal Hessian under pure ML converges and covers, while paired
  coverage holds at 17 of 19 with the same two seeds missing, both of
  them at `sigma_log_r` around 0.2 where the penalty is weak. The fits
  that were already calibrated take a small systematic shift (`mu_log_r`
  by -0.006, p = 0.001; SEs a median 4.2% narrower). Hence the `NULL`
  default: it is a lever for a fit whose dispersion variance came back
  near zero or that failed outright, rather than a correction expected
  to hold across fits.

  When both this and `omega_sigma_prior` are set they must be equal: the
  engine applies one Penalized-Complexity prior across every block it
  regularizes. Ignored for Poisson / ZIP.

- verbose:

  Unused (kept for backward compatibility); the engine is silent.

## Value

A list of class `nmix_re_fit`: `mu_lambda`, `mu_p` (community means),
`vcov` (their joint covariance from the AGHQ marginal Hessian,
`(p_lambda + p_p)` square; marginalizes the community-covariance
uncertainty rather than plugging in `Sigma`), `Sigma_lambda`, `Sigma_p`
(community covariances), `b_lambda`, `b_p` (per-species BLUP deviations,
`n_species` rows), `log_lik` (AGHQ marginal), `converged`, `K_max`,
`n_quad`, `lkj_eta`, and (when `mixture = "NB"`) the dispersion
summaries: `mu_log_r` (community-mean log-dispersion; its SE is the
trailing `vcov` diagonal), `sigma_log_r` (the per-species log-dispersion
SD), `b_logr` (per-species deviations), `r_s` (per-species sizes
\\\exp(\mu\_{\log r} + b^{\log r}\_s)\\), and `r` equal to
\\\exp(\mu\_{\log r})\\ (the community-mean size, the LogNormal median).

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Doser, J. et
al. (2023). spAbundance. `msNMix()`.

## See also

[`nmix_laplace()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace.md)
(single species),
[`nmix_site_marginal()`](https://gillescolling.com/tulpaObs/reference/nmix_site_marginal.md)
(the per-species marginal primitive),
[`tulpa_re_aghq()`](https://gillescolling.com/tulpa/reference/tulpa_re_aghq.html)
(the shared random-effect integrator).
