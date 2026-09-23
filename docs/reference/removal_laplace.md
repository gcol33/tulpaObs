# Laplace fit of the removal-sampling abundance model

Maximum-likelihood fit (non-spatial, fixed effects only) of the removal
/ sequential-depletion abundance model with a Poisson or
negative-binomial abundance mixing distribution. Pass `k` removes
`y_{ik} ~ Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the individuals
still present, and the latent `N_i ~ Poisson(lambda_i)` (or
`NegBin(lambda_i, r)`) is summed out in closed form (truncation
`K_max`). Optimisation is the same analytic inner-Newton /
observed-Fisher path the N-mixture uses; the per-site detection arm sees
the depleting available count rather than the full `N`.

## Usage

``` r
removal_laplace(
  y,
  site_idx,
  X_lambda,
  X_p,
  mixture = c("P", "NB"),
  beta_lambda_init = NULL,
  beta_p_init = NULL,
  log_r_init = NULL,
  r_max = 1e+05,
  K_max = NULL,
  max_iter = 100L,
  tol = 1e-06,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer vector of per-pass removals (long form), passes in order.

- site_idx:

  Integer vector, same length as `y`, 1-based site index; a site's
  passes must appear in pass order (the depletion order).

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_p:

  Numeric matrix `[n_obs x p_p]` of per-pass detection covariates.

- mixture:

  `"P"` (Poisson, default) or `"NB"` (negative binomial).

- beta_lambda_init, beta_p_init, log_r_init:

  Optional warm starts.

- r_max:

  Upper bound on the NB size `r` (NB only, default `1e5`).

- K_max:

  Marginal-sum truncation; defaults to `max(site removal total) + 100`.
  Must be `>=` the largest per-site removal total.

- max_iter:

  Newton iteration budget (default 100).

- tol:

  Gradient-norm convergence tolerance (default 1e-6).

- verbose:

  Print per-iteration progress.

## Value

A list of class `nmix_fit` (the shared count-marginal fit shape):
`beta_lambda`, `beta_p`, `mixture`, `log_r` / `r` (NB), `log_lik`,
`vcov`, `H_obs`, `n_iter`, `converged`, `grad_norm`, per-site `mean_N` /
`var_N` / `boundary_weight`.

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Dorazio, R.
M., Jelks, H. L., Jordan, F. (2005). Improving removal-based estimates
of abundance. *Biometrics* 61, 1093-1101.
