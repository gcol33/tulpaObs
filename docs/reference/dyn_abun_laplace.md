# Maximum-likelihood fit of the Dail-Madsen open N-mixture

Fits the Dail-Madsen (2011) open-population N-mixture (Poisson initial
abundance, binomial survival, Poisson recruitment, binomial detection)
by maximising the exact HMM forward marginal likelihood with an analytic
gradient (BFGS). The observed-information covariance is the inverse of
the negative finite-difference Jacobian of the analytic gradient at the
mode. All four arms (initial abundance `lambda`, detection `p`, survival
`omega`, recruitment `gamma`) are site-level.

## Usage

``` r
dyn_abun_laplace(
  y_flat,
  n_sites,
  T,
  J,
  K_max,
  X_lambda,
  X_p,
  X_omega,
  X_gamma,
  mixture = "poisson",
  log_r_init = NULL,
  max_iter = 300L,
  tol = 1e-08,
  verbose = FALSE
)
```

## Arguments

- y_flat:

  Integer vector, flattened detection counts (site-major, then season,
  then visit; `-1` marks a missing visit).

- n_sites, T, J:

  Sites, primary seasons, secondary visits per season.

- K_max:

  Abundance-state truncation (states `0..K_max`).

- X_lambda, X_p, X_omega, X_gamma:

  Numeric `[n_sites x p_arm]` design matrices.

- mixture:

  Initial-abundance distribution: `"poisson"` (default) or `"negbin"`
  (negative-binomial).

- log_r_init:

  Optional starting value (log scale) for the negative-binomial size
  parameter `r`; used only when `mixture = "negbin"`. Default `NULL`
  starts from `log(2)`.

- max_iter:

  BFGS iteration budget (default 300).

- tol:

  Convergence tolerance (`optim` `reltol`, default 1e-8).

- verbose:

  Print convergence status.

## Value

A list of class `dyn_abun_fit` with `beta_lambda`, `beta_p`,
`beta_omega`, `beta_gamma`, `log_lik`, `vcov`, `H_obs`, per-site
`mean_N1`, and `converged`.

## References

Dail, D., Madsen, L. (2011). Models for estimating abundance from
repeated counts of an open metapopulation. *Biometrics* 67, 577-587.
