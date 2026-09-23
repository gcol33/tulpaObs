# Laplace fit of the binned distance-sampling abundance model

Maximum-likelihood fit (non-spatial, fixed effects only) of a binned
distance-sampling abundance model with a Poisson or negative-binomial
abundance mixing distribution and a half-normal or hazard-rate detection
key. Latent abundance `N_i ~ Poisson(lambda_i)` (or
`NegBin(lambda_i, r)`) is summed out in closed form (truncation
`K_max`); the detection integrals over distance bins are evaluated by
Gauss-Legendre quadrature.

## Usage

``` r
distance_laplace(
  y,
  X_lambda,
  X_sigma,
  cutpoints,
  key = c("halfnorm", "hazard"),
  transect = c("line", "point"),
  mixture = c("P", "NB"),
  beta_lambda_init = NULL,
  beta_sigma_init = NULL,
  eta_b_init = NULL,
  log_r_init = NULL,
  r_max = 1e+05,
  K_max = NULL,
  headroom = NULL,
  quad_order = 64L,
  max_iter = 100L,
  tol = 1e-06,
  verbose = FALSE,
  quad_xptr = NULL
)
```

## Arguments

- y:

  Integer matrix `[n_sites x n_bins]` of per-bin detected counts.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_sigma:

  Numeric matrix `[n_sites x p_sigma]` of detection-scale covariates
  (the `log sigma` linear predictor).

- cutpoints:

  Numeric bin edges, length `n_bins + 1` (`0 = c_0 < c_1 < ... < c_B`).

- key:

  `"halfnorm"` (default) or `"hazard"` detection key.

- transect:

  `"line"` (default, distances uniform) or `"point"` (radial, density
  proportional to distance).

- mixture:

  `"P"` (Poisson, default) or `"NB"` (negative binomial).

- beta_lambda_init, beta_sigma_init, eta_b_init, log_r_init:

  Optional warm starts.

- r_max:

  Upper bound on the NB size `r` (NB only, default `1e5`).

- K_max:

  Marginal-sum truncation; defaults to `3 * max(site total) + 100`
  (applied per site, see `headroom`).

- headroom:

  Latent-N states summed above each site's own detected total. `NULL`
  (default) derives it from `K_max`: an unset `K_max` caps each site at
  its own total plus `2 * max(site total) + 100`, an explicit `K_max`
  truncates globally and uncapped. A negative value disables the
  per-site cap.

- quad_order:

  Gauss-Legendre nodes per bin (default 64).

- max_iter:

  Newton iteration budget (default 100).

- tol:

  Gradient-norm convergence tolerance (default 1e-6).

- verbose:

  Print per-iteration progress.

- quad_xptr:

  Optional pre-built quadrature (from `cpp_distance_build_quad()`) for a
  caller fitting several species/starts against the same `cutpoints` /
  `transect` / `quad_order`; built once and reused instead of every call
  rebuilding it. Internal; defaults to building it fresh.

## Value

A list of class `distance_fit` with `beta_lambda`, `beta_sigma`, `shape`
/ `eta_b` (hazard-rate), `log_r` / `r` (NB), `log_lik`, `vcov`, `H_obs`,
per-site `mean_N` / `var_N` / `p_det` / `boundary_weight`.

## References

Buckland, S. T., et al. (2001). Introduction to Distance Sampling.
Oxford. Royle, J. A., Dawson, D. K., Bates, S. (2004). Modeling
abundance effects in distance sampling. *Ecology* 85, 1591-1597.
