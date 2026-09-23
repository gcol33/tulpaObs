# Spatial Royle (2004) N-mixture model via nested Laplace

Nested-Laplace fit of the spatial N-mixture model \$\$N_i \sim
\mathrm{Poisson}(\lambda_i), \qquad y\_{ij} \| N_i \sim
\mathrm{Binomial}(N_i, p\_{ij}),\$\$ with abundance linear predictor
\\\log \lambda_i = X\_\lambda^{(i)} \beta\_\lambda + z\_{u(i)}\\ where
\\z \sim \mathrm{ICAR}(\tau)\\ is an intrinsic conditional
autoregressive field on the user-supplied adjacency graph, and detection
linear predictor \\\mathrm{logit}\\ p\_{ij} = X_p^{(ij)} \beta_p\\.

The hyperparameter \\\tau\\ (ICAR precision) is integrated out by an
outer grid: at each \\\tau_k\\ the inner Newton finds the joint mode of
\\(\beta\_\lambda, \beta_p, z)\\ and the Laplace log-marginal \\\log p(y
\mid \tau_k)\\ is accumulated. Posterior weights over \\\tau\\ normalise
the grid.

Inner Newton uses the marginal observed Fisher information matrix for
curvature (with the Var\\\[N \mid y_i\]\\ rank-1 correction encoding
cross-arm coupling) and falls back to the complete-data Fisher block
when the observed-info matrix is not PSD. A small diagonal ridge keeps
the Cholesky stable through the (intercept, constant-\\z\\) structural
null direction; \\z\\ is centered to sum zero after every step.

## Usage

``` r
nmix_laplace_icar(
  y,
  site_idx,
  map_site_to_unit,
  X_lambda,
  X_p,
  adj_row_ptr,
  adj_col_idx,
  n_neighbors,
  n_spatial,
  tau_grid = NULL,
  mixture = c("P", "NB"),
  r_grid = NULL,
  beta_lambda_init = NULL,
  beta_p_init = NULL,
  z_init = NULL,
  K_max = NULL,
  max_iter = 100L,
  tol = 1e-06,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer vector of observed counts (long form, one entry per visit).

- site_idx:

  Integer vector, 1-based site index for each visit, same length as `y`.

- map_site_to_unit:

  Integer vector of length `n_sites`, 1-based spatial unit index for
  each site. Sites can share spatial units; the data contribution to
  `z[u]` aggregates across all sites that map to `u`.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_p:

  Numeric matrix `[n_obs x p_p]` of detection covariates.

- adj_row_ptr, adj_col_idx, n_neighbors:

  CSR adjacency for the ICAR graph on the `n_spatial` units.
  `adj_row_ptr` has length `n_spatial + 1`, `adj_col_idx` lists 0-based
  neighbours, and `n_neighbors[s]` is the row degree of unit s.

- n_spatial:

  Number of spatial units.

- tau_grid:

  Optional numeric vector of \\\tau\\ grid points. Defaults to
  `exp(seq(log(0.3), log(30), length.out = 9))`.

- mixture:

  Abundance mixing distribution: `"P"` (Poisson, default) or `"NB"`
  (negative binomial). Under `"NB"` the NB size \\r\\ is integrated as
  an additional outer grid dimension (alongside \\\tau\\); the posterior
  `r_mean` / `r_sd` are reported from the grid weights.

- r_grid:

  Optional numeric vector of NB size grid points (NB only). Defaults to
  `exp(seq(log(0.5), log(40), length.out = 6))`. Ignored under Poisson.

- beta_lambda_init:

  Optional warm start; default `c(log(mean(y)+0.1), 0, ...)`.

- beta_p_init:

  Optional warm start; default `rep(0, p_p)`.

- z_init:

  Optional warm start for the spatial field; default zeros.

- K_max:

  Truncation for the per-site marginal sum over N. Defaults to
  `max(y) + 100`. Returned `boundary_max` flags any grid point whose
  worst site puts non-trivial mass on `K_max` – raise `K_max` if it
  exceeds 1e-4.

- max_iter, tol:

  Inner Newton iteration budget and gradient-norm tol.

- verbose:

  Print per-iteration and per-grid-point progress.

## Value

A list of class `nmix_spatial_fit`:

- `tau_grid` – input grid

- `log_marginal` – log marginal at each tau (up to a tau-independent
  constant)

- `weights` – normalised grid weights (sum to 1)

- `tau_mean`, `tau_sd` – posterior moments of tau

- `modes` – `[n_grid x (p_lambda + p_p + n_spatial)]` matrix of inner
  modes

- `beta_lambda_mean`, `beta_p_mean` – weighted-mean coefficient
  estimates

- `z_mean` – weighted-mean spatial field

- `n_iter`, `converged`, `grad_norm`, `log_lik`, `boundary_max` –
  per-grid diagnostics

- `p_lambda`, `p_p`, `n_spatial`, `K_max` – echoed dimensions

- `call` – matched call

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Besag, J.,
York, J., Mollie, A. (1991). Bayesian image restoration with two
applications in spatial statistics. *Ann. Inst. Statist. Math.* 43,
1-20. Rue, H., Martino, S., Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *JRSS-B* 71, 319-392.
