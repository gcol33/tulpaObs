# Proper CAR Royle (2004) N-mixture model via nested Laplace

Nested-Laplace fit of the spatial N-mixture model with a proper
conditional autoregressive prior on the abundance-arm spatial offset:
\$\$N_i \sim \mathrm{Poisson}(\lambda_i), \qquad y\_{ij} \| N_i \sim
\mathrm{Binomial}(N_i, p\_{ij}),\$\$ \\\log \lambda_i = X\_\lambda^{(i)}
\beta\_\lambda + z\_{u(i)}\\, where \\z \mid \tau, \rho \sim
\mathrm{N}(0, \[\tau (D - \rho W)\]^{-1})\\. Both hyperparameters are
integrated over an outer 2D grid; the inner Newton step shares the
kernel with
[`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md)
(ICAR is the \\\rho = 1\\ limit).

Unlike the ICAR fit, no sum-to-zero centering is applied – \\Q(\rho)\\
is full rank for \\\rho \< 1\\. The per-rho \\\log \|Q(\rho)\|\\ is
precomputed once via a dense Cholesky on the \\n\_{\mathrm{spatial}}
\times n\_{\mathrm{spatial}}\\ precision matrix.

## Usage

``` r
nmix_laplace_car_proper(
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
  rho_grid = NULL,
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

  Optional numeric vector of \\\tau\\ grid points (defaults to
  `exp(seq(log(0.3), log(30), length.out = 7L))`).

- rho_grid:

  Optional numeric vector of \\\rho\\ grid points in the valid
  eigenvalue interval. Defaults to a 5-point grid in \\(0, 1)\\ –
  callers that want eigenvalue-derived bounds should compute them via
  tulpa::spatial_car_proper(adjacency)`$`rho_bounds and pass an explicit
  grid in that interval.

- mixture:

  Abundance mixing distribution `"P"` (default) or `"NB"`; under `"NB"`
  the NB size \\r\\ is integrated as an additional outer grid axis.

- r_grid:

  Optional NB size grid (NB only); see
  [`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md).

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

- `theta_grid` – `[n_grid x 2]` matrix of (tau, rho) per grid point

- `tau_grid`, `rho_grid` – input axes

- `log_det_Q_rho` – precomputed log determinants per rho

- `log_marginal`, `weights`, `tau_mean`, `tau_sd`, `rho_mean`, `rho_sd`
  – as in
  [`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md)
  but with the rho marginal added

- other diagnostic fields and named coefficient means as in ICAR

## References

Cressie, N. (1993). Statistics for Spatial Data. Wiley. Rue, H., Held,
L. (2005). Gaussian Markov Random Fields. CRC.
