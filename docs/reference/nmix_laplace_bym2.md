# BYM2 Royle (2004) N-mixture model via nested Laplace

Nested-Laplace fit of the spatial N-mixture model with a BYM2 (Riebler
et al. 2016) prior on the abundance-arm spatial offset. The offset
decomposes as \$\$\phi_u = \sigma \left(\sqrt{\rho / s} \\ v_u +
\sqrt{1 - \rho} \\ w_u\right),\$\$ with \\v \sim \mathrm{ICAR}\\
(unscaled, sum-to-zero) and \\w \sim \mathrm{N}(0, I)\\ iid. \\s\\ is
the Riebler scaling factor, the geometric mean of the marginal variances
\\\mathrm{diag}(Q^{+})\\ of the intrinsic ICAR field; \\\sigma\\ is then
the joint marginal standard deviation of \\\phi\\, and \\\rho \in \[0,
1\]\\ is the spatial fraction of variance.

The inner Newton works in the joint state \\x = (\beta\_\lambda,
\beta_p, v, w)\\ (dimension \\p\_\lambda + p_p + 2
n\_{\mathrm{spatial}}\\). At the converged mode the Laplace log-marginal
is accumulated; the outer 2D grid integrates over \\(\sigma, \rho)\\.

## Usage

``` r
nmix_laplace_bym2(
  y,
  site_idx,
  map_site_to_unit,
  X_lambda,
  X_p,
  adj_row_ptr,
  adj_col_idx,
  n_neighbors,
  n_spatial,
  sigma_grid = NULL,
  rho_grid = NULL,
  mixture = c("P", "NB"),
  r_grid = NULL,
  scale_factor = NULL,
  beta_lambda_init = NULL,
  beta_p_init = NULL,
  v_init = NULL,
  w_init = NULL,
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

- sigma_grid:

  Numeric vector of \\\sigma\\ (joint sd) grid points. Defaults to
  `exp(seq(log(0.2), log(3), length.out = 5L))`.

- rho_grid:

  Numeric vector of spatial-fraction grid points in \\\[0, 1\]\\.
  Defaults to `c(0.05, 0.3, 0.5, 0.7, 0.95)`.

- mixture:

  Abundance mixing distribution `"P"` (default) or `"NB"`; under `"NB"`
  the NB size \\r\\ is integrated as an additional outer grid axis.

- r_grid:

  Optional NB size grid (NB only); see
  [`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md).

- scale_factor:

  Optional scalar Riebler scaling factor \\s\\, the geometric mean of
  the marginal variances \\\mathrm{diag}(Q^{+})\\ of the intrinsic ICAR
  field. If `NULL`, it is computed from the adjacency via dense
  eigendecomposition.

- beta_lambda_init:

  Optional warm start; default `c(log(mean(y)+0.1), 0, ...)`.

- beta_p_init:

  Optional warm start; default `rep(0, p_p)`.

- v_init, w_init:

  Optional warm starts (each length `n_spatial`) for the two BYM2 latent
  components – the structured (spatial) and unstructured (iid) effects.
  `NULL` (default) starts both at zero.

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

- `theta_grid` – `[n_grid x 2]` matrix of (sigma, rho) per grid point

- `sigma_grid`, `rho_grid` – input axes

- `scale_factor` – the Riebler scaling factor

- `log_marginal`, `weights`, `sigma_mean`, `sigma_sd`, `rho_mean`,
  `rho_sd` – posterior summaries of the joint hyperparameters

- `modes` – `[n_grid x (p_lambda + p_p + 2 n_spatial)]` per-grid modes;
  columns `(p_lambda + p_p + 1) .. (p_lambda + p_p + n_spatial)` are
  `v`, the next `n_spatial` columns are `w`

- `beta_lambda_mean`, `beta_p_mean` – weighted-mean coefficients

- `v_mean`, `w_mean` – weighted-mean ICAR and iid components

- `phi_mean` – weighted-mean total offset \\\phi\\

- other diagnostic fields as in
  [`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md)

## References

Riebler, A., Sorbye, S. H., Simpson, D., Rue, H. (2016). An intuitive
Bayesian spatial model for disease mapping that accounts for scaling.
*Statistical Methods in Medical Research* 25, 1145-1165. Morris, M.,
Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., DiMaggio, C.
(2019). Bayesian hierarchical spatial models: Implementing the BYM2
model in Stan. *Spatial and Spatio-temporal Epidemiology* 31.
