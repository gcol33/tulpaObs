# Continuous-field (SPDE) Royle (2004) N-mixture model via nested Laplace

Nested-Laplace fit of the spatial N-mixture model with a continuous
Matern (SPDE) field on the abundance arm: \$\$N_i \sim
\mathrm{Poisson}(\lambda_i), \qquad y\_{ij} \| N_i \sim
\mathrm{Binomial}(N_i, p\_{ij}),\$\$ \\\log \lambda_i = X\_\lambda^{(i)}
\beta\_\lambda + (A u)\_i\\, where \\u \sim \mathrm{N}(0,
Q(\mathrm{range}, \sigma)^{-1})\\ is a Gaussian Markov random field on
the FEM mesh with the proper Matern precision \\Q\\, and \\A\\
(\\n\_{\mathrm{sites}} \times n\_{\mathrm{mesh}}\\) projects mesh nodes
onto sites. The detection arm is \\\mathrm{logit}\\ p\_{ij} = X_p^{(ij)}
\beta_p\\.

The SPDE hyperparameters \\(\mathrm{range}, \sigma)\\ are integrated out
by the SAME outer grid the areal path
([`nmix_laplace_icar()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_icar.md)
/
[`nmix_laplace_bym2()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_bym2.md)
/
[`nmix_laplace_car_proper()`](https://gillescolling.com/tulpaObs/reference/nmix_laplace_car_proper.md))
uses for \\(\tau\[, \rho\], \sigma)\\: at each grid point the inner
Newton finds the joint mode of \\(\beta\_\lambda, \beta_p, u)\\ and the
Laplace log-marginal \\\log p(y \mid \mathrm{range}\_k, \sigma_k)\\ is
accumulated, then the grid is normalised into posterior weights. The
Matern PC priors on \\(\mathrm{range}, \sigma)\\ (from the
[`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term) enter each grid point's log-marginal, mirroring `fit_spde()`. The
precision \\Q\\ (and its \\\log\|Q\|\\) is built once per grid point on
the R side via the same FEM assembly the occupancy SPDE path uses
([`tulpa::tulpa_spde_precision_Q`](https://gillescolling.com/tulpa/reference/tulpa_spde_precision_Q.html)),
so the C++ kernel stays agnostic to the precision parameterisation.

Unlike the intrinsic ICAR field, \\Q\\ is full rank, so the (intercept,
field-mean) direction is identified by \\Q\\ itself: no sum-to-zero
centering and no constrained-covariance projection (the proper CAR path
is the areal analogue).

## Usage

``` r
nmix_laplace_spde(
  y,
  site_idx,
  X_lambda,
  X_p,
  spatial,
  mixture = c("P", "NB"),
  r_grid = NULL,
  range_grid = NULL,
  sigma_grid = NULL,
  beta_lambda_init = NULL,
  beta_p_init = NULL,
  u_init = NULL,
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

  Integer vector, 1-based site index for each visit.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_p:

  Numeric matrix `[n_obs x p_p]` of detection covariates.

- spatial:

  A `tobs_spatial` SPDE term (`type == "spde"`) carrying the FEM
  `tulpa_spec` (`A`, `C0_diag`, `G`, `n_mesh`, `nu`, `prior_range`,
  `prior_sigma`). The projection `A` must have `n_sites` rows.

- mixture:

  Abundance mixing distribution: `"P"` (Poisson, default) or `"NB"`
  (negative binomial). Under `"NB"` the NB size \\r\\ is integrated as
  an additional outer grid dimension; the posterior `r_mean` / `r_sd`
  are reported from the grid weights.

- r_grid:

  Optional numeric vector of NB size grid points (NB only). Defaults to
  `exp(seq(log(0.5), log(40), length.out = 6))`.

- range_grid, sigma_grid:

  Optional numeric vectors of SPDE range / sigma grid points. Defaults
  centre a log-spaced grid on the PC-prior medians.

- beta_lambda_init:

  Optional warm start; default `c(log(mean(y)+0.1), 0, ...)`.

- beta_p_init:

  Optional warm start; default `rep(0, p_p)`.

- u_init:

  Optional warm start for the mesh field; default zeros.

- K_max:

  Truncation for the per-site marginal sum over N. Defaults to
  `max(y) + 100`.

- max_iter, tol:

  Inner Newton iteration budget and gradient-norm tol.

- verbose:

  Print per-iteration / per-grid-point progress.

## Value

A list of class `nmix_spatial_fit` with the same shape as the areal
fitters plus `u_mean` (the posterior-mean mesh field) and `range_mean` /
`range_sd` / `sigma_mean` / `sigma_sd`.

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Lindgren,
F., Rue, H., Lindstrom, J. (2011). An explicit link between Gaussian
fields and Gaussian Markov random fields: the SPDE approach. *JRSS-B*
73, 423-498. Fuglstad, G.-A., Simpson, D., Lindgren, F., Rue, H. (2019).
Constructing priors that penalize the complexity of Gaussian random
fields. *JASA* 114.
