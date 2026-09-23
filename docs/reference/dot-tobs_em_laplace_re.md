# Fit single-season occupancy with formula random effects via Laplace + a variance-component EM (internal).

Fit single-season occupancy with formula random effects via Laplace + a
variance-component EM (internal).

## Usage

``` r
.tobs_em_laplace_re(
  model,
  re,
  priors = NULL,
  max_iter = 100L,
  tol = 1e-05,
  damping = 0.3,
  aghq = TRUE,
  n_quad = 9L,
  lkj_eta = 1.5,
  verbose = TRUE
)
```

## Arguments

- model:

  A single-season `tobs_model`.

- re:

  A list of `tobs_re` specs. Each enters either the occupancy or the
  detection predictor (its `$shared = c(occ, det)` membership); the iid
  intercept, uncorrelated-slope, and correlated-slope forms are
  supported on both arms. A term shared across both predictors is
  rejected (use NUTS).

- priors:

  Prior spec; not applied on this path (warns when active).

- max_iter, tol, damping:

  EM controls.

- aghq:

  Logical; run the adaptive Gauss-Hermite debias pass on the variance
  components after the EM converges (default `TRUE`). See `R/re_aghq.R`.

- n_quad:

  Quadrature points per random-effect dimension for the AGHQ debias
  (default 9).

- lkj_eta:

  LKJ shape for the RE-correlation regularization in the AGHQ debias
  (default 1.5; `1` disables it). See `R/re_aghq.R`.

- verbose:

  Print per-iteration progress.
