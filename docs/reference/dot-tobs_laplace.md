# Fit a tobs model with Laplace approximation (internal)

Uses tulpa's generic EM+Laplace engine
([`tulpa::tulpa_em_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_em_laplace.html))
with occupancy-specific E-step and M-step encoding, augmented with a
weakly-informative quadratic prior on the fixed-effect coefficients (see
[`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md)).
The prior is attached to each M-step block as a per-block `beta_prior`
(see `.attach_priors_to_blocks()` in `R/occu_priors.R`), which tulpa
threads through every phase – the EM iterations and the MI / Gibbs
correction refits alike. Supports all built-in model types. Called from
[`.tobs_fit_model()`](https://gillescolling.com/tulpaObs/reference/dot-tobs_fit_model.md);
not user-facing.

## Usage

``` r
.tobs_laplace(
  model,
  spatial = NULL,
  re = NULL,
  priors = NULL,
  sigma_beta = 10,
  max_iter = 50L,
  tol = 1e-04,
  damping = 0.3,
  correction = c("auto", "mi", "gibbs", "none"),
  n_imputations = 20L,
  n_gibbs = 10L,
  seed = NULL,
  approx = c("gaussian_laplace", "simplified_laplace"),
  re_aghq = TRUE,
  n_quad = 9L,
  lkj_eta = 1.5,
  latent_prior = NULL,
  heldout_state = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A `tobs_model` from
  [`.tobs_build_model()`](https://gillescolling.com/tulpaObs/reference/dot-tobs_build_model.md).

- spatial:

  Optional `tobs_spatial` spec (NULL for non-spatial).

- priors:

  Optional prior spec from
  [`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md).
  `NULL` -\> use the package defaults. Pass `FALSE` (or `"none"`) to
  disable the prior and recover the historical unpenalised MAP behavior
  – the penalised objective is
  `Q(beta) = -log L(beta) + sum (beta_j - mu_j)^2 / (2 sd_j^2)`, so
  `sd_j = Inf` yields a zero penalty term.

- sigma_beta:

  Reserved for future use (NUTS-side beta prior); ignored by the
  EM-Laplace path.

- max_iter, tol, damping:

  EM controls.

- correction:

  Post-EM correction (`"none"`, `"mi"`, `"gibbs"`). MI / Gibbs run
  tulpa's post-EM Rubin-pooled correction; the fixed-effect prior (when
  active) threads into the correction refits, so the corrected fit is
  penalised the same way as the EM point estimate.

- n_imputations:

  Number of MI draws when `correction = "mi"`.

- verbose:

  Print per-iteration progress.
