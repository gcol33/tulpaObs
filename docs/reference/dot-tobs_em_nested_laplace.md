# Fit an occupancy tobs model via nested-Laplace

Internal driver: assembles a multi-block latent prior from `spatial`,
`temporal`, and `re`, then routes the state M-step block through
[`tulpa::tulpa_nested_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_nested_laplace.html)
via `.tobs_laplace(latent_prior = )`. Supports single-season,
integrated, community, and dynamic occupancy (the state block is `"occ"`
in every callback set).

## Usage

``` r
.tobs_em_nested_laplace(
  model,
  spatial = NULL,
  temporal = NULL,
  re = NULL,
  priors = NULL,
  sigma_beta = 10,
  max_iter = 25L,
  tol = 0.001,
  damping = 0.3,
  heldout_state = NULL,
  grids = list(),
  verbose = TRUE
)
```

## Arguments

- heldout_state:

  Optional integer indices of state-block rows to treat as held-out
  (INLA NA-response prediction targets): they are dropped from the
  likelihood (`n_trials = 0`) but kept in the design so their latent
  value is informed by the prior. `NULL` (default) fits with no held-out
  rows.
