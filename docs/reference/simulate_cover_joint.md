# Simulate cover-hurdle data with a shared BYM2 spatial field

Generates synthetic data matching the joint nested-Laplace cover-hurdle
parameterisation (`sigma_occ * z` on the occurrence arm,
`alpha * sigma * z` on the cover arm) where `z` is a BYM2 latent field
on the supplied region adjacency. The cover-arm field amplitude
`sigma_pos = alpha * sigma` is the derived ratio surfaced as the
`"alpha"` column of `fit$joint$theta_grid`.

## Usage

``` r
simulate_cover_joint(
  N = 300L,
  adj,
  beta_occ = c(-0.3, 0.7),
  beta_pos = c(0.4, -0.5),
  sigma = 0.6,
  rho = 0.7,
  alpha = 1,
  positive = c("beta", "lognormal", "gaussian"),
  phi = 30,
  sigma_pos_resid = 0.4,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 300).

- adj:

  `n_s x n_s` integer adjacency matrix for the BYM2 field.

- beta_occ:

  Length-2 occurrence coefficients (intercept + slope).

- beta_pos:

  Length-2 cover-arm coefficients (intercept + slope).

- sigma:

  Marginal spatial-field amplitude (default 0.6); the occurrence-arm
  linear predictor adds `sigma * z[region]`.

- rho:

  BYM2 mixing parameter, in `[0, 1]` (default 0.7). `rho = 1` is pure
  ICAR; `rho = 0` is pure IID.

- alpha:

  Cover-arm scaling: cover-arm linear predictor adds
  `alpha * sigma * z[region]` (so `sigma_pos = alpha * sigma`).

- positive:

  Likelihood for the positive arm: `"beta"` or `"lognormal"`.

- phi:

  Beta precision when `positive = "beta"` (default 30).

- sigma_pos_resid:

  Lognormal residual SD when `positive = "lognormal"` (default 0.4);
  independent of the spatial `sigma`.

- seed:

  Optional integer seed.

## Value

A list with:

- data:

  Data frame with `x` and `region` (factor).

- y:

  Length-N numeric vector of cover values (zeros where `occur = 0`).

- adj:

  The adjacency matrix passed in (for downstream
  [`tulpa::spatial_bym2()`](https://gillescolling.com/tulpa/reference/spatial_bym2.html)).

- truth:

  Named list with `beta_occ`, `beta_pos`, `sigma`, `rho`, `alpha`,
  `sigma_pos = alpha * sigma`, `positive`, plus the simulated `phi_f`,
  `theta_f`, `w_s`, `region`.

## Details

Both BYM2 sub-blocks (`phi`, `theta`) are drawn IID Normal **and
demeaned to mean zero** before scaling, matching the sum-to-zero
constraint applied to both sub-blocks by
[`tulpa::spatial_bym2()`](https://gillescolling.com/tulpa/reference/spatial_bym2.html)
inside the nested-Laplace engine (see `.joint_inner_var()`). Without the
demean each seed carries `mean(w_s) ~ N(0, sigma^2 / n_s)` and the
constrained intercept identified by the engine targets
`beta_pos_0_truth + alpha * mean(w_s)` rather than `beta_pos_0_truth`;
coverage of the *population* truth then collapses with alpha while the
engine-reported posterior is correctly calibrated for the constrained
parameter. See `example/validation/SUMMARY.md` in INLAabun (Demo 3) for
the diagnosis.

## Examples

``` r
adj <- matrix(0L, 10, 10)
for (s in 1:10) for (j in setdiff(c(s - 1L, s + 1L), c(0L, 11L)))
  adj[s, j] <- 1L
sim <- simulate_cover_joint(N = 100, adj = adj, alpha = 1.0, seed = 1)
head(sim$data)
```
