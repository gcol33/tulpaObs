# Simulate community (multispecies) joint occupancy-cover data

Per-species joint occupancy-detection + cover model with Gaussian
community hyperpriors on the per-species coefficients of all three arms:
`beta_occ_s ~ N(mu_occ, diag(sd_occ^2))`,
`beta_p_s ~ N(mu_p, diag(sd_p^2))`,
`beta_pos_s ~ N(mu_pos, diag(sd_pos^2))`, then per cell
`z_{s,i} ~ Bernoulli(psi_{s,i})`, per visit
`y_{s,i,j} | z = 1 ~ Bernoulli(p_{s,i,j})`, and per detected visit
`c_{s,i,j} ~ f_pos(eta_pos_{s,i,j}, disp)`. The site and visit
covariates are shared across species (community covariates). The
returned `y` / `y_pos` are 3D arrays `[n_sites x J x n_species]`
suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)
(`y_pos` is `NA` where not detected).

## Usage

``` r
simulate_ms_occu_cover(
  n_species = 12,
  N = 120,
  J = 5,
  n_occ_covs = 1,
  n_det_covs = 1,
  n_pos_covs = 1,
  mu_occ = NULL,
  mu_p = NULL,
  mu_pos = NULL,
  sd_occ = 0.5,
  sd_p = 0.4,
  sd_pos = 0.4,
  positive = c("lognormal", "beta", "gaussian"),
  phi = 30,
  sigma_pos = 0.4,
  seed = NULL
)
```

## Arguments

- n_species:

  Number of species (default 12).

- N:

  Number of sites / cells (default 120).

- J:

  Number of replicate visits (default 5).

- n_occ_covs, n_det_covs, n_pos_covs:

  Number of covariates on each arm (drawn IID standard normal, shared
  across species).

- mu_occ, mu_p, mu_pos:

  Community-mean coefficient vectors `c(intercept, slopes...)` on each
  arm's link scale. Defaults pick weakly-informative values.

- sd_occ, sd_p, sd_pos:

  Per-coefficient community SD on each arm (length 1, recycled, or one
  per coefficient). Default 0.5 / 0.4 / 0.4.

- positive:

  `"lognormal"` (default) or `"beta"`.

- phi:

  Beta precision when `positive = "beta"` (default 30).

- sigma_pos:

  Lognormal residual SD when `positive = "lognormal"` (default 0.4).

- seed:

  Optional integer seed.

## Value

A list with `y` (3D detection array), `y_pos` (3D cover array, `NA`
where not detected), `data` (per-cell covariate frame), `visit_data`
(per-visit covariate frame, `N*J` rows in site-major order), `species`
(species names), and `truth` (community means / SDs, per-species
coefficients, the dispersion, and the latent state).

## Examples

``` r
sim <- simulate_ms_occu_cover(n_species = 4, N = 30, J = 3, seed = 1)
dim(sim$y)
```
