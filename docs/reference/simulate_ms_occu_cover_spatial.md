# Simulate a reduced-rank spatial-factor community occu_cover data set (K = 1)

Generates occupancy / detection / cover data for `n_species` species
sharing one latent ICAR spatial factor `w` on the cell graph, with
per-species loadings `L_s` on the occupancy state predictor and Gaussian
community priors on the per-species arm coefficients. This is the
ground-truth generator for the Stage-1 reduced-rank spatial JSDM.

## Usage

``` r
simulate_ms_occu_cover_spatial(
  adj,
  n_species = 8L,
  J = 4L,
  K = 1L,
  n_occ_covs = 1L,
  n_det_covs = 1L,
  n_pos_covs = 1L,
  mu_occ = NULL,
  mu_p = NULL,
  mu_pos = NULL,
  sd_occ = 0.4,
  sd_p = 0.4,
  sd_pos = 0.3,
  mean_load = 0,
  sd_load = 1,
  cover_factor = FALSE,
  mean_load_pos = 0,
  sd_load_pos = 1,
  field = c("icar", "car_proper", "bym2"),
  rho = 0.9,
  phi = 0.7,
  sigma_pos = 0.4,
  positive = c("lognormal", "beta"),
  seed = NULL
)
```

## Arguments

- adj:

  N x N 0/1 adjacency matrix of the cell graph (required); `N` cells.

- n_species:

  Number of species.

- J:

  Number of detection visits per cell.

- K:

  Number of shared latent spatial factors (`1 <= K <= n_species`).
  `K = 1` is the Stage-1 single-field case (loading vector, field
  vector); `K > 1` draws `K` ICAR fields with lower-triangular,
  positive-diagonal loadings and returns the `S x K` loading matrix /
  `N x K` field matrix.

- n_occ_covs, n_det_covs, n_pos_covs:

  Number of (Gaussian) covariates on the occupancy, detection, and cover
  arms; each arm also has an intercept.

- mu_occ, mu_p, mu_pos:

  Community mean coefficient vectors (intercept first). `NULL` picks
  sensible defaults of the right length.

- sd_occ, sd_p, sd_pos:

  Community RE SDs (diagonal `Sigma_.`); scalar (recycled) or
  per-coefficient.

- mean_load, sd_load:

  Mean and SD of the per-species loadings `L_s`.

- cover_factor:

  Logical; when `TRUE` the same shared fields `W` also load on the cover
  (positive) predictor through a free `S x K` loading matrix `L_pos`
  (the cover-arm factor). The cover-factor draws are gated, so `FALSE`
  (the default) reproduces the no-factor RNG stream exactly.
  `truth$L_pos` carries the generating cover loadings.

- mean_load_pos, sd_load_pos:

  Mean and SD of the cover-arm loadings `L_pos` (used only when
  `cover_factor = TRUE`).

- field:

  Areal structure of the shared latent factors: `"icar"` (improper
  intrinsic CAR, the default), `"car_proper"` (proper CAR with a
  correlation `rho`, field precision `tau (D - rho A)`), or `"bym2"`
  (the Riebler 2016 reparameterised convolution with a spatial-variance
  fraction `phi`). `"icar"` reproduces the Stage 1-2 RNG stream byte for
  byte; `truth$field` / `truth$rho` / `truth$phi` carry the choice.

- rho:

  Proper-CAR correlation in `[0, 1)`, used only when
  `field = "car_proper"` (the `rho -> 1` limit is the ICAR field).

- phi:

  BYM2 spatial-variance fraction in `[0, 1]`, used only when
  `field = "bym2"` (`phi = 1` is the pure ICAR field, `phi = 0` pure
  iid).

- sigma_pos:

  Lognormal cover residual SD (on the log scale).

- positive:

  Cover family; only `"lognormal"` in Stage 1.

- seed:

  Optional RNG seed.

## Value

A list with `y` (N x J x n_species 0/1 detections), `y_pos` (N x J x
n_species cover, `NA` off detected visits), `data` (cell-level covariate
frame), `species`, and `truth` (the canonical-form generating
parameters: `mu_*`, `sd_*`, per-species `b_*`, loadings `L`, field `w`,
`psi`, `z`).

## Examples

``` r
adj <- matrix(0, 10, 10)
adj[cbind(1:9, 2:10)] <- 1
adj <- adj + t(adj)
sim <- simulate_ms_occu_cover_spatial(adj, n_species = 4, J = 3, seed = 1)
dim(sim$y)
```
