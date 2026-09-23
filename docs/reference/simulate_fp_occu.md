# Simulate multistate false-positive occupancy data

Latent occupancy `z_i ~ Bernoulli(psi_i)` with
`logit psi = X_psi beta_psi`, observed through `J` replicate visits
emitting a detection state `y in {0,1,2}`: at an occupied site a visit
detects with probability `p11` and the detection is certain (state 2)
with probability `b` else ambiguous (state 1); at an unoccupied site a
visit yields a false-positive ambiguous detection (state 1) with
probability `p10`. Returns an `N x J` integer matrix in `{0, 1, 2}`
suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md).

## Usage

``` r
simulate_fp_occu(
  N = 300,
  J = 5,
  n_occ_covs = 1,
  beta_psi = NULL,
  p11 = 0.6,
  p10 = 0.05,
  b = 0.5,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 300).

- J:

  Number of visits (default 5).

- n_occ_covs:

  Number of occupancy covariates (default 1).

- beta_psi:

  Occupancy coefficients (logit). Default
  `c(qlogis(0.5), runif(n_occ_covs, -0.6, 0.6))`.

- p11, p10, b:

  True detection, false-positive, and certain-classification
  probabilities (scalars; defaults 0.6, 0.05, 0.5).

- seed:

  Optional random seed.

## Value

A list with `y` (N x J state matrix), `data` (occupancy covariates), and
`truth` (coefficients, per-site `psi`, scalar `p11`/`p10`/`b`, latent
`z`).

## Examples

``` r
sim <- simulate_fp_occu(N = 50, J = 4, seed = 1)
table(sim$y)
```
