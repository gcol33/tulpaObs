# Format multi-species occupancy data

Format multi-species occupancy data

## Usage

``` r
tobs_format_ms(
  y,
  occ.covs = NULL,
  det.covs = NULL,
  coords = NULL,
  species.names = NULL
)
```

## Arguments

- y:

  3D array (n_sites x max_visits x n_species) or named list of matrices.

- occ.covs:

  Data.frame of site-level covariates.

- det.covs:

  Named list of detection covariates.

- coords:

  Optional n_sites x 2 coordinate matrix.

- species.names:

  Optional character vector of species names.

## Value

An `tobs_data` object with multi-species structure.

## Examples

``` r
sim <- simulate_ms_occu(N = 40, J = 3, n.species = 4, seed = 1)
dat <- tobs_format_ms(sim$y, occ.covs = sim$data,
                      species.names = paste0("sp", 1:4))
dat$n_species
```
