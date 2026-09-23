# Format occupancy data from matrices/lists

Format occupancy data from matrices/lists

## Usage

``` r
tobs_format(y, occ.covs = NULL, det.covs = NULL, coords = NULL, species = NULL)
```

## Arguments

- y:

  Detection history matrix (n_sites x max_visits) or 3D array.

- occ.covs:

  Data.frame of site-level covariates, or named list.

- det.covs:

  Named list of detection covariates. Each element is a vector (length
  n_sites, constant across visits) or matrix (n_sites x max_visits).

- coords:

  Optional n_sites x 2 coordinate matrix.

- species:

  Optional character or integer species identifier.

## Value

An `tobs_data` object.

## Examples

``` r
sim <- simulate_occu(N = 50, J = 3, seed = 1)
dat <- tobs_format(sim$y,
                   occ.covs = sim$data[, c("occ_cov1", "occ_cov2")],
                   det.covs = list(det_cov1 = sim$data$det_cov1))
dat
```
