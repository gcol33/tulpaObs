# Suggest spatial and temporal aggregation for an identifiable occupancy model

Scans candidate spatial cell sizes and yearly clusterings of
long-format, typically single-visit plot data and scores each by how
well a single-season occupancy model separates occupancy (`psi`) from
detection (`p`). Replication is manufactured by pooling plots into a
spatial cell (assumes plots in a cell share `psi`) and years into a
contiguous block (assumes occupancy is closed across the block); an
occasion is any record in a (cell, block) bucket, so single-visit plots,
repeat-visited plots, and the no-pooling limit all flow through one code
path.

## Usage

``` r
occu_aggregation_scan(
  data,
  response,
  coords,
  year,
  plot = NULL,
  cell_sizes = NULL,
  block_lengths = NULL,
  score = c("info", "count"),
  family = c("occupancy"),
  control = list()
)
```

## Arguments

- data:

  Long data.frame, one row per plot-year record (or plot-year-visit).

- response:

  Name of the 0/1 detection column.

- coords:

  Length-2 character vector naming the x and y coordinate columns.

- year:

  Name of the integer year column.

- plot:

  Optional name of the plot-id column; enables the within-cell
  homogeneity proxy.

- cell_sizes:

  Numeric vector of candidate cell edge lengths. `NULL` auto-proposes a
  geometric ladder from the nearest-neighbour spacing to half the
  extent.

- block_lengths:

  Integer vector of candidate contiguous block lengths (in years).
  `NULL` runs mean-shift changepoint segmentation of the per-year naive
  occupancy series and uses the resulting blocks as the single temporal
  candidate.

- score:

  `"info"` (curvature-based, default) or `"count"` (structural).

- family:

  Observation family. Only `"occupancy"` is supported; the per-family
  constant-model scorer is the extension point for others.

- control:

  List of tuning knobs: `n_cell_sizes` (auto ladder length, 6),
  `nn_sample` (NN-spacing sample cap, 1500), `min_block` (min years per
  segment, 1), `cp_penalty` (explicit changepoint penalty, else auto),
  `cp_penalty_mult` (auto-penalty scale, 1), `max_condition`
  (conditioning cap for the identifiable flag, 1e6), `optim_control`
  (passed to `optim`).

## Value

A `tobs_aggregation_scan` object: `candidates` (scored data.frame),
`recommended` (chosen row or `NULL`), `segmentation` (auto changepoint
blocks, if used), and scan metadata.

## Details

Two scoring modes:

- `"info"` fits the constant (intercept-only) occupancy model per
  candidate and reports the posterior SEs of `psi`/`p` and the smallest
  eigenvalue and condition number of the 2x2 `(logit psi, logit p)`
  observed information. A confounding ridge shows up directly as a
  near-zero eigenvalue.

- `"count"` checks only the structural necessary conditions (units with
  two or more occasions, detected units, a non-degenerate naive
  detection rate).

The recommended candidate is the **least-pooling** identifiable one
(smallest mean occasions per unit, ties broken by largest information
eigenvalue): pool only as much as identifiability requires. When no
candidate is identifiable the data cannot support a standard occupancy
model and a Royle-Nichols or count model is the alternative.

## Examples

``` r
set.seed(1)
d <- data.frame(east = runif(200), north = runif(200),
                year = sample(2001:2006, 200, replace = TRUE),
                detected = rbinom(200, 1, 0.3))
scan <- occu_aggregation_scan(d, response = "detected",
                              coords = c("east", "north"), year = "year",
                              cell_sizes = c(0.2, 0.4),
                              block_lengths = c(2, 3))
scan
```
