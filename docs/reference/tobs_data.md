# Convert long-format data to a site x visit observation object

Convert long-format data to a site x visit observation object

## Usage

``` r
tobs_data(
  df,
  y,
  site,
  visit,
  type = c("occurrence", "abundance", "cover", "positive"),
  occ.covs = NULL,
  det.covs = NULL,
  coords = NULL,
  sites = NULL,
  visits = NULL,
  cover.floor = 0,
  compact = FALSE
)
```

## Arguments

- df:

  Data.frame in long format (one row per site-visit).

- y:

  Character, name of the response column. Its meaning is set by `type`:
  a 0/1 detection (`"occurrence"`), an integer count (`"abundance"`), or
  a continuous cover proportion in `[0, 1]` (`"cover"`).

- site:

  Character, name of site identifier column.

- visit:

  Character, name of visit/replicate column.

- type:

  Response kind. `"occurrence"` (default) and `"abundance"` build an
  integer site x visit matrix; `"cover"` and `"positive"` build a double
  matrix and preserve continuous values (no coercion to integer).
  `"cover"` is a proportion in `[0, 1]` (the beta cover arm);
  `"positive"` is a positive real `(0, Inf)` (the lognormal / gamma
  cover arm), validated as non-negative with no upper bound. Both apply
  the `cover.floor` absence policy below.

- occ.covs:

  Character vector of site-level covariate names.

- det.covs:

  Character vector of visit-level covariate names. A named column that
  is a factor or character is preserved as categorical: the downstream
  detection / positive design expands it to k - 1 dummies for a k-level
  factor, with the factor's first level (for a character column, the
  first level after sorting the unique values) as the reference. Numeric
  columns are kept as a continuous covariate.

- coords:

  Character vector of length 2 for coordinate columns.

- sites:

  Optional site level set. When supplied, the site x visit grid uses
  these site identifiers in this order (rather than `df`'s
  first-appearance order), so a subset of `df` pivots onto a fixed,
  externally-defined site grid; site values in `df` outside the set
  error.

- visits:

  Optional visit level set, analogous to `sites` for the visit axis
  (default the sorted unique visits in `df`).

- cover.floor:

  For `type = "cover"` / `"positive"`, the threshold at or below which a
  value is stored as `NA` rather than as a positive observation (default
  `0`). The positive arm of a hurdle is positive-only, so a `0` is an
  absence handled by the occurrence arm, not a positive observation;
  sending it to `NA` keeps a sampled-absent or unsampled cell from
  entering the positive arm as a fabricated zero (which, padded across a
  grid, flattens the spatial field). Set `cover.floor = -Inf` to keep
  every value verbatim.

- compact:

  Logical (default `FALSE`). When `TRUE`, return a compact (ragged)
  `tobs_data`: the response is stored as one row per valid site-visit (a
  `tobs_ragged` carrier) rather than as a padded
  `[n_sites x max_visits]` matrix, and each detection covariate as a
  length-V vector in the same canonical `order(site, visit)`. The
  compact layout has no per-site visit cap (its memory is the number of
  observations, not the padded grid) and is consumed by the joint
  nested-Laplace
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  engine, which works one valid visit at a time. Two compact calls on
  the same `df` / `site` / `visit` align row-for-row, so an occurrence
  response and a cover response can be paired directly.

## Value

An `tobs_data` object. With `compact = TRUE` its `$y` is a `tobs_ragged`
carrier and its `$det.covs` are length-V vectors.

## Examples

``` r
sim <- simulate_occu(N = 30, J = 3, seed = 1)
long <- data.frame(site = rep(seq_len(30), times = 3),
                   visit = rep(1:3, each = 30),
                   det = as.vector(sim$y),
                   occ_cov1 = sim$data$occ_cov1,
                   effort = runif(90))
dat <- tobs_data(long, y = "det", site = "site", visit = "visit",
                 occ.covs = "occ_cov1", det.covs = "effort")
dat
```
