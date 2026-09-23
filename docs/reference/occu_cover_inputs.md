# Build occu_cover() model inputs from a long, plot-level frame

Assemble the paired occurrence / cover response arms, the visit-level
covariate grid, and the site-level design frame that an
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
fit consumes, from one long (one row per site-visit / plot) data frame.
This is the same construction
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) runs
internally when a single
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
fit is handed a long frame; call it directly when you want to inspect
the arms before fitting, or pass them on as `y = `, `y_pos = `,
`visits = `, `data = ` to
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

The occurrence and cover arms are pivoted onto one shared
`(sites, visits)` grid, so they align row-for-row, and the site-level
design is the first row per site. The cover arm is positive-only
(meaningful where the species is detected): in the dense grid absent /
unsampled cover cells are set to `0`, and in the compact (ragged)
carrier they are `NA` in the values and never read by the joint engine.

## Usage

``` r
occu_cover_inputs(
  data,
  site,
  visit,
  response,
  y_pos,
  occ.covs = NULL,
  det.covs = NULL,
  coords = NULL,
  compact = TRUE,
  positive = "beta"
)
```

## Arguments

- data:

  A long / plot-level data frame: one row per site-visit.

- site, visit, response, y_pos:

  Column names: the site identifier, the visit / replicate, the 0/1
  detection response, and the continuous cover response (used only where
  `response == 1`). The accepted range of `y_pos` follows `positive`.

- occ.covs:

  Optional character vector of site-level covariate columns. The default
  (`NULL`) keeps the full first-row-per-site frame and lets the
  occurrence / detection / positive formulas select what they reference.

- det.covs:

  Optional character vector of visit-level covariate columns.

- coords:

  Optional length-2 character vector of coordinate columns.

- compact:

  Logical (default `TRUE`). Build the compact (ragged) arms the joint
  nested-Laplace
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
  engine consumes (no per-site visit cap), or the dense
  `[n_sites x max_visits]` grid.

- positive:

  The cover-arm distribution, matching
  [`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md):
  `"beta"` (default) stores `y_pos` as a proportion in `[0, 1]`
  ([`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md)
  type `"cover"`); `"lognormal"` / `"gamma"` store it as a positive real
  `(0, Inf)` (type `"positive"`, validated non-negative with no upper
  bound). Pick the value you pass to `occu_cover(response = )`.

## Value

A list with `y`, `y_pos`, `visits`, `site_data`, `coords`, `n_sites`,
`max_visits`, `n_visits`, and `compact`.

## See also

[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md),
[`tobs_data()`](https://gillescolling.com/tulpaObs/reference/tobs_data.md).

## Examples

``` r
sim <- simulate_occu_cover(N = 20, J = 3, seed = 1)
long <- data.frame(site = rep(1:20, each = 3), visit = rep(1:3, 20),
                   occur = as.vector(t(sim$y)),
                   cover = as.vector(t(ifelse(is.na(sim$y_pos), 0,
                                              sim$y_pos))),
                   det_cov1 = sim$visit_data$det_cov1,
                   occ_cov1 = rep(sim$data$occ_cov1, each = 3))
inp <- occu_cover_inputs(long, site = "site", visit = "visit",
                         response = "occur", y_pos = "cover",
                         det.covs = "det_cov1", compact = FALSE,
                         positive = "lognormal")
dim(inp$y)
```
