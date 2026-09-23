# Simulate a multi-season integrated occupancy data set

Draws from the
[`dyn_int_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_int_occu.md)
model: a dynamic (multi-season) occupancy process (`psi1`, colonization
`gamma`, extinction `eps`) observed by `S` detection sources with
per-source detection probability.

## Usage

``` r
simulate_dyn_int_occu(
  N = 200,
  T_seasons = 4,
  S = 2,
  J = 3,
  psi1 = 0.5,
  gamma = 0.3,
  eps = 0.2,
  p = c(0.4, 0.6),
  field = NULL,
  trend = NULL,
  source_seasons = NULL,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- T_seasons:

  Number of seasons (default 4).

- S:

  Number of detection sources (default 2).

- J:

  Visits per source (scalar or length-`S`; default 3).

- psi1, gamma, eps:

  Season-1 occupancy, colonization, extinction (defaults 0.5 / 0.3 /
  0.2).

- p:

  Per-source detection probabilities (length-`S`; default 0.4 / 0.6).

- field:

  Optional per-site shared areal field (length `N`) added to the
  first-season occupancy logit – the shared field of the multi-season
  integrated spatial model (stIntPGOcc). Default `NULL` (no field).

- trend:

  Optional per-site varying-coefficient (SVC) areal field (length `N`);
  with a covariate `w ~ N(0,1)` stored in `data`, the first-season logit
  gains `w * trend` on top of `field` – the svcTIntPGOcc surface.
  Default `NULL`. When set, `data` carries a `cell` node index (`1..N`)
  and the covariate `w` for the bar `spatial(~ 1 + w || cell, graph)`.

- source_seasons:

  Optional length-`S` list; `source_seasons[[s]]` is the integer vector
  of seasons source `s` observes (partial season overlap). The seasons a
  source does not cover are set to `NA` in its array – the staggered
  survey where sources rarely share the full season grid. Default `NULL`
  (every source observes every season).

- seed:

  Optional random seed.

## Value

A list with `y` (a length-`S` list of `[N x J x T]` arrays), `data`,
`sources`, and `truth`.

## Examples

``` r
sim <- simulate_dyn_int_occu(N = 40, T_seasons = 3, seed = 1)
lapply(sim$y, dim)
```
