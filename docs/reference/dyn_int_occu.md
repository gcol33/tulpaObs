# Multi-season integrated occupancy family

Dynamic (multi-season) occupancy observed by several detection sources
(spOccupancy `tIntPGOcc`): the product of a dynamic occupancy HMM
(season-1 occupancy `psi1`, colonization `gamma`, extinction `eps`) and
integrated occupancy (a per-season emission that pools multiple
detection sources). Pooling sources across seasons is how colonization /
extinction estimates come out of individually-sparse opportunistic data.
The latent occupancy sequence integrates out by the two-state HMM
forward recursion; the exact marginal is maximised with an
observed-information vcov.

## Usage

``` r
dyn_int_occu()
```

## Value

A `tobs_family` object for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Details

`y` is a length-`S` list of `[sites x visits x seasons]` detection
arrays (0/1/NA). The state `formula` models `logit psi1`;
`colonization = ~ ...` and `extinction = ~ ...` the site-level
transitions (required, as in
[`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md));
`detection` the shared per-source detection design (each source carries
its own coefficients). v1: every source covers all sites and the same
season grid, constant transitions, site-level detection.

## Examples

``` r
# \donttest{
sim <- simulate_dyn_int_occu(N = 200, T_seasons = 4, S = 2, seed = 1)
fit <- tobs(~ 1, data = sim$data, family = dyn_int_occu(),
            detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
            y = sim$y, sources = sim$sources, control = list(verbose = FALSE))
coef(fit)
# }
```
