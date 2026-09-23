# Multi-season occupancy with a temporal (AR1) year effect (tPGOcc)

The spOccupancy `tPGOcc` model: per-`(site, season)` Bernoulli occupancy
with a shared AR1 year random effect on the occupancy logit and NO
colonization / extinction dynamics (that is
[`dyn_occu()`](https://gillescolling.com/tulpaObs/reference/dyn_occu.md)).
Use it for an occupancy trend over years with a temporal random effect.
Fit with `method = "pg_gibbs"` (the Polya-Gamma Gibbs sampler
spOccupancy uses); `y` is a 3D array
`[n_sites x n_seasons x max_visits]` (or a list of per-season
`[n_sites x max_visits]` matrices), the occupancy `formula` and
`detection` formula are site-level.

## Usage

``` r
t_occu()
```

## Value

A `tobs_family` object.

## Examples

``` r
f <- t_occu()
f
```
