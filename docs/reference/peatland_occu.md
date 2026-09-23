# Peatland occupancy survey (synthetic)

A single-season occupancy detection-history dataset for 120 wetland
sites surveyed on 4 visits, drawn from a known generative model
(MacKenzie et al. 2002). Occupancy increases with site `wetness` and
decreases with `elevation`; per-visit detection increases with survey
`effort`. Twelve detection cells are set to `NA` to mimic incomplete
histories. Synthetic, not field data; the generative coefficients are in
`truth`.

## Usage

``` r
peatland_occu
```

## Format

A list with components:

- y:

  120 x 4 integer matrix of detections (0/1/`NA`).

- occ.covs:

  Data frame of site covariates: `elevation`, `wetness` (both
  standardised).

- det.covs:

  Named list with one 120 x 4 visit covariate matrix, `effort`
  (standardised).

- coords:

  120 x 2 matrix of `x` / `y` coordinates in the unit square.

- truth:

  List of generative quantities: `beta_psi`, `beta_p`, the realised
  latent occupancy `z`, and per-site `psi`.

## See also

[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md),
[`simulate_occu()`](https://gillescolling.com/tulpaObs/reference/simulate_occu.md)

## Examples

``` r
data(peatland_occu)
str(peatland_occu$y)
# \donttest{
fit <- tobs(~ elevation + wetness, data = peatland_occu$occ.covs,
            family = occu(), detection = ~ effort,
            y = peatland_occu$y, visits = peatland_occu$det.covs)
coef(fit)
# }
```
