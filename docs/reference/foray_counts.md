# Repeated point-count abundance survey (synthetic)

A binomial N-mixture dataset (Royle 2004) of repeated counts at 100
sites over 3 visits, drawn from a known generative model. Latent
abundance increases with `shrub` cover and decreases with `elevation`;
per-visit detection increases with survey `effort`. Synthetic, not field
data; the generative coefficients and the realised latent abundances are
in `truth`.

## Usage

``` r
foray_counts
```

## Format

A list with components:

- y:

  100 x 3 integer matrix of counts.

- occ.covs:

  Data frame of site covariates: `elevation`, `shrub` (both
  standardised).

- det.covs:

  Named list with one 100 x 3 visit covariate matrix, `effort`
  (standardised).

- truth:

  List of generative quantities: `beta_lambda`, `beta_p`, the realised
  latent abundance `N`, and per-site `lambda`.

## See also

[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md),
[`simulate_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_abun.md)

## Examples

``` r
data(foray_counts)
summary(as.vector(foray_counts$y))
# \donttest{
fit <- tobs(~ elevation + shrub, data = foray_counts$occ.covs,
            family = abun(), detection = ~ effort,
            y = foray_counts$y, visits = foray_counts$det.covs)
coef(fit)
# }
```
