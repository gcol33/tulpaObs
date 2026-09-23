# Time-to-detection occupancy family

Occupancy where a survey records the TIME to first detection rather than
a 0/1 outcome (Garrard et al. 2008; unmarked `occuTTD`). At an occupied
site the time-to-detection is exponential with rate `lambda` (constant
hazard); a survey of length `surveyLength` that reaches its end without
a detection is censored. An unoccupied site never detects. The state
`formula` models `logit psi` (occupancy); `detection` models
`log lambda` (the site-level detection rate). The latent occupancy state
integrates out in closed form (two states), so the fit maximises the
exact marginal with an observed-information vcov.

## Usage

``` r
occu_ttd(surveyLength = 1)
```

## Arguments

- surveyLength:

  Survey length `Tmax` (the censoring time): a scalar, a length-`N`
  vector, or an `N x J` matrix. Default 1.

## Value

A `tobs_family` object for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Details

`y` is an `N x J` matrix of detection times: a value in
`(0, surveyLength)` is a detection; a value `>= surveyLength` is a
non-detection (censored); `NA` is a survey not conducted.

## Examples

``` r
# \donttest{
sim <- simulate_occu_ttd(N = 200, J = 4, Tmax = 3, seed = 1)
fit <- tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 3),
            detection = ~ rate_cov1, y = sim$y, control = list(verbose = FALSE))
coef(fit)
# }
```
