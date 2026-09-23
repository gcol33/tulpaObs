# Royle-Nichols occupancy family

Occupancy with abundance-induced detection heterogeneity (Royle &
Nichols 2003; unmarked `occuRN`). Latent abundance
`N_i ~ Poisson(lambda_i)` drives the per-visit detection probability
`1 - (1 - r_ij)^{N_i}`, where `r_ij` is the per-individual detection
probability. The state `formula` models `log lambda` (abundance);
`detection` models `logit r`. Detection is site-level by default;
passing `visits` makes it visit-varying (`logit r_ij` gains the
visit-level covariates), exactly as for the occupancy / N-mixture front
doors. The latent `N` marginalises in closed form (a Poisson sum to
`K_max`), so the fit maximises the exact marginal with an
observed-information vcov.

## Usage

``` r
royle_nichols(K_max = NULL)
```

## Arguments

- K_max:

  Upper summation bound for the latent abundance (default: a data-driven
  Poisson-tail guess).

## Value

A `tobs_family` object for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Examples

``` r
# \donttest{
sim <- simulate_royle_nichols(N = 150, J = 5, seed = 1)
fit <- tobs(~ x, data = sim$data, family = royle_nichols(),
            detection = ~ 1, y = sim$y, control = list(verbose = FALSE))
coef(fit)

# Visit-varying detection via `visits`:
sv <- simulate_royle_nichols(N = 150, J = 5, beta_r_visit = 0.8, seed = 1)
fv <- tobs(~ x, data = sv$data, family = royle_nichols(),
           detection = ~ w, y = sv$y, visits = sv$visits,
           control = list(verbose = FALSE))
coef(fv)
# }
```
