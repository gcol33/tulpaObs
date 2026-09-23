# N-mixture abundance family

Latent Poisson (or NB) abundance with binomial detection per visit
(Royle 2004).

## Usage

``` r
abun(K_max = NULL, mixture = c("poisson", "negbin", "zip", "zinb"))
```

## Arguments

- K_max:

  upper bound for the latent-abundance marginal sum (the exact
  integration over `N` is truncated at `K_max`). `NULL` (default) lets
  the engine pick `max(y) + 100` (matching
  [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html));
  raise it if a fit warns that the posterior over `N` puts mass on the
  boundary.

- mixture:

  latent-abundance distribution, both fitted via tulpa's closed-form
  marginal Laplace engine. `"poisson"` is the default; `"negbin"`
  (negative binomial, `Var(N) = lambda + lambda^2 / r`) adds an
  overdispersion parameter `r`. Non-spatially the log size `log_r` is
  estimated jointly with the coefficients and reported with an SE; on
  the areal-spatial path (`method = "nested_laplace"`) `r` is integrated
  over the outer hyperparameter grid and reported as a posterior mean /
  sd. Named `mixture` (after
  [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html))
  to avoid collision with the model-type `family` argument of
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md).

## Value

A `tobs_family` object.

## Details

Royle's (2004) N-mixture model: latent abundance
`N_i ~ Poisson(lambda_i)` with `log lambda_i = X_lambda beta_lambda`,
and counts `y_ij | N_i ~ Binomial(N_i, p_ij)` with
`logit p_ij = X_p beta_p`. The abundance formula is the
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula`; the per-visit detection formula is `detection`. The marginal
likelihood integrates `N` out exactly, so the fit is a direct Laplace
approximation (no EM): `method = "laplace"` (fixed effects) or
`method = "nested_laplace"` (an areal
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
`car()` offset on the abundance arm).

## Examples

``` r
# \donttest{
sim <- simulate_abun(N = 120, J = 4, n_abund_covs = 1, n_det_covs = 1, seed = 1)
fit <- tobs(~ abund_cov1, data = sim$data, family = abun(),
            detection = ~ det_cov1, y = sim$y, method = "laplace")
summary(fit)
# }
```
