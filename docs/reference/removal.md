# Removal-sampling family

Sequential-depletion removal sampling: latent abundance
`N_i ~ Poisson(lambda_i)` (or negative binomial) observed through `K`
ordered removal passes, where pass `k` removes
`Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the individuals still
present. The declining catch sequence identifies detection `p` and
abundance `N`; the latent `N` is summed out in closed form (truncation
`K.max`), so the fit is a direct Laplace approximation (no EM), with a
NUTS path over the same marginal.

## Usage

``` r
removal(K.max = NULL, mixture = c("poisson", "negbin"))
```

## Arguments

- K.max:

  upper bound for the latent-abundance marginal sum (the exact
  integration over `N` is truncated at `K.max`). `NULL` (default) lets
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

The [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula` is the abundance (`log lambda`) model; `detection` is the
per-pass detection (`logit p`) model. The response `y` is an
`n_sites x K` integer matrix of per-pass removals with the passes in
column order; complete pass sequences are required (no `NA`).

## Engines

- `method = "laplace"` – a direct Laplace approximation over the
  closed-form marginal (Poisson or negative binomial); a grouped random
  effect is integrated by adaptive Gauss-Hermite quadrature.

- `method = "nested_laplace"` – an areal
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  field on the abundance arm, optionally composed with a
  [`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  field, through the shared areal-BFGS driver.

- `method = "nuts"` – the in-tree C++ NUTS sampler over the same
  marginal, with a fixed-hyper areal or temporal field block on the
  abundance arm (one of the two, not both).

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Dorazio, R.
M., Jelks, H. L., Jordan, F. (2005). Improving removal-based estimates
of abundance. *Biometrics* 61, 1093-1101.

## Examples

``` r
f <- removal(K.max = 100)
f
```
