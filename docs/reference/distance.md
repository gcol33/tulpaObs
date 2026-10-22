# Binned distance-sampling family

Latent abundance `N_i ~ Poisson(lambda_i)` (or negative binomial) in a
covered region, observed through a half-normal or hazard-rate detection
function over distance bins. With `B` bins the detected counts are
multinomial over `(bin 1, ..., bin B, undetected)` with cell
probabilities `pi_b = integral_bin g(x; sigma) f(x) dx` and
`1 - sum_b pi_b`, where `f(x)` is the distance density (uniform for a
line transect, proportional to distance for a point transect). The
latent `N` is summed out in closed form (truncation `K.max`), so the fit
is a direct Laplace approximation (no EM), with a NUTS path over the
same marginal; the per-bin detection integrals are evaluated by
Gauss-Legendre quadrature.

## Usage

``` r
distance(
  key = c("halfnorm", "hazard"),
  transect = c("line", "point"),
  cutpoints = NULL,
  K.max = NULL,
  mixture = c("poisson", "negbin")
)
```

## Arguments

- key:

  detection-function key. `"halfnorm"` (default) or `"hazard"` (the
  hazard-rate shape `b` is estimated as a scalar, reported as
  `log_shape`).

- transect:

  `"line"` (default; perpendicular distances uniform) or `"point"`
  (radial distances, density proportional to distance).

- cutpoints:

  numeric distance-bin edges, length `n_bins + 1`
  (`0 = c_0 < c_1 < ... < c_B`). Required.

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
site-level detection-scale (`log sigma`) model. The response `y` is an
`n_sites x n_bins` integer matrix of per-bin detected counts. The bin
edges and transect geometry travel with the family:
`distance(cutpoints = ...)`.

## Engines

- `method = "laplace"` – a direct Laplace approximation over the
  closed-form marginal (Poisson or negative binomial); a grouped random
  effect on the abundance arm is integrated by adaptive Gauss-Hermite
  quadrature.

- `method = "nested_laplace"` – an areal
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  field on the abundance arm, optionally composed with a
  [`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  field, through the shared areal-BFGS driver; a hazard-rate fit
  recovers the shape alongside the field.

- `method = "nuts"` – the in-tree C++ NUTS sampler over the same
  marginal, with a fixed-hyper areal or temporal field block on the
  abundance arm (one of the two, not both); the hazard-rate log-shape
  rides alongside.

## References

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L.,
Borchers, D. L., Thomas, L. (2001). Introduction to Distance Sampling.
Oxford. Royle, J. A., Dawson, D. K., Bates, S. (2004). Modeling
abundance effects in distance sampling. *Ecology* 85, 1591-1597.

## Examples

``` r
# \donttest{
sim <- simulate_distance(N = 200, key = "halfnorm", transect = "line",
                         n.abund.covs = 1, n.sigma.covs = 1, seed = 1)
fit <- tobs(~ abund_cov1, data = sim$data,
            family = distance(key = "halfnorm", transect = "line",
                              cutpoints = sim$cutpoints),
            detection = ~ sigma_cov1, y = sim$y, method = "laplace")
summary(fit)
# }
```
