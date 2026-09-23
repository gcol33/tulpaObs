# Binned distance-sampling family

Latent abundance `N_i ~ Poisson(lambda_i)` (or negative binomial) in a
covered region, observed through a half-normal or hazard-rate detection
function over distance bins. With `B` bins the detected counts are
multinomial over `(bin 1, ..., bin B, undetected)` with cell
probabilities `pi_b = integral_bin g(x; sigma) f(x) dx` and
`1 - sum_b pi_b`, where `f(x)` is the distance density (uniform for a
line transect, proportional to distance for a point transect). The
latent `N` is summed out in closed form (truncation `K_max`), so the fit
is a direct Laplace approximation (no EM), with a NUTS path over the
same marginal; the per-bin detection integrals are evaluated by
Gauss-Legendre quadrature.

## Usage

``` r
distance(
  key = c("halfnorm", "hazard"),
  transect = c("line", "point"),
  cutpoints = NULL,
  K_max = NULL,
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

The [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula` is the abundance (`log lambda`) model; `detection` is the
site-level detection-scale (`log sigma`) model. The response `y` is an
`n_sites x n_bins` integer matrix of per-bin detected counts. The bin
edges and transect geometry travel with the family:
`distance(cutpoints = ...)`.

## References

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L.,
Borchers, D. L., Thomas, L. (2001). Introduction to Distance Sampling.
Oxford. Royle, J. A., Dawson, D. K., Bates, S. (2004). Modeling
abundance effects in distance sampling. *Ecology* 85, 1591-1597.

## Examples

``` r
# \donttest{
sim <- simulate_distance(N = 200, key = "halfnorm", transect = "line",
                         n_abund_covs = 1, n_sigma_covs = 1, seed = 1)
fit <- tobs(~ abund_cov1, data = sim$data,
            family = distance(key = "halfnorm", transect = "line",
                              cutpoints = sim$cutpoints),
            detection = ~ sigma_cov1, y = sim$y, method = "laplace")
summary(fit)
# }
```
