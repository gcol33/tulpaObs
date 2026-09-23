# Comprehensive model checking

The
[`tulpa::check_model()`](https://gillescolling.com/tulpa/reference/check_model.html)
method for `tobs_fit`: the criteria, the posterior-predictive check and
the goodness-of-fit tests in one call, printed as a report and drawn as
the diagnostic panel.

## Usage

``` r
# S3 method for class 'tobs_fit'
check_model(
  object,
  coords = NULL,
  n.samples = 250,
  nsim = NULL,
  seed = NULL,
  plot = TRUE,
  ...
)
```

## Arguments

- object:

  A `tobs_fit` object.

- coords:

  Optional `n_sites x 2` coordinate matrix. Adds Moran's I on the
  occupancy residuals and the correlogram panel.

- n.samples:

  Posterior samples for simulation tests.

- nsim:

  The same budget under the name
  [`tulpa::check_model()`](https://gillescolling.com/tulpa/reference/check_model.html)'s
  default uses; overrides `n.samples` when given.

- seed:

  Optional random seed for the simulation tests; the caller's RNG stream
  is restored afterwards.

- plot:

  Draw the panel. `FALSE` prints the report alone, for a log or a
  headless run.

- ...:

  Must be empty; an unrecognised argument is an error.

## Value

Invisibly, the diagnostic results: `waic`, `dic`, `cpo`, `ppc`,
`zero_inflation`, `dispersion`, `pit`, `uniformity`, and `moran` when
`coords` is supplied.

## Details

The engine's default method builds its panel from
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) as numeric
vectors, which a latent-state fit does not have –
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns
`list(psi, p, z)` and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) one series per
process, so the response cannot be resolved. This method reads the same
quantities through the family's own doors instead:
[`pit_residuals()`](https://gillescolling.com/tulpa/reference/pit_residuals.html),
[`ppc()`](https://gillescolling.com/tulpaObs/reference/ppc.md),
[`test_dispersion()`](https://gillescolling.com/tulpa/reference/test_dispersion.html)
and
[`tulpa::moran_i()`](https://gillescolling.com/tulpa/reference/moran_i.html).
