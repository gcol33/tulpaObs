# Single-season occupancy family

Latent Bernoulli occupancy state with binomial detection per visit. The
MacKenzie et al. (2002) single-season model.

## Usage

``` r
occu()
```

## Value

A `tobs_family` object.

## Engines

- `method = "laplace"` – penalized EM + Laplace with Gaussian marginals,
  the default.

- `method = "laplace_sla"` – the Laplace fit with simplified-Laplace
  skew-corrected marginals.

- `method = "laplace_gibbs"` / `method = "laplace_mi"` – the Laplace fit
  with a post-EM Gibbs / multiple-imputation variance correction
  (Rubin-pooled refits, not an MCMC chain).

- `method = "pg_gibbs"` – a Polya-Gamma Gibbs sampler over the exact
  posterior (the spOccupancy `PGOcc` engine): site-level detection, an
  areal
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  field is accepted; a
  [`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  / [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`svc()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  term errors with a pointer.

- `method = "nested_laplace"` – multi-block nested Laplace: areal,
  temporal and random-effect hyperparameters are integrated on an outer
  grid, and [`predict()`](https://rdrr.io/r/stats/predict.html)
  interpolates all-NA detection sites from the field.

- `method = "nested_laplace_sla"` – the nested route with skew-corrected
  marginals.

- `method = "nuts"` – the unified C++ NUTS sampler over the structured
  terms (spatial, temporal,
  [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
  [`svc()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)),
  reporting Rhat / ESS.

## Examples

``` r
f <- occu()
f
```
