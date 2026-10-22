# Joint species distribution family (no detection process)

Multivariate occurrence with shared latent factors. No observation
replication – treats observed presence/absence as the response.

## Usage

``` r
jsdm()
```

## Value

A `tobs_family` object.

## Engines

- `method = "laplace"` – the community Laplace-EM shared with
  [`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
  (a logit-link community GLMM);
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors alone (lfJSDM) run here through the block-coordinate driver.

- `method = "nuts"` – samples the exact joint community posterior
  (community means, per-species deviations, community covariance) over
  the Bernoulli response, warm-started at the Laplace-EM mode;
  non-spatial, an areal or
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  term errors with a pointer.

- `method = "pg_gibbs"` – a per-species Polya-Gamma conjugate Gibbs with
  community mean and Inverse-Gamma variance updates, giving a calibrated
  community-variance posterior; non-spatial.

- `method = "nested_laplace"` – a shared areal field
  ([`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)),
  with or without
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors (sfJSDM), by block coordinate ascent; the areal term is
  required on this route.

## Examples

``` r
f <- jsdm()
f
```
