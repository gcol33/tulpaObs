# Multispecies (community) relative-abundance GLMM family

The community analogue of
[`count()`](https://gillescolling.com/tulpaObs/reference/count.md): a
per-species GLMM on an observed count or continuous response with
Gaussian community hyperpriors on the per-species coefficients (the
spAbundance `msAbund` model). No detection and no latent state – the
abundance counterpart of
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)
without the occupancy layer. Poisson / negative-binomial (log link) or
Gaussian (identity). `y` is an `n_sites x n_species` matrix (or a named
list of `n_species` count vectors), one value per site per species;
`species` names the columns.

## Usage

``` r
ms_count(response = c("poisson", "negbin", "gaussian", "binomial"))
```

## Arguments

- response:

  One of `"poisson"`, `"negbin"`, `"gaussian"`, or `"binomial"`. The
  binomial response is the community `k`-of-`n` GLMM (community
  `svcPGBinom`): supply the per-site (or per-`site x species`) trial
  count as `trials =` on
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
  (default 1, i.e. Bernoulli, which is the
  [`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md)
  response). With `trials > 1` the community-mean intercept carries a
  small first-order-Laplace bias of order `1 / n_species` (a few
  hundredths on the logit scale at 20 species, shrinking with more
  species; the slope and the `trials = 1` case are unbiased) – the same
  character as the negative-binomial slope attenuation noted above.

## Value

A `tobs_family` object.

## Details

As with the other community Laplace-EM families
([`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
[`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)),
the community means are recovered essentially unbiased for the Gaussian
(identity link, exact Laplace) and Poisson responses; the
negative-binomial slope carries a mild first-order-Laplace (PQL)
attenuation of order the community variance (a few percent, shrinking
with the number of species and observations per species). The
community-mean Wald intervals are calibrated to the package rubric
(pooled coverage at least 0.85).

`method = "nuts"` samples the exact joint posterior (community means,
per-species deviations, and the community covariance) for all three
responses – the negative binomial carrying a per-species dispersion
random effect, the Gaussian a per-species free residual variance – which
removes the Laplace-EM's negative-binomial attenuation and returns
calibrated, non-Gaussian community intervals. A shared areal field or
latent factors are available through `method = "nested_laplace"` / the
[`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term; see the package overview for the spatial and factor variants.

## Engines

- `method = "laplace"` – the shared community Laplace-EM, the default;
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors alone run here through the block-coordinate driver.

- `method = "nested_laplace"` – a shared areal field
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  on the abundance formula (the spAbundance `sfMsAbund` model), with or
  without
  [`latent()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  factors, by block coordinate ascent; the field is required on this
  route.

- `method = "nuts"` – the exact joint community posterior sampled by the
  in-tree C++ target, warm-started at the Laplace-EM mode, as described
  above; non-spatial.

- `method = "pg_gibbs"` – a per-species Polya-Gamma conjugate Gibbs with
  community mean and Inverse-Gamma variance updates, for the logistic
  `"binomial"` response only (another response errors with a pointer);
  non-spatial.

## See also

[`count()`](https://gillescolling.com/tulpaObs/reference/count.md)
(single species),
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md),
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)

## Examples

``` r
# \donttest{
sim <- simulate_ms_count(N = 120, n.species = 8, seed = 1)
fit <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
            species = colnames(sim$y), method = "laplace")
summary(fit)
# }
```
