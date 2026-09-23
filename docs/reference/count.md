# Count / relative-abundance GLMM family (no detection process)

A generalized linear model on an observed count (or continuous) response
directly, with no detection sub-model and no latent abundance to
marginalize. The abundance analogue of
[`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md): where
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md) fits
Royle's N-mixture (latent `N` plus imperfect detection), `count()` fits
the observed response as a plain GLMM. This is the "relative abundance"
model of `spAbundance` (`abund`): one value per site, `log`-link Poisson
/ negative-binomial counts or an identity-link Gaussian response.

## Usage

``` r
count(response = c("poisson", "negbin", "gaussian", "binomial"))
```

## Arguments

- response:

  The response distribution: `"poisson"` (log link), `"negbin"`
  (negative binomial, log link, an estimated size / dispersion),
  `"gaussian"` (identity link, an estimated residual variance), or
  `"binomial"` (logit link, `k` successes out of `n` trials per site).
  The binomial response is the detection-free binomial GLMM of
  `spOccupancy` (`svcPGBinom`): supply the per-site trial count as
  `trials =` on
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
  (default 1, i.e. a Bernoulli response). Unlike the count / Gaussian
  areal fit, a binomial areal field is identified (the variance is
  pinned by `n`), so an
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  /
  [`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  term composes with any `trials`.

## Value

A `tobs_family` object.

## Details

The response is a numeric vector (one value per site), supplied via
`y =` or on a two-sided `formula` left-hand side
(`count.value ~ predictors`).

## See also

[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md)
(N-mixture: latent abundance + detection),
[`jsdm()`](https://gillescolling.com/tulpaObs/reference/jsdm.md)
(occurrence GLMM, no detection).

## Examples

``` r
f <- count("poisson")
f
```
