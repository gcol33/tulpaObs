# Community (multispecies) binned distance-sampling family

The community analogue of
[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md):
per-species binned distance sampling with Gaussian community hyperpriors
on the per-species abundance and detection-scale coefficients, so rare
species borrow strength from common ones through the shared community
means and covariances – the same pooling that stabilises
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md),
on the distance-sampling response.

## Usage

``` r
ms_distance(
  key = c("halfnorm", "hazard"),
  transect = c("line", "point"),
  cutpoints = NULL,
  K_max = NULL,
  mixture = c("poisson", "negbin")
)
```

## Arguments

- key:

  Detection function: `"halfnorm"` (default) or `"hazard"`. Under the
  hazard-rate key the scalar log-shape is shared across species.

- transect:

  Transect geometry: `"line"` (default) or `"point"`.

- cutpoints:

  Distance-bin edges, length `dim(y)[2] + 1`, strictly increasing and
  starting at `>= 0`.

- K_max:

  Truncation for the latent abundance sum. Defaults to
  `3 * max(rowSums(y)) + 100`.

- mixture:

  Abundance mixing distribution. `"poisson"` only; the negative-binomial
  size is not yet carried as a per-species random effect.

## Value

A `tobs_family` object.

## Details

Per species `s`, site `i`, distance bin `b`:

    N_{s,i}        ~ Poisson(lambda_{s,i})
    y_{s,i,.} | N  ~ Multinomial(N_{s,i}; pi_{s,i,1..B}, 1 - p_det)
    log lambda_{s,i} = X_i     . (mu_lambda + b_lambda_s)
    log sigma_{s,i}  = X_sig_i . (mu_sigma  + b_sigma_s)
    b_lambda_s ~ N(0, Sigma_lambda),  b_sigma_s ~ N(0, Sigma_sigma)

with `pi_b` the integral of the detection function over bin `b`. The
latent `N_{s,i}` integrates out per species-site in closed form, exactly
as for
[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md);
the per-species coefficient deviations are the random effects,
integrated by the shared community Laplace-EM.

Adding `latent(n)` to the abundance formula gives residual species
co-occurrence through `n` per-site latent factors with per-species
loadings; adding an areal term
([`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md))
or
[`spde()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
alongside gives a shared spatial field. Both route through the same
block-coordinate engine as the other community families.

## Inputs

`y` is a 3D array `[n_sites x n_bins x n_species]` or a named list of
`n_sites x n_bins` per-bin count matrices. `formula` is the abundance
(`log lambda`) predictor, `detection` the detection-scale (`log sigma`)
predictor, and `species` the species labels.

## See also

[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md)
(single species),
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
(community N-mixture).

## Examples

``` r
f <- ms_distance(key = "halfnorm", transect = "line",
                 cutpoints = c(0, 25, 50, 75, 100))
f
```
