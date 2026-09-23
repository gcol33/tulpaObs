# Joint distance + removal sampling family

The `unmarked` `gdistremoval` model (Amundson et al. 2014): a
**single-season** point-count design in which the detected individuals
are recorded two ways at once – their distance band (distance sampling)
and the removal period of first detection (removal sampling). This is
NOT the open-population distance model (`distsampOpen`); the abundance
is static.

## Usage

``` r
gdistremoval(transect = c("line", "point"), cutpoints = NULL)
```

## Arguments

- transect:

  Transect geometry: `"line"` (default) or `"point"`.

- cutpoints:

  Distance-bin edges, length `ncol(y) + 1`, strictly increasing and
  starting at `>= 0`.

## Value

A `tobs_family` object.

## Details

Site abundance `N_i ~ Poisson(lambda_i)`; the detected birds are
cross-classified by a distance band and a removal period. Writing
`pdist_i` for the overall distance detection, `prem_i` for the overall
removal detection, the total detected is a binomial thinning of `N_i`,
and Poisson is closed under binomial thinning, so the marginal is
closed-form: \$\$\mathrm{ysum}\_i \sim \mathrm{Poisson}(\lambda_i\\
p^{dist}\_i\\ p^{rem}\_i),\$\$ with the distance-band counts and
removal-period counts as two conditional multinomials – the
[`double_observer()`](https://gillescolling.com/tulpaObs/reference/double_observer.md)
Poisson-multinomial pattern, here with a distance multinomial
(half-normal key band integrals) and a depleting-removal multinomial
`pi_k = r (1 - r)^{k-1}`.

Three site-level arms: log abundance `lambda` (the
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula`), log distance scale `sigma` (`detection`), and logit
per-period removal capture `r` (`removal = ~ ...`, default
intercept-only).

## Inputs

`y` is an `n_sites x n_bins` integer matrix of per-distance-band counts;
`y_rem` an `n_sites x n_periods` integer matrix of per-removal-period
counts. The per-site row totals must match (the same detected birds
cross-classified).

## References

Amundson, C. L., Royle, J. A., Handel, C. M. (2014). A hierarchical
model combining distance sampling and time removal to estimate detection
probability during avian point counts. *The Auk* 131, 476-494.

## See also

[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md),
[`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md),
[`double_observer()`](https://gillescolling.com/tulpaObs/reference/double_observer.md).

## Examples

``` r
f <- gdistremoval(transect = "point", cutpoints = c(0, 10, 20, 30, 40))
f
```
