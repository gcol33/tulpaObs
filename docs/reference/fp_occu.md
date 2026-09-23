# Multistate false-positive occupancy family

Occupancy with both false negatives and false positives in the detection
process (e.g. acoustic classifiers, citizen-science misidentification),
using the Miller et al. (2011) confirmed-detection design that makes the
model robustly identifiable. Each visit yields a detection state
`y in {0, 1, 2}`: `0` = no detection, `1` = ambiguous detection (a true
detection OR a false positive), `2` = certain / confirmed detection
(only possible when the site is truly occupied). The latent occupancy
`z` marginalises in closed form (two states), so the fit maximises the
exact marginal likelihood directly (analytic-gradient BFGS,
observed-information covariance) with a NUTS path over the same
marginal.

## Usage

``` r
fp_occu()
```

## Value

A `tobs_family` object.

## Details

Four site-level logit arms: occupancy `psi` (the
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula`), true detection `p11` (`detection`), false-positive rate
`p10`, and the probability a true detection is certain `b`. The `p10`
and `b` predictors default to intercept-only and are set with the
`p10 = ~ ...` and `certainty = ~ ...` arguments to
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
(`certainty` is the `b` arm). The response `y` is an `n_sites x J`
integer matrix in `{0, 1, 2}` (NA visits dropped).

## References

Miller, D. A. W., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur. *Ecology* 92, 1422-1428. Royle,
J. A., Link, W. A. (2006). Generalized site occupancy models allowing
for false positive and false negative errors. *Ecology* 87, 835-841.

## Examples

``` r
f <- fp_occu()
f
```
