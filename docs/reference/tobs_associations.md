# Residual species-association matrices from a spatial-factor community fit

For a reduced-rank spatial-factor community occupancy-cover fit (a
[`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md)
model carrying an
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md),
or
[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
field on the occupancy arm), return the residual species associations
the shared latent fields imply – the spatial-JSDM / HMSC output. The K
unit-scale fields induce a residual occupancy covariance `L L'` across
species; the reported matrix is its correlation form, identified
(invariant to the factor rotation, sign, and column scale). When the fit
also carries a cover-arm factor, the cover association
`corr(L_pos L_pos')` and the joint cross-arm association (standardized
`L_occ L_pos'`, species occupancy vs species cover) are available too.

## Usage

``` r
tobs_associations(
  object,
  type = c("occupancy", "cover", "cross"),
  summary = NULL
)
```

## Arguments

- object:

  A fitted `tobs_fit` from a spatial-factor
  [`ms_occu_cover()`](https://gillescolling.com/tulpaObs/reference/ms_occu_cover.md).

- type:

  Which association to return: `"occupancy"` (default), `"cover"`, or
  `"cross"` (occupancy vs cover). `"cover"` and `"cross"` require a
  cover-arm factor (an
  [`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)/[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)/[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
  term on the cover formula).

- summary:

  One of `"median"`, `"estimate"` (rotation-invariant point estimate at
  the posterior mode), `"lower"`, `"upper"`; or `NULL` (default) to
  return the full list of all four `S x S` matrices plus the interval
  `prob` and draw count.

## Value

An `S x S` correlation matrix (when `summary` is given) or a list of the
point estimate, posterior median, and interval bounds.

## Details

Each matrix is summarised from draws of the loading posterior (not the
plug-in mode), so a central interval accompanies the point estimate.

## Examples

``` r
# \donttest{
adj <- matrix(0, 16, 16)
adj[cbind(1:15, 2:16)] <- 1
adj <- adj + t(adj)
sim <- simulate_ms_occu_cover_spatial(adj, n_species = 6, K = 2, J = 4,
                                      seed = 1)
fit <- tobs(~ occ_cov1 + icar(graph = adj), data = sim$data,
            family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
            positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
            species = sim$species, method = "laplace",
            control = list(n.factors = 2, max.iter = 30, verbose = FALSE))
round(tobs_associations(fit, summary = "median"), 2)
# }
```
