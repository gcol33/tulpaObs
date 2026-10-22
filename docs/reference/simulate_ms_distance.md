# Simulate community (multispecies) binned distance-sampling data

Per-species binned distance sampling with Gaussian community
hyperpriors: `beta_lambda_s ~ N(mu_lambda, diag(sd_lambda^2))`,
`beta_sigma_s ~ N(mu_sigma, diag(sd_sigma^2))`, then
`N_{s,i} ~ Poisson(lambda_{s,i})` and the detected individuals allocated
to distance bins by the half-normal detection function. The returned `y`
is a 3D array `[n_sites x n_bins x n_species]` suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`ms_distance()`](https://gillescolling.com/tulpaObs/reference/ms_distance.md).

## Usage

``` r
simulate_ms_distance(
  n.species = 10,
  N = 100,
  cutpoints = c(0, 25, 50, 75, 100),
  transect = c("line", "point"),
  key = c("halfnorm", "hazard"),
  shape = 0,
  n.abund.covs = 1,
  n.det.covs = 0,
  mu.lambda = NULL,
  mu.sigma = NULL,
  sd.lambda = 0.4,
  sd.sigma = 0.2,
  n.factors = 0,
  load.sd = 0.5,
  field = NULL,
  quad.order = 64L,
  seed = NULL
)
```

## Arguments

- n.species:

  Number of species (default 10).

- N:

  Number of sites (default 100).

- cutpoints:

  Distance-bin edges (default `c(0, 25, 50, 75, 100)`).

- transect:

  Transect geometry: `"line"` (default) or `"point"`.

- key:

  Detection function: `"halfnorm"` (default) or `"hazard"`.

- shape:

  Hazard-rate log-shape, shared across species (ignored under the
  half-normal key). Default 0.

- n.abund.covs, n.det.covs:

  Number of abundance / detection-scale covariates (default 1 and 0).

- mu.lambda:

  Community-mean abundance coefficients on the log scale. Default
  `c(log(30), rep(0.4, n.abund.covs))`.

- mu.sigma:

  Community-mean detection-scale coefficients on the log scale. Default
  `c(log(40), rep(0, n.det.covs))`.

- sd.lambda, sd.sigma:

  Per-coefficient community SDs. Default 0.4 and 0.2.

- n.factors:

  If `> 0`, add `n.factors` per-site latent factors with per-species
  loadings to `log lambda` (the lfMsDS truth). Default 0.

- load.sd:

  SD of the factor loadings (default 0.5).

- field:

  Optional length-`N` shared spatial field added to every species'
  `log lambda`. When given alongside `n.factors > 0` the loadings are
  centred across species, so the field owns the shared spatial mean.

- quad.order:

  Gauss-Legendre nodes per bin used to integrate the per-bin detection
  probabilities (default 64, matching
  [`ms_distance()`](https://gillescolling.com/tulpaObs/reference/ms_distance.md)).
  Set it to the `quad.order` the model will be fit at: the rule is
  `(cutpoints, transect, quad.order)`, so a different order integrates a
  different pi and the data would come from a model the fit does not
  use.

- seed:

  Optional random seed.

## Value

A list with `y`, `data`, `species`, `cutpoints`, and `truth` (community
means / SDs, per-species coefficients, `lambda`, `sigma`, and – with
factors – the `loadings`, `factors` and the implied residual correlation
`cor_res`). The latent `N` is drawn inside the shared C++ simulator and
is not returned.

## Examples

``` r
sim <- simulate_ms_distance(n.species = 4, N = 30, seed = 1)
dim(sim$y)
```
