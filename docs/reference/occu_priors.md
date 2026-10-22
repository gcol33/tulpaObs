# Weakly-informative priors for occupancy Laplace fits

Constructs a prior specification consumed by
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) when a
Laplace method (`method = "laplace"` etc.) is used for the occupancy
family. The penalty is applied as a quadratic term
`+ sum((beta_j - mu_j)^2 / (2 * sd_j^2))` added to the
negative-log-posterior in each M-step submodel block by
[`tulpa::tulpa_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_laplace.html)
(via the per-block `beta_prior`).

## Usage

``` r
occu_priors(
  p.intercept = list(mean = 0, sd = 1.5),
  p.slope = list(mean = 0, sd = 2.5),
  beta.occ.intercept = list(mean = 0, sd = 2),
  beta.occ.slope = list(mean = 0, sd = 5)
)
```

## Arguments

- p.intercept:

  Prior on the detection intercept (logit scale). A list
  `list(mean = numeric(1), sd = numeric(1))`. Default
  `list(mean = 0, sd = 1.5)`.

- p.slope:

  Prior on detection slopes (every non-intercept detection coefficient).
  Default `list(mean = 0, sd = 2.5)`.

- beta.occ.intercept:

  Prior on the occupancy (psi / psi1) intercept. Default
  `list(mean = 0, sd = 2)`.

- beta.occ.slope:

  Prior on occupancy slopes. Default `list(mean = 0, sd = 5)`.

## Value

An `occu_priors` object, ready to pass to `tobs(..., priors = ...)`.

## Details

The defaults are weakly informative: they pull the detection intercept
toward `p = 0.5` strongly enough to break the psi-p identifiability
ridge at small `J`, but stay wide enough not to bias estimates at N \>=
600 with informative covariates. Setting any `sd` to `Inf` disables that
component of the prior (same penalised objective, just `1/Inf^2 = 0`).

## Examples

``` r
# disable the detection-slope penalty
priors <- occu_priors(p.slope = list(mean = 0, sd = Inf))
priors

# \donttest{
sim <- simulate_occu(N = 100, J = 3, n.occ.covs = 1, n.det.covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            priors = occu_priors(),
            control = list(verbose = FALSE, progress = FALSE))
coef(fit)
# }
```
