# Internal engine entry point

Dispatches to Laplace (default) or NUTS for a built `tobs_model`. Not
user-facing; called from
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) via the
per-family `.dispatch_*` helpers. Spatial / temporal / random-effect /
SVC / latent structure is read from the structured terms the formula
carried (`model$structured_terms`), not from arguments – there is a
single user-facing specification path.

## Usage

``` r
.tobs_fit_model(
  model,
  method = c("laplace", "nested_laplace", "nuts", "pg_gibbs"),
  priors = NULL,
  sigma.beta = NULL,
  max.iter = 100L,
  tol = 1e-04,
  damping = 0.7,
  n.iter = NULL,
  n.warmup = NULL,
  n.thin = NULL,
  n.chains = NULL,
  n.threads = NULL,
  sigma.logr = NULL,
  max.treedepth = NULL,
  adapt.delta = NULL,
  seed = NULL,
  approx = c("gaussian_laplace", "simplified_laplace"),
  correction = "none",
  n.gibbs = 10L,
  n.imputations = 20L,
  re.aghq = TRUE,
  n.quad = NULL,
  re.lkj = NULL,
  K.max = NULL,
  mixture = "poisson",
  integration = c("grid", "ccd"),
  sigma.grid = NULL,
  rho.grid = NULL,
  tau.grid = NULL,
  range.grid = NULL,
  verbose = TRUE,
  ...
)
```
