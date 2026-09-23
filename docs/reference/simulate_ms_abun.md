# Simulate community (multispecies) N-mixture abundance data

Per-species Royle (2004) N-mixture with Gaussian community hyperpriors:
`beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`,
`beta_p_s ~ N(mu_p, Sigma_p)`, then `N_{s,i} ~ Poisson(lambda_{s,i})`
(or `NegBin(mu = lambda, size = r_s)` with a per-species size
`r_s = exp(mu_log_r + b_logr_s)`) and
`y_{s,i,j} ~ Binomial(N_{s,i}, p_{s,i,j})`. The returned `y` is a 3D
array `[n_sites x J x n_species]` suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md).

## Usage

``` r
simulate_ms_abun(
  n_species = 12,
  N = 80,
  J = 4,
  n_abund_covs = 1,
  n_det_covs = 1,
  mu_lambda = NULL,
  mu_p = NULL,
  sd_lambda = 0.5,
  sd_p = 0.4,
  mixture = c("poisson", "negbin", "zip", "zinb"),
  size = 3,
  sigma_logr = 0.4,
  omega = 0.3,
  sigma_omega = 0.4,
  graph = NULL,
  sigma.field = 0.6,
  seed = NULL
)
```

## Arguments

- n_species:

  Number of species (default 12).

- N:

  Number of sites (default 80).

- J:

  Number of replicate visits (default 4).

- n_abund_covs:

  Number of abundance covariates (default 1).

- n_det_covs:

  Number of detection covariates (default 1).

- mu_lambda:

  Community-mean abundance coefficients `c(intercept, slopes...)` on the
  log scale. Default `c(log(3), rep(0.4, n_abund_covs))`.

- mu_p:

  Community-mean detection coefficients on the logit scale. Default
  `c(0.3, rep(-0.3, n_det_covs))`.

- sd_lambda:

  Per-coefficient community SD for the abundance arm (the sqrt-diagonal
  of `Sigma_lambda`). Length 1 (recycled) or `1 + n_abund_covs`. Default
  0.5.

- sd_p:

  Per-coefficient community SD for the detection arm. Default 0.4.

- mixture:

  Abundance mixing distribution: `"poisson"` (default), `"negbin"`, or
  their zero-inflated counterparts `"zip"` / `"zinb"` (a per-species
  structural-zero share).

- size:

  Community-mean negative-binomial size, equal to \\\exp(\mu\_{\log
  r})\\ (ignored under Poisson). Default 3. The per-species sizes are
  \\r_s = \exp(\mu\_{\log r} + b^{\log r}\_s)\\ with \\b^{\log r}\_s
  \sim N(0, \sigma\_{\log r}^2)\\.

- sigma_logr:

  Standard deviation of the per-species log-dispersion random effect
  `log_r_s` (used only under `mixture = "negbin"` / `"zinb"`). Default
  0.4. Set to 0 for a shared (single-`r`) community.

- omega:

  Community-mean structural-zero probability, equal to
  \\\mathrm{plogis}(\mu\_\omega)\\ (used only under `mixture = "zip"` /
  `"zinb"`). Default 0.3. The per-species structural-zero probabilities
  are \\\omega_s = \mathrm{plogis}(\mu\_\omega + b^\omega_s)\\ with
  \\b^\omega_s \sim N(0, \sigma\_\omega^2)\\.

- sigma_omega:

  Standard deviation of the per-species structural-zero-logit random
  effect `logit_omega_s` (used only under zero-inflation). Default 0.4.
  Set to 0 for a shared (single-`omega`) community.

- graph:

  Optional `N x N` 0/1 adjacency matrix. When supplied, a single shared
  spatial field `f` (one value per site) is drawn from a proper GMRF on
  the graph, centered and scaled to standard deviation `sigma.field`,
  and added to every species' abundance log-linear predictor
  (`log lambda_{s,i} = X_lambda_i . beta_lambda_s + f_i`). `N` is taken
  from `nrow(graph)`. The truth field is returned in `truth$field`.

- sigma.field:

  Standard deviation of the shared spatial field (used only when `graph`
  is supplied). Default 0.6.

- seed:

  Optional random seed.

## Value

A list with `y` (3D count array), `data` (site covariate frame),
`species` (species names), and `truth` (community means / SDs, the
per-species coefficients, lambda, p, latent N, the shared `field` when
`graph` is given, and – under `negbin` – `mu_log_r`, `sigma_log_r`, the
per-species `b_logr` / `r_s`).

`truth` carries each community mean twice. `mu_lambda`, `mu_p`,
`mu_log_r` and `mu_omega` are the POPULATION constants the per-species
values were drawn around. `mu_lambda_real`, `mu_p_real`, `mu_log_r_real`
and `mu_omega_real` are the mean of the values this seed actually drew,
which sits `sd / sqrt(n_species)` from the constant – 0.12 on `mu_log_r`
at `sigma_logr = 0.5` and 18 species.

Which one to score against follows from the question. An interval for a
community mean targets the POPULATION constant and its width includes
that `sd^2 / n_species` term, so a coverage measurement scores against
`mu_log_r`. A point-recovery or interval-SCALE measurement scores
against `mu_log_r_real`: against the constant, the draw carries about
two thirds of the spread across seeds, and a seed block that happens to
draw wide reads as a miscalibrated interval. That is what the 18-species
fixture measured – its seed block drew the species means 18-21% wider
than `sigma_logr / sqrt(18)`, and both arms calibrate once the draw is
put at its expectation. The `_real` entries are supplied so the split
does not have to be rebuilt at each call site;
`tests/testthat/helper-community-mean.R` is the assertion that consumes
them.

## Examples

``` r
sim <- simulate_ms_abun(n_species = 4, N = 30, J = 3, seed = 1)
dim(sim$y)
```
