# Simulate Dail-Madsen open N-mixture data

Latent `N_1 ~ Poisson(lambda)`; for `t >= 2`,
`N_t = Binomial(N_(t-1), omega_{t-1}) + Poisson(gamma_{t-1})`; observed
via `Binomial(N_t, p)` over `J` secondary visits in each of `T` primary
seasons. Returns a 3D array `N x J x T` suitable for
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
[`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md).

## Usage

``` r
simulate_dyn_abun(
  N = 150,
  T = 4,
  J = 3,
  n_abund_covs = 1,
  beta_lambda = NULL,
  p = 0.5,
  omega = 0.6,
  gamma = 1,
  beta_omega = NULL,
  beta_gamma = NULL,
  mixture = c("poisson", "negbin"),
  r = 2,
  zi = 0,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 150).

- T:

  Number of primary seasons (default 4).

- J:

  Number of secondary visits per season (default 3).

- n_abund_covs:

  Number of initial-abundance covariates (default 1).

- beta_lambda:

  Initial-abundance coefficients (log). Default
  `c(log(5), runif(n_abund_covs, -0.4, 0.4))`.

- p, omega, gamma:

  Detection, apparent-survival, and recruitment parameters (scalars;
  defaults 0.5, 0.6, 1.0). `omega` / `gamma` set the constant rate used
  when `beta_omega` / `beta_gamma` are `NULL`.

- beta_omega, beta_gamma:

  Optional length-2 coefficients `c(intercept, slope)` for a
  season-varying survival (logit link) or recruitment (log link) rate
  driven by a per-(site, interval) season covariate. `NULL` (default)
  keeps the constant `omega` / `gamma`.

- mixture:

  Initial-abundance distribution: `"poisson"` (default) or `"negbin"`
  (negative-binomial `N_1 ~ NB(mean = lambda, size = r)`).

- r:

  Negative-binomial size for `mixture = "negbin"` (default 2).

- zi:

  Structural-zero share for a zero-inflated model (default 0). With
  probability `zi` a site is never occupied (`N_t = 0` for all seasons),
  so all its counts are zero; fit such data with
  `dyn_abun(mixture = "zip")` / `"zinb"`.

- seed:

  Optional random seed.

## Value

A list with `y` (N x J x T count array), `data` (covariates, including
the `[N x (T-1)]` `season_cov` matrix column when a season-varying rate
is used), and `truth` (coefficients, per-site `lambda`, the realised
`omega` / `gamma` rates, and `r` under `"negbin"`).

## Details

Survival and recruitment are constant across seasons by default.
Supplying `beta_omega` (logit-scale) or `beta_gamma` (log-scale) instead
makes the rate depend on a per-(site, interval) season covariate `z`:
`omega_{i,t} = plogis(beta_omega[1] + beta_omega[2] * z_{i,t})`,
`gamma_{i,t} = exp(beta_gamma[1] + beta_gamma[2] * z_{i,t})`. The
covariate is returned as a `[N x (T-1)]` matrix column `season_cov` of
`data`, ready for `omega = ~ season_cov`.

## Examples

``` r
sim <- simulate_dyn_abun(N = 40, T = 3, J = 3, seed = 1)
dim(sim$y)
```
