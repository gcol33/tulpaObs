# Simulate an open-population distance-sampling data set

Draws from the
[`distsamp_open()`](https://gillescolling.com/tulpaObs/reference/distsamp_open.md)
model: an open metapopulation (`N_1 ~ Poisson / NB(lambda)`,
`N_t = Binomial(N_{t-1}, omega) + Poisson(gamma)`) observed by distance
sampling (half-normal key, scale `sigma`) at each primary period.

## Usage

``` r
simulate_distsamp_open(
  N = 200,
  cutpoints = c(0, 10, 20, 30, 40),
  n_seasons = 4L,
  transect = "line",
  n_abund_covs = 1,
  n_det_covs = 1,
  beta_lambda = NULL,
  beta_sigma = NULL,
  omega = 0.7,
  gamma = 2.5,
  mixture = c("poisson", "negbin", "zip", "zinb"),
  size = 3,
  zi = 0.3,
  dynamics = c("constant", "notrend", "trend", "autoreg", "ricker", "gompertz"),
  K = NULL,
  r = 0.3,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites (default 200).

- cutpoints:

  Distance-bin edges `0 = c_0 < ... < c_B`.

- n_seasons:

  Number of primary periods (default 4).

- transect:

  `"line"` (default) or `"point"`.

- n_abund_covs, n_det_covs:

  Number of abundance / distance covariates.

- beta_lambda, beta_sigma:

  Coefficients on the log-abundance and log-scale arms. Defaults give
  moderate abundance / detection.

- omega, gamma:

  Apparent survival probability and recruitment rate (intercept-only
  defaults 0.7 / 2.5).

- mixture:

  Initial-abundance distribution: `"poisson"` (default), `"negbin"` /
  `"zinb"` (draws `N_1` negative-binomial with size `size`), or `"zip"`
  / `"zinb"` (a `zi` fraction of sites are structural zeros).

- size:

  Negative-binomial size (overdispersion `r`) when `mixture` draws an NB
  initial abundance (default 3).

- zi:

  Structural-zero probability when `mixture` is `"zip"` / `"zinb"`
  (default 0.3).

- dynamics:

  Population-dynamics form (see
  [`distsamp_open()`](https://gillescolling.com/tulpaObs/reference/distsamp_open.md)):
  `"constant"` (default), `"notrend"`, `"trend"`, `"autoreg"`,
  `"ricker"`, or `"gompertz"`. The alternative dynamics draw the
  abundance sequence from the matching transition (Poisson initial
  abundance only).

- K:

  Carrying capacity for `"ricker"` / `"gompertz"` (default
  `4 * lambda0`); `r` the intrinsic growth rate for those forms (default
  0.3). For `"trend"` / `"autoreg"` the recruitment multiplier is
  `gamma`.

- r:

  Intrinsic growth rate for `"ricker"` / `"gompertz"` (default 0.3).

- seed:

  Optional random seed.

## Value

A list with `y` (`[n_sites x n_bins x n_seasons]` distance-band counts),
`data`, and `truth`.

## Examples

``` r
sim <- simulate_distsamp_open(N = 40, n_seasons = 3, seed = 1)
dim(sim$y)
```
