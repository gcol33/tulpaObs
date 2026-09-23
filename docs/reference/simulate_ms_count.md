# Simulate community relative-abundance (count / continuous) data

Per-species GLMM with Gaussian community hyperpriors on the coefficients
(the
[`ms_count()`](https://gillescolling.com/tulpaObs/reference/ms_count.md)
/ spAbundance `msAbund` model): `g(mu_{s,i}) = X_i beta_s`,
`beta_s ~ N(beta_comm_mean, diag(beta_comm_sd^2))`, `y_{s,i}` drawn
Poisson / negative-binomial (log link) or Gaussian (identity). No
detection.

## Usage

``` r
simulate_ms_count(
  N = 120,
  n_species = 10,
  beta_comm_mean = c(1, 0.5),
  beta_comm_sd = c(0.4, 0.3),
  response = c("poisson", "negbin", "gaussian", "binomial"),
  size = 2,
  size.log.sd = 0.3,
  sd = 1,
  trials = 10,
  seed = NULL
)
```

## Arguments

- N:

  Number of sites.

- n_species:

  Number of species.

- beta_comm_mean, beta_comm_sd:

  Community mean and SD of the per-species coefficients (intercept
  first). Length sets the number of covariates.

- response:

  One of `"poisson"`, `"negbin"`, `"gaussian"`, `"binomial"`.

- size:

  Negative-binomial community size (mean of the per-species `log_r`);
  `size.log.sd` is its across-species SD.

- size.log.sd:

  Across-species SD of `log(size)` (negbin only).

- sd:

  Gaussian residual SD.

- trials:

  Binomial trial count (`response = "binomial"`): a scalar (shared
  across sites and species) or a length-`N` per-site vector. Default 10.

- seed:

  Optional RNG seed.

## Value

A list with `y` (an `N x n_species` matrix), `data`, and `truth`.

## Examples

``` r
sim <- simulate_ms_count(N = 50, n_species = 4, seed = 1)
head(sim$y)
```
