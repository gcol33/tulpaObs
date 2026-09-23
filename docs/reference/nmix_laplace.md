# Laplace fit of the Royle (2004) N-mixture model

Maximum-likelihood fit (non-spatial, fixed effects only) of the Royle
(2004) N-mixture model with a Poisson or negative-binomial abundance
mixing distribution: \$\$N_i \sim \mathrm{Poisson}(\lambda_i)
\quad\text{or}\quad N_i \sim \mathrm{NegBin}(\mathrm{mean}=\lambda_i,
\mathrm{size}=r), \qquad y\_{ij} \| N_i \sim \mathrm{Binomial}(N_i,
p\_{ij}),\$\$ with abundance linear predictor \\\log \lambda_i =
X\_\lambda^{(i)} \beta\_\lambda\\ and detection linear predictor
\\\mathrm{logit}\\ p\_{ij} = X_p^{(ij)} \beta_p\\. The NB uses the
`neg_binomial_2` convention (size \\r\\, variance \\\lambda +
\lambda^2/r\\); Poisson is the \\r \to \infty\\ limit.

Optimisation uses inner Newton on \\(\beta\_\lambda, \beta_p)\\ with the
marginal observed Fisher information matrix as the curvature
(complete-data Fisher fallback at iterates where the observed-info
Hessian is not PSD). Under `mixture = "NB"` the dispersion \\\theta =
\log r\\ is a single global scalar, profiled outside the beta-Newton by
block coordinate ascent on its analytic profile score. All gradients and
Hessians are analytical (no numerical differentiation), so the fit
converges in a small number of iterations and is typically 20-50x faster
than
[`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html)'s
BFGS-with-numerical-derivatives path on equivalent problems.

This is the non-spatial entry point; a spatial nested-Laplace path that
places a latent prior block on the abundance arm is a separate entry, to
be added (see `tulpa_nested_laplace_nmix()` once shipped).

## Usage

``` r
nmix_laplace(
  y,
  site_idx,
  X_lambda,
  X_p,
  mixture = c("P", "NB"),
  beta_lambda_init = NULL,
  beta_p_init = NULL,
  log_r_init = NULL,
  r_max = 1e+05,
  K_max = NULL,
  headroom = NULL,
  max_iter = 100L,
  tol = 1e-06,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer vector of observed counts, one entry per visit (long form).

- site_idx:

  Integer vector, same length as `y`, 1-based site index indicating
  which site each visit belongs to.

- X_lambda:

  Numeric matrix `[n_sites x p_lambda]` of abundance covariates.

- X_p:

  Numeric matrix `[n_obs x p_p]` of detection covariates (long form, row
  order matches `y`).

- mixture:

  Abundance mixing distribution: `"P"` (Poisson, default) or `"NB"`
  (negative binomial). Under `"NB"` an extra dispersion parameter
  `log_r` (log NB size) is estimated jointly and reported with its
  standard error.

- beta_lambda_init:

  Optional numeric `[p_lambda]` warm start. Defaults to
  `c(log(mean(y) + 0.1), 0, 0, ...)`.

- beta_p_init:

  Optional numeric `[p_p]` warm start. Defaults to zeros.

- log_r_init:

  Optional scalar warm start for `log_r` (NB only). Defaults to a
  method-of-moments estimate from the site count totals, clamped to a
  sensible range.

- r_max:

  Upper bound on the NB size `r` (NB only, default `1e5`). If the
  optimiser pins `log_r` at `log(r_max)` the data are consistent with
  Poisson and a warning recommends `mixture = "P"`.

- K_max:

  Marginal-sum truncation. Defaults to `max(y) + 100` (matches
  [`unmarked::pcount`](https://ecoverseR.github.io/unmarked/reference/pcount.html)),
  applied per site (see `headroom`). The returned `boundary_weight` per
  site flags any site whose posterior over N puts non-trivial mass on
  its truncation; raise `K_max` if any such weight exceeds ~1e-4.

- headroom:

  Latent-N states summed above each site's own `max(y_i)`. `NULL`
  (default) derives it from `K_max`: an unset `K_max` caps each site at
  `max(y_i) + 100`, an explicit one truncates globally and uncapped. A
  negative value disables the per-site cap.

- max_iter:

  Newton iteration budget (default 100).

- tol:

  Gradient-norm convergence tolerance (default 1e-6).

- verbose:

  Print per-iteration `(log_lik, grad_norm, boundary_max)`.

## Value

A list of class `nmix_fit`:

- `beta_lambda`, `beta_p` – MLE coefficient vectors

- `mixture` – `"P"` or `"NB"`

- `log_r`, `r` – estimated log NB size and the NB size `exp(log_r)` (NB
  only; `NA` under Poisson). `vcov` carries `log_r` as its last
  coordinate, so `sqrt(diag(vcov))["log_r"]` is its standard error.

- `log_lik` – marginal log-likelihood at the mode

- `vcov` – variance-covariance matrix (marginal observed Fisher
  inverse), over `(beta_lambda, beta_p)` and, under NB, `log_r`

- `vcov_ok` – whether the final observed-info Cholesky succeeded

- `H_obs` – final marginal observed Fisher information matrix

- `n_iter` – (outer) iterations consumed

- `converged` – whether the convergence tolerance was reached

- `grad_norm` – gradient norm at termination

- `mean_N`, `var_N` – per-site posterior mean and variance of N \| y

- `boundary_weight` – per-site posterior weight on `N = K_max`

- `call` – matched call

## References

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics* 60, 108-115. Dennis, E.
B., Morgan, B. J. T., Ridout, M. S. (2015). Computational aspects of
N-mixture models. *Biometrics* 71, 237-246.
