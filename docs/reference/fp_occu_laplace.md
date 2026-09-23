# Maximum-likelihood fit of the multistate false-positive occupancy model

Fits the Miller et al. (2011) false-positive occupancy model with
confirmed detections (`y` in `{0, 1, 2}`) by maximising the exact
two-state marginal likelihood with an analytic gradient (BFGS). The
observed-information covariance is the inverse of the negative
finite-difference Jacobian of the analytic gradient at the mode. All
four arms (occupancy `psi`, true detection `p11`, false-positive `p10`,
certain-classification `b`) are site-level logit predictors.

## Usage

``` r
fp_occu_laplace(
  y,
  site_idx,
  X_psi,
  X_p11,
  X_p10,
  X_b,
  sigma_beta = NULL,
  max_iter = 200L,
  tol = 1e-08,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer vector of detection states (`0`/`1`/`2`), long form (valid
  visits only).

- site_idx:

  Integer vector, same length as `y`, 1-based site index.

- X_psi, X_p11, X_p10, X_b:

  Numeric `[n_sites x p_arm]` design matrices for the four arms.

- sigma_beta:

  Optional Gaussian prior SD on the coefficients (a mild ridge for
  stability); `NULL` (default) is the unpenalised MLE.

- max_iter:

  BFGS iteration budget (default 200).

- tol:

  Convergence tolerance (`optim` `reltol`, default 1e-8).

- verbose:

  Print convergence status.

## Value

A list of class `fp_occu_fit` with `beta_psi`, `beta_p11`, `beta_p10`,
`beta_b`, `log_lik`, `vcov`, `H_obs`, per-site posterior occupancy `w1`,
and `converged`.

## References

Miller, D. A. W., et al. (2011). Improving occupancy estimation when two
types of observational error occur. *Ecology* 92, 1422-1428.
