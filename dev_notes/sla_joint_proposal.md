# SLA on the joint nested-Laplace path — design proposal

**Goal.** Wire simplified-Laplace (RMC 2009 §3.2) marginals through the
`tulpa_nested_laplace_joint` path so `tobs(family = cover(...), engine =
"nested_laplace", approx = "simplified_laplace")` returns SN-corrected
marginals on `beta_occ` and `beta_pos` instead of the current
Gaussian-at-the-mode fallback (`message("simplified_laplace is currently
wired only for the single-Laplace cover path; falling back to
gaussian_laplace marginals on the joint nested-Laplace fit.")` — see
`R/family_cover_hurdle.R:48-52`).

## What the joint engine already exposes (no tulpa-side change needed)

`tulpa_nested_laplace_joint(..., store_Q = TRUE)` returns:

- `modes`             — `[n_grid x n_x]` per-grid joint mode of the inner latent.
- `weights`           — outer-grid posterior weights (sum to 1).
- `theta_grid`        — outer hyperparameter values (`sigma`, `rho`, `sigma_pos`, `phi_pos`).
- `arm_layout`        — `beta_start[k]` / `re_start[k]` / `phi_start` / `theta_start` / `n_x`.
- `Q_csc_*_per_grid`  — sparse lower-tri CSC of the per-grid joint precision.

`tulpaObs::.joint_inner_var(fit, beta_idx)` already inverts each `Q_k` and
applies the sum-to-zero constraint correction on `phi` / `theta` blocks,
returning per-grid `Var(beta_j | data, theta_k)`. That gives the
**Mean-of-Var** half of the marginal variance; **Var-of-means** comes from
the per-grid mode shifts. Both are wired in `family_cover_hurdle.R:717-739`.

What's still missing for SLA: the **inner third cumulant** per grid point
and the **third-moment mixture formula** for combining them across the
outer grid.

## Math: from per-grid skewness to marginal skewness

At each grid point `k`, the inner posterior over `x = (beta, w_field)` is
approximated by the multivariate Gaussian `N(mode_k, Q_k^{-1})`. SLA
replaces this with a skew-normal moment-matched to `(mode, Hessian-inverse,
3rd-cumulant-from-l3-at-mode)`. The per-grid inner third cumulant of
`beta_j` along its Sigma column is

    kappa_3[beta_j | theta_k]
      = d^3 / dh^3  L(x_hat_k + h * Sigma_k[, beta_j_col])  evaluated at h = 0
      = Sum_{a,b,c}  l'''_{abc}  *  v_j[a]  *  v_j[b]  *  v_j[c],
      where v_j = Q_k^{-1} e_{beta_j_col}  (with sum-to-zero correction).

`L` is the *joint inner log-likelihood at theta_k*: sum of arm
log-likelihoods + log-prior on the field at the fixed `theta_k`. The prior
is Gaussian in `x`, so its third derivative is zero — only the arm
likelihoods contribute to `l'''`. The constraint correction on the spatial
block applies to `v_j` (project out the constrained directions) so the
displacement does not move the field along its null space.

Marginal skewness combines per-grid mean, variance, and skewness:

    M3_j = Sum_k w_k * [ (mu_kj - mu_j)^3
                       + 3 * (mu_kj - mu_j) * sigma_kj^2
                       + gamma_kj * sigma_kj^3 ]
    gamma_j = M3_j / sigma_j^3

where `mu_j = Sum_k w_k * mu_kj`, `sigma_j^2 = var_of_means + mean_of_var`
(both already computed), and `gamma_kj` is the standardised per-grid
skewness `kappa_3[beta_j | theta_k] / sigma_kj^3`.

The marginal Gaussian-vs-skew check is then exactly the existing
`.sla_replace_draws(draws, means, sds, gamma)` path: resample each beta
column from a moment-matched skew-normal, with `cap = 0.5` and
`SN_GAMMA_MAX` fallback. No changes to that helper.

## Implementation

### New file: `R/sla_cover_hurdle_joint.R` (~220 LOC)

```r
# Joint inner log-likelihood at fixed grid point k.
# Inputs: x (length n_x latent), k (grid index), fit, enc, positive.
# Output: scalar log-lik (arm Bernoulli + arm Beta/Gaussian; prior term
#         dropped because the prior is Gaussian and its 3rd deriv is 0,
#         so it cancels under FD along any direction).
.loglik_cover_joint_at_grid <- function(x, k, fit, enc, positive) {
  L <- fit$arm_layout
  beta_occ <- x[L$beta_start[1] + seq_len(L$p[1])]
  beta_pos <- x[L$beta_start[2] + seq_len(L$p[2])]
  field    <- .joint_field_at_all_units(x, k, fit)        # length n_s

  # Donor (occ) arm sees field directly with amplitude sigma_occ_k.
  # Copy   (pos) arm sees field scaled by sigma_pos_k / sigma_occ_k.
  s_occ_k <- fit$theta_grid[k, "sigma"]
  s_pos_k <- fit$theta_grid[k, "sigma_pos"]
  spi_occ <- enc$..spi_full
  spi_pos <- enc$..spi_pos

  eta_occ <- as.numeric(enc$occ_data$X %*% beta_occ) +
             s_occ_k * field[spi_occ]
  eta_pos <- as.numeric(enc$pos_data$X %*% beta_pos) +
             s_pos_k * field[spi_pos]

  ll_occ <- sum(enc$occ_data$y * eta_occ - log1p(exp(eta_occ)))
  ll_pos <- if (positive == "beta") {
    phi_k <- fit$theta_grid[k, "phi_pos"]
    .ll_beta_arm(eta_pos, enc$pos_data$y, phi_k)
  } else {
    sig_k <- fit$theta_grid[k, "sigma_pos_noise"] %||% 1   # lognormal noise
    sum(dnorm(enc$pos_data$y, eta_pos, sig_k, log = TRUE))
  }
  ll_occ + ll_pos
}
```

(Note: the field at the **all-units** vector is needed because each arm
gathers via its own `spatial_idx`. The mode stores phi+theta per unit, so
`field[s] = sqrt(rho) * phi_k[s] + sqrt(1-rho) * theta_k[s]` for BYM2 or
`field[s] = phi_k[s]` for ICAR / CAR_proper, before per-arm sigma scaling.
Reuse the existing `.joint_field_at_obs_copy` helper, generalised to
return per-unit values rather than per-obs.)

### Per-grid SLA gamma via FD

```r
.sla_inner_gamma_joint <- function(k, beta_idx_arm, fit, enc, positive,
                                   beta_starts) {
  n_x  <- fit$arm_layout$n_x
  Qk   <- .restore_Qk_from_csc(fit, k)                   # n_x x n_x sparse
  # Constraint-corrected solve: V = Q^{-1} - Q^{-1} A^T (A Q^{-1} A^T)^{-1} A Q^{-1}
  # We need V[, beta_idx_arm] (n_x x p_arm).
  V <- .joint_constrained_solve_columns(Qk, fit$arm_layout, beta_idx_arm)
  beta_hat <- fit$modes[k, ]

  gamma <- numeric(length(beta_idx_arm))
  for (j in seq_along(beta_idx_arm)) {
    v_j     <- V[, j]
    sigma_j <- sqrt(max(v_j[beta_idx_arm[j]], 0))
    if (sigma_j <= 0) next
    eps_h <- .Machine$double.eps^(1/5)
    h     <- eps_h * sigma_j / max(sqrt(sum(v_j^2)), .Machine$double.eps)
    L_p2 <- .loglik_cover_joint_at_grid(beta_hat + 2*h*v_j, k, fit, enc, positive)
    L_p1 <- .loglik_cover_joint_at_grid(beta_hat +   h*v_j, k, fit, enc, positive)
    L_m1 <- .loglik_cover_joint_at_grid(beta_hat -   h*v_j, k, fit, enc, positive)
    L_m2 <- .loglik_cover_joint_at_grid(beta_hat - 2*h*v_j, k, fit, enc, positive)
    d3   <- (L_p2 - 2*L_p1 + 2*L_m1 - L_m2) / (2 * h^3)
    gamma[j] <- d3 / sigma_j^3
  }
  gamma
}
```

### Marginal third-moment mixture

```r
.sla_combine_grid_skewness <- function(modes_j, weights, sigma2_kj,
                                       gamma_kj, mu_j, sigma2_j) {
  # mu_j     : marginal mean (already computed)
  # sigma2_j : marginal variance (already computed)
  # All inputs length n_grid (modes_j) or scalar (mu_j / sigma2_j).
  dmu <- modes_j - mu_j
  M3  <- sum(weights * (dmu^3 + 3 * dmu * sigma2_kj +
                        gamma_kj * sqrt(sigma2_kj)^3))
  sigma_j <- sqrt(pmax(sigma2_j, 0))
  if (sigma_j <= 0) return(0)
  M3 / sigma_j^3
}
```

### Orchestrator and decode wiring

```r
.sla_compute_cover_hurdle_joint <- function(fit, enc, positive) {
  L         <- fit$arm_layout
  bocc_idx  <- L$beta_start[1] + seq_len(L$p[1])
  bpos_idx  <- L$beta_start[2] + seq_len(L$p[2])
  weights   <- fit$weights
  n_grid    <- length(weights)

  inner_var <- .joint_inner_var(fit, c(bocc_idx, bpos_idx))         # already exists
  gamma_grid <- matrix(NA_real_, n_grid,
                       length(bocc_idx) + length(bpos_idx))
  for (k in seq_len(n_grid)) {
    gamma_grid[k, seq_along(bocc_idx)]  <- .sla_inner_gamma_joint(
      k, bocc_idx, fit, enc, positive)
    gamma_grid[k, length(bocc_idx) + seq_along(bpos_idx)] <- .sla_inner_gamma_joint(
      k, bpos_idx, fit, enc, positive)
  }

  # Marginal means / variances reused from the existing decode.
  modes_occ <- fit$modes[, bocc_idx, drop = FALSE]
  modes_pos <- fit$modes[, bpos_idx, drop = FALSE]
  mu_occ    <- as.numeric(crossprod(weights, modes_occ))
  mu_pos    <- as.numeric(crossprod(weights, modes_pos))
  vom_occ   <- as.numeric(crossprod(weights, modes_occ^2)) - mu_occ^2
  vom_pos   <- as.numeric(crossprod(weights, modes_pos^2)) - mu_pos^2
  mov_occ   <- as.numeric(crossprod(weights, inner_var[, seq_along(bocc_idx)]))
  mov_pos   <- as.numeric(crossprod(weights,
                                    inner_var[, length(bocc_idx) + seq_along(bpos_idx)]))

  sigma2_occ <- vom_occ + mov_occ
  sigma2_pos <- vom_pos + mov_pos

  gamma_occ <- vapply(seq_along(bocc_idx), function(j) {
    .sla_combine_grid_skewness(
      modes_occ[, j], weights,
      sigma2_kj = inner_var[, j],
      gamma_kj  = gamma_grid[, j],
      mu_j      = mu_occ[j],
      sigma2_j  = sigma2_occ[j])
  }, numeric(1))
  gamma_pos <- vapply(seq_along(bpos_idx), function(j) {
    .sla_combine_grid_skewness(
      modes_pos[, j], weights,
      sigma2_kj = inner_var[, length(bocc_idx) + j],
      gamma_kj  = gamma_grid[, length(bocc_idx) + j],
      mu_j      = mu_pos[j],
      sigma2_j  = sigma2_pos[j])
  }, numeric(1))

  list(gamma_occ = gamma_occ, gamma_pos = gamma_pos, valid = TRUE,
       reason = "ok")
}
```

### Decode integration

In `decode_cover_hurdle_joint()` (currently lines 870-914), add an
`approx` argument and pass through from `.dispatch_cover()`. When
`approx = "simplified_laplace"`:

```r
sla_res <- .sla_compute_cover_hurdle_joint(fits$joint, enc, fits$positive)
sla_draws <- .sla_build_cover_hurdle_draws(
  beta_occ, se_occ, beta_pos, se_pos, sla_res)
draws_occ <- sla_draws$draws_occ
draws_pos <- sla_draws$draws_pos
sla_status <- sla_draws$sla_status
skew_occ <- if (isTRUE(sla_res$valid)) sla_res$gamma_occ else NULL
skew_pos <- if (isTRUE(sla_res$valid)) sla_res$gamma_pos else NULL
```

`.sla_build_cover_hurdle_draws()` is already grid-agnostic — reuses the
moment-matched SN per coefficient. No changes there.

### Dispatcher

Remove the early `message("simplified_laplace is currently wired only ...")`
and forward `approx` into `decode_cover_hurdle_joint()`.

## LOC budget

| Component                                | LOC |
|------------------------------------------|----:|
| Joint inner log-lik (occ + Beta + LN)    |  50 |
| `.joint_constrained_solve_columns`       |  35 |
| `.sla_inner_gamma_joint` (FD per j)      |  45 |
| `.sla_combine_grid_skewness`             |  20 |
| `.sla_compute_cover_hurdle_joint`        |  55 |
| Decode wiring + dispatcher edit          |  20 |
| **Total**                                | ~225 |

## Verification plan

1. **Gaussian limit**: when all `gamma_kj ≈ 0` (large-data regime), marginal
   `gamma_j` should reduce to the mixture-skewness from per-grid mode shifts
   alone — i.e., the formula degenerates to the Gaussian-mixture skewness
   already implied by `var_of_means`. Sanity test: simulate a high-N cover
   hurdle, verify `gamma_j → 0`.

2. **Cross-check vs separate-Laplace SLA**: in the special case where the
   spatial field has effectively zero amplitude (e.g. `sigma = 0.01`), the
   joint engine should reduce to the two-arm separate Laplace, and the
   SLA gamma should match `.sla_compute_cover_hurdle()` to 3 decimals.
   Add a test in `tests/testthat/test-sla-cover-joint.R`.

3. **Coverage improvement on small-n_pos**: rerun INLAabun D7 Cell B (the
   sparse-positive joint, ~35 positive obs) under `approx =
   "simplified_laplace"`. Expect intercept coverage to move from ~0.80
   (Gaussian) toward 0.93+ — the SLA gain at small n_pos is what motivated
   shipping SLA on the standalone path in the first place.

## Open questions

- The FD step size `h = eps^(1/5) * sigma_j / ||v_j||` is the standard
  choice for the 5-point central rule. With `n_x ≈ 50` and the field
  contributing most of the displacement, `||v_j||` could be large enough
  that `h` becomes tiny and cancellation dominates. Falls under "test it
  on the sparse-positive case first; if `gamma` becomes noisy, scale `h`
  per-direction by `1/sqrt(n_x)` or switch to a higher-order rule."
- For `engine = "nested_laplace"` with non-cover families (occu, abun
  when their joint paths land), the orchestrator generalises but the
  per-arm log-likelihood evaluator changes. Keep the abstraction local
  to cover-hurdle for now; lift it when the second consumer arrives.

## Out of scope

- Cross-arm correlations in the SN draws. Single-Laplace SLA already
  treats arms as independent; the joint path inherits this. Joint
  correlations within an arm are not preserved either — matches the
  single-arm SLA spec from `dev_notes/upstream_tulpa_sla_spec.md §2.1`.
- The phi marginal under `cover('beta', engine = 'nested_laplace')`.
  Same Gaussian-mixture-over-grid treatment as currently; SLA on a 1-D
  outer hyperparameter would be a separate small extension.
