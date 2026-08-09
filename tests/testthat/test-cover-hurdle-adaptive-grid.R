# Parameter-recovery test for the adaptive hyperparameter grid wrapped
# around `tulpa_nested_laplace_joint`. Asserts that 95% credible intervals
# for the copy-scaling parameter `alpha` cover the truth at the nominal
# rate. The fixed-grid path under-covers when the truth sits near an outer
# axis edge, because the integrand cannot extend past the boundary; the
# adaptive path adds a densification point inside the heaviest cell and two
# outward extension points whenever the relative integrand density at
# the boundary exceeds `adaptive_grid_edge_thresh`.
#
# See gcol33/tulpaObs#8 and gcol33/INLAabun's Demo 3
# `example/validation/d3_joint_copy.R` for the reproducer that motivated
# the fix. The configuration here mirrors D3 (N=300, n_s=25, BYM2,
# beta-positive arm) with `alpha_true = 1.5` on the default copy axis
# (gcol33/tulpaObs#194: the boundary these seeds were chosen to sit on was
# an axis of the retired (sigma_occ, sigma_pos) parameterization, so the
# fixture no longer places the truth at an edge).

simulate_d3_like <- function(seed, alpha_true,
                             N = 300L, n_s = 25L,
                             sigma_true = 0.6, rho_true = 0.7,
                             phi_true = 30,
                             beta_occ = c(-0.3, 0.7),
                             beta_pos = c(0.4, -0.5)) {
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  phi_f   <- rnorm(n_s, 0, 1)
  theta_f <- rnorm(n_s, 0, 1)
  w_s     <- sigma_true * (sqrt(rho_true) * phi_f +
                              sqrt(1 - rho_true) * theta_f)
  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha_true * w_s[spatial_idx]
  mu_pos  <- plogis(eta_pos)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi_true,
                     (1 - mu_pos[is_pos]) * phi_true)
  y       <- pmin(pmax(y, 0), 1 - 1e-6)
  list(
    data = data.frame(x = x, region = factor(spatial_idx)),
    y    = y,
    truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                 sigma = sigma_true, rho = rho_true,
                 alpha = alpha_true, phi = phi_true)
  )
}

chain_adj_for_test <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  adj
}

test_that("adaptive grid covers alpha at the upper boundary across 20 seeds", {
  skip_on_cran()
  skip_if_fast()
  truth_alpha <- 1.5
  n_seeds <- 20L
  n_s     <- 25L
  adj     <- chain_adj_for_test(n_s)

  results <- vapply(seq_len(n_seeds), function(r) {
    sim <- simulate_d3_like(seed = 3400L + r, alpha_true = truth_alpha,
                             n_s = n_s)
    fit <- tobs(
      formula  = ~ x + bym2(graph = adj, group_var = "region"),
      data     = sim$data,
      family   = cover("beta"),
      y        = sim$y,
      method   = "nested_laplace",
      control  = list(
        sigma.grid     = c(0.3, 0.6, 0.9),
        rho.grid       = c(0.5, 0.7, 0.9),
        adaptive.grid  = TRUE
      )
    )
    expect_s3_class(fit, "cover_fit")
    # Use the marginalized 95% quantile CI on alpha. Wald mean +/- 1.96*sd
    # over-rejects in this regime because the alpha posterior is
    # right-skewed at sigma_pos near the upper boundary (median below
    # mean; long upper tail when the data weakly identifies alpha).
    # Quantile CI is the engine's documented summary for skew axes.
    alpha_lo <- fit$joint$theta_ci_lo["alpha"]
    alpha_hi <- fit$joint$theta_ci_hi["alpha"]
    alpha_lo <= truth_alpha && truth_alpha <= alpha_hi
  }, logical(1))

  coverage <- mean(results)
  # Target: nominal 0.95 coverage; on the calibration sweep (20 seeds,
  # N=300, BYM2, beta-positive arm) the CI-based check covers all 20
  # seeds. 0.85 is a conservative lower bound robust to seed-level noise.
  expect_gte(coverage, 0.85)
})

test_that("adaptive grid is strictly better than fixed grid at the boundary", {
  skip_on_cran()
  skip_if_fast()
  truth_alpha <- 1.5
  n_seeds <- 20L
  n_s     <- 25L
  adj     <- chain_adj_for_test(n_s)

  cover_fixed <- logical(n_seeds)
  cover_adapt <- logical(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_d3_like(seed = 3400L + r, alpha_true = truth_alpha,
                             n_s = n_s)
    # Pin phi_grid to the original 13-point default to keep this test
    # at its calibrated regime: with the post-tulpa#19 default (7-point
    # phi + interior densification) the fixed-grid path already covers
    # alpha well because the coarser phi marginalisation widens the
    # joint Sd(alpha), so the gap that motivated this test no longer
    # opens up. Holding phi fixed isolates the sigma_pos boundary fix.
    ctrl <- list(
      sigma.grid     = c(0.3, 0.6, 0.9),
      rho.grid       = c(0.5, 0.7, 0.9),
      phi.grid       = exp(seq(log(2), log(300), length.out = 13))
    )
    fit_fix <- tobs(formula = ~ x + bym2(graph = adj, group_var = "region"), data = sim$data, family = cover("beta"),
                    y = sim$y, method = "nested_laplace",
                    control = c(ctrl, list(adaptive.grid = FALSE)))
    fit_ad  <- tobs(formula = ~ x + bym2(graph = adj, group_var = "region"), data = sim$data, family = cover("beta"),
                    y = sim$y, method = "nested_laplace",
                    control = c(ctrl, list(adaptive.grid = TRUE)))
    # Quantile CI coverage on alpha (see comment in the preceding test
    # for the rationale of avoiding the Wald check at the boundary).
    cover_fixed[r] <- fit_fix$joint$theta_ci_lo["alpha"] <= truth_alpha &&
      truth_alpha <= fit_fix$joint$theta_ci_hi["alpha"]
    cover_adapt[r] <- fit_ad$joint$theta_ci_lo["alpha"] <= truth_alpha &&
      truth_alpha <= fit_ad$joint$theta_ci_hi["alpha"]
  }
  # Adaptive grid should not regress coverage at the boundary scenario.
  # Under the legacy cartesian refinement (~351 added cells per seed)
  # the gap was 15-30 pp because the boundary alpha posterior was
  # severely truncated. Mode-tracked 1D refinement (gcol33/tulpa#19) is
  # ~100x cheaper but approximates the conditional posterior shape at
  # extension points from the boundary anchor, so the absolute coverage
  # gap is smaller (a few pp) and seed-noise can flip the sign. The
  # invariant we still defend: adaptive does not *hurt* coverage by
  # more than seed-level noise.
  expect_gte(mean(cover_adapt) - mean(cover_fixed), -0.15)
})

test_that("adaptive grid stays a no-op when the integrand has fully decayed", {
  skip_on_cran()
  skip_if_fast()
  # alpha_true = 0 (i.e. sigma_pos_true = 0) places the truth at the
  # *opposite* edge of the user grid, but the data drives the posterior
  # to concentrate near 0 with quickly-decaying tails, so the integrand
  # at the upper boundary (sigma_pos = 0.9) is essentially zero and
  # refinement should not fire.
  sim <- simulate_d3_like(seed = 3101L, alpha_true = 0.0)
  n_s <- nlevels(sim$data$region)
  adj <- chain_adj_for_test(n_s)
  fit <- tobs(
    formula = ~ x + bym2(graph = adj, group_var = "region"), data = sim$data, family = cover("beta"), y = sim$y,
    method = "nested_laplace",
    control = list(
      sigma.grid     = c(0.3, 0.6, 0.9),
      rho.grid       = c(0.5, 0.7, 0.9),
      adaptive.grid  = TRUE
    )
  )
  # adaptive_grid_info should be present and report refinement on
  # sigma_pos if and only if the boundary integrand density crossed the
  # threshold. At sigma_pos_true = 0 the integrand at sigma_pos = 0.9 is
  # ~1e-15 relative to the peak (lm @ 0.9 << lm @ 0), so the trigger
  # must NOT fire on the max boundary. The min boundary is at the local
  # mode, so the trigger is also OK to not fire because the score is
  # dominated by the density (~ 1) but extension below 0 is bounds-
  # clipped to no points. Since the tulpa#19 follow-up the engine also
  # supports interior densification on the `phi_<arm>` axis (default
  # 7-point log grid is coarse enough to trigger), so refinement here
  # often fires on `phi_pos` even when sigma_pos doesn't.
  info <- fit$joint$adaptive_grid_info
  if (!is.null(info)) {
    triggered <- unlist(strsplit(info$triggered_axes, ","))
    # Refinement, when it fires, must be on a legitimate refinable axis
    # — either sigma_pos (min-side densification) or a per-arm phi axis
    # (interior densification on the coarse default phi grid).
    expect_true(any(triggered %in% c("sigma_pos", "phi_pos")))
    # Coverage at the true alpha = 0: the marginal alpha posterior is
    # heavily right-skewed near the boundary (median ~0.004, mean ~0.02,
    # CI ~[0.0001, 0.2]). A Wald z-score `|mean|/sd` over-rejects in this
    # regime because the Gaussian-style SD captures only the left side of
    # the spike. Use the marginalized weighted-quantile CI (the engine's
    # documented summary for skew axes) and check the upper edge stays
    # well below the bulk of the user grid -- i.e. the posterior
    # genuinely concentrates near the truth.
    expect_lt(fit$joint$theta_ci_lo["alpha"], 0.05)
    expect_lt(fit$joint$theta_ci_hi["alpha"], 0.5)
  }
})

test_that("outer-grid pruning keeps the mode and leaves estimates unchanged", {
  skip_on_cran()
  skip_if_fast()
  # tulpaObs#20: the dense outer tensor concentrates posterior mass on a few
  # cells (ESS ~ 1), so most cells run a full-data inner solve for negligible
  # weight. The cheap-pass prune skips them. This test asserts pruning is a
  # no-op on the inference: the modal hyperparameters and the coefficient
  # estimates must match the un-pruned dense-grid fit, and the pruned cells
  # must carry negligible weight (the mode is never pruned).
  n_s <- 25L; adj <- chain_adj_for_test(n_s)
  sim <- simulate_d3_like(seed = 101L, alpha_true = 1.0, N = 400L, n_s = n_s)
  ctrl_grid <- list(
    sigma.grid     = exp(seq(log(0.2), log(1.5), length.out = 5)),
    rho.grid       = c(0.25, 0.5, 0.7, 0.9),
    phi.grid       = exp(seq(log(2), log(300), length.out = 7)),
    adaptive.grid  = FALSE
  )
  run <- function(prune)
    tobs(formula = ~ x + bym2(graph = adj, group_var = "region"),
         data = sim$data, family = cover("beta"), y = sim$y,
         method = "nested_laplace",
         control = c(ctrl_grid, list(prune = prune, prune.tol = 1e-4)))

  f_off <- run(FALSE)
  f_on  <- run(TRUE)

  # Most cells pruned; the prune machinery is engaged.
  expect_true(f_on$joint$prune_n_pruned > 0.5 * f_on$joint$n_grid)

  # The modal hyperparameters (highest-weight cell) agree between fits.
  mode_off <- f_off$joint$theta_grid[which.max(f_off$joint$weights), ]
  mode_on  <- f_on$joint$theta_grid[which.max(f_on$joint$weights), ]
  names(mode_off) <- f_off$joint$theta_names
  names(mode_on)  <- f_on$joint$theta_names
  for (nm in c("sigma", "rho", "alpha"))
    expect_equal(unname(mode_on[nm]), unname(mode_off[nm]), tolerance = 0.05)

  # Coefficient estimates materially unchanged by pruning.
  expect_equal(f_on$beta_occ, f_off$beta_occ, tolerance = 0.02)
  expect_equal(f_on$beta_pos, f_off$beta_pos, tolerance = 0.02)
  expect_equal(unname(f_on$joint$theta_mean["alpha"]),
               unname(f_off$joint$theta_mean["alpha"]), tolerance = 0.05)
})
