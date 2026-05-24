# Coverage test for beta_pos_0 in the BYM2 joint nested-Laplace engine.
#
# Closes the prior-mismatch loop diagnosed in INLAabun's d3 validation
# sweep: without demeaning the BYM2 sub-blocks `phi_f` and `theta_f`
# in the simulator, mean(w_s) per seed has SD ~ sigma / sqrt(n_s),
# the constrained cover-arm intercept identified by the engine targets
# `beta_pos_0_truth + alpha * mean(w_s_sim)` rather than the population
# truth, and 95% CI coverage of the population truth collapses to
# 0.43-0.47 at alpha_true >= 1.0. With demeaning (which is what
# `simulate_cover_joint()` does by construction), coverage returns to
# nominal.
#
# This test fixes:
#   (1) regressions in `.joint_inner_var`'s constraint correction; and
#   (2) regressions in `simulate_cover_joint()` that drop the demean.
#
# See `R/sim_cover_hurdle.R` (`simulate_cover_joint`) and
# `R/family_cover_hurdle.R` (`.joint_inner_var`) for the docstrings.

test_that("joint nested_laplace beta_pos_0 covers nominally at alpha=1 (BYM2)", {
  skip_on_cran()
  n_seeds      <- 20L
  n_s          <- 25L
  N            <- 200L
  alpha_true   <- 1.0
  beta_pos_0_truth <- 0.4

  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }

  in_ci <- logical(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_cover_joint(
      N        = N,
      adj      = adj,
      alpha    = alpha_true,
      positive = "beta",
      beta_pos = c(beta_pos_0_truth, -0.5),
      seed     = 7000L + r
    )
    fit <- tobs(
      formula = ~ x + bym2(graph = adj, group_var = "region"),
      data    = sim$data,
      family  = cover("beta"),
      y       = sim$y,
      method  = "nested_laplace",
      control = list(
        sigma.grid     = c(0.3, 0.6, 0.9),
        rho.grid       = c(0.5, 0.7, 0.9),
        sigma.pos.grid = c(0.0, 0.3, 0.6, 0.9, 1.2, 1.5),
        adaptive.grid  = FALSE
      )
    )
    expect_s3_class(fit, "cover_fit")
    expect_true(is.finite(fit$beta_pos[1]) && is.finite(fit$se_pos[1]))
    lo <- fit$beta_pos[1] - 1.96 * fit$se_pos[1]
    hi <- fit$beta_pos[1] + 1.96 * fit$se_pos[1]
    in_ci[r] <- (beta_pos_0_truth >= lo && beta_pos_0_truth <= hi)
  }

  # Nominal 95%; allow margin for 20-seed MC noise. The pre-fix
  # regime (un-demeaned simulator at alpha=1) sat at 0.43-0.47, so
  # 0.80 is well clear of that floor while leaving headroom for
  # honest seed-to-seed noise.
  cov_hat <- mean(in_ci)
  expect_gte(cov_hat, 0.80)
})


test_that("simulate_cover_joint() demeans both BYM2 sub-blocks", {
  # Direct unit test: the demeaning is the load-bearing detail of the
  # helper, so check it explicitly rather than infer from coverage.
  adj <- matrix(0L, 10L, 10L)
  for (s in seq_len(10L)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, 11L))) adj[s, j] <- 1L
  }
  sim <- simulate_cover_joint(N = 80L, adj = adj, alpha = 1.0, seed = 99)
  expect_lt(abs(mean(sim$truth$phi_f)),   1e-12)
  expect_lt(abs(mean(sim$truth$theta_f)), 1e-12)
  # Composite field is a linear combination of two zero-mean inputs,
  # so it's also zero-mean by construction.
  expect_lt(abs(mean(sim$truth$w_s)),     1e-12)
})
