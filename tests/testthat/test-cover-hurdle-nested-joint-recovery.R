# Multi-seed parameter-recovery tests for cover() dispersion scalars
# (`sigma_pos` for the lognormal arm, `phi_pos` for the beta arm).
#
# These complement the shape / single-seed tests in
# `test-cover-hurdle-nested-joint.R` and `test-cover-hurdle-beta.R` by
# asserting `|estimate - truth| < tol` across 10 seeds on the paths
# where the engine recovers cleanly.
#
# Triggered by INLAabun's validation harness, which previously dropped
# `sigma_pos` / `phi_pos` from its summary tables (a reporting bug) and
# in passing surfaced a recovery problem in the joint-engine beta path
# (phi profile is downward-biased; see below).
#
# Path coverage:
#   * joint  nested_laplace + lognormal -> sigma_pos       (passes here)
#   * joint  nested_laplace + beta      -> phi_pos         (known biased,
#                                                          test skipped
#                                                          with reason)
#   * separate-hurdle       + beta      -> phi_pos         (passes here)
#
# The separate-hurdle lognormal multi-seed sigma_pos recovery already
# lives in `test-cover-hurdle-lognormal.R::"repeat fits recover truth
# in aggregate"` and is not duplicated here.

simulate_joint_lognormal_for_recovery <- function(N = 400, n_s = 30,
                                                  sigma = 0.6, rho = 0.7,
                                                  alpha = 1.0,
                                                  beta_occ = c(0.2, 0.7),
                                                  beta_pos = c(-1.5, 0.3),
                                                  sigma_pos_true = 0.4,
                                                  seed = 11) {
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  phi_f   <- rnorm(n_s, 0, 1)
  theta_f <- rnorm(n_s, 0, 1)
  w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)

  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))

  is_pos  <- occur == 1L
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
  log_y   <- rnorm(N, eta_pos, sigma_pos_true)
  y       <- ifelse(is_pos, exp(log_y), 0)
  # cover() validates y in [0, 1]; the cfg above keeps the lognormal median
  # well below 1. Cap catches rare outliers.
  y <- pmin(y, 1 - 1e-6)

  list(data = data.frame(x = x, region = factor(spatial_idx)), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                    sigma_pos = sigma_pos_true, alpha = alpha))
}

simulate_separate_beta_for_recovery <- function(N = 400,
                                                beta_occ = c(-0.4, 0.8),
                                                beta_pos = c(0.4, -0.5),
                                                phi = 30, seed = 5000) {
  set.seed(seed)
  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  occur   <- rbinom(N, 1L, plogis(eta_occ))
  mu_pos  <- plogis(eta_pos)
  y       <- numeric(N); is_pos <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi,
                     (1 - mu_pos[is_pos]) * phi)
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos, phi = phi))
}

chain_adj_for_test <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  adj
}

test_that("joint nested_laplace recovers sigma_pos (lognormal) across 10 seeds", {
  skip_on_cran()
  truth_sigma <- 0.4
  n_seeds <- 10L
  n_s     <- 30L
  adj     <- chain_adj_for_test(n_s)
  sigma_hats <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_joint_lognormal_for_recovery(
      N = 400, n_s = n_s, sigma_pos_true = truth_sigma, seed = 1000L + r
    )
    spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")
    fit <- tobs(
      formula  = ~ x,
      data     = sim$data,
      family   = cover("lognormal"),
      y        = sim$y,
      spatial  = spatial,
      engine   = "nested_laplace",
      control  = list(
        sigma_grid = c(0.3, 0.6, 0.9),
        rho_grid   = c(0.5, 0.7, 0.9),
        alpha_grid = c(0.5, 1.0, 1.5)
      )
    )
    expect_s3_class(fit, "cover_fit")
    expect_true(is.finite(fit$sigma_pos) && fit$sigma_pos > 0)
    sigma_hats[r] <- fit$sigma_pos
  }
  # Recovery: mean estimate within 15% of truth and worst-seed bias < 25%.
  # Tolerances reflect what the probe in dev_notes/probe_joint_recovery.R
  # measured (mean rel err ~ 2%, max rel err ~ 9%) plus a buffer.
  rel_err <- abs(sigma_hats - truth_sigma) / truth_sigma
  expect_lt(abs(mean(sigma_hats) - truth_sigma) / truth_sigma, 0.15)
  expect_lt(max(rel_err), 0.25)
})

test_that("separate-hurdle beta recovers phi_pos across 10 seeds", {
  skip_on_cran()
  truth_phi <- 30
  n_seeds <- 10L
  phi_hats <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_separate_beta_for_recovery(
      N = 400, phi = truth_phi, seed = 5000L + r
    )
    fit <- tobs(
      formula = ~ x,
      data    = sim$data,
      family  = cover(positive = "beta"),
      y       = sim$y
    )
    expect_s3_class(fit, "cover_fit")
    expect_true(is.finite(fit$phi_pos) && fit$phi_pos > 0)
    phi_hats[r] <- fit$phi_pos
  }
  # Recovery: mean within 10% of truth, worst seed within 25%.
  # Probe in dev_notes/probe_separate_beta.R measured mean rel err ~4%,
  # max rel err ~19% at n_pos median ~165 with truth_phi = 30.
  rel_err <- abs(phi_hats - truth_phi) / truth_phi
  expect_lt(abs(mean(phi_hats) - truth_phi) / truth_phi, 0.10)
  expect_lt(max(rel_err), 0.25)
})

test_that("joint nested_laplace beta phi_pos: KNOWN biased (skipped)", {
  # The joint engine pre-fits beta phi on the positive subset alone via
  # tulpa_laplace_beta() (no spatial term). The between-region positive-arm
  # variation that would carry the field signal is absorbed into residual
  # precision, biasing phi downward. INLAabun validation Demo 7 Cell B
  # measured ~70% downward bias at n_pos ~ 45; dev_notes/probe_joint_recovery2.R
  # measured ~55% bias at n_pos ~ 320. Filed as gcol33/tulpaObs#5 -
  # a passing recovery test requires either re-profiling phi on the joint
  # objective (folds back into Phase 3 phi posterior integration) or
  # subtracting alpha_hat * w_s_hat from the positive arm before the
  # phi pre-fit (mirror the lognormal fix in #4).
  skip("blocked on gcol33/tulpaObs#5")
})
