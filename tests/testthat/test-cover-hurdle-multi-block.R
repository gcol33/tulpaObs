# Phase J-D wire-through: cover_hurdle with multi-block latent prior
# (BYM2 spatial + AR1 temporal + IID observer) via tulpa's multi-block
# joint nested-Laplace engine.
#
# Verifies:
#   * tobs(cover, ~ x + bym2() + temporal() + re(), engine =
#     "nested_laplace") fits end-to-end and returns a `cover_fit`.
#   * Beta posterior-weighted estimates land on the correct side of zero.
#   * The fit's underlying joint object is the multi-block class.
#   * Multi-block hyperparameter summary is attached (block_moments
#     populated, alpha finite).
#
# Uses positive = "beta" so phi is integrated on the outer grid. (The
# lognormal arm now integrates its noise SD the same way.)

simulate_cover_multi_block <- function(N = 400, n_s = 16L, n_years = 6L,
                                       n_obs = 12L,
                                       sigma = 0.6, rho = 0.7, alpha = 1.2,
                                       sigma_year = 0.3, sigma_obs = 0.25,
                                       rho_ar = 0.6,
                                       beta_occ = c(0.2, 0.7),
                                       beta_pos = c(-0.5, -0.3),
                                       phi_b = 30,
                                       seed = 7001) {
  set.seed(seed)
  # 4 x 4 chain-of-chains adjacency for BYM2.
  grid_n <- as.integer(sqrt(n_s))
  if (grid_n * grid_n != n_s) {
    stop("simulate_cover_multi_block(): n_s must be a perfect square")
  }
  adj <- matrix(0L, n_s, n_s)
  for (i in seq_len(grid_n)) for (j in seq_len(grid_n)) {
    s <- (i - 1L) * grid_n + j
    if (i > 1L)      adj[s, (i - 2L) * grid_n + j] <- 1L
    if (i < grid_n)  adj[s, i * grid_n + j]        <- 1L
    if (j > 1L)      adj[s, (i - 1L) * grid_n + (j - 1L)] <- 1L
    if (j < grid_n)  adj[s, (i - 1L) * grid_n + (j + 1L)] <- 1L
  }

  s_idx <- sample.int(n_s, N, replace = TRUE)
  t_idx <- sample.int(n_years, N, replace = TRUE)
  o_idx <- sample.int(n_obs, N, replace = TRUE)

  phi_f   <- rnorm(n_s); phi_f   <- phi_f   - mean(phi_f)
  theta_f <- rnorm(n_s); theta_f <- theta_f - mean(theta_f)
  w_unit  <- sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f
  ar_t    <- numeric(n_years)
  ar_t[1L] <- rnorm(1, 0, sigma_year)
  for (t in 2:n_years) {
    ar_t[t] <- rho_ar * ar_t[t - 1L] +
      rnorm(1, 0, sigma_year * sqrt(1 - rho_ar^2))
  }
  iota_o  <- rnorm(n_obs, 0, sigma_obs)

  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x +
             sigma * w_unit[s_idx] + ar_t[t_idx] + iota_o[o_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))

  eta_pos <- beta_pos[1] + beta_pos[2] * x +
             (alpha * sigma) * w_unit[s_idx] +
             ar_t[t_idx] + iota_o[o_idx]
  mu      <- plogis(eta_pos)
  cov_pos <- rbeta(N, mu * phi_b, (1 - mu) * phi_b)
  cov_pos <- pmin(pmax(cov_pos, 1e-6), 1 - 1e-6)
  y       <- ifelse(occur == 1L, cov_pos, 0)

  list(
    data = data.frame(
      x      = x,
      region = factor(s_idx, levels = seq_len(n_s)),
      year   = t_idx,
      obs    = o_idx
    ),
    y     = y,
    adj   = adj,
    truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                 sigma = sigma, alpha = alpha, rho = rho,
                 sigma_year = sigma_year, sigma_obs = sigma_obs)
  )
}

test_that("cover(beta) with spatial + temporal + RE fits via multi-block", {
  sim <- simulate_cover_multi_block(N = 400, seed = 7001)
  adj <- sim$adj

  # 6 spatial x 2 temporal x 2 RE x 4 phi = 96 cells (above the engine's
  # 50-cell warn threshold but well under the 1000-cell hard cap).
  fit <- suppressWarnings(tobs(
    formula  = ~ x + bym2(graph = adj, group_var = "region") +
                 temporal(year, type = "ar1") + re(obs, type = "iid"),
    data     = sim$data,
    family   = cover("beta"),
    y        = sim$y,
    method   = "nested_laplace",
    control  = list(
      sigma_grid         = c(0.3, 0.6, 1.0),
      rho_grid           = c(0.5, 0.85),
      sigma_pos_grid     = c(0.4, 0.8, 1.2),
      tau_temporal_grid  = c(4, 16),
      rho_temporal_grid  = c(0.3, 0.7),
      sigma_re_grid      = c(0.15, 0.4),
      phi_grid           = exp(seq(log(5), log(80), length.out = 4)),
      adaptive_grid      = FALSE
    )
  ))

  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "beta")
  expect_true(fit$converged)
  expect_true(inherits(fit$joint, "tulpa_nested_laplace_joint_multi"))

  # Beta point estimates land on the correct side of zero.
  expect_true(all(is.finite(fit$beta_occ)))
  expect_true(all(is.finite(fit$beta_pos)))
  expect_gt(fit$beta_occ[2], 0)   # truth 0.7
  expect_lt(fit$beta_pos[2], 0)   # truth -0.3

  # phi_pos finite and inside the user grid.
  expect_true(is.finite(fit$phi_pos))
  expect_gt(fit$phi_pos, 1)
  expect_lt(fit$phi_pos, 200)

  # Multi-block hyperparameter summary: spatial + temporal + RE blocks all
  # report sensible posterior means inside their grids.
  # The R-facing outer grid for the BYM2 copy block lives in (sigma,
  # alpha) space; the C++ kernel sees (sigma_occ, sigma_pos) but those
  # names are not exposed at the R block_moments layer (see tulpa's
  # nested_laplace_joint.R "API contract" comment around L672).
  bm <- fit$joint$block_moments
  expect_length(bm, 3L)
  expect_named(bm[[1L]]$mean, c("sigma", "alpha", "rho"))
  expect_named(bm[[2L]]$mean, c("tau", "rho"))
  expect_named(bm[[3L]]$mean, "sigma")
  expect_true(bm[[1L]]$mean[["sigma"]] > 0.1 &&
              bm[[1L]]$mean[["sigma"]] < 1.5)
  expect_true(bm[[1L]]$mean[["alpha"]] > 0.3 &&
              bm[[1L]]$mean[["alpha"]] < 3.0)
  expect_true(bm[[2L]]$mean[["tau"]] >= 4 &&
              bm[[2L]]$mean[["tau"]] <= 16)
  expect_true(bm[[3L]]$mean[["sigma"]] >= 0.15 &&
              bm[[3L]]$mean[["sigma"]] <= 0.4)

  # alpha = sigma_pos / sigma_occ exposed via theta_mean on the copy
  # block. In the multi-block path joint_grid columns are prefixed
  # `b<N>.` (the copy block here is block 1, so `b1.alpha`).
  expect_true("b1.alpha" %in% names(fit$joint$theta_mean))
  expect_true(is.finite(fit$joint$theta_mean[["b1.alpha"]]))
  expect_gt(fit$joint$theta_mean[["b1.alpha"]], 0)
})


test_that("cover(): multi-block rejects method = 'laplace'", {
  sim <- simulate_cover_multi_block(N = 200, seed = 7002)
  adj <- sim$adj
  expect_error(
    tobs(
      formula  = ~ x + bym2(graph = adj, group_var = "region") +
                   temporal(year, type = "ar1"),
      data     = sim$data,
      family   = cover("beta"),
      y        = sim$y,
      method   = "laplace"
    ),
    regexp = "require method = 'nested_laplace'"
  )
})


test_that("cover(): multi-block resolves character group / time columns", {
  # Smoke test that the temporal()/re() terms resolve both factor (region,
  # obs) and plain integer (year) columns to index codes. The full fit runs
  # end-to-end in the first test; this checks the lightweight code path.
  sim <- simulate_cover_multi_block(N = 200, seed = 7004)
  sim$data$year   <- as.integer(sim$data$year)
  sim$data$obs    <- factor(sim$data$obs)
  adj <- sim$adj

  # Grid trips tulpa's >50-cell soft warning (the copy block inserts midpoint
  # cells); grid size is irrelevant to this column-resolution smoke test, so
  # suppress it, as the recovery test above does for the same reason.
  fit <- suppressWarnings(tobs(
    formula  = ~ x + bym2(graph = adj, group_var = "region") +
                 temporal(year, type = "iid") + re(obs, type = "iid"),
    data     = sim$data,
    family   = cover("beta"),
    y        = sim$y,
    method   = "nested_laplace",
    control  = list(
      sigma_grid         = c(0.4, 0.8),
      rho_grid           = c(0.5, 0.85),
      sigma_pos_grid     = c(0.6, 1.0),
      sigma_temporal_grid = c(0.2, 0.4),
      sigma_re_grid      = c(0.2, 0.4),
      phi_grid           = c(10, 40),
      adaptive_grid      = FALSE
    )
  ))
  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)
  expect_true(inherits(fit$joint, "tulpa_nested_laplace_joint_multi"))
})
