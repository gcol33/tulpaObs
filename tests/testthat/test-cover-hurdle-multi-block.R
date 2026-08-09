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
#   * Each block's posterior mean follows its own simulated truth across a
#     pair of fits on one grid.
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

# gcol33/tulpaObs#199. A posterior-weighted mean over a pinned axis is a convex
# combination of that axis's nodes, so a band drawn at or outside the node span
# holds by construction. Five of this test's bands were drawn that way, two of
# them exactly at the pin. What replaces them is a paired fit: one grid, two
# simulated truths, and one ordering per block, so each assertion fails when
# that block's own truth stops moving.
#
# Grid (every axis pinned, `adaptive.grid = FALSE`): 3 sigma x 3 alpha x 2 rho x
# 4 tau x 1 rho_temporal x 3 sigma_re x 4 phi = 864 base cells. Each truth below
# sits inside its axis's span and on none of its nodes. The block axes are
# asserted below to be exactly these; the phi axis can pick up a few extra nodes
# from the engine's var-of-means consistency pass, which is a separate mechanism
# from the adaptive grid this switches off.
.mb_grid <- list(
  sigma.grid         = c(0.2, 0.45, 0.9),
  rho.grid           = c(0.5, 0.85),
  alpha.grid         = c(0.4, 1.0, 2.5),
  tau.temporal.grid  = c(1, 4, 16, 64),
  rho.temporal.grid  = 0.6,
  sigma.re.grid      = c(0.06, 0.2, 0.7),
  phi.grid           = c(6, 15, 38, 95),
  adaptive.grid      = FALSE
)

.mb_fit <- function(sim) {
  suppressWarnings(tobs(
    formula  = ~ x + bym2(graph = sim$adj, group_var = "region") +
                 temporal(year, type = "ar1") + re(obs, type = "iid"),
    data     = sim$data,
    family   = cover("beta"),
    y        = sim$y,
    method   = "nested_laplace",
    control  = .mb_grid
  ))
}

test_that("cover(beta) with spatial + temporal + RE fits via multi-block", {
  skip_if_fast()
  sim <- simulate_cover_multi_block(N = 400, seed = 7001)
  fit <- .mb_fit(sim)

  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "beta")
  expect_true(fit$converged)
  expect_true(inherits(fit$joint, "tulpa_nested_laplace_joint_multi"))

  # Beta point estimates land on the correct side of zero.
  expect_true(all(is.finite(fit$beta_occ)))
  expect_true(all(is.finite(fit$beta_pos)))
  expect_gt(fit$beta_occ[2], 0)   # truth 0.7
  expect_lt(fit$beta_pos[2], 0)   # truth -0.3
  expect_true(is.finite(fit$phi_pos))

  # Multi-block hyperparameter summary: three blocks, named moments.
  # The R-facing outer grid for the BYM2 copy block lives in (sigma, alpha)
  # space; the C++ kernel sees `b1.sigma_occ` / `b1.sigma_pos`, materialized
  # from (sigma, alpha) at the kernel-call boundary and not exposed at the R
  # block_moments layer (tulpa's `.joint_call_kernel_via_multi()` in
  # R/nested_laplace_joint_backends.R carries the per-type axis table).
  bm <- fit$joint$block_moments
  expect_length(bm, 3L)
  expect_named(bm[[1L]]$mean, c("sigma", "alpha", "rho"))
  expect_named(bm[[2L]]$mean, c("tau", "rho"))
  expect_named(bm[[3L]]$mean, "sigma")

  # Every block axis the driver placed is the one asked for, per block.
  tg <- fit$joint$theta_grid
  expect_equal(sort(unique(tg[, "b1.sigma"])), .mb_grid$sigma.grid)
  expect_equal(sort(unique(tg[, "b1.alpha"])), .mb_grid$alpha.grid)
  expect_equal(sort(unique(tg[, "b2.tau"])),   .mb_grid$tau.temporal.grid)
  expect_equal(sort(unique(tg[, "b3.sigma"])), .mb_grid$sigma.re.grid)

  expect_true("b1.alpha" %in% names(fit$joint$theta_mean))
  expect_true(is.finite(fit$joint$theta_mean[["b1.alpha"]]))
})


test_that("cover(): each multi-block hyperparameter follows its own truth", {
  skip_if_fast()
  # Same grid, two truths. A is the fixture's own configuration; B raises the
  # copy coefficient and both non-spatial SDs, lowers the field SD, and makes
  # the beta arm far more disperse. Measured on seed 7001 (tulpa 0.0.163):
  #
  #   quantity       A        B        ratio  margin asserted
  #   b1.sigma     0.8933   0.4480     1.99x   1.5x
  #   b1.alpha     1.0013   1.5717     1.57x   1.3x
  #   b2.tau      51.4964   8.7709     5.87x   2.0x
  #   b3.sigma     0.2395   0.7000     2.92x   2.0x
  #   phi_pos     38.0000   6.0092     6.32x   3.0x
  #
  # Every ordering holds on all five of seeds 7001 / 7011-7014, and each rule
  # fails when its own component of B is put back to A's value (the sigma_year
  # reversal takes down `b2.tau` AND `b1.sigma`, the year effect and the field
  # being partly confounded at this fixture size).
  fit_a <- .mb_fit(simulate_cover_multi_block(N = 400, seed = 7001))
  fit_b <- .mb_fit(simulate_cover_multi_block(N = 400, seed = 7001,
                                              sigma      = 0.35,
                                              alpha      = 2.2,
                                              sigma_year = 0.8,
                                              sigma_obs  = 0.6,
                                              phi_b      = 8))
  a <- fit_a$joint$block_moments
  b <- fit_b$joint$block_moments

  # Spatial field SD: truth 0.6 -> 0.35.
  expect_gt(a[[1L]]$mean[["sigma"]], 1.5 * b[[1L]]$mean[["sigma"]])
  # Copy coefficient: truth 1.2 -> 2.2.
  expect_gt(b[[1L]]$mean[["alpha"]], 1.3 * a[[1L]]$mean[["alpha"]])
  # AR1 precision: truth SD 0.3 -> 0.8, so tau falls.
  expect_gt(a[[2L]]$mean[["tau"]], 2.0 * b[[2L]]$mean[["tau"]])
  # Observer RE SD: truth 0.25 -> 0.6.
  expect_gt(b[[3L]]$mean[["sigma"]], 2.0 * a[[3L]]$mean[["sigma"]])
  # Beta precision: truth 30 -> 8.
  expect_gt(fit_a$phi_pos, 3.0 * fit_b$phi_pos)
})


# gcol33/tulpaObs#192. The multi-block copy spec used to carry the copy axis on
# `sigma_pos_grid`, a field tulpa's copy resolver does not read, so the axis
# integrated was the engine's own default whatever the caller asked for -- two
# grids an order of magnitude apart gave a bit-identical `log_marginal`. The
# assertion is therefore that the knob CHANGES the fit: a band the default also
# satisfies is what let this survive.
test_that("cover(): control$alpha.grid places the multi-block copy axis", {
  skip_if_fast()
  sim <- simulate_cover_multi_block(N = 300, seed = 7005)
  adj <- sim$adj

  fit_at <- function(alpha_grid) {
    suppressWarnings(tobs(
      formula  = ~ x + bym2(graph = adj, group_var = "region") +
                   temporal(year, type = "ar1"),
      data     = sim$data,
      family   = cover("beta"),
      y        = sim$y,
      method   = "nested_laplace",
      control  = list(
        sigma.grid        = c(0.4, 0.8),
        rho.grid          = 0.85,
        alpha.grid        = alpha_grid,
        tau.temporal.grid = 9,
        rho.temporal.grid = 0.5,
        phi.grid          = c(12, 40),
        adaptive.grid     = FALSE
      )
    ))
  }

  lo <- fit_at(c(0.2, 0.6))
  hi <- fit_at(c(1.4, 2.2))

  # The axis the driver placed is the one asked for, not the engine default.
  expect_equal(sort(unique(lo$joint$theta_grid[, "b1.alpha"])), c(0.2, 0.6))
  expect_equal(sort(unique(hi$joint$theta_grid[, "b1.alpha"])), c(1.4, 2.2))

  # And it reaches the likelihood: two disjoint axes cannot integrate to the
  # same marginal.
  expect_false(isTRUE(all.equal(sum(exp(lo$joint$log_marginal -
                                        max(lo$joint$log_marginal))),
                                sum(exp(hi$joint$log_marginal -
                                        max(hi$joint$log_marginal))))))
  expect_false(isTRUE(all.equal(max(lo$joint$log_marginal),
                                max(hi$joint$log_marginal))))
  expect_false(isTRUE(all.equal(unname(lo$beta_pos), unname(hi$beta_pos))))
})

test_that("cover(): the retired sigma.pos.grid knob is refused, not ignored", {
  sim <- simulate_cover_multi_block(N = 120, seed = 7006)
  expect_error(
    tobs(
      formula  = ~ x + bym2(graph = sim$adj, group_var = "region") +
                   temporal(year, type = "ar1"),
      data     = sim$data,
      family   = cover("beta"),
      y        = sim$y,
      method   = "nested_laplace",
      control  = list(sigma.pos.grid = c(0.4, 0.8, 1.2))
    ),
    "control$alpha.grid", fixed = TRUE)
})


test_that("cover(): multi-block rejects method = 'laplace'", {
  skip_if_fast()
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
  skip_if_fast()
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
      sigma.grid         = c(0.4, 0.8),
      rho.grid           = c(0.5, 0.85),
      sigma.temporal.grid = c(0.2, 0.4),
      sigma.re.grid      = c(0.2, 0.4),
      phi.grid           = c(10, 40),
      adaptive.grid      = FALSE
    )
  ))
  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)
  expect_true(inherits(fit$joint, "tulpa_nested_laplace_joint_multi"))
})
