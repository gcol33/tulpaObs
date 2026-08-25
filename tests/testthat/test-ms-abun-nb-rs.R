# Community / multispecies N-mixture with a PER-SPECIES negative-binomial
# dispersion random effect:
#
#   log_r_s ~ N(mu_log_r, sigma_log_r)   (r_s ~ LogNormal)
#
# so (mu_log_r, sigma_log_r) are community hyperparameters alongside
# (mu_lambda, mu_p, Sigma_lambda, Sigma_p). The oracle widens the per-species RE
# vector to (b_lambda_s, b_p_s, b_logr_s); the third covariance block is the
# scalar sigma_log_r^2.
#
# Point recovery and CI coverage for this fixture are in
# test-ms-abun-nb-rs-recovery.R and test-ms-abun-nb-rs-coverage.R.

test_that("ms_abun(negbin) floors the log_r block at the converged AGHQ order", {
  skip_on_cran()
  skip_if_fast()
  # . Two Gauss-Hermite nodes leave no freedom to represent curvature, so where
  # the 1-D log_r posterior is not near-Gaussian the marginal can come out
  # arbitrarily sharp -- a 17x collapse of fit$sds[["log_r"]] with the fit
  # converged, mu_log_r ordinary and nothing warned. The rule is converged at
  # three nodes, so the scalar order is a floor and a request for two does not
  # get two. Asserted on the per-block grid the fit reports rather than on an SE,
  # because the collapse is data-dependent: it needs a fixture whose realised
  # dispersion spread is wide, and only some seeds are.
  sim <- simulate_ms_abun(n_species = 3, N = 12, J = 2,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(3), 0.2), mu_p = c(0.5, -0.2),
                          sd_lambda = 0.3, sd_p = 0.3,
                          mixture = "negbin", size = 5, sigma_logr = 0.4,
                          seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, n.quad = 3L, n.quad.scalar = 2L,
                             max.iter = 5L))

  # Blocks are ordered [lambda | p | log_r]. The coefficient blocks take the
  # requested n.quad; the scalar block is floored, not clamped down to 2.
  grid <- fit$ms_community$n_quad_grid
  expect_length(grid, 3L)
  expect_identical(as.integer(grid[1:2]), c(3L, 3L))
  expect_identical(as.integer(grid[[3L]]), .TOBS_MIN_SCALAR_NQUAD)
  expect_gte(as.integer(grid[[3L]]), 3L)
})
