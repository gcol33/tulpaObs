# test-count-gof.R - dispersion / zero-inflation / outlier goodness-of-fit for
# the count families (abun / removal / distance / dyn_abun). Previously these
# families errored via .tobs_gof_require_single (single-season only); now they
# score the per-site TOTAL count against its posterior predictive distribution.

expect_gof <- function(r) {
  expect_true(is.finite(r$ratio))
  expect_true(is.finite(r$p.value) && r$p.value >= 0 && r$p.value <= 1)
}

test_that("count GOF tests run on abun and recover ~1 dispersion when well-fit", {
  skip_on_cran()
  sim <- simulate_abun(N = 200, J = 4, n_abund_covs = 2, n_det_covs = 1, seed = 7)
  fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data, family = abun(K_max = 60),
              detection = ~ det_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  d <- test_dispersion(fit, n.samples = 100)
  z <- test_zero_inflation(fit, n.samples = 100)
  o <- test_outliers(fit, n.samples = 100)
  expect_gof(d); expect_gof(z); expect_gof(o)
  # Data generated from the model: dispersion ratio near 1, not extreme p.
  expect_gt(d$ratio, 0.4); expect_lt(d$ratio, 2.5)
  expect_gt(d$p.value, 0.02); expect_lt(d$p.value, 0.98)
})

test_that("count GOF tests run on removal / distance / dyn_abun", {
  skip_on_cran()
  skip_if_fast()
  simr <- simulate_removal(N = 150, K = 5, n_abund_covs = 2, n_det_covs = 1, seed = 8)
  fitr <- tobs(~ abund_cov1 + abund_cov2, data = simr$data, family = removal(K_max = 60),
               detection = ~ det_cov1, y = simr$y, method = "laplace",
               control = list(verbose = FALSE))
  expect_gof(test_dispersion(fitr, n.samples = 80))
  expect_gof(test_zero_inflation(fitr, n.samples = 80))

  cuts5 <- c(0, 10, 20, 30, 40, 50)
  simd <- simulate_distance(N = 200, cutpoints = cuts5, key = "halfnorm",
                            transect = "line", n_abund_covs = 2, n_sigma_covs = 1,
                            seed = 11)
  fitd <- tobs(~ abund_cov1 + abund_cov2, data = simd$data,
               family = distance(key = "halfnorm", transect = "line",
                                 cutpoints = simd$cutpoints),
               detection = ~ sigma_cov1, y = simd$y, method = "laplace",
               control = list(verbose = FALSE))
  expect_gof(test_dispersion(fitd, n.samples = 80))
  expect_gof(test_outliers(fitd, n.samples = 80))

  simda <- simulate_dyn_abun(N = 100, T = 4, J = 3, n_abund_covs = 1, seed = 9)
  fitda <- tobs(~ abund_cov1, data = simda$data, family = dyn_abun(K_max = 35),
                detection = ~ 1, y = simda$y, method = "laplace",
                control = list(verbose = FALSE))
  expect_gof(test_dispersion(fitda, n.samples = 50))
  expect_gof(test_zero_inflation(fitda, n.samples = 50))
})
