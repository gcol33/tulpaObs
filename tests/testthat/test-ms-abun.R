# Community / multispecies N-mixture (ms_abun): parameter recovery, community-
# mean CI coverage, and S3 method coverage. The fit is tulpa's C++ Laplace-EM
# (per-species coefficient RE with Gaussian community covariances).

test_that("ms_abun recovers community means and per-species coefficients", {
  skip_on_cran()
  set.seed(11)
  sim <- simulate_ms_abun(n_species = 14, N = 90, J = 4,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.4, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4, seed = 11)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(), detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means: within ~2.5 SE of truth on every coordinate.
  est <- fit$means
  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  z <- abs(est - truth) / fit$sds
  expect_true(all(z < 2.5))

  # Per-species coefficients: high correlation with the simulated truth.
  cm <- fit$ms_community
  expect_gt(min(diag(cor(cm$coef_lambda, sim$truth$beta_lambda))), 0.9)
  expect_gt(min(diag(cor(cm$coef_p,      sim$truth$beta_p))),      0.85)

  # Community SDs track the realized per-species spread (not biased to ~0).
  emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))
  est_sd <- c(cm$sd_lambda, cm$sd_p)
  expect_true(all(est_sd > 0.5 * emp_sd & est_sd < 1.6 * emp_sd))
})

test_that("ms_abun community-mean 95% CIs cover at the nominal rate", {
  skip_on_cran()
  n_seed <- 20L
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 12, N = 60, J = 4,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.4),
                            sd_lambda = 0.5, sd_p = 0.4, seed = 100 + s)
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
                family = ms_abun(), detection = ~ det_cov1,
                species = sim$species, method = "laplace",
                control = list(verbose = FALSE))
    truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered <- c(covered, truth >= lo & truth <= hi)
  }
  # Nominal 95%; allow Monte-Carlo slack on 16 x 4 = 64 intervals.
  expect_gt(mean(covered), 0.85)
})

test_that("ms_abun S3 methods work", {
  skip_on_cran()
  set.seed(3)
  sim <- simulate_ms_abun(n_species = 8, N = 40, J = 3, seed = 3)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(), detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))

  cf <- coef(fit)
  expect_true(is.list(cf) && all(c("lambda", "p") %in% names(cf)))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 8L * (2L + 2L))   # n_species x (p_lambda + p_p)
  expect_true(all(c("species", "arm", "term", "estimate") %in% names(re)))

  fv <- fitted(fit)
  expect_equal(dim(fv$lambda), c(40L, 8L))
  expect_true(all(fv$lambda > 0))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(40L, 3L, 8L))

  expect_type(nobs(fit), "integer")
})
