# test-royle-nichols.R - Royle-Nichols occupancy (abundance-induced detection
# heterogeneity; unmarked occuRN). The latent N ~ Poisson marginalises in closed
# form; the exact marginal is maximised with an observed-information vcov.

test_that("royle_nichols() constructor returns a working tobs_family", {
  f <- royle_nichols()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "royle_nichols")
  expect_identical(f$status, "working")
  expect_identical(f$default_engine, "laplace")
})

test_that("royle_nichols() recovers the generative truth", {
  skip_on_cran()
  skip_if_fast()
  truth <- c("lambda_(Intercept)" = 0.3, "lambda_x" = 0.5, "r_(Intercept)" = -0.8)
  S <- 15
  est <- matrix(NA_real_, S, 3)
  cov95 <- matrix(NA, S, 3)
  for (s in seq_len(S)) {
    sim <- simulate_royle_nichols(N = 300, J = 6, beta_lambda = c(0.3, 0.5),
                                  beta_r = -0.8, seed = 200 + s)
    fit <- tobs(~ x, data = sim$data, family = royle_nichols(),
                detection = ~ 1, y = sim$y, control = list(verbose = FALSE))
    expect_true(isTRUE(fit$convergence$converged))
    m <- fit$means; se <- fit$sds
    est[s, ] <- m[names(truth)]
    cov95[s, ] <- truth >= m[names(truth)] - 1.96 * se[names(truth)] &
                  truth <= m[names(truth)] + 1.96 * se[names(truth)]
  }
  # Point recovery: mean estimate within 12% of truth on each parameter.
  for (j in 1:3) {
    expect_lt(abs(mean(est[, j]) - truth[j]) / abs(truth[j]), 0.12)
  }
  # 95% Wald CI coverage near nominal over the seeds (a 15-seed estimate is
  # noisy; the 25-seed probe in dev_notes shows 0.92-1.00). Guard against gross
  # under-coverage rather than the exact rate.
  expect_gte(min(colMeans(cov95)), 0.75)
})

test_that("royle_nichols() S3 surface works", {
  skip_on_cran()
  sim <- simulate_royle_nichols(N = 150, J = 5, seed = 1)
  fit <- tobs(~ x, data = sim$data, family = royle_nichols(),
              detection = ~ 1, y = sim$y, control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")

  fv <- fitted(fit)
  expect_named(fv, c("lambda", "r", "p"))
  expect_length(fv$lambda, 150L)
  expect_true(all(fv$p >= 0 & fv$p <= 1))

  # predict(): abundance (default) on the intensity scale, detection in [0,1].
  expect_true(all(predict(fit) > 0))
  expect_true(all(predict(fit, type = "detection") >= 0 &
                  predict(fit, type = "detection") <= 1))
  expect_equal(length(predict(fit, newdata = data.frame(x = c(-1, 0, 1)))), 3L)

  rr <- residuals(fit)
  expect_length(rr$occ, 150L)

  # simulate(): posterior replicate detection histories in the input shape.
  sims <- simulate(fit, nsim = 3)
  expect_length(sims, 3L)
  expect_equal(dim(sims[[1]]), c(150L, 5L))

  # WAIC / DIC / CPO score the exact per-site marginal.
  expect_true(is.finite(tobs_waic(fit)$waic))
  expect_true(is.finite(tobs_dic(fit)$dic))
  expect_true(is.finite(tobs_cpo(fit)$lpml))

  co <- coef(fit)
  expect_true(is.list(co) || is.numeric(co))
})

test_that("royle_nichols() rejects unsupported methods", {
  sim <- simulate_royle_nichols(N = 60, J = 4, seed = 2)
  expect_error(
    tobs(~ x, data = sim$data, family = royle_nichols(), detection = ~ 1,
         y = sim$y, method = "nuts", control = list(verbose = FALSE)),
    "royle_nichols")
})
