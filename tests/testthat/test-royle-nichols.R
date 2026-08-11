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
  S <- 25
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
  # Point recovery: mean estimate within 12% of truth on each parameter
  # (measured max rel-bias over these seeds is ~0.02).
  for (j in 1:3) {
    expect_lt(abs(mean(est[, j]) - truth[j]) / abs(truth[j]), 0.12)
  }
  # 95% Wald CI coverage at the graduation floor (>= 0.85 over >= 20 seeds).
  # Measured min coverage over these fixed seeds is 0.88 (deterministic).
  expect_gte(min(colMeans(cov95)), 0.85)
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
  expect_true(is.finite(waic(fit)$waic))
  expect_true(is.finite(dic(fit)$dic))
  expect_true(is.finite(cpo(fit)$lpml))

  co <- coef(fit)
  expect_true(is.list(co) || is.numeric(co))
})

test_that("royle_nichols() visit marginal reduces to the site marginal", {
  # With detection constant within a site the full per-visit product must equal
  # the (k, n) sufficient-statistic form exactly.
  set.seed(1)
  S <- 25L; J <- 6L; K <- 45L
  lambda   <- runif(S, 0.5, 3)
  r_site   <- runif(S, 0.1, 0.6)
  y        <- matrix(rbinom(S * J, 1L, 0.3), S, J)
  site_ll  <- tulpaObs:::.rn_site_loglik(lambda, r_site, rowSums(y == 1L),
                                         rep(J, S), K)
  site_idx <- rep(seq_len(S), each = J)
  visit_ll <- tulpaObs:::.rn_visit_loglik(lambda, r_site[site_idx],
                                          as.vector(t(y)), site_idx, K, S)
  expect_lt(max(abs(site_ll - visit_ll)), 1e-8)
})

test_that("royle_nichols() recovers visit-varying detection", {
  skip_on_cran()
  skip_if_fast()
  truth <- c("lambda_(Intercept)" = 0.3, "lambda_x" = 0.5,
             "r_(Intercept)" = -0.8, "r_w" = 0.8)
  S <- 20
  est <- matrix(NA_real_, S, 4)
  cov95 <- matrix(NA, S, 4)
  for (s in seq_len(S)) {
    sim <- simulate_royle_nichols(N = 300, J = 6, beta_lambda = c(0.3, 0.5),
                                  beta_r = -0.8, beta_r_visit = 0.8,
                                  seed = 400 + s)
    fit <- tobs(~ x, data = sim$data, family = royle_nichols(),
                detection = ~ w, y = sim$y, visits = sim$visits,
                control = list(verbose = FALSE))
    expect_true(isTRUE(fit$convergence$converged))
    m <- fit$means[names(truth)]; se <- fit$sds[names(truth)]
    est[s, ] <- m
    cov95[s, ] <- truth >= m - 1.96 * se & truth <= m + 1.96 * se
  }
  # Point recovery within 12% of truth (measured max rel-bias ~0.04).
  for (j in 1:4) {
    expect_lt(abs(mean(est[, j]) - truth[j]) / abs(truth[j]), 0.12)
  }
  # 95% Wald CI coverage at the graduation floor (measured min ~0.95).
  expect_gte(min(colMeans(cov95)), 0.85)
})

test_that("royle_nichols() visit-varying S3 surface works", {
  skip_on_cran()
  sim <- simulate_royle_nichols(N = 150, J = 5, beta_r_visit = 0.8, seed = 7)
  fit <- tobs(~ x, data = sim$data, family = royle_nichols(), detection = ~ w,
              y = sim$y, visits = sim$visits, control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true("r_w" %in% names(fit$means))

  fv <- fitted(fit)
  expect_equal(dim(fv$r), c(150L, 5L))          # per-visit detection grid
  expect_equal(dim(fv$p), c(150L, 5L))
  expect_true(all(fv$p >= 0 & fv$p <= 1, na.rm = TRUE))

  # predict(detection) returns the fitted grid; new visit covariates error.
  expect_equal(dim(predict(fit, type = "detection")), c(150L, 5L))
  expect_error(predict(fit, type = "detection",
                       newdata = data.frame(w = 0)), "visit-level")

  expect_length(residuals(fit)$occ, 150L)
  expect_equal(dim(simulate(fit)), c(150L, 5L))
  expect_true(is.finite(waic(fit)$waic))
  expect_true(is.finite(dic(fit)$dic))
})

test_that("royle_nichols() rejects unsupported methods", {
  sim <- simulate_royle_nichols(N = 60, J = 4, seed = 2)
  expect_error(
    tobs(~ x, data = sim$data, family = royle_nichols(), detection = ~ 1,
         y = sim$y, method = "nuts", control = list(verbose = FALSE)),
    "royle_nichols")
})
