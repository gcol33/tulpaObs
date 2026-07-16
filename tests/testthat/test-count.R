# count() -- GLMM on the observed count / continuous response directly (no
# detection, no latent state); the relative-abundance model of spAbundance
# (abund). Poisson / negative-binomial (log link) and Gaussian (identity).
#
# Recovery-grade tests (per the "statistical code needs recovery tests" rule):
# point recovery against simulated truth + 95% CI coverage across seeds, plus
# dispersion recovery for negbin / gaussian. The structural tests check the
# family wiring, method registry, dispatch guards, and S3 surface.

test_that("count() family is wired and reports its supported methods", {
  f <- count()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "count")
  expect_identical(f$status, "working")
  expect_identical(f$params$response, "poisson")
  expect_identical(f$response, "vector")
  expect_identical(count("negbin")$params$response, "negbin")
  expect_identical(count("gaussian")$params$response, "gaussian")
  expect_identical(count("binomial")$params$response, "binomial")

  # the method registry is the single source of truth for supported backends
  expect_true("laplace" %in% tulpaObs:::.tobs_family_methods$count)
})

test_that("count() dispatch guards reject unsupported inputs", {
  sim <- simulate_count(N = 40, beta = c(0.5, 0.3), seed = 1)

  # a detection formula is meaningless (no detection process)
  expect_error(
    tobs(~ x, data = sim$data, family = count(), y = sim$y, detection = ~ 1),
    "no detection process")

  # a structured term is not yet wired (must error, not silently drop)
  expect_error(
    tobs(~ x + re(g), data = cbind(sim$data, g = gl(4, 10)),
         family = count(), y = sim$y),
    "not yet wired")

  # Poisson / negbin require non-negative integer counts
  expect_error(
    tobs(~ x, data = sim$data, family = count("poisson"),
         y = sim$y + 0.5),
    "non-negative integer")
})

test_that("count() accepts the response on a two-sided formula LHS", {
  skip_if_fast()
  sim <- simulate_count(N = 200, beta = c(0.7, 0.5), response = "poisson",
                        seed = 3)
  df  <- cbind(sim$data, yy = sim$y)
  fit <- tobs(yy ~ x, data = df, family = count("poisson"),
              control = list(progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_equal(unname(unlist(coef(fit))), sim$truth$beta, tolerance = 0.15)
  # the response given twice must error
  expect_error(
    tobs(yy ~ x, data = df, family = count("poisson"), y = sim$y),
    "given twice")
})

test_that("count() S3 surface works (fitted / predict / residuals / WAIC)", {
  skip_if_fast()
  sim <- simulate_count(N = 250, beta = c(0.8, 0.6), response = "poisson",
                        seed = 5)
  fit <- tobs(~ x, data = sim$data, family = count("poisson"), y = sim$y,
              control = list(progress = FALSE))

  expect_length(fitted(fit)$mu, 250)
  expect_true(all(fitted(fit)$mu > 0))
  for (ty in c("deviance", "pearson", "response")) {
    r <- residuals(fit, type = ty)$mu
    expect_length(r, 250)
    expect_true(all(is.finite(r)))
  }
  # predict() uses the point-estimate linear predictor, so it must equal
  # fitted() in-sample and exp(X . means) on newdata (the package convention:
  # fitted / predict plug in `means`, matching jsdm and the occupancy families).
  expect_equal(unname(predict(fit, newdata = sim$data)), unname(fitted(fit)$mu),
               tolerance = 1e-8)
  nd <- data.frame(x = c(-1, 0, 1))
  pr <- predict(fit, newdata = nd)
  b  <- fit$means[seq_len(2)]
  expect_equal(pr, exp(b[1] + b[2] * c(-1, 0, 1)), tolerance = 1e-6)
  expect_true(all(diff(pr) > 0))   # monotone increasing in x (positive slope)
  expect_identical(nobs(fit), 250L)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
  expect_true(w$p_waic > 0 && w$p_waic < 6)   # ~2 fixed effects
})

test_that("Poisson count fit recovers truth + 95% CI coverage", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(0.8, 0.6, -0.4)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_count(N = 400, beta = beta, response = "poisson",
                          seed = 100 + s)
    fit <- tobs(~ x + x2, data = sim$data, family = count("poisson"),
                y = sim$y, control = list(progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.05)
  expect_true(all(colMeans(cover) >= 0.85))
})

test_that("Negative-binomial count fit recovers coefficients + size", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(1.0, 0.5)
  size <- 2
  n_seed <- 20L
  est <- matrix(NA_real_, n_seed, length(beta))
  siz <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_count(N = 500, beta = beta, response = "negbin",
                          size = size, seed = 200 + s)
    fit <- tobs(~ x, data = sim$data, family = count("negbin"), y = sim$y,
                control = list(progress = FALSE))
    est[s, ] <- unname(unlist(coef(fit)))
    siz[s]   <- fit$count_dispersion$phi
  }
  expect_equal(colMeans(est), beta, tolerance = 0.06)
  # the size is harder to pin down; the median should be within ~40%
  expect_equal(stats::median(siz), size, tolerance = 0.8)
})

test_that("Gaussian count fit recovers coefficients + residual variance", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(2.0, -0.7)
  sd_true <- 1.5
  n_seed <- 20L
  est <- matrix(NA_real_, n_seed, length(beta))
  vr  <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_count(N = 400, beta = beta, response = "gaussian",
                          sd = sd_true, seed = 300 + s)
    fit <- tobs(~ x, data = sim$data, family = count("gaussian"), y = sim$y,
                control = list(progress = FALSE))
    est[s, ] <- unname(unlist(coef(fit)))
    vr[s]    <- fit$count_dispersion$phi
  }
  expect_equal(colMeans(est), beta, tolerance = 0.05)
  expect_equal(mean(vr), sd_true^2, tolerance = 0.3)
})


# --- binomial GLMM (spOccupancy svcPGBinom family, no replicates) ----------

test_that("count(response = 'binomial') validates trials and the response", {
  sim <- simulate_count(N = 60, beta = c(0.3, 0.5), response = "binomial",
                        trials = 8, seed = 1)

  # `trials` only applies to the binomial response
  expect_error(
    tobs(~ x, data = sim$data, family = count("poisson"),
         y = round(exp(sim$data$x)), trials = 3),
    "trials.*only.*binomial|binomial")

  # successes cannot exceed trials
  expect_error(
    tobs(~ x, data = sim$data, family = count("binomial"),
         y = sim$y, trials = 1),
    "<=|exceed|k <= n")

  # trials length must be scalar or one per site
  expect_error(
    tobs(~ x, data = sim$data, family = count("binomial"),
         y = sim$y, trials = c(8, 8)),
    "scalar or one value per site|length")

  # default trials = 1 gives a Bernoulli response
  s2 <- simulate_count(N = 60, beta = c(0, 0.5), response = "binomial",
                       trials = 1, seed = 2)
  fit <- tobs(~ x, data = s2$data, family = count("binomial"), y = s2$y,
              control = list(progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
})

test_that("count() S3 surface works for the binomial response", {
  skip_if_fast()
  sim <- simulate_count(N = 250, beta = c(0.4, 0.7), response = "binomial",
                        trials = 12, seed = 5)
  fit <- tobs(~ x, data = sim$data, family = count("binomial"), y = sim$y,
              trials = 12, control = list(progress = FALSE))
  # fitted = expected successes n * p, in [0, n]
  fv <- fitted(fit)$mu
  expect_length(fv, 250)
  expect_true(all(fv >= 0 & fv <= 12))
  for (ty in c("deviance", "pearson", "response")) {
    r <- residuals(fit, type = ty)$mu
    expect_length(r, 250)
    expect_true(all(is.finite(r)))
  }
  # predict(newdata) returns the per-trial probability in (0, 1)
  pr <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
  expect_true(all(pr > 0 & pr < 1))
  expect_true(all(diff(pr) > 0))              # monotone in x (positive slope)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
})

test_that("binomial count fit recovers truth + 95% CI coverage (trials > 1)", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(0.3, 0.8, -0.5)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_count(N = 400, beta = beta, response = "binomial",
                          trials = 10, seed = 400 + s)
    fit <- tobs(~ x + x2, data = sim$data, family = count("binomial"),
                y = sim$y, trials = 10, control = list(progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.05)
  expect_true(all(colMeans(cover) >= 0.85))
})

test_that("Bernoulli count fit recovers truth + coverage (trials = 1)", {
  skip_if_fast()
  skip_on_cran()
  # trials = 1 is svcPGBinom's setting; a single Bernoulli per site needs more
  # sites than the trials > 1 case for the same precision.
  beta <- c(-0.2, 0.9)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_count(N = 800, beta = beta, response = "binomial",
                          trials = 1, seed = 500 + s)
    fit <- tobs(~ x, data = sim$data, family = count("binomial"),
                y = sim$y, trials = 1, control = list(progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.06)
  expect_true(all(colMeans(cover) >= 0.85))
})
