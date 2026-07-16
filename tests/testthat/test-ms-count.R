# ms_count() -- community / multispecies relative-abundance GLMM (the spAbundance
# msAbund model): per-species GLMM on an observed count / continuous response with
# Gaussian community hyperpriors on the coefficients, no detection and no latent
# state. The community analogue of count(), fit by the shared community Laplace-EM.
#
# Recovery-grade (per the "statistical code needs recovery tests" rule):
# community-mean recovery + 95% coverage of the population community mean across
# seeds (the community-mean covariance must propagate the between-species
# variance), plus dispersion recovery for negbin / gaussian. Structural tests
# cover the family wiring, the method registry, and the dispatch gates.

test_that("ms_count() family is wired and reports its supported methods", {
  f <- ms_count()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "ms_count")
  expect_identical(f$status, "working")
  expect_identical(f$params$response, "poisson")
  expect_identical(ms_count("negbin")$params$response, "negbin")
  expect_identical(ms_count("gaussian")$params$response, "gaussian")
  expect_identical(ms_count("binomial")$params$response, "binomial")
  expect_true("laplace" %in% tulpaObs:::.tobs_family_methods$ms_count)
})

test_that("ms_count() dispatch guards reject unsupported inputs", {
  sim <- simulate_ms_count(N = 40, n_species = 5, seed = 1)

  # a detection formula is meaningless (no detection process)
  expect_error(
    tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
         species = colnames(sim$y), detection = ~ 1),
    "no detection process")

  # a structured term is not yet wired (must error, not silently drop)
  expect_error(
    tobs(~ x + re(g), data = cbind(sim$data, g = gl(4, 10)),
         family = ms_count(), y = sim$y, species = colnames(sim$y)),
    "not yet wired")

  # nested_laplace on the community count family needs a shared areal field
  expect_error(
    tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
         species = colnames(sim$y), method = "nested_laplace"),
    "shared areal field")

  # Poisson / negbin require non-negative integer counts
  bad <- sim; bad$y[1, 1] <- 0.5
  expect_error(
    tobs(~ x, data = bad$data, family = ms_count("poisson"), y = bad$y,
         species = colnames(bad$y)),
    "non-negative integer")
})

test_that("ms_count() S3 surface works (coef / ranef / fitted / simulate / WAIC)", {
  skip_on_cran()
  sim <- simulate_ms_count(N = 120, n_species = 8, response = "poisson", seed = 2)
  fit <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
              species = colnames(sim$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_length(unlist(coef(fit)), 2L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_true(all(is.finite(confint(fit))))
  rf <- ranef(fit)
  expect_equal(nrow(rf), 8L * 2L)              # 8 species x 2 coefs
  ft <- fitted(fit)$mu
  expect_equal(dim(ft), c(120L, 8L))
  # per-species fitted correlates with the observed counts
  expect_gt(stats::cor(as.numeric(ft), as.numeric(sim$y)), 0.4)
  sm <- simulate(fit)
  expect_equal(dim(sm), c(120L, 8L))
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
})

test_that("ms_count() Laplace fit accepts missing (NA) site x species entries", {
  skip_on_cran()
  sim <- simulate_ms_count(N = 120, n_species = 8, response = "poisson", seed = 3)
  y <- sim$y
  set.seed(3)
  y[sample(length(y), floor(0.15 * length(y)))] <- NA   # ~15% missing at random
  y[10:110, 2] <- NA                                     # ragged: species 2 sparse
  fit <- tobs(~ x, data = sim$data, family = ms_count(), y = y,
              species = colnames(y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_equal(fit$N, sum(!is.na(y)))                    # N counts observed entries
  expect_true(all(is.finite(unlist(coef(fit)))))
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("Poisson community count recovers community means with ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(1, 0.5)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_count(N = 150, n_species = 20, beta_comm_mean = beta,
                             beta_comm_sd = c(0.4, 0.3), response = "poisson",
                             seed = s)
    fit <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
                species = colnames(sim$y), method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
  }
  # community-mean recovery: the mean over seeds centres on the population mean
  expect_equal(colMeans(est), beta, tolerance = 0.05)
  # pooled coverage of the POPULATION community mean (Vf propagates the
  # between-species variance); the package rubric floor is 0.85 pooled.
  expect_gte(mean(cover), 0.85)
})

test_that("binomial community count recovers means with ~95% coverage (#125)", {
  skip_if_fast()
  skip_on_cran()
  # The binomial (trials > 1) community MEAN carries a small O(1/n_species)
  # first-order-Laplace intercept bias (measured: ~0.06 at S = 20, ~0.036 at
  # S = 40, absent at trials = 1 / the jsdm bernoulli case); the slope is
  # unbiased. Same order and character as the documented negbin-slope
  # attenuation in this family. Recover at S = 30 (bias ~0.04) and assert the
  # coverage, which is the calibration that matters.
  beta <- c(0.2, 0.6)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_count(N = 150, n_species = 30, beta_comm_mean = beta,
                             beta_comm_sd = c(0.4, 0.3), response = "binomial",
                             trials = 10, seed = 800 + s)
    fit <- tobs(~ x, data = sim$data, family = ms_count("binomial"), y = sim$y,
                species = colnames(sim$y), trials = 10, method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
  }
  # Absolute bias bound (the logit-scale intercept sits near 0, so a relative
  # tolerance would be far tighter than the documented O(1/S) bias warrants).
  expect_lt(max(abs(colMeans(est) - beta)), 0.07)
  expect_gte(mean(cover), 0.85)
})

test_that("Gaussian + negbin community count recover means + dispersion", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(1, 0.5)
  n_seed <- 20L

  # Gaussian (identity link, exact Laplace): unbiased community means + residual
  # variance recovery.
  est_g <- matrix(NA_real_, n_seed, length(beta)); vr <- numeric(n_seed)
  cov_g <- matrix(FALSE, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_count(N = 150, n_species = 20, beta_comm_mean = beta,
                             beta_comm_sd = c(0.4, 0.3), response = "gaussian",
                             sd = 0.8, seed = s)
    fit <- tobs(~ x, data = sim$data, family = ms_count("gaussian"), y = sim$y,
                species = colnames(sim$y), method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    b  <- unname(unlist(coef(fit))); se <- sqrt(diag(vcov(fit)))
    est_g[s, ] <- b
    cov_g[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
    vr[s]      <- mean(fit$ms_dispersion$variance)
  }
  expect_equal(colMeans(est_g), beta, tolerance = 0.05)
  expect_gte(mean(cov_g), 0.85)
  expect_equal(mean(vr), 0.8^2, tolerance = 0.1)

  # Negbin (per-species dispersion RE): community-mean recovery + community
  # log-size recovery. The slope carries a mild Laplace-EM attenuation (~10 %),
  # so the tolerance is a touch wider and the pooled coverage floor is 0.80.
  est_n <- matrix(NA_real_, n_seed, length(beta)); mlr <- numeric(n_seed)
  cov_n <- matrix(FALSE, n_seed, length(beta))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_count(N = 150, n_species = 20, beta_comm_mean = beta,
                             beta_comm_sd = c(0.4, 0.3), response = "negbin",
                             size = 3, size.log.sd = 0.3, seed = s)
    fit <- tobs(~ x, data = sim$data, family = ms_count("negbin"), y = sim$y,
                species = colnames(sim$y), method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    b  <- unname(unlist(coef(fit))); se <- sqrt(diag(vcov(fit)))
    est_n[s, ] <- b
    cov_n[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
    mlr[s]     <- fit$ms_dispersion$mu_log_r
  }
  expect_equal(colMeans(est_n), beta, tolerance = 0.07)
  expect_gte(mean(cov_n), 0.80)
  expect_equal(stats::median(mlr), log(3), tolerance = 0.3)
})
