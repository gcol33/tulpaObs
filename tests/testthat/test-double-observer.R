# =============================================================================
# test-double-observer.R - double-observer abundance (double_observer();
# unmarked multinomPois with the double-observer pi-function).
#
# N ~ Poisson(lambda), two independent observers with detection p1 / p2; the
# three observable cell counts (obs1-only, obs2-only, both) are, by
# Poisson-multinomial thinning, independent Poissons with means lambda*pi_c, so
# the marginal is closed form (no latent-N sum). Fit maximises it (optim BFGS)
# with an observed-information vcov. Non-spatial laplace.
# =============================================================================

test_that("double_observer() constructor + gates", {
  f <- double_observer()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "double_observer")
  sim <- simulate_double_observer(N = 40, seed = 1)
  expect_error(
    tobs(~ abund_cov1, data = sim$data, family = double_observer(),
         detection = ~ det_cov1, y = sim$y, method = "nuts"),
    "laplace")
  # y must be N x 3.
  expect_error(
    tobs(~ abund_cov1, data = sim$data, family = double_observer(),
         detection = ~ det_cov1, y = sim$y[, 1:2]),
    "N x 3")
})

test_that("double_observer() fits + full S3 surface", {
  sim <- simulate_double_observer(N = 200, beta_lambda = c(log(8), 0.4),
                                  beta_p1 = c(stats::qlogis(0.5), 0.2),
                                  beta_p2 = c(stats::qlogis(0.45), -0.1),
                                  seed = 3)
  fit <- tobs(~ abund_cov1, data = sim$data, family = double_observer(),
              detection = ~ det_cov1, y = sim$y,
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$convergence$converged))
  expect_true(all(c("lambda_(Intercept)", "p1_(Intercept)", "p2_(Intercept)")
                  %in% names(fit$means)))
  fv <- fitted(fit)
  expect_named(fv, c("lambda", "p1", "p2", "cell10", "cell01", "cell11"))
  expect_true(all(fv$lambda > 0))
  expect_true(all(fv$p1 > 0 & fv$p1 < 1) && all(fv$p2 > 0 & fv$p2 < 1))
  expect_length(predict(fit, type = "abundance"), 200L)
  expect_equal(ncol(predict(fit, type = "detection")), 2L)
  w <- tobs_waic(fit, n.draws = 200L)
  expect_true(is.finite(w$waic))
  s2 <- simulate(fit, nsim = 1)
  expect_equal(dim(s2), c(200L, 3L))
  expect_length(residuals(fit)$occ, 200L)
})

test_that("double_observer() recovers lambda + per-observer detection (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 15L
  bl0 <- log(9); bl1 <- 0.4; bp1 <- stats::qlogis(0.55); bp2 <- stats::qlogis(0.4)
  li <- ls <- q1 <- q2 <- rep(NA_real_, n_seed)
  hit <- tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_double_observer(
      N = 250, n_abund_covs = 1, n_det_covs = 1,
      beta_lambda = c(bl0, bl1), beta_p1 = c(bp1, 0), beta_p2 = c(bp2, 0),
      seed = 900 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, family = double_observer(),
           detection = ~ det_cov1, y = sim$y,
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    m <- fit$means
    li[s] <- m[["lambda_(Intercept)"]]; ls[s] <- m[["lambda_abund_cov1"]]
    q1[s] <- stats::plogis(m[["p1_(Intercept)"]])
    q2[s] <- stats::plogis(m[["p2_(Intercept)"]])
    tot <- tot + 1L
    if (abs(m[["lambda_abund_cov1"]] - bl1) <= 1.96 * fit$sds[["lambda_abund_cov1"]])
      hit <- hit + 1L
  }
  expect_lt(abs(mean(li, na.rm = TRUE) - bl0), 0.08)
  expect_lt(abs(mean(ls, na.rm = TRUE) - bl1), 0.06)
  expect_lt(abs(mean(q1, na.rm = TRUE) - 0.55), 0.05)
  expect_lt(abs(mean(q2, na.rm = TRUE) - 0.40), 0.05)
  expect_gte(hit / tot, 0.85)
})
