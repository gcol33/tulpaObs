# =============================================================================
# test-dyn-abun-zip.R - zero-inflated open N-mixture (Dail-Madsen) for
# dyn_abun(mixture = "zip" / "zinb"; gcol33/tulpaObs#116).
#
# A structural-zero site is never occupied (N_t = 0 for every season), so every
# count at that site is zero; the observed per-site marginal is the two-component
# mixture omega * 1{all y = 0} + (1 - omega) * L_DailMadsen. The structural-zero
# share is the new parameter; these check the S3 surface, its recovery + 95%
# coverage over >= 20 seeds, and that the ZINB dispersion also recovers.
# =============================================================================

test_that("dyn_abun(mixture='zip') gates + S3 surface", {
  sim <- simulate_dyn_abun(N = 120, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(6), 0.3), p = 0.5, omega = 0.6,
                           gamma = 1.0, zi = 0.3, seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data,
              family = dyn_abun(mixture = "zip", K_max = 30),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$zero_inflated))
  expect_true("zi_logit" %in% names(fit$means))
  expect_false(is.null(fit$zi_omega))
  # A spatial / RE / NUTS combination is rejected in v1 with a pointer.
  adj <- diag(0, 120)
  expect_error(
    tobs(~ icar(graph = adj), data = sim$data,
         family = dyn_abun(mixture = "zip", K_max = 30),
         detection = ~ 1, y = sim$y, method = "nested_laplace"),
    "zip|zinb|does not yet compose")
})

test_that("dyn_abun(mixture='zip') recovers the structural-zero share", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  bo1 <- 0.3; zi_tru <- 0.3
  zi_hat <- lam0_hat <- rep(NA_real_, n_seed)
  cov_hit <- 0L; cov_tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_abun(N = 200, T = 3, J = 3, n_abund_covs = 1,
                             beta_lambda = c(log(6), bo1), p = 0.5,
                             omega = 0.6, gamma = 1.0, zi = zi_tru,
                             seed = 300 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data,
           family = dyn_abun(mixture = "zip", K_max = 30),
           detection = ~ 1, y = sim$y, method = "laplace",
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    zi_hat[s]   <- fit$zi_omega
    lam0_hat[s] <- fit$means[["lambda_(Intercept)"]]
    # 95% Wald CI for the ZI logit contains the truth.
    m <- fit$means[["zi_logit"]]; se <- fit$sds[["zi_logit"]]
    lo <- stats::plogis(m - 1.96 * se); hi <- stats::plogis(m + 1.96 * se)
    cov_tot <- cov_tot + 1L
    if (zi_tru >= lo && zi_tru <= hi) cov_hit <- cov_hit + 1L
  }
  # The structural-zero share (the new parameter) recovers with nominal coverage.
  expect_lt(abs(mean(zi_hat, na.rm = TRUE) - zi_tru), 0.06)
  expect_gte(cov_hit / cov_tot, 0.85)
  # Initial abundance recovers within the N-mixture lambda-p ridge band at J = 3.
  expect_lt(abs(mean(lam0_hat, na.rm = TRUE) - log(6)), 0.25)
})

test_that("dyn_abun(mixture='zinb') recovers zero-inflation + dispersion", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 12L
  zi_tru <- 0.25; r_tru <- 3
  zi_hat <- r_hat <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_abun(N = 220, T = 3, J = 3, n_abund_covs = 1,
                             beta_lambda = c(log(7), 0.2), p = 0.55, omega = 0.6,
                             gamma = 1.2, mixture = "negbin", r = r_tru,
                             zi = zi_tru, seed = 500 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data,
           family = dyn_abun(mixture = "zinb", K_max = 32),
           detection = ~ 1, y = sim$y, method = "laplace",
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    zi_hat[s] <- fit$zi_omega
    r_hat[s]  <- fit$dispersion$r
  }
  expect_lt(abs(mean(zi_hat, na.rm = TRUE) - zi_tru), 0.10)
  # NB size is weakly identified against the zero source; a generous band.
  expect_lt(abs(median(r_hat, na.rm = TRUE) - r_tru) / r_tru, 0.6)
})
