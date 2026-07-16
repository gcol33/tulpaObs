# =============================================================================
# test-occu-ttd.R - time-to-detection occupancy (occu_ttd(); unmarked occuTTD).
#
# Two-state occupancy marginal with a censored-exponential time-to-detection
# emission (gcol33/tulpaObs#116): z ~ Bernoulli(psi), t | z=1 ~ Exp(rate lambda)
# censored at surveyLength. Latent z integrates out in closed form; the exact
# marginal is maximised (optim BFGS) with an observed-information vcov -- the
# royle_nichols() recipe with a continuous emission. Non-spatial laplace only.
# =============================================================================

test_that("occu_ttd() constructor + gates", {
  f <- occu_ttd(surveyLength = 3)
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "occu_ttd")
  expect_equal(f$params$surveyLength, 3)
  sim <- simulate_occu_ttd(N = 40, J = 3, Tmax = 2, seed = 1)
  # NUTS is not a supported engine (v1 laplace only).
  expect_error(
    tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 2),
         detection = ~ rate_cov1, y = sim$y, method = "nuts"),
    "laplace")
  # Visit-level rate covariates rejected (rate is site-level in v1).
  vd <- data.frame(z = rnorm(40 * 3))
  expect_error(
    tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 2),
         detection = ~ rate_cov1, y = sim$y, visits = vd),
    "site-level")
})

test_that("occu_ttd() fits + full S3 surface", {
  sim <- simulate_occu_ttd(N = 200, J = 5, beta_psi = c(qlogis(0.6), 0.7),
                           beta_rate = c(log(0.7), -0.4), Tmax = 3, seed = 3)
  fit <- tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 3),
              detection = ~ rate_cov1, y = sim$y,
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$convergence$converged))
  expect_true(all(c("psi_(Intercept)", "rate_(Intercept)") %in% names(fit$means)))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "rate", "p"))
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$rate > 0))
  expect_true(all(fv$p >= 0 & fv$p <= 1))
  expect_equal(predict(fit, type = "state"), fv$psi)
  expect_equal(predict(fit, type = "detection"), fv$rate)

  w <- tobs_waic(fit, n.draws = 200L)
  expect_true(is.finite(w$waic) && w$p_waic > 0)
  s2 <- simulate(fit, nsim = 1)
  expect_true(all(s2[!is.na(s2)] >= 0 & s2[!is.na(s2)] <= 3))   # in [0, Tmax]
  expect_length(residuals(fit)$occ, fit$model$n_sites)
})

test_that("occu_ttd() recovers psi + rate coefficients (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 15L
  bpsi0 <- qlogis(0.6); bpsi1 <- 0.7; brate0 <- log(0.7); brate1 <- -0.4
  pi_ <- ps <- ri <- rs <- rep(NA_real_, n_seed)
  hit <- tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu_ttd(N = 300, J = 5, beta_psi = c(bpsi0, bpsi1),
                             beta_rate = c(brate0, brate1), Tmax = 3,
                             seed = 700 + s)
    fit <- tryCatch(
      tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 3),
           detection = ~ rate_cov1, y = sim$y,
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    m <- fit$means
    pi_[s] <- m[["psi_(Intercept)"]]; ps[s] <- m[["psi_psi_cov1"]]
    ri[s]  <- m[["rate_(Intercept)"]]; rs[s] <- m[["rate_rate_cov1"]]
    tot <- tot + 1L
    if (abs(m[["rate_rate_cov1"]] - brate1) <= 1.96 * fit$sds[["rate_rate_cov1"]])
      hit <- hit + 1L
  }
  expect_lt(abs(mean(pi_, na.rm = TRUE) - bpsi0),  0.12)
  expect_lt(abs(mean(ps,  na.rm = TRUE) - bpsi1),  0.10)
  expect_lt(abs(mean(ri,  na.rm = TRUE) - brate0), 0.08)
  expect_lt(abs(mean(rs,  na.rm = TRUE) - brate1), 0.06)
  expect_gte(hit / tot, 0.85)
})
