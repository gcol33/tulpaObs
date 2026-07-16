# =============================================================================
# test-abun-zip.R - zero-inflated N-mixture (ZIP / ZINB) for abun(mixture =).
#
# Structural-zero mixture over the Royle marginal (gcol33/tulpaObs#116): a share
# `omega` of sites are structural zeros (N = 0 regardless of lambda), the rest
# Pois/NB(lambda). A pure-R additive layer over the shared per-site Royle pieces
# (nmix_site_marginal); intercept-only structural-zero probability, non-spatial
# laplace. These check the constructor + gates, the S3 surface, and multi-seed
# parameter recovery (betas + omega [+ r]) against simulated truth.
# =============================================================================

test_that("abun(mixture = 'zip' / 'zinb') constructor + gates", {
  expect_equal(abun(mixture = "zip")$params$mixture, "zip")
  expect_equal(abun(mixture = "zinb")$params$mixture, "zinb")
  # ZI path is non-spatial laplace only (v1); a structured term / other engine
  # errors with a pointer rather than silently dropping it.
  sim <- simulate_abun(N = 40, J = 3, mixture = "zip", omega = 0.3, seed = 1)
  expect_error(
    tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
         family = abun(mixture = "zip"), method = "nuts"),
    "laplace")
})

test_that("abun(mixture = 'zip') fits + exposes the structural-zero logit", {
  sim <- simulate_abun(N = 150, J = 5, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(6), 0.5), beta_p = c(0.3, -0.3),
                       mixture = "zip", omega = 0.35, seed = 2)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = abun(mixture = "zip", K_max = 60L), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$mixture, "zip")
  expect_true(isTRUE(fit$zero_inflated))
  expect_true(is.finite(fit$zi_omega) && fit$zi_omega > 0 && fit$zi_omega < 1)
  # The structural-zero logit is a model coordinate: in means / vcov / sds with a
  # name (surfaced by summary), like the negbin log_r -- not in the per-process
  # coef() list.
  expect_true("logit_omega" %in% names(fit$means))
  expect_true("logit_omega" %in% rownames(vcov(fit)))
  expect_true(is.finite(fit$sds[["logit_omega"]]))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("abun(mixture = 'zip') recovers betas + omega (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 15L
  bl0 <- log(6); bl1 <- 0.5; bp0 <- 0.3; bp1 <- -0.3; om_t <- 0.35
  li <- ls <- pi_ <- ps <- om <- rep(NA_real_, n_seed)
  hit <- tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_abun(N = 250, J = 5, n_abund_covs = 1, n_det_covs = 1,
                         beta_lambda = c(bl0, bl1), beta_p = c(bp0, bp1),
                         mixture = "zip", omega = om_t, seed = 300 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
           family = abun(mixture = "zip", K_max = 60L), method = "laplace",
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    m <- fit$means
    li[s] <- m[["lambda_(Intercept)"]]; ls[s] <- m[["lambda_abund_cov1"]]
    pi_[s] <- m[["p_(Intercept)"]];      ps[s] <- m[["p_det_cov1"]]
    om[s] <- stats::plogis(m[["logit_omega"]])
    # 95% Wald coverage on the abundance slope (truth bl1).
    tot <- tot + 1L
    if (abs(m[["lambda_abund_cov1"]] - bl1) <= 1.96 * fit$sds[["lambda_abund_cov1"]])
      hit <- hit + 1L
  }
  expect_lt(abs(mean(li, na.rm = TRUE) - bl0), 0.10)
  expect_lt(abs(mean(ls, na.rm = TRUE) - bl1), 0.08)
  expect_lt(abs(mean(pi_, na.rm = TRUE) - bp0), 0.10)
  expect_lt(abs(mean(ps, na.rm = TRUE) - bp1), 0.08)
  expect_lt(abs(mean(om, na.rm = TRUE) - om_t), 0.05)
  expect_gte(hit / tot, 0.85)
})

test_that("abun(mixture = 'zinb') recovers betas + omega + size (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  # ZINB has a known zero-source confounding (structural omega vs NB
  # overdispersion r); test in a well-identified regime -- higher counts + more
  # visits, so the count-distribution shape separates the two mechanisms.
  n_seed <- 12L
  bl0 <- log(10); bl1 <- 0.4; om_t <- 0.25; size_t <- 6
  li <- ls <- om <- rr <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_abun(N = 400, J = 8, n_abund_covs = 1, n_det_covs = 1,
                         beta_lambda = c(bl0, bl1), beta_p = c(0.5, -0.2),
                         mixture = "zinb", size = size_t, omega = om_t,
                         seed = 500 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
           family = abun(mixture = "zinb", K_max = 150L), method = "laplace",
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    m <- fit$means
    li[s] <- m[["lambda_(Intercept)"]]; ls[s] <- m[["lambda_abund_cov1"]]
    om[s] <- stats::plogis(m[["logit_omega"]]); rr[s] <- exp(m[["log_r"]])
  }
  expect_lt(abs(median(li, na.rm = TRUE) - bl0), 0.12)
  expect_lt(abs(median(ls, na.rm = TRUE) - bl1), 0.10)
  expect_lt(abs(median(om, na.rm = TRUE) - om_t), 0.07)
  # Size r is weakly identified (the zero-source ridge); a lenient bound on the
  # central tendency, not per-seed.
  expect_gt(median(rr, na.rm = TRUE), 3)
})
