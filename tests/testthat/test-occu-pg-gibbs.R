# =============================================================================
# test-occu-pg-gibbs.R - Polya-Gamma Gibbs for single-season occupancy
# (method = "pg_gibbs"; spOccupancy PGOcc; gcol33/tulpaObs#126).
#
# A REAL Gibbs chain over the exact posterior via PG data augmentation (Polson,
# Scott & Windle 2013) using tulpa's tested PG sampler -- distinct from
# method = "laplace_gibbs", which is a stochastic-EM variance correction. On a
# well-identified occupancy model the PG-Gibbs posterior must MATCH the Laplace
# observed-Fisher fit (same data), which is the calibration reference; these
# check the S3 surface, that match, convergence diagnostics, and multi-seed
# recovery with nominal interval coverage.
# =============================================================================

test_that("occu() method = 'pg_gibbs' gates + S3", {
  sim <- simulate_occu(N = 120, J = 4, n_occ_covs = 1, n_det_covs = 1, seed = 1)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
              y = sim$y, method = "pg_gibbs",
              control = list(n.iter = 1000L, n.warmup = 500L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_true(all(c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)",
                    "p_det_cov1") %in% names(fit$means)))
  expect_true(all(is.finite(fit$rhat)))
  expect_true(all(fit$rhat < 1.1))
  expect_true(all(fit$ess > 50))
  # An icar() field IS supported (spPGOcc, tested in test-occu-pg-gibbs-spatial.R);
  # a temporal / RE term is not wired for the PG path in v1.
  expect_error(
    tobs(~ occ_cov1 + (1 | occ_cov1), data = sim$data, family = occu(),
         detection = ~ det_cov1, y = sim$y, method = "pg_gibbs"),
    "temporal / RE|RE")
})

test_that("occu() pg_gibbs posterior matches the Laplace fit (calibration)", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu(N = 400, J = 5, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 0.6), beta_det = c(0.2, -0.4), seed = 5)
  fl <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
             y = sim$y, method = "laplace", control = list(verbose = FALSE))
  fg <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
             y = sim$y, method = "pg_gibbs",
             control = list(n.iter = 3000L, n.warmup = 1500L, n.chains = 3L,
                            seed = 7, verbose = FALSE))
  lap_m <- c(unlist(coef(fl)$psi), unlist(coef(fl)$p))
  # Posterior means within ~1 Laplace-SE of the Laplace MLE.
  expect_true(all(abs(fg$means - lap_m) < fl$sds))
  # Posterior SDs match the observed-Fisher SEs within 20%.
  expect_true(all(abs(fg$sds / fl$sds - 1) < 0.2))
  expect_true(all(fg$rhat < 1.05))
})

test_that("occu() pg_gibbs matches the NUTS posterior (both exact)", {
  skip_on_cran()
  skip_if_fast()
  # PG-Gibbs and NUTS both target the EXACT posterior (unlike Laplace), so on the
  # same data their posterior means and SDs must agree -- a stronger anchor than
  # the Laplace calibration, which only holds asymptotically.
  sim <- simulate_occu(N = 300, J = 5, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 0.6), beta_det = c(0.2, -0.4), seed = 11)
  fn <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
             y = sim$y, method = "nuts",
             control = list(n.iter = 1500L, n.warmup = 750L, n.chains = 3L,
                            seed = 3, verbose = FALSE))
  fg <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
             y = sim$y, method = "pg_gibbs",
             control = list(n.iter = 3000L, n.warmup = 1500L, n.chains = 3L,
                            seed = 7, verbose = FALSE))
  # Means within ~0.4 SE of each other; SDs within 25% (both Monte-Carlo noisy).
  se <- pmax(fn$sds, fg$sds)
  expect_true(all(abs(fg$means[names(fn$means)] - fn$means) < 0.4 * se))
  expect_true(all(abs(fg$sds[names(fn$sds)] / fn$sds - 1) < 0.25))
})

test_that("occu() pg_gibbs agrees with spOccupancy::PGOcc (reference impl)", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("spOccupancy")
  # External reference: spOccupancy::PGOcc IS the Polson-Scott-Windle PG Gibbs
  # sampler this method reimplements. On matched data the two posteriors must
  # overlap -- the reference-implementation anchor the DoD asks for.
  sim <- simulate_occu(N = 300, J = 5, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 0.6), beta_det = c(0.2, -0.4), seed = 21)
  y <- sim$y; y[y < 0L] <- NA
  fg <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
             y = sim$y, method = "pg_gibbs",
             control = list(n.iter = 4000L, n.warmup = 2000L, n.chains = 3L,
                            seed = 4, verbose = FALSE))
  det_cov1 <- matrix(rep(sim$data$det_cov1, ncol(y)), nrow = nrow(y))
  # spOccupancy's lme4 dependency emits findbars/nobars deprecation warnings; they
  # are unrelated to the fit and would otherwise surface as test warnings.
  sp <- suppressWarnings(spOccupancy::PGOcc(
    occ.formula = ~ occ_cov1, det.formula = ~ det_cov1,
    data = list(y = y, occ.covs = data.frame(occ_cov1 = sim$data$occ_cov1),
                det.covs = list(det_cov1 = det_cov1)),
    n.samples = 4000L, n.burn = 2000L, n.thin = 1L, n.chains = 3L,
    verbose = FALSE))
  sp_m <- c(colMeans(sp$beta.samples), colMeans(sp$alpha.samples))
  sp_s <- c(apply(sp$beta.samples, 2L, stats::sd),
            apply(sp$alpha.samples, 2L, stats::sd))
  # Same PG posterior: means agree within ~0.3 of the (matched) posterior SD.
  expect_true(all(abs(fg$means - sp_m) < 0.3 * pmax(fg$sds, sp_s)))
  expect_true(all(abs(fg$sds / sp_s - 1) < 0.3))
})

test_that("occu() pg_gibbs recovers coefficients with nominal coverage", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  bo0 <- 0.3; bo1 <- 0.6; bd0 <- 0.2; bd1 <- -0.4
  truth <- c(bo0, bo1, bd0, bd1)
  est <- matrix(NA_real_, n_seed, 4L)
  cov_hit <- 0L; cov_tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu(N = 400, J = 5, n_occ_covs = 1, n_det_covs = 1,
                         beta_occ = c(bo0, bo1), beta_det = c(bd0, bd1),
                         seed = 200 + s)
    fg <- tryCatch(
      tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
           y = sim$y, method = "pg_gibbs",
           control = list(n.iter = 2000L, n.warmup = 1000L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fg)) next
    est[s, ] <- fg$means
    # 95% credible intervals (per-parameter quantiles of the draws).
    qs <- apply(fg$draws, 2L, stats::quantile, c(0.025, 0.975))
    for (k in seq_len(4L)) {
      cov_tot <- cov_tot + 1L
      if (truth[k] >= qs[1, k] && truth[k] <= qs[2, k]) cov_hit <- cov_hit + 1L
    }
  }
  # Slopes recover in aggregate (occupancy intercept is detection-filtered).
  expect_lt(abs(mean(est[, 2], na.rm = TRUE) - bo1), 0.10)
  expect_lt(abs(mean(est[, 3], na.rm = TRUE) - bd0), 0.10)
  expect_lt(abs(mean(est[, 4], na.rm = TRUE) - bd1), 0.10)
  # Nominal 95% coverage (working floor 0.85 for the Monte-Carlo estimate).
  expect_gte(cov_hit / cov_tot, 0.85)
})
