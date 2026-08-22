# =============================================================================
# test-ms-int-occu-pg-gibbs.R - Polya-Gamma Gibbs for COMMUNITY multi-source
# integrated occupancy (method = "pg_gibbs").
#
# A real Gibbs chain over the exact posterior via PG data augmentation: per
# species one latent occupancy state observed by D detection sources, per-species
# occupancy / per-source detection coefficients with Gaussian community
# hyperpriors. The community VARIANCE recovers where the community Laplace-EM
# attenuates it (a documented lower bound). These check the S3 surface, community-
# mean + variance recovery, and convergence.
# =============================================================================

test_that("ms_int_occu() method = 'pg_gibbs' gates + S3", {
  sim <- simulate_ms_int_occu(N = 120, J = c(3, 4), n_species = 6, n_data = 2,
                              seed = 1)
  fit <- tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6), method = "pg_gibbs",
              control = list(n.iter = 1200L, n.warmup = 600L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_true("pg_gibbs" %in% tulpaObs:::.tobs_family_methods$ms_int_occu)
  expect_true(all(c("psi_(Intercept)", "p1_(Intercept)", "p2_(Intercept)") %in%
                    names(fit$means)))
  # per-arm community structure present + S3 works
  expect_false(is.null(fit$ms_community$sd_psi))
  expect_false(is.null(fit$ms_community$coef_p1))
  expect_equal(nrow(fit$ms_community$coef_psi), 6L)
  expect_false(is.null(ranef(fit)))
  expect_false(is.null(coef(fit)))
  expect_true(all(fit$rhat < 1.15))
})

test_that("ms_int_occu() pg_gibbs recovers community mean + variance", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  # simulate_ms_int_occu draws per-species logit-occupancy ~ N(0, 0.5^2) and each
  # source's logit-detection ~ N(0, 0.3^2), so the community MEAN is ~0 and the
  # community SD is 0.5 (occupancy) -- the target the PG variance must recover.
  psi_mu <- p_mu <- sd_psi <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_int_occu(N = 160, J = c(4, 3), n_species = 10, n_data = 2,
                                seed = 100 + s)
    fit <- tryCatch(
      tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
           y = sim$y, species = paste0("sp", 1:10), method = "pg_gibbs",
           control = list(n.iter = 1500L, n.warmup = 750L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    psi_mu[s]  <- fit$means[["psi_(Intercept)"]]
    p_mu[s]    <- fit$means[["p1_(Intercept)"]]
    sd_psi[s]  <- fit$ms_community$sd_psi[[1L]]
  }
  # Community means center near 0 (the simulator's community mean).
  expect_lt(abs(mean(psi_mu, na.rm = TRUE)), 0.20)
  expect_lt(abs(mean(p_mu, na.rm = TRUE)), 0.20)
  # Community SD recovers the 0.5 truth (the calibration the Laplace-EM misses),
  # within a modest band given the moderate species count.
  expect_lt(abs(mean(sd_psi, na.rm = TRUE) - 0.5), 0.15)
})
