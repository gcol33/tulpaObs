# =============================================================================
# test-ms-occu.R - community (multispecies) single-season occupancy (ms_occu()).
#
# Per-species occupancy + detection coefficient RE with INDEPENDENT per-arm
# Gaussian community covariances (the spOccupancy msPGOcc model), fit by the
# shared community Laplace-EM (R/community_em.R). Replaces the legacy generic-
# engine community path, which fit a pooled GLM on the Laplace route and forced
# one shared species intercept across psi and p on NUTS (gcol33/tulpaObs#30).
# status = "working" but community-variance recovery for binary detection
# carries the documented Laplace small-cluster attenuation, so per-species
# detection recovery is checked loosely.
# =============================================================================


test_that("ms_occu() constructor returns a tobs_family", {
  f <- ms_occu()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_occu")
  expect_equal(f$default_engine, "laplace")
})


test_that("ms_occu() recovers community means + per-species coefs", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 130, J = 4, n_species = 16,
                          beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                          alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                          seed = 41)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(),
              detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(16)),
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means within ~2.5 SE of truth.
  truth <- c("psi_(Intercept)" = 0, "psi_x" = 0.6, "p_(Intercept)" = 0.2)
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi[, 1], sim$truth$beta_species[, 1]), 0.60)  # occ int
  expect_gt(cor(cm$coef_psi[, 2], sim$truth$beta_species[, 2]), 0.50)  # occ slope
  expect_gt(cor(cm$coef_p[, 1],   sim$truth$alpha_species[, 1]), 0.45) # det int

  expect_true(all(cm$sd_psi > 0 & cm$sd_p > 0))
})


test_that("ms_occu() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  n_seed <- 12L
  covered <- logical(0)
  truth <- c("psi_(Intercept)" = 0, "psi_x" = 0.6, "p_(Intercept)" = 0.2)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_occu(N = 120, J = 4, n_species = 14,
                            beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                            alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                            seed = 700 + s)
    fit <- tryCatch(
      tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
           y = sim$y, species = paste0("sp", seq_len(14)),
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  expect_gt(mean(covered), 0.80)
})


test_that("ms_occu() S3 methods work, incl. richness", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 60, J = 3, n_species = 8,
                          beta_comm_mean = c(0, 0.5), alpha_comm_mean = c(0.2),
                          seed = 5)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(8)),
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_setequal(names(cf), c("psi", "p"))
  expect_setequal(names(cf$psi), c("(Intercept)", "x"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 8L * 3L)               # 8 species x (psi int + psi x + p int)
  expect_setequal(unique(re$arm), c("psi", "p"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi), c(60L, 8L))
  expect_equal(dim(fv$p),   c(60L, 8L))
  expect_equal(dim(fv$z),   c(60L, 8L))
  expect_true(all(fv$psi > 0 & fv$psi < 1))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(60L, 3L, 8L))

  rich <- tobs_richness(fit)
  expect_s3_class(rich, "data.frame")
  expect_equal(nrow(rich), 60L)
  expect_true(all(rich$mean >= 0 & rich$mean <= 8))

  expect_type(nobs(fit), "integer")
})


test_that("ms_occu() capability gates: laplace only, species required", {
  sim <- simulate_ms_occu(N = 30, J = 3, n_species = 4, seed = 1)
  # nested_laplace / nuts no longer offered (community NUTS needs independent
  # per-arm RE blocks in the sampler; gcol33/tulpaObs#30).
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)), method = "nuts"),
    "not available")
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nested_laplace"),
    "not available")
})
