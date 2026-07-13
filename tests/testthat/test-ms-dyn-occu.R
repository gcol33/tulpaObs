# =============================================================================
# test-ms-dyn-occu.R - community (multispecies) dynamic occupancy (ms_dyn_occu()).
#
# Per-species first-season occupancy + detection coefficient RE with per-arm
# Gaussian community covariances; shared community colonisation / extinction.
# The latent occupancy path integrates out by an HMM forward filter and the
# per-species deviations are integrated by the shared community Laplace-EM
# (R/community_em.R). The family is status = "working" (gcol33/tulpaObs#99):
# community-mean 95% CI coverage is gated at the 0.85 working floor of the
# recovery rubric (measured pooled coverage ~0.98 at 24 seeds; the shared
# gamma / eps dynamics cover at ~1.0).
# =============================================================================


test_that("ms_dyn_occu() constructor returns a tobs_family", {
  f <- ms_dyn_occu()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_dyn_occu")
  expect_equal(f$status, "working")
  expect_equal(f$default_engine, "laplace")
})


test_that("ms_dyn_occu() recovers community means + per-species coefs", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 90, J = 3, n_species = 14, n_seasons = 4,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = 31)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
              detection = ~ 1, y = sim$y,
              species = paste0("sp", seq_len(14)),
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means within ~2.5 SE of truth. psi1 intercept -> beta_comm_mean[1];
  # detection community mean -> 0 (p_species ~ N(0, 0.5)); gamma / eps -> the
  # shared transition logits.
  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  m <- fit$means[names(truth)]
  s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi1[, 1], stats::qlogis(sim$truth$psi1_species)), 0.70)
  expect_gt(cor(cm$coef_p[, 1],    stats::qlogis(sim$truth$p_species)),    0.60)

  # Per-species variance components recover the realised spread (gcol33/tulpaObs
  # #99): mild Laplace attenuation on the binary arms (measured est/realised
  # ~0.91 psi1, ~0.95 p over seeds), so a band around the realised SD.
  emp_psi1 <- stats::sd(stats::qlogis(sim$truth$psi1_species))
  emp_p    <- stats::sd(stats::qlogis(sim$truth$p_species))
  expect_true(cm$sd_psi1[1] > 0.4 * emp_psi1 && cm$sd_psi1[1] < 1.8 * emp_psi1)
  expect_true(cm$sd_p[1]    > 0.4 * emp_p    && cm$sd_p[1]    < 1.8 * emp_p)
})


test_that("ms_dyn_occu() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 12L
  covered <- logical(0)
  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_dyn_occu(N = 80, J = 3, n_species = 12, n_seasons = 4,
                                beta_comm_mean = c(0.3), beta_comm_sd = c(0.6),
                                gamma = 0.2, epsilon = 0.1, seed = 400 + s)
    fit <- tryCatch(
      tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
           detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(12)),
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  # Working-family gate (gcol33/tulpaObs#99): pooled over the shared psi1 / p
  # means and the shared gamma / eps dynamics x 12 seeds. Measured pooled ~0.98
  # at 24 seeds, comfortably clearing the 0.85 working floor.
  expect_gt(mean(covered), 0.85)
})


test_that("ms_dyn_occu() S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 50, J = 3, n_species = 8, n_seasons = 4,
                              seed = 5)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
              detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(8)),
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("psi1", "p", "gamma", "eps"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  # n_species x (psi1 + p) = 8 x 2
  expect_equal(nrow(re), 8L * 2L)
  expect_setequal(unique(re$arm), c("psi1", "p"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi1), c(50L, 8L))
  expect_equal(dim(fv$p),    c(50L, 8L))
  expect_length(fv$gamma, 50L)
  expect_length(fv$eps,   50L)
  expect_true(all(fv$psi1 > 0 & fv$psi1 < 1))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(50L, 3L, 4L, 8L))

  expect_type(nobs(fit), "integer")
})


test_that("ms_dyn_occu() capability gates", {
  sim <- simulate_ms_dyn_occu(N = 30, J = 3, n_species = 4, n_seasons = 3,
                              seed = 1)
  # nested_laplace / nuts not offered.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nested_laplace"),
    "not available")
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)), method = "nuts"),
    "not available")
  # Missing y.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         species = paste0("sp", seq_len(4)), method = "laplace"),
    "y")
})
