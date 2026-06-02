# =============================================================================
# test-ms-int-occu.R - community (multispecies) integrated occupancy
# (ms_int_occu()).
#
# Multiple detection sources share one latent occupancy state per species; per-
# species occupancy + per-source detection coefficient RE with per-arm Gaussian
# community covariances. The latent state marginalises out per species-site (a
# multi-source two-state mixture, analytic gradient) and the per-species
# deviations are integrated by the shared community Laplace-EM (R/community_em.R).
# status = "experimental" -> coverage floor >= 0.80.
# =============================================================================


test_that("ms_int_occu() constructor returns a tobs_family", {
  f <- ms_int_occu()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_int_occu")
  expect_equal(f$status, "experimental")
  expect_equal(f$default_engine, "laplace")
})


test_that("ms_int_occu() recovers community means + per-species coefs", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_int_ms_occu(N = 140, J = c(3, 4), n_species = 14,
                              n_data = 2, seed = 23)
  fit <- tobs(~ 1, data = sim$data, family = ms_int_occu(),
              detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(14)),
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means -> 0 (psi_species ~ N(0, 0.5); p_det[[d]] ~ N(0, 0.3)).
  truth <- c("psi_(Intercept)" = 0, "p1_(Intercept)" = 0, "p2_(Intercept)" = 0)
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species coefficients track the simulated truth on every arm.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi[, 1], stats::qlogis(sim$truth$psi_species)), 0.70)
  expect_gt(cor(cm$coef_p1[, 1],  stats::qlogis(sim$truth$p_det[[1]])),  0.60)
  expect_gt(cor(cm$coef_p2[, 1],  stats::qlogis(sim$truth$p_det[[2]])),  0.60)

  expect_true(all(cm$sd_psi > 0 & cm$sd_p1 > 0 & cm$sd_p2 > 0))
})


test_that("ms_int_occu() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 12L
  covered <- logical(0)
  truth <- c("psi_(Intercept)" = 0, "p1_(Intercept)" = 0, "p2_(Intercept)" = 0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_int_ms_occu(N = 130, J = c(3, 4), n_species = 12,
                                n_data = 2, seed = 600 + s)
    fit <- tryCatch(
      tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
           y = sim$y, species = paste0("sp", seq_len(12)),
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  expect_gt(mean(covered), 0.80)
})


test_that("ms_int_occu() S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_int_ms_occu(N = 80, J = c(3, 4), n_species = 8,
                              n_data = 2, seed = 5)
  fit <- tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(8)),
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("psi", "p1", "p2"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  # n_species x (psi + p1 + p2) = 8 x 3
  expect_equal(nrow(re), 8L * 3L)
  expect_setequal(unique(re$arm), c("psi", "p1", "p2"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi), c(80L, 8L))
  expect_equal(length(fv$p), 2L)
  expect_equal(dim(fv$p$p1), c(80L, 8L))

  ys <- simulate(fit, nsim = 1)
  expect_length(ys, 2L)
  expect_equal(dim(ys[[1]]), c(80L, 3L, 8L))
  expect_equal(dim(ys[[2]]), c(80L, 4L, 8L))

  expect_type(nobs(fit), "integer")
})


test_that("ms_int_occu() capability gates", {
  sim <- simulate_int_ms_occu(N = 40, J = c(3, 3), n_species = 4,
                              n_data = 2, seed = 1)
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nested_laplace"),
    "not available")
  # Missing y.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
         species = paste0("sp", seq_len(4)), method = "laplace"),
    "y")
  # Partial site overlap is not supported (full overlap only).
  bad <- sim$y
  bad[[2]] <- bad[[2]][1:30, , , drop = FALSE]
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
         y = bad, species = paste0("sp", seq_len(4)), method = "laplace"),
    "overlap")
})
