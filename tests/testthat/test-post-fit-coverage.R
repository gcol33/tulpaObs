# test-post-fit-coverage.R - link-aware marginal effects and community-wide
# species richness. tobs_marginal_effect() / predict_terms() previously
# hardcoded the logit link, silently returning logit-of-log-lambda on count
# fits; tobs_richness() was restricted to ms_occu().

test_that("tobs_marginal_effect returns the intensity scale on an abun fit", {
  skip_on_cran()
  sim <- simulate_abun(N = 150, J = 4, n_abund_covs = 2, n_det_covs = 1, seed = 7)
  fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data, family = abun(K_max = 60),
              detection = ~ det_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  me <- tobs_marginal_effect(fit, "abund_cov1", process = "abundance",
                             n_points = 20)
  # Intensity (lambda), not a probability: positive and able to exceed 1.
  expect_true(all(me$estimate > 0))
  expect_gt(max(me$estimate), 1)
  # Detection arm stays on the probability scale.
  med <- tobs_marginal_effect(fit, "det_cov1", process = "detection",
                              n_points = 20)
  expect_true(all(med$estimate >= 0 & med$estimate <= 1))
})

test_that("tobs_richness works on all three community occupancy families", {
  skip_on_cran()
  so <- simulate_ms_occu(N = 40, J = 3, n_species = 6, seed = 1)
  fo <- tobs(~ x, data = so$data, family = ms_occu(), detection = ~ 1, y = so$y,
             species = paste0("sp", 1:6), method = "laplace",
             control = list(verbose = FALSE))
  ro <- tobs_richness(fo)
  expect_equal(nrow(ro), 40L)
  expect_true(all(ro$mean >= 0 & ro$mean <= 6))

  skip_if_fast()
  si <- simulate_ms_int_occu(N = 50, J = c(3, 4), n_species = 6, seed = 3)
  fi <- tobs(~ 1, data = si$data, family = ms_int_occu(), detection = ~ 1,
             y = si$y, species = paste0("sp", 1:6), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(nrow(tobs_richness(fi)), 50L)

  sd_ <- simulate_ms_dyn_occu(N = 40, J = 3, n_species = 6, n_seasons = 3,
                              gamma = 0.2, epsilon = 0.1, seed = 2)
  fd <- tobs(~ 1, data = sd_$data, family = ms_dyn_occu(), detection = ~ 1,
             y = sd_$y, species = paste0("sp", 1:6), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(nrow(tobs_richness(fd)), 40L)
})
