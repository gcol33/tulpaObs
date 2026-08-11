# gcol33/tulpaObs#222: pit_residuals() on single-season occu() used to pin every
# detected site's PIT at 1 and every all-zero site's PIT at its CDF value plus a
# tiny 1/n_draws jitter, so test_uniformity() rejected on a correctly specified
# model. A single fixture's KS p-value is not evidence of calibration -- that is
# what let the degeneracy survive -- so this asserts the KS p-values THEMSELVES
# look roughly uniform across many independent fits on self-generated data.

test_that("pit_residuals() on occu() is not degenerate at 1 for detected sites", {
  sim <- simulate_occu(N = 300, J = 6, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 1.0), beta_det = c(0.7, 0.6), seed = 1)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
              y = sim$y, method = "laplace", control = list(verbose = FALSE))
  pit <- pit_residuals(fit, n.samples = 250)

  never_detected <- mean(rowSums(sim$y > 0, na.rm = TRUE) == 0)
  # Before the fix, ~all detected sites piled up at exactly 1.
  expect_lt(mean(pit > 0.999), never_detected + 0.05)
  expect_true(all(pit >= 0 & pit <= 1))
})

test_that("pit_residuals() KS p-values are not degenerate over seeds (calibration)", {
  skip_if_fast()
  n_seeds <- 20L
  pvals <- vapply(seq_len(n_seeds), function(seed) {
    sim <- simulate_occu(N = 150, J = 5, n_occ_covs = 1, n_det_covs = 1,
                         beta_occ = c(0.2, 0.8), beta_det = c(0.5, 0.5), seed = seed)
    fit <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
                y = sim$y, method = "laplace", control = list(verbose = FALSE))
    pit <- pit_residuals(fit, n.samples = 250)
    unname(test_uniformity(pit)$p.value)
  }, numeric(1))

  # A calibrated PIT check rejects at alpha only ~alpha of the time; the
  # degenerate kernel rejected essentially always (p < 1e-10 territory, driven
  # by the pile-up at exactly 1, not sampling variation).
  expect_lt(mean(pvals < 0.05), 0.40)
  expect_gt(min(pvals), 1e-4)
})
