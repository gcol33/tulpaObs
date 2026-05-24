# Parameter-recovery + reproducibility tests for the newly-reachable
# stochastic Laplace correction routes (method = "laplace_gibbs" /
# "laplace_mi"). These run the unpenalised EM with a post-EM Rubin-pooled
# correction (tulpa's tulpa_em_laplace); tobs() seeds the R-side draws so the
# pooled estimate reproduces. The fixed-effect priors are auto-disabled on
# these routes (gcol33/tulpa#27), so this is a clean parameter-recovery check.

sim_occu_fixed <- function(seed = 41, N = 400L, J = 5L,
                           beta_occ = c(0.4, -0.8), beta_det = c(0.0, 0.4)) {
  simulate_occu(N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L,
                beta_occ = beta_occ, beta_det = beta_det, seed = seed)
}

test_that("method = 'laplace_gibbs' recovers occupancy/detection fixed effects", {
  s <- sim_occu_fixed(seed = 41)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_gibbs",
              control = list(seed = 123, verbose = FALSE))

  expect_identical(fit$method, "laplace_gibbs")
  # The Gibbs route runs the unpenalised EM (no fixed-effect prior).
  expect_null(fit$priors)
  # The seed used for the stochastic correction is recorded for reproducibility.
  expect_identical(fit$seed, 123L)

  cf_psi <- coef(fit)$psi
  cf_p   <- coef(fit)$p
  expect_lt(abs(cf_psi[["(Intercept)"]] - s$truth$beta_occ[1]), 0.35)
  expect_lt(abs(cf_psi[["occ_cov1"]]    - s$truth$beta_occ[2]), 0.35)
  expect_lt(abs(cf_p[["(Intercept)"]]   - s$truth$beta_det[1]), 0.35)
  expect_lt(abs(cf_p[["det_cov1"]]      - s$truth$beta_det[2]), 0.35)
})

test_that("a fixed seed makes method = 'laplace_gibbs' reproducible", {
  s <- sim_occu_fixed(seed = 41)
  f <- function() tobs(~ occ_cov1, data = s$data, y = s$y,
                       detection = ~ det_cov1, family = occu(),
                       method = "laplace_gibbs",
                       control = list(seed = 7, verbose = FALSE))
  a <- f(); b <- f()
  expect_equal(a$means, b$means)
})

test_that("method = 'laplace_mi' recovers fixed effects and records its seed", {
  s <- sim_occu_fixed(seed = 42)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_mi",
              control = list(seed = 99, n.imputations = 25L, verbose = FALSE))

  expect_identical(fit$method, "laplace_mi")
  expect_identical(fit$seed, 99L)
  cf_psi <- coef(fit)$psi
  expect_lt(abs(cf_psi[["(Intercept)"]] - s$truth$beta_occ[1]), 0.35)
  expect_lt(abs(cf_psi[["occ_cov1"]]    - s$truth$beta_occ[2]), 0.35)
})
