# Parameter-recovery + reproducibility tests for the stochastic Laplace
# correction routes (method = "laplace_gibbs" / "laplace_mi"). These run a
# post-EM Rubin-pooled correction on top of tulpa's EM+Laplace. Since
# gcol33/tulpa#27, the weakly-informative fixed-effect prior threads into the
# correction refits, so these routes apply the same default prior as
# method = "laplace" (pass priors = FALSE to recover the unpenalised
# correction). tobs() seeds the R-side hard-z draws so the pooled estimate
# reproduces.

sim_occu_fixed <- function(seed = 41, N = 400L, J = 5L,
                           beta_occ = c(0.4, -0.8), beta_det = c(0.0, 0.4)) {
  simulate_occu(N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L,
                beta_occ = beta_occ, beta_det = beta_det, seed = seed)
}

test_that("method = 'laplace_gibbs' recovers occupancy/detection fixed effects", {
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_gibbs",
              control = list(seed = 123, verbose = FALSE))

  expect_identical(fit$method, "laplace_gibbs")
  # The Gibbs route now carries the default weakly-informative prior, threaded
  # through the correction refits (gcol33/tulpa#27).
  expect_s3_class(fit$priors, "occu_priors")
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
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  f <- function() tobs(~ occ_cov1, data = s$data, y = s$y,
                       detection = ~ det_cov1, family = occu(),
                       method = "laplace_gibbs",
                       control = list(seed = 7, verbose = FALSE))
  a <- f(); b <- f()
  expect_equal(a$means, b$means)
})

test_that("method = 'laplace_mi' recovers fixed effects and records its seed", {
  skip_if_fast()
  s <- sim_occu_fixed(seed = 42)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_mi",
              control = list(seed = 99, n.imputations = 25L, verbose = FALSE))

  expect_identical(fit$method, "laplace_mi")
  expect_identical(fit$seed, 99L)
  expect_s3_class(fit$priors, "occu_priors")
  cf_psi <- coef(fit)$psi
  expect_lt(abs(cf_psi[["(Intercept)"]] - s$truth$beta_occ[1]), 0.35)
  expect_lt(abs(cf_psi[["occ_cov1"]]    - s$truth$beta_occ[2]), 0.35)
})

test_that("priors = FALSE recovers the unpenalised Gibbs correction", {
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_gibbs", priors = FALSE,
              control = list(seed = 123, verbose = FALSE))
  expect_null(fit$priors)
  cf_psi <- coef(fit)$psi
  expect_lt(abs(cf_psi[["occ_cov1"]] - s$truth$beta_occ[2]), 0.35)
})

test_that("the prior threads into the Gibbs correction at the small-J ridge", {
  # At small N, J the psi-p ridge biases the detection slope toward zero unless
  # the prior breaks it. The penalised Gibbs correction should track the truth
  # better than the unpenalised one -- evidence the prior actually flows into
  # the correction refits, not just the EM point estimate.
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 12L
  truth_slope <- 0.8
  N <- 200L; J <- 4L
  est_pen <- numeric(0); est_unp <- numeric(0)

  for (s in seq_len(n_seeds)) {
    set.seed(5000L + s)
    x_occ <- rnorm(N); x_det <- rnorm(N)
    psi   <- plogis(0.5 + 0.5 * x_occ)
    p_tru <- plogis(0.0 + truth_slope * x_det)
    z     <- rbinom(N, 1L, psi)
    y     <- matrix(0L, N, J)
    for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
    dat <- data.frame(x_occ = x_occ, x_det = x_det)

    fp <- tryCatch(tobs(~ x_occ, data = dat, family = occu(), detection = ~ x_det,
                        y = y, method = "laplace_gibbs",
                        control = list(seed = 11, verbose = FALSE)),
                   error = function(e) NULL)
    fu <- tryCatch(tobs(~ x_occ, data = dat, family = occu(), detection = ~ x_det,
                        y = y, method = "laplace_gibbs", priors = FALSE,
                        control = list(seed = 11, verbose = FALSE)),
                   error = function(e) NULL)
    if (!is.null(fp)) est_pen <- c(est_pen, unname(fp$means[["p_x_det"]]))
    if (!is.null(fu)) est_unp <- c(est_unp, unname(fu$means[["p_x_det"]]))
  }

  expect_true(length(est_pen) >= floor(0.8 * n_seeds))
  expect_true(length(est_unp) >= floor(0.8 * n_seeds))
  # Penalised correction bias must be materially smaller than unpenalised.
  expect_lt(abs(mean(est_pen) - truth_slope), abs(mean(est_unp) - truth_slope))
})
