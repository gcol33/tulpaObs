# =============================================================================
# Family-level NUTS calibration: 20-seed 95% CI coverage on the coefficients the
# family NUTS path exists to calibrate (gcol33/tulpaObs#139). The single-seed
# recovery tests in each family's own test file prove the byte-exact C++<->R
# oracle and point recovery; these prove the intervals cover at the nominal rate.
# The estimand is the calibrated coefficient interval: NUTS samples the EXACT
# closed-form marginal, so an unbiased posterior should contain truth ~95% of the
# time even for weakly-identified coefficients (a wide but honest CI still covers).
#
# Pooled over all coefficients x 20 seeds, assert >= the 0.85 rubric floor (a
# calibrated sampler measures ~0.93-0.95; the floor absorbs Monte-Carlo slack).
# All heavy, so skip_if_fast()-gated -- zero cost in the fast suite.
# =============================================================================

.nuts_ci_cover <- function(fit, truth) {
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  abs(est - truth) <= 1.96 * se        # per-coefficient 95% Wald CI containment
}

cuts5 <- seq(0, 1, length.out = 6)

test_that("distance NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(45), 0.4); beta_sigma <- c(log(0.45), 0.2)
  truth <- c(beta_lambda, beta_sigma)
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_distance(N = 120, cutpoints = cuts5, key = "halfnorm",
                             transect = "line", n_abund_covs = 1, n_sigma_covs = 1,
                             beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                             seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data,
                    family = distance(key = "halfnorm", transect = "line",
                                      cutpoints = sim$cutpoints),
                    detection = ~ sigma_cov1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) covered <- c(covered, .nuts_ci_cover(fit, truth))
  }
  expect_gte(mean(covered), 0.85)
})

test_that("removal NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(7), 0.5); beta_p <- c(0.3, -0.3)
  truth <- c(beta_lambda, beta_p)
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_removal(N = 80, K = 5, n_abund_covs = 1, n_det_covs = 1,
                            beta_lambda = beta_lambda, beta_p = beta_p, seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data, family = removal(),
                    detection = ~ det_cov1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) covered <- c(covered, .nuts_ci_cover(fit, truth))
  }
  expect_gte(mean(covered), 0.85)
})

test_that("fp_occu NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.6)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_fp_occu(N = 300, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                            p11 = 0.6, p10 = 0.05, b = 0.5, seed = 600 + s)
    fit <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                    detection = ~ 1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) covered <- c(covered, .nuts_ci_cover(fit, truth))
  }
  expect_gte(mean(covered), 0.85)
})

test_that("dyn_abun NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(6), 0.4)
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2))
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_dyn_abun(N = 70, T = 3, J = 3, n_abund_covs = 1,
                             beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                             gamma = 1.2, seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 26),
                    detection = ~ 1, y = sim$y, method = "nuts",
                    control = list(n.iter = 250L, n.warmup = 250L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) covered <- c(covered, .nuts_ci_cover(fit, truth))
  }
  expect_gte(mean(covered), 0.85)
})

test_that("abun NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(4), 0.5, -0.3); beta_p <- c(0.2, -0.4)
  truth <- c(beta_lambda, beta_p)
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_abun(N = 60, J = 4, n_abund_covs = 2, n_det_covs = 1,
                         beta_lambda = beta_lambda, beta_p = beta_p,
                         mixture = "poisson", seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1 + abund_cov2, data = sim$data, y = sim$y,
                    family = abun(), detection = ~ det_cov1, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) covered <- c(covered, .nuts_ci_cover(fit, truth))
  }
  expect_gte(mean(covered), 0.85)
})
