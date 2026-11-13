# Parameter-recovery + reproducibility tests for the stochastic Laplace
# correction routes. method = "laplace_gibbs" runs tulpa's exact
# MH-within-Gibbs chain over p(theta, Z | y) from the EM fit and reports that
# chain; method = "laplace_mi" runs Rubin-pooled complete-data refits. Both
# carry the same weakly-informative fixed-effect prior as method = "laplace".
# tobs() seeds the R-side latent draws so both reproduce.

sim_occu_fixed <- function(seed = 41, N = 400L, J = 5L,
                           beta_occ = c(0.4, -0.8), beta_det = c(0.0, 0.4)) {
  simulate_occu(N = N, J = J, n.occ.covs = 1L, n.det.covs = 1L,
                beta.occ = beta_occ, beta.det = beta_det, seed = seed)
}

test_that("method = 'laplace_gibbs' recovers occupancy/detection fixed effects", {
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_gibbs",
              control = list(seed = 123, verbose = FALSE))

  expect_identical(fit$method, "laplace_gibbs")
  # The Gibbs route now carries the default weakly-informative prior, threaded
  # through the correction refits.
  expect_s3_class(fit$priors, "occu_priors")
  # The seed used for the stochastic correction is recorded for reproducibility.
  expect_identical(fit$seed, 123L)

  cf_psi <- coef(fit, arm = "psi")
  cf_p   <- coef(fit, arm = "p")
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
  cf_psi <- coef(fit, arm = "psi")
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
  cf_psi <- coef(fit, arm = "psi")
  expect_lt(abs(cf_psi[["occ_cov1"]] - s$truth$beta_occ[2]), 0.35)
})

test_that("method = 'laplace_gibbs' reports the exact chain, not a Gaussian", {
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  fit <- tobs(~ occ_cov1, data = s$data, y = s$y, detection = ~ det_cov1,
              family = occu(), method = "laplace_gibbs",
              control = list(seed = 3, n.gibbs = 400L, verbose = FALSE))

  # 400 sweeps, half warmup: 200 kept rows, each paired with its latent draw.
  expect_identical(nrow(fit$draws), 200L)
  expect_identical(fit$draws_kind, "chain")
  expect_identical(fit$inference_tier, 1L)
  expect_length(fit$latent_draws, 200L)
  expect_identical(fit$gibbs_control$warmup, 200L)

  # The reported moments are the chain's own, and the chain is diagnosed.
  expect_equal(unname(fit$means), unname(colMeans(fit$draws)))
  expect_equal(unname(fit$sds), unname(apply(fit$draws, 2L, stats::sd)))
  expect_true(all(is.finite(fit$convergence$rhat)))
  expect_setequal(fit$convergence$parameter, names(fit$means))
  expect_gt(fit$mean_accept, 0.2)
})

test_that("the fixed-effect prior enters the Gibbs chain's target", {
  # A prior pinned at 0 with SD 0.02 on the detection slope (data SE ~0.065)
  # must hold the chain near 0: the precision-weighted mean is ~0.03. Under
  # the default prior the slope sits near its truth 0.4.
  skip_if_fast()
  s <- sim_occu_fixed(seed = 41)
  f <- function(pr) tobs(~ occ_cov1, data = s$data, y = s$y,
                         detection = ~ det_cov1, family = occu(),
                         method = "laplace_gibbs", priors = pr,
                         control = list(seed = 5, verbose = FALSE))
  pinned <- f(occu_priors(p.slope = list(mean = 0, sd = 0.02)))
  free   <- f(NULL)
  expect_lt(abs(pinned$means[["p_det_cov1"]]), 0.06)
  expect_gt(free$means[["p_det_cov1"]], 0.2)
})

test_that("method = 'laplace_gibbs' refuses a spatial term", {
  sim <- simulate_occu(N = 40L, J = 3L, seed = 1L)
  sim$data$x <- stats::runif(40L); sim$data$yc <- stats::runif(40L)
  expect_error(
    tobs(~ spde(x, yc), data = sim$data, family = occu(), detection = ~ 1,
         y = sim$y, method = "laplace_gibbs", control = list(verbose = FALSE)),
    "non-spatial models only")
})

# The corrections hand the M-step a hard latent draw. That draw is complete data,
# so each refit's Hessian is the complete-data information Rubin's rules pool as
# within-imputation variance. Encoding it through the EM's M = 1000
# pseudo-binomial shrinks the state arm's within-variance 1000-fold, leaving the
# pooled SE almost pure between-draw variance, about a third of the exact-marginal
# SE on the occu fixture below. The reference is the plain Laplace route, whose
# SEs come from the exact-marginal Hessian. The mean tolerance is in units of that
# SE, and the band on the SE ratio leaves room for the Monte-Carlo noise of the
# default draw count.
.gibbs_vs_laplace <- function(fits, info) {
  ref <- fits$laplace
  for (m in c("laplace_gibbs", "laplace_mi")) {
    f <- fits[[m]]
    expect_identical(names(f$means), names(ref$means), info = paste(info, m))
    ratio <- f$sds / ref$sds
    expect_true(all(ratio > 0.7 & ratio < 1.4),
                info = sprintf("%s %s SE ratio: %s", info, m,
                               paste(round(ratio, 2), collapse = " ")))
    shift <- abs(f$means - ref$means) / ref$sds
    expect_true(all(shift < 0.75),
                info = sprintf("%s %s mean shift / SE: %s", info, m,
                               paste(round(shift, 2), collapse = " ")))
  }
}

.fit_three_routes <- function(fit_one) {
  ms <- c("laplace", "laplace_gibbs", "laplace_mi")
  stats::setNames(lapply(ms, function(m) fit_one(m, if (m == "laplace")
    list(verbose = FALSE) else list(verbose = FALSE, seed = 7L))), ms)
}

test_that("the corrections pool the exact-marginal SEs on single-season occu", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_occu(N = 150L, J = 3L, seed = 4L)
  fits <- .fit_three_routes(function(m, ctl)
    tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
         y = sim$y, method = m, control = ctl))
  .gibbs_vs_laplace(fits, "occu")
})

test_that("dyn_occu fits under laplace_gibbs / laplace_mi and pools calibrated SEs", {
  skip_if_fast()
  skip_on_cran()
  # A hard draw samples each site's whole occupancy path (forward-filter
  # backward-sample) and is encoded by counting its transitions. A site with no
  # interval starting in the origin state contributes no transition trial; a
  # padded one-trial row per such site pulls gamma about 2.4 SE low here.
  sim <- simulate_dyn_occu(N = 120L, J = 4L, n.seasons = 4L, seed = 3L)
  fits <- .fit_three_routes(function(m, ctl)
    tobs(~ 1, data = sim$data, family = dyn_occu(), detection = ~ 1,
         colonization = ~ 1, extinction = ~ 1, y = sim$y, method = m,
         control = ctl))
  .gibbs_vs_laplace(fits, "dyn_occu")

  sv <- simulate_dyn_occu(N = 150L, J = 3L, n.seasons = 5L,
                          beta.gamma = c(-1, 0.8), seed = 1L)
  fits <- .fit_three_routes(function(m, ctl)
    tobs(~ 1, data = sv$data, family = dyn_occu(), detection = ~ 1,
         colonization = ~ gamma_cov, extinction = ~ 1, y = sv$y, method = m,
         control = ctl))
  .gibbs_vs_laplace(fits, "dyn_occu season-varying")
})

test_that("the corrections pool the exact-marginal SEs on int_occu", {
  skip_if_fast()
  skip_on_cran()
  N <- 150L
  set.seed(22L)
  x_cov <- rnorm(N); det_cov <- rnorm(N)
  z <- rbinom(N, 1, plogis(0.2 + 0.7 * x_cov))
  mk <- function(J, p0) {
    p <- plogis(p0 + 0.4 * det_cov)
    y <- matrix(0L, N, J)
    for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
    y
  }
  dat <- data.frame(occ_cov = x_cov, det_cov = det_cov)
  yy <- list(src1 = mk(4L, -0.2), src2 = mk(3L, -0.5))
  fits <- .fit_three_routes(function(m, ctl)
    tobs(~ occ_cov, data = dat, family = int_occu(), detection = ~ det_cov,
         y = yy, method = m, control = ctl))
  .gibbs_vs_laplace(fits, "int_occu")
})
