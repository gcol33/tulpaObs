# Parameter-recovery + reproducibility tests for the stochastic Laplace
# correction routes (method = "laplace_gibbs" / "laplace_mi"). These run a
# post-EM Rubin-pooled correction on top of tulpa's EM+Laplace. Since, the
# weakly-informative fixed-effect prior threads into the correction refits,
# so these routes apply the same default prior as method = "laplace" (pass
# priors = FALSE to recover the unpenalised correction). tobs() seeds the
# R-side hard-z draws so the pooled estimate reproduces.

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
  # through the correction refits.
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
  sim <- simulate_dyn_occu(N = 120L, J = 4L, n_seasons = 4L, seed = 3L)
  fits <- .fit_three_routes(function(m, ctl)
    tobs(~ 1, data = sim$data, family = dyn_occu(), detection = ~ 1,
         colonization = ~ 1, extinction = ~ 1, y = sim$y, method = m,
         control = ctl))
  .gibbs_vs_laplace(fits, "dyn_occu")

  sv <- simulate_dyn_occu(N = 150L, J = 3L, n_seasons = 5L,
                          beta_gamma = c(-1, 0.8), seed = 1L)
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
