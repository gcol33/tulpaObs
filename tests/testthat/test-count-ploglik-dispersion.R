# =============================================================================
# test-count-ploglik-dispersion.R -- the count families' pointwise
# log-likelihood reads its trailing coordinates by name and scores the
# mixture the fit maximised.
#
# build_nmix_fit() appends up to three blocks after the (lambda, p) coefficients
# -- the structural-zero logit, the random-effect variance / BLUP block, and
# log_r -- each conditionally, so the coordinate one past the coefficients is
# log_r only when the other two are absent. dyn_abun estimates log_r as a named
# column and carries its own zero-inflation layer. The reference each block
# asserts against is the marginal the fit itself maximised: scored at the mode,
# the pointwise log-likelihood must sum to `fit$log_lik`, and scoring the wrong
# dispersion / dropping the mixture must move that total materially.
# =============================================================================

# The draw matrix at the fitted mode, in the layout the kernels read.
at_mode <- function(fit) {
  matrix(fit$means[colnames(fit$draws)], nrow = 1L,
         dimnames = list(NULL, colnames(fit$draws)))
}
ctl <- list(verbose = FALSE, progress = FALSE)

test_that("abun(mixture = 'zip') scores the structural-zero mixture", {
  skip_on_cran()
  sim <- simulate_abun(N = 120, J = 4, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(6), 0.5), beta_p = c(0.3, -0.3),
                       mixture = "zip", omega = 0.35, seed = 2)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = abun(mixture = "zip"), method = "laplace", control = ctl)
  th <- at_mode(fit)
  expect_true("logit_omega" %in% colnames(th))

  # The pointwise kernel reproduces the marginal the ZIP fitter maximised.
  expect_lt(abs(sum(.tobs_ploglik_from_draws(fit$model, th)) - fit$log_lik), 1e-6)

  # Discriminates: without the structural-zero layer (the former behaviour, which
  # also read logit_omega positionally as log_r) the total is far away.
  no_zi <- th[, setdiff(colnames(th), "logit_omega"), drop = FALSE]
  expect_gt(abs(sum(.tobs_ploglik_from_draws(fit$model, no_zi)) - fit$log_lik), 1)

  expect_true(is.finite(waic(fit)$elpd_waic))
})

test_that("abun() + (1 | g) scores a finite pointwise log-likelihood", {
  skip_on_cran()
  # The RE variance components carry deliberately NA draws, so reading the
  # coordinate past the coefficients as log_r gave exp(NA) for every site.
  set.seed(7); N <- 90L; ngrp <- 6L
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = 0.6)
  dat <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  lam <- exp(as.numeric(stats::model.matrix(~ x1, dat) %*% c(log(5), 0.4)) + b[grp])
  Nlat <- stats::rpois(N, lam)
  y <- matrix(NA_integer_, N, 4L)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(4L, Nlat[i], stats::plogis(0.5))

  fit <- tobs(~ x1 + (1 | g), data = dat, detection = ~ 1, y = y,
              family = abun(), method = "laplace", control = ctl)
  expect_true(any(is.na(fit$draws[, grep("^sigma_", colnames(fit$draws))])))
  expect_true(all(is.finite(.tobs_ploglik_from_draws(fit$model, fit$draws))))
  expect_true(is.finite(waic(fit)$elpd_waic))
})

test_that("an areal negbin abun() is scored as negbin, not Poisson", {
  skip_on_cran()
  # The areal fit integrates the size over the outer grid and carries no log_r
  # draw coordinate, so the size travels to the kernels on the model.
  chain <- function(n) { A <- matrix(0L, n, n)
    for (i in seq_len(n - 1L)) { A[i, i + 1L] <- 1L; A[i + 1L, i] <- 1L }; A }
  set.seed(11); n <- 40L
  dat <- data.frame(x1 = stats::rnorm(n))
  lam <- exp(log(7) + 0.4 * dat$x1)
  Nlat <- stats::rnbinom(n, size = 3, mu = lam)
  y <- matrix(NA_integer_, n, 4L)
  for (i in seq_len(n)) y[i, ] <- stats::rbinom(4L, Nlat[i], stats::plogis(0.6))

  fit <- tobs(~ x1 + icar(graph = chain(n)), data = dat, detection = ~ 1, y = y,
              family = abun(mixture = "negbin"), method = "nested_laplace",
              control = ctl)
  expect_false("log_r" %in% colnames(fit$draws))
  model <- .tobs_model_with_nb_size(fit)
  expect_equal(model$nb_r, as.numeric(fit$nmix_dispersion$r))
  expect_true(all(is.finite(.tobs_count_nb_size(model, fit$draws))))

  # Discriminates: dropping the size falls back to the Poisson marginal.
  th <- at_mode(fit)
  pois <- model; pois$nb_r <- NULL
  expect_gt(abs(sum(.tobs_ploglik_from_draws(model, th)) -
                sum(.tobs_ploglik_from_draws(pois, th))), 1)
})

test_that("dyn_abun(mixture = 'negbin') is scored at the estimated log_r", {
  skip_on_cran()
  sim <- simulate_dyn_abun(N = 50, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(8), 0), p = 0.6, omega = 0.7,
                           gamma = 1.2, mixture = "negbin", r = 3, seed = 4)
  fit <- tobs(~ 1, data = sim$data, detection = ~ 1, y = sim$y,
              family = dyn_abun(mixture = "negbin"), method = "laplace",
              control = ctl)
  th <- at_mode(fit)
  expect_true("log_r" %in% colnames(th))
  expect_lt(abs(sum(.tobs_ploglik_from_draws(fit$model, th)) - fit$log_lik), 1e-6)

  # Discriminates: the kernel used to be handed eta_logr = 0, i.e. size r = 1.
  at_r1 <- th; at_r1[1L, "log_r"] <- 0
  expect_gt(abs(sum(.tobs_ploglik_from_draws(fit$model, at_r1)) - fit$log_lik), 1)
})

test_that("dyn_abun(mixture = 'zinb') scores both the size and the ZI layer", {
  skip_on_cran()
  skip_if_fast()   # 68 s: the ZINB open-population forward is the file's one
                   # block too slow for the smoke tier (measured 2026-08-22).
  sim <- simulate_dyn_abun(N = 60, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(8), 0), p = 0.6, omega = 0.7,
                           gamma = 1.2, mixture = "negbin", r = 4, zi = 0.3,
                           seed = 5)
  fit <- tobs(~ 1, data = sim$data, detection = ~ 1, y = sim$y,
              family = dyn_abun(mixture = "zinb"), method = "laplace",
              control = ctl)
  th <- at_mode(fit)
  # zi_logit is the structural-zero coordinate; omega_* is the SURVIVAL arm.
  expect_true(all(c("zi_logit", "log_r") %in% colnames(th)))
  expect_lt(abs(sum(.tobs_ploglik_from_draws(fit$model, th)) - fit$log_lik), 1e-6)

  # Discriminates on each layer separately: "zinb" used to miss `use_nb`, and the
  # structural-zero component was absent from the kernel entirely.
  no_zi <- th[, setdiff(colnames(th), "zi_logit"), drop = FALSE]
  expect_gt(abs(sum(.tobs_ploglik_from_draws(fit$model, no_zi)) - fit$log_lik), 1)
  at_r1 <- th; at_r1[1L, "log_r"] <- 0
  expect_gt(abs(sum(.tobs_ploglik_from_draws(fit$model, at_r1)) - fit$log_lik), 1)

  expect_true(is.finite(waic(fit)$elpd_waic))
})

test_that("the posterior-mean row carries the draw names the kernels read", {
  skip_on_cran()
  # .tobs_loglik_at_mean() feeds DIC; an unnamed mean row would score the plain
  # Poisson marginal for a negbin fit.
  sim <- simulate_abun(N = 100, J = 4, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(6), 0.4), beta_p = c(0.3, -0.2),
                       mixture = "negbin", size = 3, seed = 9)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = abun(mixture = "negbin"), method = "laplace",
              control = ctl)
  mean_row <- matrix(colMeans(fit$draws), nrow = 1L,
                     dimnames = list(NULL, colnames(fit$draws)))
  ref <- as.numeric(.tobs_ploglik_from_draws(fit$model, mean_row))
  expect_equal(.tobs_loglik_at_mean(fit), ref)

  at_r1 <- mean_row; at_r1[1L, "log_r"] <- 0
  expect_gt(abs(sum(ref) - sum(.tobs_ploglik_from_draws(fit$model, at_r1))), 1)
  expect_true(is.finite(dic(fit)$dic))
})
