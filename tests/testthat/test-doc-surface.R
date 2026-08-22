# =============================================================================
# test-doc-surface.R -- the three behaviour halves of the documentation audit
# (gcol33 issue on API.md and the diagnostics doors).
#
# Each was a promise the code did not keep: a documented return field that
# existed on one dispatch branch only, arguments the uniformity test forwards
# that fell into `...`, and a pre-fit check that reported a clean bill of health
# on every family it never ran on.
# =============================================================================

ctl <- list(verbose = FALSE, progress = FALSE)

test_that("test_outliers() returns the same field names on every family", {
  skip_if_fast()
  skip_on_cran()
  so <- simulate_occu(N = 50, J = 3, seed = 11)
  fo <- tobs(~ occ_cov1, data = so$data, family = occu(), detection = ~ 1,
             y = so$y, method = "laplace", control = ctl)
  sa <- simulate_abun(N = 50, J = 3, n_abund_covs = 1, n_det_covs = 1,
                      beta_lambda = c(log(5), 0.3), beta_p = c(0.4, -0.2),
                      seed = 11)
  fa <- tobs(~ abund_cov1, data = sa$data, family = abun(),
             detection = ~ det_cov1, y = sa$y, method = "laplace",
             control = ctl)
  ro <- test_outliers(fo, n.samples = 40L)
  ra <- test_outliers(fa, n.samples = 40L)
  # The occupancy branch used to return `n_outliers` and no `ratio`, so
  # `$observed` was NULL on one family and `$n_outliers` NULL on the other and
  # no caller could read the statistic by name across both.
  for (nm in c("observed", "expected", "ratio", "p.value")) {
    expect_false(is.null(ro[[nm]]), info = nm)
    expect_false(is.null(ra[[nm]]), info = nm)
  }
  expect_null(ro$n_outliers)
  expect_true(is.numeric(ro$observed) && length(ro$observed) == 1L)
})

test_that("pit_residuals() honours the arguments test_uniformity forwards", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_occu(N = 60, J = 3, seed = 12)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ 1,
              y = sim$y, method = "laplace", control = ctl)
  # tulpa::test_uniformity() calls pit_residuals(object, observed=, nsim=,
  # seed=). None is a prefix of `n.samples`, so all three used to fall into
  # `...`: nsim was ignored and the draw selection ran unseeded.
  expect_identical(pit_residuals(fit, seed = 5L), pit_residuals(fit, seed = 5L))
  expect_false(identical(pit_residuals(fit, seed = 5L),
                         pit_residuals(fit, seed = 6L)))
  # A seeded call must leave the caller's stream where it found it.
  set.seed(99); before <- stats::runif(1)
  set.seed(99); invisible(pit_residuals(fit, seed = 5L)); after <- stats::runif(1)
  expect_identical(before, after)
  # Same seed, two nsim: the larger request must actually draw more.
  k1 <- tulpa::test_uniformity(fit, nsim = 40L, seed = 3L)$statistic
  k2 <- tulpa::test_uniformity(fit, nsim = 40L, seed = 3L)$statistic
  expect_identical(k1, k2)
  expect_false(identical(k1, tulpa::test_uniformity(fit, nsim = 400L,
                                                    seed = 3L)$statistic))
})

test_that("tobs_check_id() says when its pre-fit checks did not run", {
  skip_if_fast()
  skip_on_cran()
  so <- simulate_occu(N = 50, J = 3, seed = 13)
  fo <- tobs(~ occ_cov1, data = so$data, family = occu(), detection = ~ 1,
             y = so$y, method = "laplace", control = ctl)
  expect_true(suppressMessages(tobs_check_id(fo$model))$prefit_checked)

  sa <- simulate_abun(N = 50, J = 3, n_abund_covs = 1, n_det_covs = 1,
                      beta_lambda = c(log(5), 0.3), beta_p = c(0.4, -0.2),
                      seed = 13)
  fa <- tobs(~ abund_cov1, data = sa$data, family = abun(),
             detection = ~ det_cov1, y = sa$y, method = "laplace",
             control = ctl)
  # The gate used to be `%in% c("single", "community")`, and no family assigns
  # "community", so every family but single-season occupancy skipped the checks
  # and still reported "No identifiability issues detected."
  expect_true(suppressMessages(tobs_check_id(fa$model))$prefit_checked)

  sc <- simulate_count(N = 80, beta = c(0.8, 0.5), response = "poisson",
                       seed = 13)
  fc <- tobs(~ x, data = sc$data, y = sc$y, family = count(),
             method = "laplace", control = ctl)
  r <- suppressMessages(tobs_check_id(fc$model))
  expect_false(r$prefit_checked)
  expect_match(paste(utils::capture.output(
                 tobs_check_id(fc$model), type = "message"), collapse = " "),
               "did not run")
  expect_false("community" %in% .TOBS_CHECK_ID_TYPES)
})
