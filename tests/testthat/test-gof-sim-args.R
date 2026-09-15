# =============================================================================
# test-gof-sim-args.R -- the goodness-of-fit doors honour the arguments their
# generics are called with.
#
# test_dispersion / test_outliers / test_zero_inflation / check_model /
# pit_residuals are tulpa generics taking `...`, and tulpa's default methods
# spell the budget `nsim` and take `seed`. A tobs method that only knew
# `n.samples` dropped both into `...`: the budget ran at its default and the
# seed reached nothing.
# =============================================================================

ctl <- list(verbose = FALSE, progress = FALSE)

test_that("the argument resolver takes either budget name and refuses the rest", {
  expect_identical(.tobs_sim_args("f()", 250, NULL, NULL, list()), 250L)
  expect_identical(.tobs_sim_args("f()", 250, 12, NULL, list()), 12L)
  expect_error(.tobs_sim_args("f()", 250, NULL, NULL, list(nsmi = 10)),
               "unused argument `nsmi`")
  expect_error(.tobs_sim_args("f()", 250, NULL, NULL, list(1, a = 2)),
               "unused arguments `<unnamed>`, `a`")
  expect_error(.tobs_sim_args("f()", 250, NULL, c(0, 1), list()),
               "`observed` is not used")
  expect_error(.tobs_sim_args("f()", 0, NULL, NULL, list()),
               "one positive number")
})

test_that("a seed scopes the RNG and restores the caller's stream", {
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  a <- .tobs_with_seed(5L, stats::runif(3))
  after <- stats::runif(1)
  expect_identical(before, after)
  expect_identical(a, .tobs_with_seed(5L, stats::runif(3)))
  expect_false(identical(a, .tobs_with_seed(6L, stats::runif(3))))
})

test_that("the dispersion p-value follows the requested alternative", {
  sim <- c(1, 2, 3, 4)
  expect_identical(.tobs_sim_p_value(sim, 3, "greater"), 0.5)
  expect_identical(.tobs_sim_p_value(sim, 3, "less"), 0.75)
  expect_identical(.tobs_sim_p_value(sim, 3, "two.sided"), 1)
  expect_identical(.tobs_sim_p_value(sim, 4, "two.sided"), 0.5)
})

test_that("the GOF tests honour nsim, seed, and refuse unknown arguments", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_occu(N = 60L, J = 4L, seed = 11L)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, method = "laplace",
              control = ctl)

  expect_length(test_dispersion(fit, nsim = 10L)$sim, 10L)
  expect_length(test_dispersion(fit, n.samples = 12L)$sim, 12L)

  expect_identical(test_dispersion(fit, nsim = 30L, seed = 1L),
                   test_dispersion(fit, nsim = 30L, seed = 1L))
  expect_identical(test_outliers(fit, nsim = 30L, seed = 1L),
                   test_outliers(fit, nsim = 30L, seed = 1L))
  expect_identical(test_zero_inflation(fit, nsim = 30L, seed = 1L),
                   test_zero_inflation(fit, nsim = 30L, seed = 1L))

  two <- test_dispersion(fit, nsim = 30L, seed = 1L, alternative = "two.sided")
  expect_identical(two$alternative, "two.sided")
  expect_identical(two$p.value,
                   .tobs_sim_p_value(two$sim, two$observed, "two.sided"))

  expect_error(test_zero_inflation(fit, nsmi = 10L), "unused argument `nsmi`")
  expect_error(test_outliers(fit, observed = sim$y), "`observed` is not used")
  expect_error(check_model(fit, plot = FALSE, nsmi = 10L),
               "unused argument `nsmi`")

  r1 <- utils::capture.output(c1 <- check_model(fit, nsim = 20L, seed = 3L,
                                                plot = FALSE))
  r2 <- utils::capture.output(c2 <- check_model(fit, nsim = 20L, seed = 3L,
                                                plot = FALSE))
  expect_identical(c1$dispersion, c2$dispersion)
  expect_length(c1$dispersion$sim, 20L)
})
