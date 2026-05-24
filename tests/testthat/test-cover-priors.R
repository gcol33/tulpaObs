# cover_priors(): opt-in fixed-effect priors for the cover hurdle.
# Verifies construction/validation, that the penalty actually threads into the
# occurrence-arm fit (a tight slope prior shrinks the estimate toward 0), and
# that unsupported combinations error instead of silently ignoring the prior.

test_that("cover_priors() constructs and validates", {
  cp <- cover_priors()
  expect_s3_class(cp, "cover_priors")
  expect_s3_class(cp, "tobs_priors_spec")
  expect_named(cp, c("occ_intercept", "occ_slope", "pos_intercept", "pos_slope"))

  expect_error(cover_priors(occ_slope = list(mean = 0, sd = -1)), "positive")
  expect_error(cover_priors(pos_intercept = list(mean = 0)), "mean")
  # Inf sd is allowed (= no penalty on that component)
  expect_silent(cover_priors(pos_slope = list(mean = 0, sd = Inf)))
})

test_that("a tight occurrence-slope prior shrinks the occ slope (lognormal)", {
  sim <- simulate_cover(N = 150L, beta_occ = c(0, 1.4), beta_pos = c(-1, 0.3),
                        sigma_pos = 0.4, seed = 11)
  d <- sim$data

  f_unpen <- tobs(~ x, data = d, family = cover(positive = "lognormal"),
                  y = sim$y, method = "laplace")
  f_pen   <- tobs(~ x, data = d, family = cover(positive = "lognormal"),
                  y = sim$y, method = "laplace",
                  priors = cover_priors(occ_slope = list(mean = 0, sd = 0.05)))

  b_unpen <- unname(f_unpen$beta_occ["x"])
  b_pen   <- unname(f_pen$beta_occ["x"])
  expect_true(is.finite(b_unpen) && is.finite(b_pen))
  # the tight N(0, 0.05) prior pulls the slope strongly toward 0
  expect_lt(abs(b_pen), abs(b_unpen))
  expect_lt(abs(b_pen), 0.5 * abs(b_unpen))
})

test_that("priors = NULL fits unpenalised (cover priors are opt-in)", {
  sim <- simulate_cover(N = 120L, beta_occ = c(0, 1.0), beta_pos = c(-1, 0.3),
                        sigma_pos = 0.4, seed = 7)
  a <- tobs(~ x, data = sim$data, family = cover(positive = "lognormal"),
            y = sim$y, method = "laplace")
  b <- tobs(~ x, data = sim$data, family = cover(positive = "lognormal"),
            y = sim$y, method = "laplace", priors = NULL)
  expect_equal(unname(a$beta_occ["x"]), unname(b$beta_occ["x"]), tolerance = 1e-8)
})

test_that("occu_priors() is rejected for cover() with a pointer to cover_priors()", {
  sim <- simulate_cover(N = 80L, seed = 3)
  expect_error(
    tobs(~ x, data = sim$data, family = cover(positive = "lognormal"),
         y = sim$y, method = "laplace", priors = occu_priors()),
    "cover_priors"
  )
})

test_that("beta arm: a tight positive-slope prior shrinks the pos slope", {
  sim <- simulate_cover(N = 160L, beta_occ = c(0.5, 0.5), beta_pos = c(-1, 0.6),
                        sigma_pos = 0.4, seed = 5)
  d <- sim$data

  f_unpen <- tobs(~ x, data = d, family = cover(positive = "beta"),
                  y = sim$y, method = "laplace")
  # beta-arm prior now threads through tulpa_laplace_beta(beta_prior=)
  f_pen   <- tobs(~ x, data = d, family = cover(positive = "beta"),
                  y = sim$y, method = "laplace",
                  priors = cover_priors(pos_slope = list(mean = 0, sd = 0.02)))
  expect_s3_class(f_pen, "tobs_fit")

  b_unpen <- unname(f_unpen$beta_pos["x"])
  b_pen   <- unname(f_pen$beta_pos["x"])
  expect_true(is.finite(b_unpen) && is.finite(b_pen))
  expect_lt(abs(b_pen), abs(b_unpen))   # tight pos_slope prior shrinks toward 0
})

test_that("nested_laplace does not yet accept cover priors (errors, no silent drop)", {
  sim <- simulate_cover(N = 80L, seed = 9)
  expect_error(
    tobs(~ x, data = sim$data, family = cover(positive = "lognormal"),
         y = sim$y, method = "nested_laplace", priors = cover_priors()),
    "nested_laplace"
  )
})
