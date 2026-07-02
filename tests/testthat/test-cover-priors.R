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

  f_unpen <- tobs(~ x, data = d, family = cover(response = "lognormal"),
                  y = sim$y, method = "laplace")
  f_pen   <- tobs(~ x, data = d, family = cover(response = "lognormal"),
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
  a <- tobs(~ x, data = sim$data, family = cover(response = "lognormal"),
            y = sim$y, method = "laplace")
  b <- tobs(~ x, data = sim$data, family = cover(response = "lognormal"),
            y = sim$y, method = "laplace", priors = NULL)
  expect_equal(unname(a$beta_occ["x"]), unname(b$beta_occ["x"]), tolerance = 1e-8)
})

test_that("occu_priors() is rejected for cover() with a pointer to cover_priors()", {
  sim <- simulate_cover(N = 80L, seed = 3)
  expect_error(
    tobs(~ x, data = sim$data, family = cover(response = "lognormal"),
         y = sim$y, method = "laplace", priors = occu_priors()),
    "cover_priors"
  )
})

test_that("beta arm: a tight positive-slope prior shrinks the pos slope", {
  sim <- simulate_cover(N = 160L, beta_occ = c(0.5, 0.5), beta_pos = c(-1, 0.6),
                        sigma_pos = 0.4, seed = 5)
  d <- sim$data

  f_unpen <- tobs(~ x, data = d, family = cover(response = "beta"),
                  y = sim$y, method = "laplace")
  # beta-arm prior now threads through tulpa_laplace_beta(beta_prior=)
  f_pen   <- tobs(~ x, data = d, family = cover(response = "beta"),
                  y = sim$y, method = "laplace",
                  priors = cover_priors(pos_slope = list(mean = 0, sd = 0.02)))
  expect_s3_class(f_pen, "tobs_fit")

  b_unpen <- unname(f_unpen$beta_pos["x"])
  b_pen   <- unname(f_pen$beta_pos["x"])
  expect_true(is.finite(b_unpen) && is.finite(b_pen))
  expect_lt(abs(b_pen), abs(b_unpen))   # tight pos_slope prior shrinks toward 0
})

test_that("nested_laplace cover threads fixed-effect priors (#54)", {
  skip_if_fast()
  skip_on_cran()
  # Spatial cover data on a chain graph; the nested-Laplace path requires an
  # areal term, so a non-spatial formula cannot exercise the prior plumbing.
  set.seed(9)
  N <- 220L; n_s <- 25L
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  w_s <- 0.6 * rnorm(n_s)
  x   <- rnorm(N)
  eta_occ <- -0.3 + 0.7 * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- 0.4 - 0.5 * x + w_s[spatial_idx]
  y <- ifelse(occur == 1L, pmin(exp(rnorm(N, eta_pos, 0.4)), 1 - 1e-6), 0)
  d <- data.frame(x = x, region = factor(spatial_idx))

  nbr <- lapply(seq_len(n_s),
                function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

  ctrl <- list(sigma.grid = c(0.4, 0.8), sigma.pos.grid = c(0.0, 0.6),
               progress = FALSE)
  f_unpen <- tobs(~ x + icar(graph = adj, group_var = "region"), data = d,
                  family = cover("lognormal"), y = y,
                  method = "nested_laplace", control = ctrl)
  # A tight prior on the positive-arm slope must shrink it toward zero relative
  # to the unpenalised fit -- proof the prior is applied, not dropped.
  f_pen <- tobs(~ x + icar(graph = adj, group_var = "region"), data = d,
                family = cover("lognormal"), y = y,
                method = "nested_laplace", control = ctrl,
                priors = cover_priors(pos_slope = list(mean = 0, sd = 0.02)))

  expect_s3_class(f_pen, "cover_fit")
  b_unpen <- unname(f_unpen$beta_pos["x"])
  b_pen   <- unname(f_pen$beta_pos["x"])
  expect_true(is.finite(b_unpen) && is.finite(b_pen))
  expect_lt(abs(b_pen), abs(b_unpen))
})
