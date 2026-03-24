test_that("single-season occupancy model runs", {
  set.seed(123)
  n_sites <- 100
  max_visits <- 3

  # Simulate data
  x_elev <- rnorm(n_sites)
  psi_true <- plogis(-0.5 + 1.2 * x_elev)
  z <- rbinom(n_sites, 1, psi_true)

  p_true <- 0.4
  y <- matrix(NA, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    if (z[i] == 1) {
      y[i, ] <- rbinom(max_visits, 1, p_true)
    } else {
      y[i, ] <- 0L
    }
  }

  site_data <- data.frame(elevation = x_elev)

  # Build model
  mod <- occ(~ elevation, ~ 1, data = site_data, y = y)
  expect_s3_class(mod, "tulpaOcc_model")
  expect_equal(mod$n_sites, n_sites)
  expect_equal(mod$max_visits, max_visits)
  expect_equal(mod$p_occ, 2)  # intercept + elevation
  expect_equal(mod$p_det, 1)  # intercept only

  # Fit (short run for testing — uses full NUTS via tulpa)
  fit <- occ_fit(mod, iter = 500, warmup = 250, seed = 42, verbose = FALSE)
  expect_s3_class(fit, "tulpaOcc_fit")
  expect_equal(fit$n_samples, 250)
  expect_equal(ncol(fit$draws), 3)  # 2 occ + 1 det

  # Check that estimates are in reasonable range
  # True: beta_occ = c(-0.5, 1.2), logit(0.4) ~ -0.405
  expect_true(fit$means[1] > -3 && fit$means[1] < 2)  # intercept
  expect_true(fit$means[2] > -1 && fit$means[2] < 4)  # elevation effect
  expect_true(fit$mean_psi > 0.1 && fit$mean_psi < 0.9)
  expect_true(fit$mean_p > 0.1 && fit$mean_p < 0.9)
})

test_that("occ model validates inputs", {
  expect_error(occ(~ 1, ~ 1, data.frame(x = 1:5), y = matrix(0, 3, 2)),
               "y has 3 rows but data has 5 rows")
  expect_error(occ(~ 1, ~ 1, data.frame(x = 1:5), y = c(0, 1, 0)),
               "y must be a matrix")
})

test_that("NA visits handled correctly", {
  n_sites <- 20
  y <- matrix(c(1, 0, NA), nrow = n_sites, ncol = 3, byrow = TRUE)
  site_data <- data.frame(x = rnorm(n_sites))
  mod <- occ(~ 1, ~ 1, data = site_data, y = y)
  expect_true(all(mod$y[, 3] == -1L))
})
