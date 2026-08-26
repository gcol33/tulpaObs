test_that("single-season occupancy model runs", {
  set.seed(123)
  n_sites <- 100
  max_visits <- 3

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

  fit <- tobs(
    formula   = ~ elevation,
    data      = site_data,
    family    = occu(),
    detection = ~ 1,
    y         = y,
    method    = "laplace",
    control   = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$n_samples > 0)
  expect_true(ncol(fit$draws) >= 3)

  # Sanity bounds only — N = 100, J = 3, p = 0.4 leaves the EM with non-
  # trivial identifiability slack along the psi*p^J ridge, so the point
  # estimates can drift well off the true (-0.5, 1.2, qlogis(0.4)) values.
  expect_true(fit$means[1] > -3 && fit$means[1] < 4)
  expect_true(fit$means[2] > -1 && fit$means[2] < 5)
  expect_true(fit$intercepts$psi > 0.05 && fit$intercepts$psi < 0.99)
  expect_true(fit$intercepts$p > 0.05 && fit$intercepts$p < 0.95)
})

test_that("tobs validates inputs", {
  expect_error(
    tobs(~ 1, data.frame(x = 1:5), family = occu(),
         detection = ~ 1, y = matrix(0, 3, 2)),
    "y has 3 rows but data has 5 rows"
  )
  expect_error(
    tobs(~ 1, data.frame(x = 1:5), family = occu(),
         detection = ~ 1, y = c(0, 1, 0)),
    "y must be a matrix"
  )
  expect_error(
    tobs(~ 1, data.frame(x = 1:3)),
    "`family` is required"
  )
  expect_error(
    tobs(~ 1, data.frame(x = 1:3), family = "occu"),
    "must be a tobs_family"
  )
})

test_that("NA visits handled correctly", {
  n_sites <- 20
  y <- matrix(c(1, 0, NA), nrow = n_sites, ncol = 3, byrow = TRUE)
  site_data <- data.frame(x = rnorm(n_sites))
  mod <- .tobs_build_model(occ_formula = ~ 1, det_formula = ~ 1,
                           data = site_data, y = y)
  expect_true(all(mod$y[, 3] == -1L))
  expect_output(print(mod), "Single-season occupancy model")
})
