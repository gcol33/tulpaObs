# Tests for SPDE spatial spec via tobs_spde()

test_that("tobs_spde creates valid specification", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  spec <- tobs_spde(coords)

  expect_s3_class(spec, "tobs_spatial")
  expect_equal(spec$type, "spde")
  expect_true(spec$n_units >= 30)
  expect_equal(spec$shared, c(TRUE, FALSE))
})

test_that("tobs_spde creates spec from formula", {
  set.seed(42)
  df <- data.frame(x = runif(20), y = runif(20))
  spec <- tobs_spde(~ x + y, data = df)

  expect_s3_class(spec, "tobs_spatial")
  expect_equal(spec$type, "spde")
})

test_that("tobs_spde print method works", {
  set.seed(42)
  spec <- tobs_spde(cbind(runif(20), runif(20)))
  expect_output(print(spec), "spde")
  expect_output(print(spec), "Matern")
})

test_that("tobs_spde with fractional nu creates valid spec", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  spec <- tobs_spde(coords, nu = 1.5)

  expect_equal(spec$nu, 1.5)
  expect_s3_class(spec, "tobs_spatial")
})
