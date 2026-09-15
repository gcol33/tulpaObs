# =============================================================================
# test-accessors-without-model-type.R -- cover() and occu_categorical() fits
# carry no `model$model_type`, and every generic door has to read that as "no
# family handler" rather than compare against NULL.
# =============================================================================

ctl <- list(verbose = FALSE, progress = FALSE)

test_that("the S3 handler lookup declines a missing model type", {
  expect_null(.tobs_s3_handler("ranef", NULL))
  expect_null(.tobs_s3_handler("ranef", character(0)))
})

test_that("ranef() on a cover() fit returns the empty frame", {
  skip_if_fast()
  skip_on_cran()
  sc <- simulate_cover(N = 200L, beta_occ = c(-0.5, 0.8), beta_pos = c(-1, 0.3),
                       sigma_pos = 0.4, response = "lognormal", seed = 51L)
  fc <- tobs(~ x, data = sc$data, family = cover("lognormal"), y = sc$y,
             method = "laplace", control = ctl)
  re <- ranef(fc)
  expect_s3_class(re, "data.frame")
  expect_identical(nrow(re), 0L)
})

test_that("occu_categorical() fits: empty ranef(), refused ppc() / pit_residuals()", {
  skip_if_fast()
  skip_on_cran()
  sk <- simulate_occu_categorical(N = 300L, seed = 1L)
  fk <- tobs(~ x, data = sk$data, family = occu_categorical(), y = sk$y,
             control = ctl)
  re <- ranef(fk)
  expect_s3_class(re, "data.frame")
  expect_identical(nrow(re), 0L)
  expect_error(ppc(fk), "supports single-season occupancy")
  expect_error(pit_residuals(fk), "supports single-season occupancy")
})
