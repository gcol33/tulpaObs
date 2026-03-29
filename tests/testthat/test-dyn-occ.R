test_that("dynamic occupancy model runs", {
  set.seed(42)
  n_sites <- 30
  n_seasons <- 3
  max_visits <- 3

  y_array <- array(0L, dim = c(n_sites, max_visits, n_seasons))
  z <- matrix(0, n_sites, n_seasons)
  z[, 1] <- rbinom(n_sites, 1, 0.6)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t-1] * (1 - rbinom(n_sites, 1, 0.1)) +
              (1 - z[, t-1]) * rbinom(n_sites, 1, 0.2)
  }
  for (i in seq_len(n_sites))
    for (t in seq_len(n_seasons))
      if (z[i, t] == 1) y_array[i, , t] <- rbinom(max_visits, 1, 0.5)

  mod <- occu(~ 1, ~ 1, data.frame(x = rnorm(n_sites)), y_array,
              col_formula = ~ 1, ext_formula = ~ 1)
  expect_s3_class(mod, "tulpaOcc")
  expect_equal(mod$model_type, "dynamic")
  expect_equal(mod$n_sites, n_sites)
  expect_equal(mod$n_seasons, n_seasons)

  fit <- occu_fit(mod, iter = 100, warmup = 50, seed = 1, verbose = FALSE)
  expect_s3_class(fit, "tulpaOcc_fit")
  expect_equal(fit$n_params, 4)
  expect_true(fit$intercepts$psi1 > 0 && fit$intercepts$psi1 < 1)
})

test_that("dynamic occu validates inputs", {
  expect_error(occu(~ 1, ~ 1, data.frame(x = 1:5), y = matrix(0, 3, 2),
                    col_formula = ~ 1),
               "3D array")
})
