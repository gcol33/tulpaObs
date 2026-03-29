test_that("JSDM model constructs and prints", {
  set.seed(42)
  n <- 20; sp <- 3
  y <- matrix(rbinom(n * sp, 1, 0.4), n, sp)
  d <- data.frame(x = rnorm(n))

  mod <- occu(~ x, data = d, y = y, jsdm = TRUE, species = TRUE)
  expect_s3_class(mod, "tulpaOcc")
  expect_equal(mod$model_type, "jsdm")
  expect_equal(mod$n_species, sp)
  expect_output(print(mod), "Joint species distribution")
})

test_that("integrated model constructs", {
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  y1 <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  y2 <- matrix(rbinom(n * 4, 1, 0.4), n, 4)

  mod <- occu(~ x, ~ 1, d, y = list(s1 = y1, s2 = y2), integrated = TRUE)
  expect_s3_class(mod, "tulpaOcc")
  expect_equal(mod$model_type, "integrated")
  expect_equal(mod$n_sources, 2)
  expect_output(print(mod), "Integrated")
})

test_that("simulation functions produce correct dimensions", {
  sim_int <- simIntOcc(N_total = 30, n_data = 2, J = c(3, 4), seed = 42)
  expect_length(sim_int$y, 2)
  expect_equal(ncol(sim_int$y[[1]]), 3)
  expect_equal(ncol(sim_int$y[[2]]), 4)

  sim_tms <- simTMsOcc(N = 10, J = 3, n_species = 3, n_seasons = 4, seed = 42)
  expect_equal(dim(sim_tms$y), c(10, 3, 4, 3))

  sim_ims <- simIntMsOcc(N = 10, J = c(3, 4), n_species = 3, seed = 42)
  expect_length(sim_ims$y, 2)
})

test_that("checkIdentifiability works", {
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  y <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  mod <- occu(~ x, ~ 1, d, y)

  result <- checkIdentifiability(mod)
  expect_type(result, "list")
  expect_true("identifiable" %in% names(result))
})



test_that("pitResiduals returns uniform-ish values", {
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  mod <- occu(~ 1, ~ 1, d, y)
  fit <- occu_fit(mod, iter = 200, warmup = 100, seed = 42, verbose = FALSE)

  pit <- pitResiduals(fit, n.samples = 50)
  expect_length(pit, n)
  expect_true(all(pit >= 0 & pit <= 1))

  # KS test
  ks <- testUniformity(pit)
  expect_s3_class(ks, "htest")
})

test_that("testDispersion returns sensible output", {
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- occu_fit(occu(~ 1, ~ 1, d, y), iter = 200, warmup = 100, seed = 42, verbose = FALSE)

  disp <- testDispersion(fit, n.samples = 20)
  expect_true(is.finite(disp$ratio))
  expect_true(disp$p.value >= 0 && disp$p.value <= 1)
})

test_that("testZeroInflation returns sensible output", {
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- occu_fit(occu(~ 1, ~ 1, d, y), iter = 200, warmup = 100, seed = 42, verbose = FALSE)

  zi <- testZeroInflation(fit, n.samples = 20)
  expect_true(is.finite(zi$ratio))
})

test_that("tulpa generic diagnostics work via inheritance", {
  set.seed(42)
  n <- 30; coords <- cbind(runif(n), runif(n))
  d <- data.frame(x = rnorm(n)); z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- occu_fit(occu(~ 1, ~ 1, d, y), verbose = FALSE)

  r <- residuals(fit)$occ
  mi <- tulpa::moranI(r, coords, k = 3)
  expect_true(is.finite(mi$I))

  dw <- tulpa::durbinWatson(r)
  expect_true(dw$DW > 0 && dw$DW < 4)

  vg <- tulpa::variogram(r, coords, n_bins = 5)
  expect_true(nrow(vg) > 0)
})

test_that("update works", {
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- occu_fit(occu(~ 1, ~ 1, d, y), verbose = FALSE)

  fit2 <- update(fit, verbose = FALSE)
  expect_s3_class(fit2, "tulpaOcc_fit")
  expect_true(fit2$n_samples > 0)
})



test_that("checkModel runs without error", {
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- occu_fit(occu(~ 1, ~ 1, d, y), iter = 200, warmup = 100, seed = 42, verbose = FALSE)

  expect_output(checkModel(fit), "tulpaOcc Model Diagnostics")
})
