.simple_fit <- function(formula = ~ 1, det = ~ 1, n = 30, seed = 42,
                       engine = "nuts", iter = 200, warmup = 100) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  ctrl <- if (engine == "nuts")
    list(iter = iter, warmup = warmup, seed = seed, verbose = FALSE)
  else
    list(verbose = FALSE)
  fit <- tobs(formula, data = d, family = occu(),
              detection = det, y = y, engine = engine,
              control = ctrl)
  list(fit = fit, y = y, d = d, n = n)
}

test_that("JSDM model fits", {
  set.seed(42)
  n <- 20; sp <- 3
  y <- matrix(rbinom(n * sp, 1, 0.4), n, sp)
  d <- data.frame(x = rnorm(n))

  fit <- tobs(~ x, data = d, family = jsdm(), y = y, species = TRUE,
              engine = "nuts",
              control = list(iter = 200, warmup = 100, seed = 42, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
})

test_that("integrated model fits", {
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  y1 <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  y2 <- matrix(rbinom(n * 4, 1, 0.4), n, 4)

  fit <- tobs(~ x, data = d, family = int_occu(),
              detection = ~ 1, y = list(s1 = y1, s2 = y2),
              engine = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
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
  mod <- .tobs_build_model(~ x, ~ 1, d, y)

  result <- checkIdentifiability(mod)
  expect_type(result, "list")
  expect_true("identifiable" %in% names(result))
})

test_that("pitResiduals returns uniform-ish values", {
  res <- .simple_fit(engine = "nuts")
  pit <- pitResiduals(res$fit, n.samples = 50)
  expect_length(pit, res$n)
  expect_true(all(pit >= 0 & pit <= 1))
  ks <- testUniformity(pit)
  expect_s3_class(ks, "htest")
})

test_that("testDispersion returns sensible output", {
  res <- .simple_fit(engine = "nuts")
  disp <- testDispersion(res$fit, n.samples = 20)
  expect_true(is.finite(disp$ratio))
  expect_true(disp$p.value >= 0 && disp$p.value <= 1)
})

test_that("testZeroInflation returns sensible output", {
  res <- .simple_fit(engine = "nuts")
  zi <- testZeroInflation(res$fit, n.samples = 20)
  expect_true(is.finite(zi$ratio))
})

test_that("tulpa generic diagnostics work via inheritance", {
  set.seed(42)
  n <- 30; coords <- cbind(runif(n), runif(n))
  res <- .simple_fit(engine = "laplace", n = n)
  r <- residuals(res$fit)$occ
  mi <- tulpa::moran_i(r, coords, k = 3)
  expect_true(is.finite(unname(mi$statistic)))
  dw <- tulpa::durbin_watson(r)
  dw_stat <- unname(dw$statistic)
  expect_true(dw_stat > 0 && dw_stat < 4)
  vg <- tulpa::tulpa_variogram(r, coords, n_bins = 5)
  expect_true(nrow(vg) > 0)
})

test_that("update works", {
  res <- .simple_fit(engine = "laplace")
  fit2 <- update(res$fit, verbose = FALSE)
  expect_s3_class(fit2, "tobs_fit")
  expect_true(fit2$n_samples > 0)
})

test_that("checkModel runs without error", {
  res <- .simple_fit(engine = "nuts")
  expect_output(checkModel(res$fit), "tobs Model Diagnostics")
})
