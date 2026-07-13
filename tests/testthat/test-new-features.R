.simple_fit <- function(formula = ~ 1, det = ~ 1, n = 30, seed = 42,
                       method = "nuts", iter = 200, warmup = 100) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  ctrl <- if (method == "nuts")
    list(n.iter = iter, n.warmup = warmup, seed = seed, verbose = FALSE)
  else
    list(verbose = FALSE)
  fit <- tobs(formula, data = d, family = occu(),
              detection = det, y = y, method = method,
              control = ctrl)
  list(fit = fit, y = y, d = d, n = n)
}

test_that("JSDM model fits", {
  skip_if_fast()
  set.seed(42)
  n <- 20; sp <- 3
  y <- matrix(rbinom(n * sp, 1, 0.4), n, sp)
  d <- data.frame(x = rnorm(n))

  fit <- tobs(~ x, data = d, family = jsdm(), y = y, species = TRUE,
              method = "nuts",
              control = list(n.iter = 200, n.warmup = 100, seed = 42, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
})

test_that("integrated model fits", {
  skip_if_fast()
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  y1 <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  y2 <- matrix(rbinom(n * 4, 1, 0.4), n, 4)

  fit <- tobs(~ x, data = d, family = int_occu(),
              detection = ~ 1, y = list(s1 = y1, s2 = y2),
              method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
})

test_that("simulation functions produce correct dimensions", {
  sim_int <- simulate_int_occu(N_total = 30, n_data = 2, J = c(3, 4), seed = 42)
  expect_length(sim_int$y, 2)
  expect_equal(ncol(sim_int$y[[1]]), 3)
  expect_equal(ncol(sim_int$y[[2]]), 4)

  sim_tms <- simulate_ms_dyn_occu(N = 10, J = 3, n_species = 3, n_seasons = 4, seed = 42)
  expect_equal(dim(sim_tms$y), c(10, 3, 4, 3))

  sim_ims <- simulate_ms_int_occu(N = 10, J = c(3, 4), n_species = 3, seed = 42)
  expect_length(sim_ims$y, 2)
})

test_that("tobs_check_id works", {
  set.seed(42)
  n <- 20
  d <- data.frame(x = rnorm(n))
  y <- matrix(rbinom(n * 3, 1, 0.3), n, 3)
  mod <- .tobs_build_model(~ x, ~ 1, d, y)

  result <- tobs_check_id(mod)
  expect_type(result, "list")
  expect_true("identifiable" %in% names(result))
})

test_that("tobs_pit_residuals returns uniform-ish values", {
  skip_if_fast()
  res <- .simple_fit(method = "nuts")
  pit <- tobs_pit_residuals(res$fit, n.samples = 50)
  expect_length(pit, res$n)
  expect_true(all(pit >= 0 & pit <= 1))
  ks <- tobs_test_uniformity(pit)
  expect_s3_class(ks, "htest")
})

test_that("tobs_test_dispersion returns sensible output", {
  skip_if_fast()
  res <- .simple_fit(method = "nuts")
  disp <- tobs_test_dispersion(res$fit, n.samples = 20)
  expect_true(is.finite(disp$ratio))
  expect_true(disp$p.value >= 0 && disp$p.value <= 1)
})

test_that("tobs_test_zero_inflation returns sensible output", {
  skip_if_fast()
  res <- .simple_fit(method = "nuts")
  zi <- tobs_test_zero_inflation(res$fit, n.samples = 20)
  expect_true(is.finite(zi$ratio))
})

test_that("tulpa generic diagnostics work via inheritance", {
  skip_if_fast()
  set.seed(42)
  n <- 30; coords <- cbind(runif(n), runif(n))
  res <- .simple_fit(method = "laplace", n = n)
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
  skip_if_fast()
  res <- .simple_fit(method = "laplace")
  fit2 <- update(res$fit, verbose = FALSE)
  expect_s3_class(fit2, "tobs_fit")
  expect_true(fit2$n_samples > 0)
})

test_that("tobs_check runs without error", {
  skip_if_fast()
  res <- .simple_fit(method = "nuts")
  expect_output(tobs_check(res$fit), "tobs Model Diagnostics")
})

test_that("tobs_check reports Moran's I when coords are supplied", {
  skip_if_fast()
  set.seed(7)
  res <- .simple_fit(method = "nuts", n = 30)
  coords <- cbind(runif(res$n), runif(res$n))
  expect_output(tobs_check(res$fit, coords = coords), "Moran's I")
})

test_that("tobs_ppc is calibrated for a correct occupancy model", {
  skip_if_fast()
  # The latent state is sampled from its full conditional given the detection
  # history, so the Bayesian p-value sits near 0.5 for a correct fit. Drawing
  # z from the prior psi instead collapsed it toward 0.
  set.seed(11)
  n <- 200
  d <- data.frame(x = rnorm(n))
  psi <- plogis(0.3 + 0.9 * d$x)
  z <- rbinom(n, 1, psi)
  y <- matrix(rbinom(n * 5, 1, z * 0.45), n, 5)
  fit <- tobs(~ x, data = d, family = occu(), detection = ~ 1, y = y,
              method = "nuts",
              control = list(n.iter = 800, n.warmup = 400, seed = 11,
                             verbose = FALSE))
  bp <- tobs_ppc(fit, n.samples = 300)$bayesian.p
  expect_gt(bp, 0.1)
  expect_lt(bp, 0.9)
})
