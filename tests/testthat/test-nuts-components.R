# Structural NUTS coverage for the four components previously flagged as
# blocked by upstream tulpa bugs: temporal(), multi-term re(), svc(),
# latent(). Smoke iter counts: these verify the end-to-end C++ NUTS path
# runs without crashing and returns a fit with the expected shape, not that
# gradients are correct or posteriors are calibrated. Recovery is a separate
# concern (see test-spde-occ.R for the recovery pattern with N=400).

ctl_nuts <- function(extra = NULL) {
  c(list(iter = 60, warmup = 30, seed = 1, verbose = FALSE), extra)
}

expect_nuts_fit <- function(fit, n_expected_cols, label = "") {
  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$n_samples > 0, label = paste(label, "n_samples > 0"))
  expect_true(is.matrix(fit$draws), label = paste(label, "draws is matrix"))
  expect_gte(ncol(fit$draws), n_expected_cols)
  expect_equal(nrow(fit$draws), fit$n_samples, label = paste(label, "draws row count"))
  expect_true(all(is.finite(fit$means)), label = paste(label, "means all finite"))
}

test_that("NUTS runs with a temporal(ar1) term attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$time <- sample.int(4, N, replace = TRUE)

  fit <- tobs(
    ~ occ_cov1 + temporal(time, type = "ar1"), data = sim$data,
    family = occu(), detection = ~ det_cov1, y = sim$y,
    method = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "temporal-ar1")
  expect_s3_class(fit$temporal, "tobs_temporal")
})

test_that("NUTS runs with two re() terms attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$grp  <- sample.int(5, N, replace = TRUE)
  sim$data$time <- sample.int(4, N, replace = TRUE)

  fit <- tobs(
    ~ occ_cov1 + re(grp) + re(time), data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y,
    method = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "multi-re")
  expect_true(is.list(fit$re))
  expect_length(fit$re, 2)
  expect_s3_class(fit$re[[1]], "tobs_re")
})

test_that("lme4 bar syntax fits identically to the equivalent re() call", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$grp <- sample.int(5, N, replace = TRUE)

  fit_bar <- tobs(
    ~ occ_cov1 + (1 | grp), data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y, method = "nuts", control = ctl_nuts()
  )
  fit_re <- tobs(
    ~ occ_cov1 + re(grp), data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y, method = "nuts", control = ctl_nuts()
  )

  expect_nuts_fit(fit_bar, n_expected_cols = 3, label = "bar (1|grp)")
  expect_s3_class(fit_bar$re[[1]], "tobs_re")
  # Same spec + same seed => bit-identical fit.
  expect_equal(fit_bar$means, fit_re$means)
})

test_that("NUTS runs with an svc() term attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$lon <- runif(N)
  sim$data$lat <- runif(N)

  fit <- tobs(
    ~ occ_cov1 + svc(lon, lat, indices = 1L, nn = 8), data = sim$data,
    family = occu(), detection = ~ det_cov1, y = sim$y,
    method = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "svc")
  expect_s3_class(fit$svc, "tobs_svc")
})

test_that("NUTS runs with a latent() term on ms_occu", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3; n_sp <- 4
  ms <- simulate_ms_occu(N = N, J = J, n_species = n_sp, seed = 1)
  sp_names <- paste0("sp", seq_len(n_sp))

  fit <- tobs(
    ~ x + latent(2), data = ms$data, family = ms_occu(),
    detection = ~ 1, y = ms$y, species = sp_names,
    method = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 2, label = "latent")
  expect_s3_class(fit$latent, "tobs_latent")
  expect_equal(fit$latent$n_factors, 2L)
})
