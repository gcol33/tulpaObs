# Structural NUTS coverage for the four components previously flagged as
# blocked by upstream tulpa bugs: tobs_temporal, multi-term tobs_re, tobs_svc,
# tobs_latent. Smoke iter counts: these verify the end-to-end C++ NUTS path
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

test_that("NUTS runs with tobs_temporal (AR1) attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$time <- sample.int(4, N, replace = TRUE)

  fit <- tobs(
    ~ occ_cov1, data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y,
    temporal = tobs_temporal(type = "ar1", time = "time"),
    engine = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "temporal-ar1")
  expect_s3_class(fit$temporal, "tobs_temporal")
})

test_that("NUTS runs with two-term tobs_re list attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  sim$data$grp  <- sample.int(5, N, replace = TRUE)
  sim$data$time <- sample.int(4, N, replace = TRUE)

  fit <- tobs(
    ~ occ_cov1, data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y,
    re = list(
      tobs_re(group = "grp",  type = "intercept"),
      tobs_re(group = "time", type = "intercept")
    ),
    engine = "nuts", control = ctl_nuts()
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "multi-re")
  expect_true(is.list(fit$re))
  expect_length(fit$re, 2)
  expect_s3_class(fit$re[[1]], "tobs_re")
})

test_that("NUTS runs with tobs_svc attached", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3
  sim <- simulate_occu(N = N, J = J, seed = 1)
  coords <- matrix(runif(N * 2), N, 2)

  fit <- tobs(
    ~ occ_cov1, data = sim$data, family = occu(),
    detection = ~ det_cov1, y = sim$y,
    engine = "nuts",
    control = ctl_nuts(list(svc = tobs_svc(indices = 1L, coords = coords, nn = 8)))
  )
  expect_nuts_fit(fit, n_expected_cols = 4, label = "svc")
  expect_s3_class(fit$svc, "tobs_svc")
})

test_that("NUTS runs with tobs_latent on ms_occu", {
  skip_on_cran()

  set.seed(1)
  N <- 40; J <- 3; n_sp <- 4
  ms <- simulate_ms_occu(N = N, J = J, n_species = n_sp, seed = 1)
  sp_names <- paste0("sp", seq_len(n_sp))

  fit <- tobs(
    ~ x, data = ms$data, family = ms_occu(),
    detection = ~ 1, y = ms$y, species = sp_names,
    engine = "nuts",
    control = ctl_nuts(list(latent = tobs_latent(n_factors = 2)))
  )
  expect_nuts_fit(fit, n_expected_cols = 2, label = "latent")
  expect_s3_class(fit$latent, "tobs_latent")
  expect_equal(fit$latent$n_factors, 2L)
})
