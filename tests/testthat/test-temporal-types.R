.tt_removal_data <- function(seed = 1L, Tt = 6L, per_t = 25L) {
  set.seed(seed)
  N <- Tt * per_t
  period <- rep(seq_len(Tt), each = per_t)
  u <- stats::rnorm(Tt, 0, 0.4); u <- u - mean(u)
  x <- stats::rnorm(N)
  lambda <- exp(log(8) + 0.5 * x + u[period])
  K <- 3L; Nn <- stats::rpois(N, lambda); y <- matrix(0L, N, K); rem <- Nn
  for (k in 1:K) { y[, k] <- stats::rbinom(N, rem, 0.45); rem <- rem - y[, k] }
  list(data = data.frame(x = x, period = period), y = y, Tt = Tt)
}

test_that(".tobs_rows keeps one column for width-one rows", {
  m <- tulpaObs:::.tobs_rows(list(c(tau = 1), c(tau = 2), c(tau = 3)), 1L)
  expect_identical(dim(m), c(3L, 1L))
  expect_equal(as.numeric(crossprod(c(0.2, 0.3, 0.5), m)), 2.3)
  m2 <- tulpaObs:::.tobs_rows(list(1:2, 3:4), 2L)
  expect_equal(m2, matrix(c(1, 2, 3, 4), 2, byrow = TRUE))
})

test_that("temporal NUTS warm start reads an absent rho as NA", {
  for (type in c("rw1", "rw2", "iid", "ar1")) {
    temporal <- list(type = type, time_idx = rep(1:5, each = 2L), n_times = 5L)
    hyper <- if (type == "ar1") c(tau = 2, rho = 0.5) else c(tau = 2)
    tw <- tulpaObs:::.tobs_nuts_temporal_warm(temporal, 10L, "removal",
      function(max_iter, tol) list(temporal_hyper = hyper))
    expect_identical(tw$n_field_units, 5L, info = type)
    expect_identical(tw$field_map, rep(1:5, each = 2L), info = type)
    expect_true(is.finite(tw$fl$tau), info = type)
    expect_identical(is.na(tw$fl$rho), type != "ar1", info = type)
  }
  expect_error(
    tulpaObs:::.tobs_nuts_temporal_warm(list(type = "rw1", time_idx = 1:3), 4L,
                                        "removal", function(...) NULL),
    "one time index per site is required for removal NUTS")
})

test_that("temporal(prior.sigma =) is validated and threaded to the engine (#437)", {
  tt <- tulpaObs:::.tobs_term_temporal(rep(1:4, 3), type = "rw1",
                                       prior.sigma = c(0.5, 0.05))
  expect_identical(tt$prior_sigma, c(0.5, 0.05))
  expect_null(tulpaObs:::.tobs_term_temporal(rep(1:4, 3))$prior_sigma)
  for (bad in list(c(-1, 0.1), c(1, 1), 1, c(1, NA)))
    expect_error(tulpaObs:::.tobs_term_temporal(rep(1:4, 3), prior.sigma = bad),
                 "prior.sigma")
  expect_error(tulpaObs:::.tobs_term_temporal(rep(1:4, 3), tau.shape = 2),
               "unused argument")
})

test_that("temporal(prior.sigma =) reaches the occupancy NUTS engine (#437)", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu(N = 60, J = 3, seed = 1)
  sim$data$time <- rep(1:6, length.out = 60)
  ctl <- list(n.iter = 300L, n.warmup = 300L, n.chains = 1L, verbose = FALSE,
              progress = FALSE, seed = 1L)
  fit_with <- function(f)
    tobs(f, data = sim$data, family = occu(), detection = ~ det_cov1,
         y = sim$y, method = "nuts", control = ctl)
  tight <- fit_with(~ occ_cov1 + temporal(time, type = "rw1",
                                          prior.sigma = c(0.01, 0.01)))
  dflt  <- fit_with(~ occ_cov1 + temporal(time, type = "rw1"))
  # The temporal log-precision is the first sampler coordinate after the fixed
  # effects; P(sigma > 0.01) = 0.01 pins the field SD near zero.
  k <- sum(grepl("^(psi|p)_", colnames(tight$draws))) + 1L
  expect_gt(mean(tight$draws[, k]), mean(dflt$draws[, k]) + 3)
})

test_that("removal() fits every documented temporal type under nested_laplace", {
  skip_on_cran()
  d <- .tt_removal_data()
  for (type in c("rw1", "rw2", "iid", "ar1")) {
    f <- stats::as.formula(sprintf("~ x + temporal(period, type = \"%s\")", type))
    fit <- tobs(f, data = d$data, family = removal(), detection = ~ 1, y = d$y,
                method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    expect_length(fit$temporal_field, d$Tt)
    expect_true(is.finite(fit$means[["lambda_x"]]), info = type)
  }
})

test_that("removal() NUTS fits every documented temporal type", {
  skip_on_cran()
  skip_if_fast()
  d <- .tt_removal_data(seed = 2L)
  for (type in c("rw1", "rw2", "iid", "ar1")) {
    f <- stats::as.formula(sprintf("~ x + temporal(period, type = \"%s\")", type))
    fit <- tobs(f, data = d$data, family = removal(), detection = ~ 1, y = d$y,
                method = "nuts",
                control = list(n.iter = 150L, n.warmup = 150L, n.chains = 1L,
                               verbose = FALSE, progress = FALSE))
    expect_true(all(is.finite(fit$means)), info = type)
    expect_lt(abs(fit$means[["lambda_x"]] - 0.5), 0.3)
  }
})
