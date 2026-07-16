# NUTS + a temporal() field on the observation families (removal / distance /
# fp_occu; gcol33/tulpaObs#114). A temporal() term on its own (no simultaneous
# areal field) is carried under NUTS by a FIXED-HYPER non-centered temporal field:
# the field precision (tau, rho) is fixed at the temporal-only nested-Laplace
# areal-BFGS estimate and the whitened raw ~ N(0, I) rides the SAME family NUTS
# field block the areal path uses (only the loading L and the per-site field map
# differ -- L is the eigen-loading of tau Q_temporal, the map is the period index).
# This is the dyn_abun NUTS+temporal recipe (already shipped) ported to the three
# remaining observation families. The recovery invariant: NUTS reproduces the
# integrated temporal field (cor high), recovers the abundance / occupancy slope,
# and samples without divergences. areal + temporal simultaneously under NUTS stays
# gated (combine them under nested_laplace); that gate is asserted too.

# One AR1 temporal field over Tt periods, per_t sites per period (demeaned).
.ont_ar1_field <- function(Tt, rho = 0.7, sig = 0.5, seed = 1L) {
  set.seed(seed); u <- numeric(Tt)
  u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
  for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
  u - mean(u)
}


test_that("removal() NUTS + temporal() recovers the AR1 field + slope, 0 divergences (#114)", {
  skip_on_cran()
  skip_if_fast()
  Tt <- 8L; per_t <- 30L; N <- Tt * per_t; n_seeds <- 6L
  fcor <- slope <- diverg <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    u <- .ont_ar1_field(Tt, seed = 200L + s)
    set.seed(200L + s); period <- rep(seq_len(Tt), each = per_t); x <- stats::rnorm(N)
    lambda <- exp(log(8) + 0.5 * x + u[period]); K <- 4L
    Nn <- stats::rpois(N, lambda); y <- matrix(0L, N, K); rem <- Nn
    for (k in 1:K) { y[, k] <- stats::rbinom(N, rem, 0.45); rem <- rem - y[, k] }
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period), family = removal(),
                         detection = ~ 1, y = y, method = "nuts",
                         control = list(n.iter = 500L, n.warmup = 500L,
                                        verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nuts")
      expect_null(fit$spatial_field)                # temporal-only: no areal field
      expect_length(fit$temporal_field, Tt)
    }
    slope[s]  <- fit$means[["lambda_x"]]
    diverg[s] <- mean(fit$nuts$divergent)
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_equal(mean(diverg[ok]), 0)                 # fixed-hyper clean geometry
  expect_lt(abs(mean(slope[ok]) - 0.5), 0.12)       # abundance slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.80)     # AR1 temporal field recovered
})


test_that("distance() NUTS + temporal() recovers the AR1 field + slope (#114)", {
  skip_on_cran()
  skip_if_fast()
  Tt <- 8L; per_t <- 30L; N <- Tt * per_t; n_seeds <- 6L
  cut <- c(0, 10, 20, 30, 40); B <- 4L; sigd <- exp(3.0)
  g_mid <- exp(-((head(cut, -1) + tail(cut, -1)) / 2)^2 / (2 * sigd^2))
  pi_b <- g_mid * diff(cut) / 40 * 0.6
  fcor <- slope <- diverg <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    u <- .ont_ar1_field(Tt, seed = 300L + s)
    set.seed(300L + s); period <- rep(seq_len(Tt), each = per_t); x <- stats::rnorm(N)
    lamd <- exp(log(40) + 0.5 * x + u[period])
    Nd <- stats::rpois(N, lamd); y <- matrix(0L, N, B)
    for (i in seq_len(N)) y[i, ] <- stats::rbinom(B, Nd[i], pi_b)
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period),
                         family = distance(key = "halfnorm", transect = "line", cutpoints = cut),
                         detection = ~ 1, y = y, method = "nuts",
                         control = list(n.iter = 500L, n.warmup = 500L,
                                        verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nuts")
      expect_length(fit$temporal_field, Tt)
    }
    slope[s]  <- fit$means[["lambda_x"]]
    diverg[s] <- mean(fit$nuts$divergent)
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_equal(mean(diverg[ok]), 0)
  expect_lt(abs(mean(slope[ok]) - 0.5), 0.12)
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.80)
})


test_that("fp_occu() NUTS + temporal() recovers the AR1 field, 0 divergences (#114)", {
  skip_on_cran()
  skip_if_fast()
  Tt <- 8L; per_t <- 40L; N <- Tt * per_t; n_seeds <- 6L
  fcor <- slope <- diverg <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    u <- .ont_ar1_field(Tt, sig = 0.7, seed = 400L + s)
    set.seed(400L + s); period <- rep(seq_len(Tt), each = per_t); x <- stats::rnorm(N)
    psi <- stats::plogis(0.3 + 0.7 * x + u[period]); z <- stats::rbinom(N, 1, psi); J <- 6L
    y <- matrix(0L, N, J)
    for (i in seq_len(N)) for (j in 1:J)
      y[i, j] <- if (z[i] == 1) sample(0:2, 1, prob = c(0.45, 0.25, 0.30))
                 else sample(0:1, 1, prob = c(0.92, 0.08))
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period), family = fp_occu(),
                         detection = ~ 1, y = y, method = "nuts",
                         control = list(n.iter = 600L, n.warmup = 600L, adapt.delta = 0.95,
                                        verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nuts")
      expect_length(fit$temporal_field, Tt)
    }
    slope[s]  <- fit$means[["psi_x"]]
    diverg[s] <- mean(fit$nuts$divergent)
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_equal(mean(diverg[ok]), 0)
  expect_lt(abs(mean(slope[ok]) - 0.7), 0.20)       # occupancy slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.70)
})


test_that("distance() hazard-key NUTS + temporal() recovers field + shape (#114)", {
  skip_on_cran()
  skip_if_fast()
  Tt <- 8L; per_t <- 30L; N <- Tt * per_t; n_seeds <- 6L
  cut <- c(0, 10, 20, 30, 40); W <- 40; B <- 4L; sigma <- 18; shape <- 3.5
  mid <- (head(cut, -1) + tail(cut, -1)) / 2
  pi_b <- (1 - exp(-(mid / sigma)^(-shape))) * diff(cut) / W
  fcor <- shp <- diverg <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    u <- .ont_ar1_field(Tt, seed = 700L + s)
    set.seed(700L + s); period <- rep(seq_len(Tt), each = per_t); x <- stats::rnorm(N)
    lamd <- exp(log(45) + 0.5 * x + u[period])
    Nd <- stats::rpois(N, lamd); y <- matrix(0L, N, B)
    for (i in seq_len(N)) y[i, ] <- stats::rbinom(B, Nd[i], pi_b)
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period),
                         family = distance(key = "hazard", transect = "line", cutpoints = cut),
                         detection = ~ 1, y = y, method = "nuts",
                         control = list(n.iter = 500L, n.warmup = 500L,
                                        verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nuts")
      expect_length(fit$temporal_field, Tt)
      expect_true("log_shape" %in% names(fit$means))   # hazard shape carried
    }
    diverg[s] <- mean(fit$nuts$divergent)
    shp[s]    <- fit$means[["log_shape"]]
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(shp)
  expect_gte(mean(ok), 0.75)
  expect_equal(mean(diverg[ok]), 0)
  expect_lt(abs(mean(shp[ok]) - log(shape)), 0.4)   # scalar shape recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.75)     # AR1 temporal field recovered
})


test_that("NUTS + areal AND temporal simultaneously stays gated (#114)", {
  # The temporal field rides the NUTS field block on its own; a simultaneous areal
  # field is not wired under NUTS (combine them under nested_laplace).
  adj <- diag(0, 6L); for (i in 1:5) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  set.seed(1); N <- 60L; cell <- rep(seq_len(6L), each = 10L)
  period <- rep(rep(1:5, each = 2L), length.out = N)
  y <- matrix(rpois(N * 3L, 3), N, 3L)
  expect_error(
    suppressWarnings(tobs(
      ~ 1 + icar(graph = adj, group_var = "cell") + temporal(period, type = "ar1"),
      data = data.frame(cell = cell, period = period), family = removal(),
      detection = ~ 1, y = y, method = "nuts",
      control = list(verbose = FALSE, progress = FALSE))),
    "temporal|areal|nested_laplace")
})
