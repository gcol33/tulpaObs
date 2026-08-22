# svc() (the continuous NNGP spatially-varying coefficient) is wired for
# single-season occupancy: sampled by the NUTS path, and fit on the Laplace
# backends through the shared areal-BFGS nested-Laplace driver (R/occu_svc.R).
# On every other family it used to be extracted and silently dropped -- the
# model fit with no spatially-varying coefficient, no error. Those must still
# error with a pointer to the recovery-tested areal-bar route. These are fast
# error paths (the guard fires before any fit), so they run in every tier.
#
# `prior_range` is supplied even though these calls are expected to fail: without
# it svc() errors in its own constructor, and that message also matches the regex
# below -- the tests would pass while never reaching the family guard.

test_that("svc() on abun() errors instead of silently dropping", {
  set.seed(1)
  y  <- matrix(stats::rpois(30 * 3, 3), 30, 3)
  df <- data.frame(lon = stats::runif(30), lat = stats::runif(30))
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)), data = df,
         family = abun(), detection = ~ 1, y = y),
    "svc|areal|spatially-varying")
})

test_that("svc() on a multi-season family still errors", {
  set.seed(3)
  n <- 24L; J <- 3L; n_seasons <- 3L
  y <- array(stats::rbinom(n * J * n_seasons, 1, 0.4), dim = c(n, J, n_seasons))
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n))
  # The SVC surfaces are wired on the single-season occupancy marginal only.
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)), data = df,
         family = dyn_occu(), detection = ~ 1, y = y,
         colonization = ~ 1, extinction = ~ 1, method = "laplace"),
    "svc|areal|spatially-varying")
})

test_that("svc() on a method without an SVC route still errors", {
  set.seed(4)
  y  <- matrix(stats::rbinom(30 * 3, 1, 0.4), 30, 3)
  df <- data.frame(lon = stats::runif(30), lat = stats::runif(30))
  # Polya-Gamma Gibbs carries no latent-field block for a continuous surface.
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)), data = df,
         family = occu(), detection = ~ 1, y = y, method = "pg_gibbs"),
    "svc|areal|spatially-varying")
})

test_that("svc() alongside another structured term errors on the Laplace route", {
  set.seed(5)
  n <- 20L
  y  <- matrix(stats::rbinom(n * 3, 1, 0.4), n, 3)
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n),
                   cell = seq_len(n))
  # The block list carries the NNGP surfaces alone; an areal field alongside
  # them is not wired, so it errors rather than being dropped.
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)) +
           icar(graph = adj, group_var = "cell"),
         data = df, family = occu(), detection = ~ 1, y = y,
         method = "laplace"),
    "svc|spatial|not wired")
})

test_that("svc() on single-season occu() FITS under laplace", {
  set.seed(2)
  n <- 40L
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n))
  z  <- stats::rbinom(n, 1, 0.5)
  y  <- matrix(stats::rbinom(n * 4, 1, 0.6 * rep(z, 4)), n, 4)
  fit <- tobs(~ svc(lon, lat, indices = 1, nn = 8, prior_range = c(0.1, 0.05)),
              data = df, family = occu(), detection = ~ 1, y = y,
              method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  # The surface is reconstructed and exposed, matching the NUTS path's naming.
  expect_length(as.numeric(fit$svc_field), n)
  expect_true(all(is.finite(as.numeric(fit$svc_field))))
  expect_s3_class(fit$svc, "tobs_svc")
})
