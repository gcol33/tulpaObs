# Nested-Laplace generalised beyond single-season occupancy: the multi-block
# latent prior is attached to the state ("occ") M-step block of the integrated,
# community, and dynamic callbacks (the same builders the single-Laplace path
# uses -- there is no `build_*_callbacks_nested` duplicate). These are
# smoke/shape tests matching the bar of test-nested-laplace-occu.R: the
# underlying multi-block engine's calibrated recovery is the same Phase D
# follow-up tracked for single-season occupancy. They assert the pipeline runs,
# attaches the prior, and recovers a site-length latent field.

chain_adj <- function(n) {
  adj <- matrix(0, n, n)
  for (i in seq_len(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  adj
}


test_that("integrated occupancy fits nested_laplace with a spatial field", {
  set.seed(1)
  n <- 40; adj <- chain_adj(n)
  z <- rbinom(n, 1, plogis(0.2 + 0.5 * (d <- rnorm(n))))
  y <- lapply(c(3, 4), function(j) {
    m <- matrix(0L, n, j)
    for (i in seq_len(n)) if (z[i]) m[i, ] <- rbinom(j, 1, 0.4)
    m
  })
  names(y) <- c("s1", "s2")

  fit <- tobs(~ elev + bym2(graph = adj), data = data.frame(elev = d),
              family = int_occu(), detection = ~ 1, y = y,
              method = "nested_laplace",
              control = list(max.iter = 6L, verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$nested_laplace$multi_prior[[1]]$type, "bym2")
  # bym2 carries 2 components per spatial unit (combined + structured).
  expect_equal(length(fit$spatial_field), 2L * n)
})


test_that("dynamic occupancy fits nested_laplace with a spatial field on psi1", {
  set.seed(3)
  n_sites <- 36; n_seasons <- 3; J <- 3; adj <- chain_adj(n_sites)
  elev <- rnorm(n_sites)
  zmat <- matrix(0L, n_sites, n_seasons)
  zmat[, 1] <- rbinom(n_sites, 1, plogis(0.2 + 0.4 * elev))
  for (t in 2:n_seasons) {
    surv <- zmat[, t - 1] * (1 - rbinom(n_sites, 1, plogis(-1)))
    col  <- (1 - zmat[, t - 1]) * rbinom(n_sites, 1, plogis(-1.5))
    zmat[, t] <- surv + col
  }
  y <- array(NA_integer_, dim = c(n_sites, J, n_seasons))
  for (i in seq_len(n_sites)) for (t in seq_len(n_seasons))
    y[i, , t] <- if (zmat[i, t]) rbinom(J, 1, 0.4) else 0L

  fit <- tobs(~ elev + icar(graph = adj), data = data.frame(elev = elev),
              family = dyn_occu(), detection = ~ 1, y = y,
              col_formula = ~ 1, ext_formula = ~ 1, method = "nested_laplace",
              control = list(max.iter = 6L, verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$nested_laplace$multi_prior[[1]]$type, "icar")
  expect_equal(length(fit$spatial_field), n_sites)
})
