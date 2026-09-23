# Draws from an outer-grid mixture (#368): the reported draw matrix of a
# grid-integrated fit is sampled from sum_k w_k N(m_k, C_k), not from the one
# Gaussian at the grid mean.

skew <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3

test_that("mixture draws keep the grid's two moments and its shape", {
  set.seed(1)
  w     <- c(0.7, 0.2, 0.1)
  modes <- rbind(c(0, 0), c(2, 1), c(5, 2))
  covs  <- list(diag(c(0.3, 0.2)), matrix(c(0.5, 0.1, 0.1, 0.4), 2),
                diag(c(1, 0.5)))
  mbar <- as.numeric(crossprod(w, modes))
  V    <- .tobs_grid_vcov(modes, w, covs, center = mbar)

  d <- .tobs_grid_mixture_draws(2e5, w, modes, covs)
  expect_equal(colMeans(d), mbar, tolerance = 0.02)
  expect_equal(stats::cov(d), V, tolerance = 0.02)
  # The right tail the third cell puts on coordinate 1 is kept; one Gaussian
  # at (mbar, V) has skew 0 by construction.
  expect_gt(skew(d[, 1]), 1)
})

test_that("trailing coordinates follow their conditional under V", {
  set.seed(2)
  w     <- c(0.6, 0.4)
  modes <- rbind(c(-1, 0), c(1.5, 1))
  covs  <- list(diag(2) * 0.2, diag(2) * 0.3)
  mbar  <- as.numeric(crossprod(w, modes))
  Vb    <- .tobs_grid_vcov(modes, w, covs, center = mbar)
  V <- matrix(0, 3, 3)
  V[1:2, 1:2] <- Vb
  V[3, 1:2] <- V[1:2, 3] <- c(0.3, -0.1)
  V[3, 3] <- 0.5
  means <- c(mbar, 2)

  d <- .tobs_grid_mixture_draws(2e5, w, modes, covs, means, V)
  expect_equal(ncol(d), 3L)
  expect_equal(colMeans(d), means, tolerance = 0.02)
  expect_equal(stats::cov(d), V, tolerance = 0.03)
})

test_that("a grid wider than `means` is cut to its leading coordinates", {
  set.seed(3)
  modes <- rbind(c(0, 0, 9), c(1, 1, 9))
  covs  <- list(diag(3) * 0.1, diag(3) * 0.1)
  d <- .tobs_grid_mixture_draws(100, c(0.5, 0.5), modes, covs, means = c(0.5, 0.5))
  expect_equal(dim(d), c(100L, 2L))
})

test_that("a cell with no within covariance is drawn at its mode", {
  set.seed(4)
  d <- .tobs_grid_mixture_draws(50, c(0, 1), rbind(0, 3), list(1, NULL))
  expect_true(all(d == 3))
})

test_that("no per-cell covariance falls to the Gaussian the fit describes", {
  set.seed(5)
  d <- .tobs_grid_mixture_draws(5e4, NULL, NULL, NULL, means = c(1, -1),
                                V = diag(c(0.5, 2)))
  expect_equal(colMeans(d), c(1, -1), tolerance = 0.03)
  expect_equal(diag(stats::cov(d)), c(0.5, 2), tolerance = 0.03)
})
