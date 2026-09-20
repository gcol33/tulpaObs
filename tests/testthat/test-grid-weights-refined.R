# A refinement pass appends slice cells to an outer grid: new levels on one axis
# placed in ONE row of the others. Such a grid is not a tensor product, so the
# weights a post-processor rebuilds from `theta_grid` must measure it cell by
# cell from `refining_axis`. The expected measure below is written out from the
# fibre geometry rather than read from the engine.
#
# Base tensor: sigma in {0.25, 0.5, 1} (log scale, even log step log 2) by
# rho in {0.2, 0.5, 0.8} (even natural step). Each base cell owns 1/9 of the
# measure. Two sigma slices at rho = 0.5, sigma in {0.35, 0.7}, re-tile that one
# row along log sigma: the fibre's interior edges are the log-midpoints of its
# five nodes and its outer edges the base tensor's, half a log step past 0.25
# and 1. Rho is not refined, so each fibre cell's share is its log-width over
# the row's total log-width 3 log 2, times the row's 1/3.

.gwr_grid <- function() {
  base <- as.matrix(expand.grid(sigma = c(0.25, 0.5, 1), rho = c(0.2, 0.5, 0.8)))
  slices <- cbind(sigma = c(0.35, 0.7), rho = 0.5)
  list(theta_grid = rbind(base, slices),
       refining_axis = c(rep("", nrow(base)), "sigma", "consistency_sigma"))
}

.gwr_expected_measure <- function(tg) {
  u <- log(c(0.25, 0.35, 0.5, 0.7, 1))
  h <- log(2) / 2
  edges <- c(u[1L] - h, (u[-1L] + u[-length(u)]) / 2, u[length(u)] + h)
  fibre <- diff(edges) / (3 * log(2)) / 3
  out <- rep(1 / 9, nrow(tg))
  row <- which(tg[, "rho"] == 0.5)
  out[row] <- fibre[match(round(log(tg[row, "sigma"]), 12), round(u, 12))]
  out
}

test_that("a refined grid's rebuilt weights follow its cell-by-cell measure", {
  g  <- .gwr_grid()
  tg <- g$theta_grid
  lm <- c(-1.2, -0.4, -0.9, -0.3, 0, -0.6, -1.5, -0.8, -1.1, -0.1, -0.2)
  fit <- c(g, list(log_marginal = lm))

  mu <- .gwr_expected_measure(tg)
  expect_equal(sum(mu), 1, tolerance = 1e-12)
  want <- exp(lm) * mu; want <- want / sum(want)

  got <- tulpaObs:::.tobs_grid_weights(fit)
  expect_equal(got, want, tolerance = 1e-10)

  # The tensor-product rule on the same cells gives every sigma level the width
  # of a whole column, which is a different answer on this grid.
  tensor <- tulpa::tulpa_normalise_weights_safe(
    lm, "outer grid", log_quad = tulpa::tulpa_grid_log_quad(tg))
  expect_gt(max(abs(tensor - want)), 1e-3)

  # A caller-supplied log-marginal is weighted by the same measure.
  lm2 <- lm + tg[, "rho"]
  want2 <- exp(lm2) * mu; want2 <- want2 / sum(want2)
  expect_equal(tulpaObs:::.tobs_grid_weights(fit, log_marginal = lm2), want2,
               tolerance = 1e-10)
})

test_that("a grid without slice cells keeps the tensor-product weights", {
  g  <- .gwr_grid()
  tg <- g$theta_grid[1:9, , drop = FALSE]
  lm <- seq(-1, 0, length.out = 9L)
  tensor <- tulpa::tulpa_normalise_weights_safe(
    lm, "outer grid", log_quad = tulpa::tulpa_grid_log_quad(tg))
  expect_equal(tulpaObs:::.tobs_grid_weights(list(theta_grid = tg,
                                                  log_marginal = lm)),
               tensor, tolerance = 1e-12)
  expect_equal(tulpaObs:::.tobs_grid_weights(list(theta_grid = tg,
                                                  log_marginal = lm,
                                                  refining_axis = rep("", 9L))),
               tensor, tolerance = 1e-12)
})

test_that("a log-marginal that does not match the grid is an error", {
  g <- .gwr_grid()
  fit <- c(g, list(log_marginal = rep(0, nrow(g$theta_grid) - 1L)))
  expect_error(tulpaObs:::.tobs_grid_weights(fit), "log_quad")
  fit2 <- list(theta_grid = g$theta_grid, log_marginal = rep(0, 11L),
               refining_axis = g$refining_axis[-1L])
  expect_error(tulpaObs:::.tobs_grid_weights(fit2), "refining")
})

test_that("a stored cell measure is the one the weights take", {
  g  <- .gwr_grid()
  lm <- seq(-1, 0, length.out = nrow(g$theta_grid))
  lq <- log(rev(seq_len(nrow(g$theta_grid))))
  fit <- c(g, list(log_marginal = lm, log_quad = lq))
  want <- exp(lm + lq); want <- want / sum(want)
  expect_equal(tulpaObs:::.tobs_grid_weights(fit), want, tolerance = 1e-12)
})

test_that("converged-cell weights rebuilt without engine weights keep the measure", {
  g  <- .gwr_grid()
  tg <- g$theta_grid
  lm <- c(-1.2, -0.4, -0.9, -0.3, 0, -0.6, -1.5, NaN, -1.1, -0.1, -0.2)
  mu <- .gwr_expected_measure(tg)
  ok <- which(is.finite(lm))
  want <- exp(lm[ok]) * mu[ok]; want <- want / sum(want)

  for (ew in list(NULL, rep(NA_real_, length(lm)), numeric(3))) {
    fit <- c(g, list(log_marginal = lm, weights = ew))
    got <- suppressWarnings(tulpaObs:::.tobs_joint_ok_cells(fit, "fixture"))
    expect_identical(got$ok_cells, ok)
    expect_equal(got$w, want, tolerance = 1e-10)
  }
})

test_that("the areal product grid is measured on its named hyperparameters", {
  blk <- function(sg, rho) list(
    cells = unlist(lapply(rho, function(r) lapply(sg, function(s)
      list(sigma = s, rho = r))), recursive = FALSE),
    to_hyper = function(cell) c(sigma = cell$sigma, rho = cell$rho))
  b1 <- blk(c(0.5, 1, 2), c(0.05, 0.3, 0.95))
  cg <- tulpaObs:::.tobs_block_cell_product(list(b1))
  M  <- tulpaObs:::.tobs_block_hyper_grid(list(b1), cg)
  expect_identical(colnames(M), c("sigma", "rho"))
  expect_equal(unname(M[, "sigma"]), rep(c(0.5, 1, 2), 3L))
  expect_equal(unname(M[, "rho"]), rep(c(0.05, 0.3, 0.95), each = 3L))

  lm <- numeric(nrow(M))
  w  <- tulpaObs:::.tobs_grid_weights(list(theta_grid = M), log_marginal = lm)
  lq <- tulpa::tulpa_grid_log_quad(M)
  expect_equal(w, exp(lq) / sum(exp(lq)), tolerance = 1e-12)
  expect_gt(max(abs(w - 1 / nrow(M))), 1e-3)

  b2 <- list(cells = list(list(tau = 1), list(tau = 4)),
             to_hyper = function(cell) c(tau = cell$tau))
  cg2 <- tulpaObs:::.tobs_block_cell_product(list(b1, b2))
  M2  <- tulpaObs:::.tobs_block_hyper_grid(list(b1, b2), cg2)
  expect_identical(colnames(M2), c("b1.sigma", "b1.rho", "b2.tau"))
  expect_identical(nrow(M2), 18L)
})

test_that("the sampled span read off a refined grid is the span its cell measure covers", {
  g <- .gwr_grid()
  warm_of <- function(tg, refining) list(joint_fit = list(
    theta_grid = tg, refining_axis = refining, axis_support = NULL))
  base <- g$theta_grid[1:9, , drop = FALSE]

  refined <- tulpaObs:::.occu_cover_nuts_hyper_bounds(
    warm_of(g$theta_grid, g$refining_axis), "bym2")
  plain <- tulpaObs:::.occu_cover_nuts_hyper_bounds(
    warm_of(base, rep("", 9L)), "bym2")

  # Slices inside the base span leave the span where the base tensor put it:
  # half a log step past 0.25 and 1 on sigma.
  expect_equal(as.numeric(refined$sigma), as.numeric(plain$sigma),
               tolerance = 1e-12)
  expect_equal(as.numeric(plain$sigma), c(0.25, 1) * c(2^-0.5, 2^0.5),
               tolerance = 1e-12)
  expect_equal(as.numeric(attr(refined$sigma, "nodes")), c(0.25, 1))
})
