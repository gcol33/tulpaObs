# INLA-style NA-response prediction on the nested-Laplace path. A single-season
# site whose detection history is all-missing (all NA) is held out of the
# likelihood (n_trials = 0) but kept in the latent field, so its occupancy is
# interpolated from neighbours. `predict(type = "state")` returns the
# marginalised per-site psi posterior with the held-out rows flagged.
#
# This is a recovery test, not a smoke test: it simulates a known smooth spatial
# field on a 2-D grid, holds out a SCATTERED set of sites (so each has observed
# neighbours -- the realistic interpolation case, not extrapolation across a
# contiguous gap), and checks the held-out predictions track the truth.

# Rook-adjacency grid graph + a smooth field.
.make_grid_occu <- function(gx = 9, gy = 9, J = 8, p_det = 0.5, seed = 11) {
  set.seed(seed)
  n <- gx * gy
  coord <- expand.grid(cx = seq_len(gx), cy = seq_len(gy))
  adj <- matrix(0, n, n)
  idx_of <- function(i, j) (j - 1) * gx + i
  for (i in seq_len(gx)) for (j in seq_len(gy)) {
    a <- idx_of(i, j)
    if (i < gx) { b <- idx_of(i + 1, j); adj[a, b] <- 1; adj[b, a] <- 1 }
    if (j < gy) { b <- idx_of(i, j + 1); adj[a, b] <- 1; adj[b, a] <- 1 }
  }
  u <- 0.9 * scale(coord$cx)[, 1] + 0.7 * scale(coord$cy)[, 1] +
       1.2 * exp(-((coord$cx - 5)^2 + (coord$cy - 5)^2) / 6)
  u <- u - mean(u)
  psi <- plogis(u)
  z <- rbinom(n, 1, psi)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i]) y[i, ] <- rbinom(J, 1, p_det)
  list(adj = adj, y = y, psi = psi, n = n)
}


test_that("all-NA sites are flagged held-out and dropped from the likelihood", {
  d <- .make_grid_occu()
  heldout <- seq(2, d$n, by = 4)
  y <- d$y; y[heldout, ] <- NA

  fit <- tobs(~ 1 + icar(graph = d$adj), data = data.frame(s = seq_len(d$n)),
              family = occu(), detection = ~ 1, y = y, method = "nested_laplace",
              control = list(max.iter = 25L, tol = 1e-4, verbose = FALSE))

  sp <- predict(fit, type = "state")
  expect_named(sp, c("row", "psi", "psi_lower", "psi_upper", "heldout"))
  expect_equal(nrow(sp), d$n)
  expect_identical(which(sp$heldout), as.integer(heldout))
  expect_true(all(sp$psi >= 0 & sp$psi <= 1))
  # The driver recorded exactly these rows as held out (n_trials = 0).
  expect_identical(fit$nested_laplace$heldout, as.integer(heldout))
})


test_that("held-out occupancy is recovered by spatial interpolation (icar)", {
  d <- .make_grid_occu()
  heldout <- seq(2, d$n, by = 4)
  y <- d$y; y[heldout, ] <- NA

  fit <- tobs(~ 1 + icar(graph = d$adj), data = data.frame(s = seq_len(d$n)),
              family = occu(), detection = ~ 1, y = y, method = "nested_laplace",
              control = list(max.iter = 30L, tol = 1e-5, verbose = FALSE))

  ho <- predict(fit, type = "state")
  ho <- ho[ho$heldout, ]
  truth <- d$psi[heldout]

  # Predictions track the held-out truth (interpolated from observed neighbours).
  expect_gt(cor(ho$psi, truth), 0.6)
  expect_lt(mean(abs(ho$psi - truth)), 0.35)
})


test_that("held-out occupancy is recovered for a bym2 field too", {
  # bym2's predictor mixes structured + unstructured components with
  # hyperparameter-dependent scales, so it cannot be reconstructed from the
  # modes alone; prediction reads the engine's per-cell fitted eta
  # (tulpa_nested_laplace()$fitted_eta), which is exact for every prior.
  d <- .make_grid_occu()
  heldout <- seq(2, d$n, by = 4)
  y <- d$y; y[heldout, ] <- NA

  fit <- tobs(~ 1 + bym2(graph = d$adj), data = data.frame(s = seq_len(d$n)),
              family = occu(), detection = ~ 1, y = y, method = "nested_laplace",
              control = list(max.iter = 30L, tol = 1e-5, verbose = FALSE))

  sp <- predict(fit, type = "state")
  expect_named(sp, c("row", "psi", "psi_lower", "psi_upper", "heldout"))
  expect_identical(which(sp$heldout), as.integer(heldout))
  ho <- sp[sp$heldout, ]
  truth <- d$psi[heldout]
  expect_gt(cor(ho$psi, truth), 0.6)
  expect_lt(mean(abs(ho$psi - truth)), 0.35)
})


test_that("held-out 95% intervals are calibrated (coverage recovery)", {
  # Calibration recovery test (not smoke): the marginalized-occupancy state pass
  # fits D_i = 1{>=1 detection} ~ Bernoulli(q_i * sigma(eta)) with the ICAR field
  # and integrates over the tau grid, so the held-out psi interval should cover
  # the true psi at ~the nominal rate. Pool over seeds for a stable estimate.
  skip_on_cran()
  covered <- logical(0)
  widths  <- numeric(0)
  for (s in 1:5) {
    d <- .make_grid_occu(seed = 40 + s)
    heldout <- seq(2, d$n, by = 4)
    y <- d$y; y[heldout, ] <- NA
    fit <- tobs(~ 1 + icar(graph = d$adj), data = data.frame(s = seq_len(d$n)),
                family = occu(), detection = ~ 1, y = y,
                method = "nested_laplace",
                control = list(max.iter = 30L, tol = 1e-5, verbose = FALSE))
    sp <- predict(fit, type = "state")
    ho <- sp[sp$heldout, ]
    truth <- d$psi[heldout]
    covered <- c(covered, truth >= ho$psi_lower & truth <= ho$psi_upper)
    widths  <- c(widths, ho$psi_upper - ho$psi_lower)
  }
  # Nominal 95%; allow slack for the Laplace + grid approximation and finite
  # pooled N. Intervals are calibrated-to-slightly-conservative by design.
  expect_gt(mean(covered), 0.85)
  # Intervals must be informative, not the trivial [0, 1].
  expect_lt(median(widths), 0.95)
})


test_that("predict(type = 'state') errors for a non-nested fit", {
  d <- .make_grid_occu(seed = 3)
  fit_lap <- tobs(~ 1, data = data.frame(s = seq_len(d$n)), family = occu(),
                  detection = ~ 1, y = d$y, method = "laplace",
                  control = list(verbose = FALSE))
  expect_error(predict(fit_lap, type = "state"), "nested_laplace")
})
