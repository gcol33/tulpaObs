# The positive arm's dispersion is pre-fit from RESIDUALS, not from the marginal
# spread of the response.
#
# A dispersion is the spread about the arm's own predictor. The marginal spread
# is not: it also carries whatever the arm's covariates and its shared field
# explain, so pre-fitting at `var(y)` sets the dispersion above the truth by
# exactly the part of the model that works. On the coupled path that value is
# PINNED -- it is the arm's `phi`, not an axis -- so the whole error lands in the
# fit. Measured on 12 seeds of `simulate_occu_cover("lognormal")` at a true
# dispersion SD of 0.4 with an ICAR field (alpha 0.6, sigma 0.8): the marginal
# spread reads 0.46 to 0.90, the residual about X + cell reads 0.34 to 0.51.

test_that("residual variance drops what the design explains", {
  set.seed(1)
  n <- 60L
  g <- rep(seq_len(20L), each = 3L)
  z <- 1 + 0.5 * rnorm(n) + rnorm(20L, 0, 0.8)[g] + rnorm(n, 0, 0.4)
  X <- cbind(1, seq_len(n) / n)

  marginal <- stats::var(z)
  v_x      <- tulpaObs:::.tobs_prefit_resid_var(z, X)
  v_xg     <- tulpaObs:::.tobs_prefit_resid_var(z, X, group = g)

  # Each nested model can only take variance away.
  expect_lt(v_x, marginal)
  expect_lt(v_xg, v_x)
  # Absorbing the group is what reaches the residual scale, since the group
  # effect is not a column of X.
  expect_lt(sqrt(v_xg), 0.6)
})

test_that("a group that cannot be absorbed falls back to the covariates", {
  set.seed(2)
  n <- 40L
  z <- rnorm(n); X <- cbind(1, rnorm(n))
  v_x <- tulpaObs:::.tobs_prefit_resid_var(z, X)

  # One row per group: the factor is saturated, absorbs the response exactly and
  # leaves nothing to estimate from.
  expect_equal(tulpaObs:::.tobs_prefit_resid_var(z, X, group = seq_len(n)), v_x)
  # A group of the wrong length is not a group.
  expect_equal(tulpaObs:::.tobs_prefit_resid_var(z, X, group = 1:3), v_x)
  expect_equal(tulpaObs:::.tobs_prefit_resid_var(z, X, group = NULL), v_x)
})

test_that("a design with nothing left to estimate from declines", {
  # `min_df` residual degrees of freedom or it is not a measurement.
  expect_true(is.na(tulpaObs:::.tobs_prefit_resid_var(c(1, 2), cbind(1, c(0, 1)))))
  expect_true(is.na(tulpaObs:::.tobs_prefit_resid_var(1, cbind(1))))
  # Non-finite rows are dropped rather than poisoning the fit.
  set.seed(3)
  z <- c(rnorm(30), NA); X <- cbind(1, rnorm(31))
  v <- tulpaObs:::.tobs_prefit_resid_var(z, X)
  expect_true(is.finite(v) && v > 0)
})

test_that("each cover family maps the residual to its own dispersion scale", {
  set.seed(4)
  n <- 60L
  g <- rep(seq_len(20L), each = 3L)
  X <- cbind(1, rnorm(n))
  eta <- as.numeric(X %*% c(0.2, 0.4)) + rnorm(20L, 0, 0.5)[g]

  # lognormal: the residual is taken on log(y) and the ENGINE carries a variance.
  y_ln <- exp(eta + rnorm(n, 0, 0.3))
  arm_ln <- list(y = y_ln, X = X, spatial_idx = g)
  phi_ln <- tulpaObs:::.occu_cover_prefit_dispersion(arm_ln, "lognormal", 99)
  expect_lt(abs(sqrt(phi_ln) - 0.3), 0.15)

  # gaussian: same residual, taken on the response itself, so a negative value
  # is legitimate and must not be logged.
  y_g <- eta + rnorm(n, 0, 0.3) - 5
  arm_g <- list(y = y_g, X = X, spatial_idx = g)
  phi_g <- tulpaObs:::.occu_cover_prefit_dispersion(arm_g, "gaussian", 99)
  expect_lt(abs(sqrt(phi_g) - 0.3), 0.15)

  # beta: the same moment match the marginal estimator used, with the RESIDUAL
  # variance in place of the marginal one -- so a precision, not a variance, and
  # larger than what the marginal spread implies.
  y_b <- plogis(eta + rnorm(n, 0, 0.3))
  arm_b <- list(y = y_b, X = X, spatial_idx = g)
  phi_b <- tulpaObs:::.occu_cover_prefit_dispersion(arm_b, "beta", 99)
  mu <- mean(y_b)
  expect_gt(phi_b, mu * (1 - mu) / stats::var(y_b) - 1)
  expect_gte(phi_b, 1)
})

test_that("the pre-fit reads only the rows the cover density scores", {
  # Per-visit, the pos arm carries EVERY valid visit; the cell-coupling spec
  # gates `f_pos` on detection, so a non-detected visit's `y_pos` is a
  # placeholder no term reads -- and for the lognormal arm that placeholder is
  # zero, whose log is not a number.
  set.seed(5)
  n <- 60L
  X <- cbind(1, rnorm(n))
  y <- exp(as.numeric(X %*% c(0.2, 0.4)) + rnorm(n, 0, 0.3))
  scored <- rep(c(TRUE, FALSE), length.out = n)
  y[!scored] <- 0                                  # the placeholder
  arm <- list(y = y, X = X, spatial_idx = rep(seq_len(20L), each = 3L))

  # Unmasked, every row is unusable and the caller's fallback stands.
  expect_identical(
    tulpaObs:::.occu_cover_prefit_dispersion(arm, "lognormal", 99), 99)
  # Masked, it recovers the dispersion.
  phi <- tulpaObs:::.occu_cover_prefit_dispersion(arm, "lognormal", 99,
                                                  scored = scored)
  expect_lt(abs(sqrt(phi) - 0.3), 0.2)

  # A mask of the wrong length is ignored rather than silently mis-aligning the
  # rows against the design.
  expect_identical(
    tulpaObs:::.occu_cover_prefit_dispersion(arm, "lognormal", 99,
                                             scored = c(TRUE, FALSE)), 99)
})

test_that("cover()'s own prefit is the same estimator, without the group", {
  # `.prefit_lognormal_sigma()` delegates, so the two paths cannot drift into
  # two different answers. It passes no group deliberately: cover() centres a
  # 7-node `phi.grid` on this value rather than pinning at it.
  set.seed(6)
  n <- 50L
  X <- cbind(1, rnorm(n))
  y <- as.numeric(X %*% c(1, 0.5)) + rnorm(n, 0, 0.7)
  enc <- list(pos_data = list(y = y, X = X))

  got <- tulpaObs:::.prefit_lognormal_sigma(enc, list())
  # The closed form it replaced, on a full-rank design.
  b   <- qr.solve(X, y)
  want <- sqrt(sum((y - X %*% b)^2) / (n - ncol(X)))
  expect_equal(got, want)

  # Degenerate input keeps the documented 1.0.
  expect_equal(tulpaObs:::.prefit_lognormal_sigma(
    list(pos_data = list(y = 1, X = cbind(1))), list()), 1.0)
})
