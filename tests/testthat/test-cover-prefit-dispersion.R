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

  # lognormal: the residual is taken on log(y), and the value comes back on the
  # FAMILY SURFACE -- an SD. The engine's variance conversion belongs to the
  # caller, because the joint fitter and the NUTS fitter pin different scales.
  y_ln <- exp(eta + rnorm(n, 0, 0.3))
  arm_ln <- list(y = y_ln, X = X, spatial_idx = g)
  phi_ln <- tulpaObs:::.occu_cover_prefit_dispersion(arm_ln, "lognormal", 99)
  expect_lt(abs(phi_ln - 0.3), 0.15)

  # gaussian: same residual, taken on the response itself, so a negative value
  # is legitimate and must not be logged.
  y_g <- eta + rnorm(n, 0, 0.3) - 5
  arm_g <- list(y = y_g, X = X, spatial_idx = g)
  phi_g <- tulpaObs:::.occu_cover_prefit_dispersion(arm_g, "gaussian", 99)
  expect_lt(abs(phi_g - 0.3), 0.15)

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
  expect_lt(abs(phi - 0.3), 0.2)

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

test_that("both engines pin the dispersion from one estimator", {
  # `.occu_cover_prefit_dispersion()` returns the FAMILY SURFACE, and the two
  # fitters convert differently: the joint path writes
  # `.cover_phi_sd_to_engine()`'s variance into `arm$phi`, the NUTS path carries
  # the surface value into `log_disp`. Returning the engine scale would hand
  # NUTS a variance where it wants an SD, and the two engines would pin
  # different dispersions for one model -- which is what
  # `test-occu-cover-spatial-nuts.R` (NUTS beta SDs against nested-Laplace SEs)
  # asserts against.
  set.seed(7)
  n <- 60L
  X <- cbind(1, rnorm(n))
  g <- rep(seq_len(20L), each = 3L)
  y <- exp(as.numeric(X %*% c(0.2, 0.4)) + rnorm(20L, 0, 0.5)[g] +
             rnorm(n, 0, 0.3))
  arm <- list(y = y, X = X, spatial_idx = g)

  surface <- tulpaObs:::.occu_cover_prefit_dispersion(arm, "lognormal", 99)
  # An SD, not a variance: a value near 0.3 rather than near 0.09.
  expect_lt(abs(surface - 0.3), 0.15)
  # The joint path's conversion squares it for the gaussian-family arm the
  # lognormal compiles to.
  engine <- tulpaObs:::.cover_phi_sd_to_engine(
    surface, tulpaObs:::.cover_pos_engine_family("lognormal"))
  expect_equal(engine, surface^2)
  # Beta states a precision on both sides, so its surface IS its engine value.
  expect_equal(tulpaObs:::.cover_phi_sd_to_engine(
    30, tulpaObs:::.cover_pos_engine_family("beta")), 30)
})


# A stated `phi.grid.pos` replaces the pre-fit. One node is a pin: the engine
# reads a scalar `phi_grid` entry as no axis and keeps the arm's `phi`, so the
# node has to be written there, and the reported dispersion is the one held.
# Several nodes are an axis whose stated nodes are all integrated, over the
# stated span.
.phi_pin_occu_cover_fit <- function(phi.grid.pos, ...) {
  N <- 20L; J <- 3L
  adj <- chain_adj(N)
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", adj = adj,
    beta_occ = c(0.2, 0.6), beta_p = c(0.4, -0.5),
    beta_pos = c(log(0.25), 0.3), sigma = 0.8, alpha = 1.0, sigma_pos = 0.4,
    seed = 707L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = adj),
    data = cbind(data.frame(site_id = seq_len(N)), sim$data),
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + share(spatial()),
    y = od$y, y_pos = y_pos, visits = od$det.covs, method = "nested_laplace",
    control = c(list(verbose = FALSE, progress = FALSE, engine = "joint",
                     phi.grid.pos = phi.grid.pos), list(...))))
}

test_that("a one-node phi.grid.pos is the dispersion occu_cover holds and reports", {
  skip_on_cran()
  held_of <- function(fit) tulpaObs:::.tobs_joint_fit(fit)$responses$pos$phi

  f4 <- .phi_pin_occu_cover_fit(0.4)
  expect_false("phi_pos" %in% colnames(tulpaObs:::.tobs_joint_fit(f4)$theta_grid))
  # The lognormal arm's engine `phi` is the residual variance.
  expect_equal(held_of(f4), 0.4^2, tolerance = 1e-12)
  expect_equal(f4$model$cover_pos_disp, 0.4, tolerance = 1e-12)
  expect_equal(f4$model$cover_pos_disp, sqrt(held_of(f4)), tolerance = 1e-12)

  # A different node is a different fit.
  f6 <- .phi_pin_occu_cover_fit(0.6)
  expect_equal(held_of(f6), 0.6^2, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(f4$means, f6$means)))
})

test_that("a multi-node phi.grid.pos is integrated on the stated nodes", {
  skip_on_cran()
  stated <- c(0.25, 0.35, 0.5)
  axis_of <- function(fit)
    sort(unique(tulpaObs:::.tobs_joint_fit(fit)$theta_grid[, "phi_pos"]))

  # With both refinement passes off the axis is exactly the stated nodes.
  fixed <- .phi_pin_occu_cover_fit(stated, adaptive.grid = FALSE,
                                   var.of.means.consistency = FALSE)
  expect_equal(axis_of(fixed), stated^2, tolerance = 1e-12)

  # An integrated dispersion is reported as its posterior, in `means`, and the
  # fit records no held value beside it.
  expect_true("phi_pos" %in% names(fixed$means))
  expect_null(fixed$model$cover_pos_disp)
  d <- tulpaObs:::.tobs_joint_draws(fixed, n = 4000L)$disp
  expect_true(all(d^2 %in% stated^2))
  expect_equal(mean(d), fixed$means[["phi_pos"]], tolerance = 0.05)

  # The refinement passes densify inside the stated span and drop no node.
  refined <- .phi_pin_occu_cover_fit(stated)
  expect_true(all(stated^2 %in% axis_of(refined)))
  expect_equal(range(axis_of(refined)), range(stated^2), tolerance = 1e-12)
})

test_that("a one-node phi.grid.pos is the dispersion occu_multiscale_cover holds", {
  skip_on_cran()
  sim <- simulate_occu_multiscale_cover(n_cells = 20L, plots_per_cell = 3L,
                                        visits_per_plot = 2L, phi = 0.4,
                                        sigma = 0.02, seed = 101L)
  fit_at <- function(v) suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot,
    positive = ~ x_cov + share(spatial(), alpha = grid(c(0, 0.5, 1, 2))),
    y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE,
                   sigma.grid = c(0.1, 0.5, 1), phi.grid.pos = v)))
  f4 <- fit_at(0.4)
  jf <- tulpaObs:::.tobs_joint_fit(f4)
  expect_false("phi_pos" %in% colnames(jf$theta_grid))
  expect_equal(jf$responses$pos$phi, 0.4^2, tolerance = 1e-12)
  expect_equal(f4$model$cover_pos_disp, 0.4, tolerance = 1e-12)

  # The held SD reaches the lognormal conditional mean exp(eta + sigma^2 / 2)
  # on both doors, fitted() and predict().
  expect_equal(tulpaObs:::.occu_mscale_cover_sigma_pos(f4), 0.4, tolerance = 1e-12)
  field <- as.numeric(f4$spatial_field)
  at_median <- tulpaObs:::.occu_mscale_cover_surface_at(
    f4$model, f4$means, field, unname(f4$means["alpha"]), 0)$cover
  expect_equal(tulpaObs:::.tobs_fitted_occu_multiscale_cover(f4)$cover,
               at_median * exp(0.4^2 / 2), tolerance = 1e-12)
  pr <- predict(f4, type = "cover")
  pr0 <- local({
    g <- f4; g$model$cover_pos_disp <- 0
    predict(g, type = "cover")
  })
  expect_equal(pr$mean, pr0$mean * exp(0.4^2 / 2), tolerance = 1e-10)

  # The default route holds its pre-fit, and reports that, not zero.
  f0 <- fit_at(NULL)
  jf0 <- tulpaObs:::.tobs_joint_fit(f0)
  expect_false("phi_pos" %in% colnames(jf0$theta_grid))
  expect_equal(tulpaObs:::.occu_mscale_cover_sigma_pos(f0),
               sqrt(jf0$responses$pos$phi), tolerance = 1e-12)
  expect_gt(tulpaObs:::.occu_mscale_cover_sigma_pos(f0), 0)
})

test_that("a one-node phi.grid is the dispersion cover() holds and reports", {
  skip_on_cran()
  set.seed(11)
  n_s <- 20L; N <- 300L
  adj <- chain_adj(n_s)
  x   <- rnorm(N)
  reg <- sample(n_s, N, replace = TRUE)
  occ <- rbinom(N, 1L, plogis(0.2 + 0.5 * x))
  y   <- ifelse(occ == 1L,
                pmin(exp(rnorm(N, log(0.2) + 0.3 * x, 0.4)), 1 - 1e-6), 0)
  d   <- data.frame(x = x, region = factor(reg, levels = seq_len(n_s)))
  fit_at <- function(v) suppressWarnings(tobs(
    formula = ~ x + icar(graph = adj, group_var = "region") +
      share(spatial(), alpha = grid(c(0.5, 1.0))),
    data = d, family = cover("lognormal"), y = y, method = "nested_laplace",
    control = list(sigma.grid = c(0.25, 0.5, 1.0), phi.grid = v)))

  f4 <- fit_at(0.4)
  jf <- tulpaObs:::.tobs_joint_fit(f4)
  expect_false("phi_pos" %in% colnames(jf$theta_grid))
  # The lognormal arm's engine `phi` is the residual variance.
  expect_equal(jf$responses$pos$phi, 0.4^2, tolerance = 1e-12)
  expect_equal(f4$sigma_pos, 0.4, tolerance = 1e-12)
  expect_equal(f4$sigma_pos_sd, 0)
  expect_equal(unique(tulpaObs:::.tobs_joint_draws(f4, n = 50L)$disp), 0.4,
               tolerance = 1e-12)

  f6 <- fit_at(0.6)
  expect_equal(tulpaObs:::.tobs_joint_fit(f6)$responses$pos$phi, 0.6^2,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(f4$beta_pos, f6$beta_pos)))
})

test_that("a one-node dispersion grid must be a positive number", {
  expect_null(tulpaObs:::.cover_phi_stated_pin(NULL))
  expect_null(tulpaObs:::.cover_phi_stated_pin(c(0.2, 0.4)))
  expect_identical(tulpaObs:::.cover_phi_stated_pin(0.4), 0.4)
  for (bad in list(0, -1, NA_real_, Inf, "x")) {
    expect_error(tulpaObs:::.cover_phi_stated_pin(bad), "one node pins")
  }
})
