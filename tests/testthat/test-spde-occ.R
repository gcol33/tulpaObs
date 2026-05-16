# Tests for SPDE spatial spec via tobs_spde()

test_that("tobs_spde creates valid specification", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  spec <- tobs_spde(coords)

  expect_s3_class(spec, "tobs_spatial")
  expect_equal(spec$type, "spde")
  expect_true(spec$n_units >= 30)
  expect_equal(spec$shared, c(TRUE, FALSE))
})

test_that("tobs_spde creates spec from formula", {
  set.seed(42)
  df <- data.frame(x = runif(20), y = runif(20))
  spec <- tobs_spde(~ x + y, data = df)

  expect_s3_class(spec, "tobs_spatial")
  expect_equal(spec$type, "spde")
})

test_that("tobs_spde print method works", {
  set.seed(42)
  spec <- tobs_spde(cbind(runif(20), runif(20)))
  expect_output(print(spec), "spde")
  expect_output(print(spec), "Matern")
})

test_that("tobs_spde with fractional nu creates valid spec", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  spec <- tobs_spde(coords, nu = 1.5)

  expect_equal(spec$nu, 1.5)
  expect_s3_class(spec, "tobs_spatial")
})

# Single-seed slim recovery at N = 400. Verifies the SPDE-Laplace EM path
# returns sensible (beta_occ, beta_det) and a positively-correlated mesh
# field. Generous tolerances — this tier catches structural regressions
# (sign flips, missing field, dropped prior).
test_that("tobs() + tobs_spde() Laplace recovers beta and the field shape", {
  skip_on_cran()
  skip_if_not_installed("tulpaMesh")

  set.seed(42)
  n_sites <- 400
  J <- 8
  coords <- cbind(runif(n_sites), runif(n_sites))

  # Smooth spatial signal on the unit square
  u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x_cov  <- rnorm(n_sites)
  beta_occ <- c(-0.5, 0.7)
  beta_det <- c(-0.3, 0.4)
  det_cov <- rnorm(n_sites)

  eta_occ <- beta_occ[1] + beta_occ[2] * x_cov + u_true
  z <- rbinom(n_sites, 1, plogis(eta_occ))
  eta_det <- beta_det[1] + beta_det[2] * det_cov
  p <- plogis(eta_det)

  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) {
    if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
  }
  dat <- data.frame(occ_cov = x_cov, det_cov = det_cov)

  sp <- tobs_spde(coords = coords, max_edge = c(0.3, 0.6), nu = 1,
                  prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))

  fit <- tobs(formula = ~ occ_cov, data = dat, family = occu(),
              detection = ~ det_cov, y = y, spatial = sp,
              engine = "laplace", control = list(verbose = FALSE))

  # Detection coefficients should be sensible
  expect_lt(abs(fit$means["p_det_cov"] - beta_det[2]), 1.0)

  # Slope on the occupancy covariate should be near truth (0.7) within
  # generous tolerance for Laplace EM at this data scale
  expect_lt(abs(fit$means["psi_occ_cov"] - beta_occ[2]), 0.5)

  # The latent SPDE mesh field should be stored on the fit
  expect_false(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field), sp$n_units)

  # And it should be positively correlated with the truth at sites.
  # Sanity threshold only — the SPDE-Laplace EM converges to a slightly
  # damped field at this n_sites / J / max_edge mix, sitting around 0.4
  # in practice and tipping above or below by a few hundredths per seed.
  field_at_sites <- as.numeric(sp$tulpa_spec$A %*% fit$spatial_field)
  expect_gt(cor(field_at_sites, u_true), 0.3)
})

# Higher-N tier: 3 seeds at N = 1500, J = 8 with a finer mesh. Asserts
# per-seed and aggregate bias bounds tighter than the N = 400 slim test.
# Calibration ran 5 seeds, max |psi_slope - 0.7| = 0.185, max |p_slope -
# 0.4| = 0.091, min field_at_sites cor = 0.43 (see
# dev_notes/_spde_highN_probe.R). Bounds below sit comfortably outside
# those observed maxima and detect regressions that the N = 400 slim
# tier is too noisy to catch.
test_that("tobs() + tobs_spde() Laplace recovery tightens at N = 1500", {
  skip_on_cran()
  skip_if_not_installed("tulpaMesh")

  n_sites <- 1500
  J <- 8
  beta_occ <- c(-0.5, 0.7)
  beta_det <- c(-0.3, 0.4)

  run_seed <- function(seed) {
    set.seed(seed)
    coords <- cbind(runif(n_sites), runif(n_sites))
    u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
    x_cov  <- rnorm(n_sites)
    det_cov <- rnorm(n_sites)
    eta_occ <- beta_occ[1] + beta_occ[2] * x_cov + u_true
    z <- rbinom(n_sites, 1, plogis(eta_occ))
    eta_det <- beta_det[1] + beta_det[2] * det_cov
    p <- plogis(eta_det)
    y <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) {
      if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
    }
    dat <- data.frame(occ_cov = x_cov, det_cov = det_cov)
    sp <- tobs_spde(coords = coords, max_edge = c(0.18, 0.45), nu = 1,
                    prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))
    fit <- tobs(formula = ~ occ_cov, data = dat, family = occu(),
                detection = ~ det_cov, y = y, spatial = sp,
                engine = "laplace", control = list(verbose = FALSE))
    field_at_sites <- as.numeric(sp$tulpa_spec$A %*% fit$spatial_field)
    list(
      psi_slope_bias = unname(fit$means["psi_occ_cov"]) - beta_occ[2],
      p_slope_bias   = unname(fit$means["p_det_cov"])   - beta_det[2],
      field_cor      = cor(field_at_sites, u_true)
    )
  }

  seeds <- c(1L, 2L, 3L)
  res <- lapply(seeds, run_seed)

  psi_bias  <- vapply(res, `[[`, numeric(1), "psi_slope_bias")
  p_bias    <- vapply(res, `[[`, numeric(1), "p_slope_bias")
  field_cor <- vapply(res, `[[`, numeric(1), "field_cor")

  # Per-seed bounds: tightened ~2x relative to the N = 400 slim test.
  for (k in seq_along(seeds)) {
    expect_lt(abs(psi_bias[k]), 0.25)
    expect_lt(abs(p_bias[k]),   0.20)
    expect_gt(field_cor[k],     0.35)
  }

  # Aggregate across seeds: average bias should be well inside the
  # per-seed band — calibration showed mean |psi_slope bias| ~= 0.086.
  expect_lt(mean(abs(psi_bias)), 0.15)
  expect_lt(mean(abs(p_bias)),   0.10)
  expect_gt(mean(field_cor),     0.40)
})
