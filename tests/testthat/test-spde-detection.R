# Parameter recovery for an SPDE mesh field on the DETECTION arm of a
# single-season occupancy model (gcol33/tulpaObs#21). The state arm is
# non-spatial; the field sits on logit(p) and is fit through the SPDE-Laplace
# EM path. Asserts (a) detection fixed-effect recovery within tolerance and
# (b) detection field-shape recovery cor(u_hat, u_truth) > 0.7.
#
# Simulator: psi non-spatial (occ_cov slope 0.7), detection logit carries a
# smooth Matern-like signal `u_true` plus a det_cov slope (0.4). The recovered
# field at sites is A %*% fit$spatial_field_det; the occupancy field is NULL
# (the field is detection-only).

`%||%` <- function(a, b) if (is.null(a)) b else a

.sim_spde_det <- function(seed, n_sites = 600, J = 10,
                          beta_occ = c(0.2, 0.7), beta_det = c(-0.2, 0.4)) {
  set.seed(seed)
  coords <- cbind(runif(n_sites), runif(n_sites))
  u_true <- 0.9 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x_cov   <- rnorm(n_sites)
  det_cov <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(beta_occ[1] + beta_occ[2] * x_cov))
  p <- plogis(beta_det[1] + beta_det[2] * det_cov + u_true)
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
  list(
    data = data.frame(occ_cov = x_cov, det_cov = det_cov,
                      lon = coords[, 1], lat = coords[, 2]),
    y = y, u_true = u_true, coords = coords
  )
}

.fit_spde_det <- function(sim) {
  tobs(
    formula = ~ occ_cov,
    data = sim$data, family = occu(),
    detection = ~ det_cov + spde(lon, lat, max_edge = c(0.25, 0.5),
                                 nu = 1, prior_range = c(0.3, 0.5),
                                 prior_sigma = c(0.8, 0.5)),
    y = sim$y, method = "laplace", control = list(verbose = FALSE)
  )
}

# Single-seed recovery at N = 600, J = 10. Calibration over 6 seeds
# (dev_notes/probe_spde_detection.R): max |p_slope - 0.4| = 0.073, min
# field_at_sites cor = 0.716, mean cor = 0.754. The detection-arm slope is
# mildly attenuated by the field; the bound below sits well outside the
# observed maximum and catches structural regressions (the field not wired,
# a sign flip, the offset dropped).
test_that("tobs() + detection-arm spde() Laplace recovers beta and the field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  sim <- .sim_spde_det(seed = 1)
  fit <- .fit_spde_det(sim)

  # Detection fixed effect recovers within a generous Laplace-EM band.
  expect_lt(abs(fit$means["p_det_cov"] - 0.4), 0.20)
  # Occupancy slope on the non-spatial state arm stays near truth.
  expect_lt(abs(fit$means["psi_occ_cov"] - 0.7), 0.30)

  # The detection field lives on its own slot; the occupancy field is NULL.
  expect_false(is.null(fit$spatial_field_det))
  expect_true(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field_det), fit$spatial$n_units)

  # Detection fixed-effect SEs are real (not NA placeholders).
  expect_true(is.finite(fit$sds["p_det_cov"]))

  # Field-shape recovery: the recovered detection field correlates with truth.
  field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field_det)
  expect_gt(cor(field_at_sites, sim$u_true), 0.7)
})

# Aggregate over 3 seeds: average detection-slope bias and field correlation
# sit comfortably inside the per-seed band (calibration mean cor = 0.754).
test_that("detection-arm spde() recovery holds across seeds", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  seeds <- c(2L, 3L, 4L)
  res <- lapply(seeds, function(s) {
    sim <- .sim_spde_det(seed = s)
    fit <- .fit_spde_det(sim)
    field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field_det)
    list(p_bias = unname(fit$means["p_det_cov"]) - 0.4,
         field_cor = cor(field_at_sites, sim$u_true))
  })

  p_bias    <- vapply(res, `[[`, numeric(1), "p_bias")
  field_cor <- vapply(res, `[[`, numeric(1), "field_cor")

  for (k in seq_along(seeds)) {
    expect_lt(abs(p_bias[k]), 0.20)
    expect_gt(field_cor[k],  0.7)
  }
  expect_lt(mean(abs(p_bias)), 0.12)
  expect_gt(mean(field_cor),   0.72)
})

# gcol33/tulpaObs#218: the field was reportable (fit$spatial_field_det) but not
# predictive -- fitted()/predict() built the detection logit from the fixed
# effects alone. fitted(fit)$p must now equal the fixed-effect predictor plus
# the projected field, and predict() with no newdata delegates to fitted().
test_that("fitted() adds the detection-arm spde() field to p", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  sim <- .sim_spde_det(seed = 1)
  fit <- .fit_spde_det(sim)

  X_det    <- fit$model$X_processes[[2]]
  beta_det <- fit$means[fit$model$process_info[[1]]$p +
                          seq_len(fit$model$process_info[[2]]$p)]
  A        <- fit$spatial$tulpa_spec$A
  eta_manual <- as.vector(X_det %*% beta_det) +
    as.numeric(A %*% fit$spatial_field_det)

  f <- fitted(fit)
  expect_equal(unname(f$p), plogis(eta_manual))
  expect_gt(cor(f$p, sim$u_true), 0.5)

  # predict() with no newdata is the in-sample fit.
  expect_identical(predict(fit)$p, f$p)

  # A fit with no detection-arm field is unaffected (the offset is a no-op).
  fit_flat <- tobs(~ occ_cov, data = sim$data, family = occu(),
                   detection = ~ det_cov, y = sim$y, method = "laplace",
                   control = list(verbose = FALSE))
  expect_null(fit_flat$spatial_field_det)
  eta_flat <- as.vector(fit_flat$model$X_processes[[2]] %*%
    fit_flat$means[fit_flat$model$process_info[[1]]$p +
                    seq_len(fit_flat$model$process_info[[2]]$p)])
  expect_equal(unname(fitted(fit_flat)$p), plogis(eta_flat))
})
