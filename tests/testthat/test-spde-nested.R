# Recovery + regression tests for a continuous SPDE field on the multi-block
# EM nested-Laplace occupancy path (method = "nested_laplace"). Closes the last
# arm of gcol33/tulpaObs#21: the non-joint multi-block engine now carries the
# SPDE mesh projection A (a many-to-one site -> mesh-node map) via the same
# make_spde_block factory + SpdeQBuilder the single-Laplace SPDE path uses
# (tulpa LatentBlock INDEXED_MULTI), so an spde() field fits end to end.
#
# Meshes are built with cutoff = 0 + the default max_edge (reliable across
# seeds; works around the tulpaMesh zero-triangle collapse, gcol33/tulpaMesh#3).

simulate_spde_occu <- function(seed, n_sites, J = 8,
                               beta_occ = c(-0.5, 0.7),
                               beta_det = c(-0.3, 0.4)) {
  set.seed(seed)
  coords <- cbind(runif(n_sites), runif(n_sites))
  # Smooth Matern-like signal on the unit square.
  u_true  <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x_cov   <- rnorm(n_sites)
  det_cov <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(beta_occ[1] + beta_occ[2] * x_cov + u_true))
  p <- plogis(beta_det[1] + beta_det[2] * det_cov)
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
  list(
    data = data.frame(occ_cov = x_cov, det_cov = det_cov,
                      lon = coords[, 1], lat = coords[, 2]),
    y = y, u_true = u_true
  )
}

fit_spde_nested <- function(d) {
  tobs(
    formula = ~ occ_cov + spde(lon, lat, cutoff = 0, nu = 1,
                               prior_range = c(0.3, 0.5),
                               prior_sigma = c(0.7, 0.5)),
    data = d$data, family = occu(),
    detection = ~ det_cov, y = d$y,
    method = "nested_laplace", control = list(verbose = FALSE)
  )
}

# A continuous spde() field on method = "nested_laplace" used to be an honest
# stop(); it now fits and stores the mesh field.
test_that("spde() + nested_laplace builds the multi-block SPDE prior and fits", {
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")
  d <- simulate_spde_occu(1L, n_sites = 120, J = 6)
  fit <- fit_spde_nested(d)

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$nested_laplace$multi_prior[[1]]$type, "spde")
  # The stored field is the n_mesh mesh realization.
  expect_false(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field), fit$spatial$tulpa_spec$n_mesh)
})

# Recovery: occupancy fixed effect + field shape via cor(A u_hat, u_truth).
# Calibration (dev_notes/probe_spde_nested.R, 4 seeds, N = 400, J = 8,
# cutoff = 0 fine mesh, mode-centred (range, sigma) grid): occupancy slope
# max |bias| ~= 0.25, detection slope max |bias| ~= 0.11, field cor
# min ~= 0.45 / mean ~= 0.49.
#
# Field-shape ceiling. With ONE latent occupancy state per site (a single
# binary outcome per site informing the field), the nested path -- which
# integrates the Matern (range, sigma) over a grid -- recovers the field at
# cor ~= 0.4-0.5. This is the documented single-latent-state-per-site ceiling
# the shipped single-Laplace state-arm tests accept (test-spde-occ.R uses
# cor > 0.3 / > 0.35). The single-Laplace SPDE path reaches cor ~= 0.65 on
# identical data because it fits the field at the FIXED prior-mode (no
# hyperparameter integration); the nested path pays a field-sharpness cost to
# integrate (range, sigma) uncertainty, so it sits at the state-arm ceiling,
# well above noise and far enough to catch a field-not-wired / sign-flip /
# dropped-projection regression. See dev_notes/spde_nested_field_ceiling.md.
test_that("spde() + nested_laplace recovers beta_occ and the field shape", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  beta_occ <- c(-0.5, 0.7)
  beta_det <- c(-0.3, 0.4)
  seeds <- 1:3

  res <- lapply(seeds, function(s) {
    d <- simulate_spde_occu(s, n_sites = 400, J = 8,
                            beta_occ = beta_occ, beta_det = beta_det)
    fit <- fit_spde_nested(d)
    field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field)
    list(
      psi_bias  = unname(fit$means["psi_occ_cov"]) - beta_occ[2],
      p_bias    = unname(fit$means["p_det_cov"])   - beta_det[2],
      field_cor = cor(field_at_sites, d$u_true)
    )
  })

  psi_bias  <- vapply(res, `[[`, numeric(1), "psi_bias")
  p_bias    <- vapply(res, `[[`, numeric(1), "p_bias")
  field_cor <- vapply(res, `[[`, numeric(1), "field_cor")

  # (a) Occupancy / detection fixed-effect recovery, per seed.
  for (k in seq_along(seeds)) {
    expect_lt(abs(psi_bias[k]), 0.30)
    expect_lt(abs(p_bias[k]),   0.20)
    # (b) Field-shape recovery at the single-latent-state-per-site ceiling.
    expect_gt(field_cor[k], 0.35)
  }

  # Aggregate bounds, comfortably inside the per-seed bands.
  expect_lt(mean(abs(psi_bias)), 0.20)
  expect_lt(mean(abs(p_bias)),   0.12)
  expect_gt(mean(field_cor),     0.40)
})

# Regression: the areal nested-Laplace path (icar / bym2 / car_proper) must
# still fit after the INDEXED_MULTI generalisation of the multi-block engine.
test_that("areal nested_laplace path still fits (regression)", {
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")
  set.seed(7)
  n_sites <- 24
  adj <- matrix(0, n_sites, n_sites)
  for (i in seq_len(n_sites - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  x <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(0.2 + 0.5 * x))
  y <- matrix(0L, n_sites, 4)
  for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(4, 1, 0.5)
  dat <- data.frame(x = x)

  fit_icar <- tobs(~ x + icar(graph = adj), data = dat, family = occu(),
                   detection = ~ 1, y = y, method = "nested_laplace",
                   control = list(max.iter = 5L, verbose = FALSE))
  expect_s3_class(fit_icar, "tobs_fit")
  expect_identical(fit_icar$nested_laplace$multi_prior[[1]]$type, "icar")

  fit_bym2 <- tobs(~ x + bym2(graph = adj), data = dat, family = occu(),
                   detection = ~ 1, y = y, method = "nested_laplace",
                   control = list(max.iter = 5L, verbose = FALSE))
  expect_s3_class(fit_bym2, "tobs_fit")
  expect_identical(fit_bym2$nested_laplace$multi_prior[[1]]$type, "bym2")
})
