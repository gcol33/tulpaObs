# Recovery tests for the continuous-mesh (SPDE) field on the N-mixture fitters
# (gcol33/tulpaObs#21). A KNOWN smooth Matern-like field is simulated on the
# abundance arm and recovered through the spde() path; the gate asserts (a)
# abundance fixed-effect recovery and (b) field-shape recovery
# cor(A u_hat, u_truth) > 0.6 at the sites.
#
# The triangulation collapses to zero triangles for some (max_edge, cutoff)
# settings at larger point counts (upstream gcol33/tulpaMesh#3); `cutoff = 0`
# with the default max_edge triangulates reliably across seeds and N, so the
# tests build the mesh that way.

# --- single-species N-mixture abundance (R/abun.R) -------------------------

test_that("abun() + spde() recovers the abundance slope and the mesh field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  beta_lambda <- c(0.6, 0.5)     # log-lambda intercept + slope
  beta_p      <- c(0.2, -0.3)    # logit-p intercept + slope

  run_seed <- function(seed, n_sites = 300, J = 6) {
    set.seed(seed)
    coords <- cbind(runif(n_sites), runif(n_sites))
    u_true <- 0.9 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
    x_cov  <- rnorm(n_sites)
    det_cov <- rnorm(n_sites)
    lambda <- exp(beta_lambda[1] + beta_lambda[2] * x_cov + u_true)
    N <- rpois(n_sites, lambda)
    p <- plogis(beta_p[1] + beta_p[2] * det_cov)
    y <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) y[i, ] <- rbinom(J, N[i], p[i])
    dat <- data.frame(abun_cov = x_cov, det_cov = det_cov,
                      lon = coords[, 1], lat = coords[, 2])
    fit <- tobs(
      formula = ~ abun_cov + spde(lon, lat, cutoff = 0, nu = 1,
                                  prior_range = c(0.3, 0.5),
                                  prior_sigma = c(0.8, 0.5)),
      data = dat, family = abun(),
      detection = ~ det_cov, y = y,
      method = "nested_laplace", control = list(verbose = FALSE))
    field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field)
    list(
      lam_slope = unname(fit$means["lambda_abun_cov"]),
      p_slope   = unname(fit$means["p_det_cov"]),
      field_cor = cor(field_at_sites, u_true),
      n_units   = length(fit$spatial_field),
      mesh_n    = fit$spatial$tulpa_spec$n_mesh)
  }

  seeds <- c(1L, 2L, 3L)
  res <- lapply(seeds, run_seed)

  lam_slope <- vapply(res, `[[`, numeric(1), "lam_slope")
  p_slope   <- vapply(res, `[[`, numeric(1), "p_slope")
  field_cor <- vapply(res, `[[`, numeric(1), "field_cor")

  # The mesh field is stored at the FEM nodes; A projects it to sites.
  for (k in seq_along(seeds)) {
    expect_equal(res[[k]]$n_units, res[[k]]$mesh_n)
    # Calibration (dev_notes/probe_spde_nmix_abun.R, 5 seeds at N = 300):
    # max |lam_slope - 0.5| = 0.082, max |p_slope + 0.3| = 0.059, min
    # field_cor = 0.927. Bounds sit outside those maxima.
    expect_lt(abs(lam_slope[k] - beta_lambda[2]), 0.20)
    expect_lt(abs(p_slope[k]   - beta_p[2]),      0.15)
    expect_gt(field_cor[k], 0.6)
  }
  expect_lt(mean(abs(lam_slope - beta_lambda[2])), 0.15)
  expect_gt(mean(field_cor), 0.85)
})


# --- community / multispecies N-mixture (R/ms_abun.R) ----------------------

test_that("ms_abun() + spde() recovers the shared mesh field across species", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  # Single seed: the per-species joint EM over the (range, sigma) grid is the
  # cost driver, so the scale (N = 120, 4 species) is kept modest. The shared
  # field is the SPDE-specific gate (cor > 0.6); the community-mean abundance
  # slope is weakly identified at small species counts, so its tolerance is
  # generous (calibration: lam_slope in 0.25-0.62, field_cor 0.96-0.98 over
  # the verified seeds; see dev_notes/probe_spde_ms_abun.R).
  set.seed(2)
  n_sites <- 120; n_species <- 4; J <- 5
  coords <- cbind(runif(n_sites), runif(n_sites))
  u_true <- 0.9 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x_cov  <- rnorm(n_sites); det_cov <- rnorm(n_sites)
  mu_lambda <- c(0.5, 0.4); mu_p <- c(0.2, -0.3)
  blam <- cbind(rnorm(n_species, 0, 0.3), rnorm(n_species, 0, 0.3))
  bp   <- cbind(rnorm(n_species, 0, 0.3), rnorm(n_species, 0, 0.3))
  y <- array(0L, dim = c(n_sites, J, n_species))
  for (s in seq_len(n_species)) {
    lambda <- exp(mu_lambda[1] + blam[s, 1] +
                  (mu_lambda[2] + blam[s, 2]) * x_cov + u_true)
    N <- rpois(n_sites, lambda)
    p <- plogis(mu_p[1] + bp[s, 1] + (mu_p[2] + bp[s, 2]) * det_cov)
    for (i in seq_len(n_sites)) y[i, , s] <- rbinom(J, N[i], p[i])
  }
  dimnames(y) <- list(NULL, NULL, paste0("sp", seq_len(n_species)))
  dat <- data.frame(abun_cov = x_cov, det_cov = det_cov,
                    lon = coords[, 1], lat = coords[, 2])

  fit <- tobs(
    formula = ~ abun_cov + spde(lon, lat, cutoff = 0, nu = 1,
                                prior_range = c(0.3, 0.5),
                                prior_sigma = c(0.8, 0.5)),
    data = dat, family = ms_abun(), species = dimnames(y)[[3]],
    detection = ~ det_cov, y = y,
    method = "nested_laplace", control = list(verbose = FALSE))

  expect_equal(length(fit$spatial_field), fit$spatial$tulpa_spec$n_mesh)
  field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field)
  expect_gt(cor(field_at_sites, u_true), 0.6)
  # Detection slope is well identified (pooled across species); abundance slope
  # is the weakly-identified community mean at small S -- a loose bound that
  # still catches a sign flip or a dropped covariate.
  expect_lt(abs(unname(fit$means["p_det_cov"])     - mu_p[2]),      0.30)
  expect_lt(abs(unname(fit$means["lambda_abun_cov"]) - mu_lambda[2]), 0.40)
})
