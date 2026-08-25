# =============================================================================
# test-ms-occu-spatial-spde.R - continuous Matern (SPDE) shared field on the
# community occupancy arm (ms_occu() + spde(); #239).
#
# Mirrors the areal recovery test in test-ms-occu-spatial.R and the continuous
# SPDE test for count() in test-count-spatial.R: the field lives on n_mesh FEM
# nodes and the barycentric projector fit$spatial$tulpa_spec$A maps it onto
# sites. This exercises the SPDE arm added to the shared community areal driver
# (community_spatial_em.h / src/ms_occu_spatial.cpp), the occupancy analogue of
# the count/N-mixture community SPDE path (cpp_nmix_community_spatial_spde).
# =============================================================================

test_that("ms_occu + spde() recovers community means + a continuous field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  set.seed(42)
  n <- 200L
  coords <- cbind(stats::runif(n), stats::runif(n))
  u <- 0.9 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x <- stats::rnorm(n)
  dat <- data.frame(x = x, lon = coords[, 1], lat = coords[, 2])
  X_occ <- stats::model.matrix(~ x, dat)

  n_species <- 10L
  mu_psi <- c(0, 0.5); sd_psi <- c(0.4, 0.3)
  mu_p <- 0.2; sd_p <- 0.3
  J <- 4L
  beta_psi <- cbind(stats::rnorm(n_species, mu_psi[1], sd_psi[1]),
                    stats::rnorm(n_species, mu_psi[2], sd_psi[2]))
  beta_p <- matrix(stats::rnorm(n_species, mu_p, sd_p), n_species, 1)
  y <- array(NA_integer_, dim = c(n, J, n_species))
  for (s in seq_len(n_species)) {
    psi <- stats::plogis(as.numeric(X_occ %*% beta_psi[s, ]) + u)
    p   <- stats::plogis(beta_p[s, 1])
    z <- stats::rbinom(n, 1, psi)
    for (i in seq_len(n)) y[i, , s] <- stats::rbinom(J, 1, z[i] * p)
  }
  species <- paste0("sp", seq_len(n_species))

  fit <- tobs(~ x + spde(lon, lat, max_edge = c(0.3, 0.6), nu = 1,
                         prior_range = c(0.3, 0.5), prior_sigma = c(0.9, 0.5)),
              data = dat, family = ms_occu(), detection = ~ 1,
              y = y, species = species, method = "nested_laplace",
              control = list(verbose = FALSE))

  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$spatial$type, "spde")
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$spatial$tulpa_spec$A))
  expect_equal(length(fit$spatial_field), fit$spatial$n_units)

  field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field)
  expect_gt(cor(field_at_sites, u), 0.6)

  truth <- c("psi_(Intercept)" = mu_psi[1], "psi_x" = mu_psi[2], "p_(Intercept)" = mu_p)
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 3))

  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi[, 1], beta_psi[, 1]), 0.4)
})

test_that("ms_occu spde() field composes with build_ms_occu_fit S3", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  set.seed(3)
  n <- 80L
  coords <- cbind(stats::runif(n), stats::runif(n))
  u <- 0.5 * sin(4 * coords[, 1])
  dat <- data.frame(lon = coords[, 1], lat = coords[, 2])
  n_species <- 5L; J <- 3L
  beta_psi <- stats::rnorm(n_species, 0, 0.3)
  beta_p <- stats::rnorm(n_species, 0, 0.3)
  y <- array(NA_integer_, dim = c(n, J, n_species))
  for (s in seq_len(n_species)) {
    psi <- stats::plogis(beta_psi[s] + u)
    p   <- stats::plogis(beta_p[s])
    z <- stats::rbinom(n, 1, psi)
    for (i in seq_len(n)) y[i, , s] <- stats::rbinom(J, 1, z[i] * p)
  }
  species <- paste0("sp", seq_len(n_species))

  fit <- tobs(~ spde(lon, lat, max_edge = c(0.3, 0.6), nu = 1,
                     prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5)),
              data = dat, family = ms_occu(), detection = ~ 1,
              y = y, species = species, method = "nested_laplace",
              control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_false(is.null(coef(fit)))
  expect_false(is.null(fit$ms_community))
})
