# =============================================================================
# test-ms-abun-zip.R - zero-inflated community N-mixture (ms_abun(mixture =
# "zip" / "zinb")).
#
# Structural absence as a per-species random effect: logit_omega_s ~ N(mu_omega,
# sigma_omega), a share omega_s of a species' sites is never occupied (N = 0).
# The per-site marginal wraps the Royle marginal in a structural-zero mixture in
# the native community AGHQ oracle (NMixCommunityOracle), integrated jointly with
# the per-species abundance / detection / (ZINB) dispersion coefficients. These
# check the constructor + gates, the S3 / fit surface, and parameter recovery
# (community means + the structural-zero share) against simulated truth.
#
# Cost note: the ZI count families have no closed-form EM, so every fit is the
# joint AGHQ optimizer (n_quad^(p_lambda + p_p + 1) nodes per species). These use
# small communities + n_quad = 2 to stay tier-3-affordable; they are
# skip_on_cran() + skip_if_fast().
# =============================================================================

test_that("ms_abun(mixture = 'zip' / 'zinb') constructor + gates", {
  expect_equal(ms_abun(mixture = "zip")$params$mixture, "zip")
  expect_equal(ms_abun(mixture = "zinb")$params$mixture, "zinb")
  # ZI is the non-spatial community fit; a shared field would silently map to
  # Poisson / negbin, so it errors with a pointer instead of dropping the
  # structural-zero share.
  sim <- simulate_ms_abun(n_species = 5, N = 30, J = 3, mixture = "zip",
                          omega = 0.3, seed = 1)
  adj <- matrix(0L, 30, 30); adj[cbind(1:29, 2:30)] <- 1L; adj <- adj + t(adj)
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj), data = sim$data, y = sim$y,
         family = ms_abun(mixture = "zip"), detection = ~ det_cov1,
         species = sim$species, method = "nested_laplace",
         control = list(verbose = FALSE)),
    "zero-inflated")
})

test_that("ms_abun(mixture = 'zip') fits, exposes + recovers the structural-zero RE", {
  skip_on_cran()
  skip_if_fast()
  set.seed(21)
  sim <- simulate_ms_abun(n_species = 6, N = 60, J = 5,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(5), 0.4), mu_p = c(0.4, -0.3),
                          sd_lambda = 0.35, sd_p = 0.3,
                          mixture = "zip", omega = 0.35, sigma_omega = 0.3,
                          seed = 21)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = ms_abun(mixture = "zip", K_max = 60L),
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE, n.quad = 2L))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$convergence$converged))
  expect_identical(fit$mixture, "zip")
  expect_true(isTRUE(fit$zero_inflated))

  # The community structural-zero logit is a model coordinate (means / vcov /
  # sds with a name), like the negbin log_r.
  expect_true("logit_omega" %in% names(fit$means))
  expect_true("logit_omega" %in% rownames(vcov(fit)))
  expect_true(is.finite(fit$sds[["logit_omega"]]))
  expect_true(all(is.finite(vcov(fit))))

  # ms_zi summary: community-mean omega + per-species omega_s.
  expect_true(is.finite(fit$zi_omega) && fit$zi_omega > 0 && fit$zi_omega < 1)
  expect_length(fit$ms_zi$omega_s, 6L)
  expect_true(all(fit$ms_zi$omega_s > 0 & fit$ms_zi$omega_s < 1))

  # ranef() carries the per-species structural-zero-logit deviation.
  re <- ranef(fit)
  expect_true("omega" %in% re$arm)

  # simulate() replicates the 3D community array (with structural zeros).
  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(60L, 5L, 6L))

  # Recovery: community abundance mean + the structural-zero share.
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(5)) /
            fit$sds[["lambda_(Intercept)"]], 3.0)
  expect_lt(abs(fit$zi_omega - 0.35), 0.12)
})

test_that("ms_abun(mixture = 'zip') community-mean CIs cover across seeds", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  om_t <- 0.3
  covered <- logical(0)
  om_hat <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 6, N = 55, J = 4,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                            sd_lambda = 0.35, sd_p = 0.3,
                            mixture = "zip", omega = om_t, sigma_omega = 0.3,
                            seed = 700 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
           family = ms_abun(mixture = "zip", K_max = 50L),
           species = sim$species, method = "laplace",
           control = list(verbose = FALSE, progress = FALSE, n.quad = 2L)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
    fx <- fit$means[setdiff(names(fit$means), "logit_omega")]
    sx <- fit$sds[setdiff(names(fit$sds), "logit_omega")]
    covered <- c(covered, truth >= fx - 1.96 * sx & truth <= fx + 1.96 * sx)
    om_hat[s] <- fit$zi_omega
  }
  # Community-mean 95% CI coverage at the graduation floor.
  expect_gt(mean(covered), 0.85)
  # Community-mean structural-zero probability recovered on average.
  expect_lt(abs(mean(om_hat, na.rm = TRUE) - om_t), 0.07)
})

test_that("ms_abun(mixture = 'zinb') wires log_r + omega through the joint path", {
  skip_on_cran()
  skip_if_fast()
  # ZINB has the known zero-source confounding (structural omega vs NB
  # overdispersion); a well-identified regime -- higher counts + more visits.
  set.seed(41)
  sim <- simulate_ms_abun(n_species = 6, N = 80, J = 6,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(9), 0.4), mu_p = c(0.5, -0.2),
                          sd_lambda = 0.35, sd_p = 0.3,
                          mixture = "zinb", size = 6, omega = 0.25,
                          sigma_omega = 0.3, seed = 41)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
              family = ms_abun(mixture = "zinb", K_max = 150L),
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE, n.quad = 2L))

  expect_equal(fit$mixture, "zinb")
  # Both trailing coordinates surfaced.
  expect_true(all(c("log_r", "logit_omega") %in% names(fit$means)))
  expect_true(is.finite(fit$sds[["log_r"]]) && fit$sds[["log_r"]] > 0)
  expect_true(is.finite(fit$sds[["logit_omega"]]))

  # Size r and structural-zero omega in the right ballpark (lenient on the
  # zero-source ridge).
  expect_true(is.finite(fit$ms_dispersion$r) && fit$ms_dispersion$r > 0)
  expect_lt(abs(log(fit$ms_dispersion$r) - log(6)), log(3.5))
  expect_lt(abs(fit$zi_omega - 0.25), 0.15)
})
