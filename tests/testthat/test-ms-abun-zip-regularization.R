# =============================================================================
# test-ms-abun-zip-regularization.R - regularization + robustness of the
# community zero-inflated N-mixture (ms_abun(mixture = "zip" / "zinb")).
#
# Two guards on the joint AGHQ ZI fit:
#   * K_max floor. The marginal sums the latent N to K_max, so a count above
#     K_max is structurally impossible -> a user-pinned K_max < max(y) made the
#     per-(species,site) marginal -Inf and the joint optimum singular. It now
#     errors early with an actionable message.
#   * Structural-zero variance prior. sigma_omega is the softest AGHQ direction
#     and at few species can collapse to the 0 boundary, flattening the marginal
#     Hessian and attenuating the recovered SD. A weak PC prior (default
#     control$omega.sigma.prior = c(1, 0.05), consuming tulpa 0.0.85's
#     tulpa_re_aghq(sigma_prior=)) adds curvature there. These check the guard,
#     the recovery of sigma_omega under the default prior, and the NULL-disable
#     threading.
# =============================================================================

test_that("ms_abun(mixture = 'zip') errors when K_max is below the largest count", {
  set.seed(11)
  sim <- simulate_ms_abun(n_species = 5, N = 40, J = 4,
                          n_abund_covs = 0, n_det_covs = 0,
                          mu_lambda = log(9), mu_p = 0.4,
                          mixture = "zip", omega = 0.25, seed = 11)
  ymax <- max(unlist(sim$y), na.rm = TRUE)
  expect_gt(ymax, 5)   # counts around lambda = 9 clear a tiny K_max
  expect_error(
    tobs(~ 1, data = sim$data, detection = ~ 1, y = sim$y,
         family = ms_abun(mixture = "zip", K_max = 5L),
         species = sim$species, method = "laplace",
         control = list(verbose = FALSE, progress = FALSE, n.quad = 2L)),
    "below the largest observed count")
})

test_that("ms_abun(mixture = 'zip') default prior recovers sigma_omega without collapse", {
  skip_on_cran()
  skip_if_fast()
  # Few species + a moderate structural-zero share: the collapse-prone regime.
  truth_sig <- 0.3
  n_seed <- 6L
  sig_prior <- rep(NA_real_, n_seed)
  sig_ml    <- rep(NA_real_, n_seed)
  om_prior  <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 6, N = 60, J = 4,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                            sd_lambda = 0.35, sd_p = 0.3,
                            mixture = "zip", omega = 0.3, sigma_omega = truth_sig,
                            seed = 500 + s)
    fit_p <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
           family = ms_abun(mixture = "zip"), species = sim$species,
           method = "laplace",
           control = list(verbose = FALSE, progress = FALSE, n.quad = 2L)),
      error = function(e) NULL)
    fit_0 <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
           family = ms_abun(mixture = "zip"), species = sim$species,
           method = "laplace",
           control = list(verbose = FALSE, progress = FALSE, n.quad = 2L,
                          omega.sigma.prior = NULL)),
      error = function(e) NULL)
    if (!is.null(fit_p) && isTRUE(fit_p$convergence$converged)) {
      sig_prior[s] <- fit_p$ms_zi$sigma_omega
      om_prior[s]  <- fit_p$zi_omega
    }
    if (!is.null(fit_0) && isTRUE(fit_0$convergence$converged))
      sig_ml[s] <- fit_0$ms_zi$sigma_omega
  }
  # The prior keeps the mean sigma_omega off the 0 boundary (no systematic
  # collapse) and near truth; the community-mean omega is unbiased on average.
  expect_gt(mean(sig_prior, na.rm = TRUE), 0.12)
  expect_lt(abs(mean(sig_prior, na.rm = TRUE) - truth_sig), 0.15)
  expect_lt(abs(mean(om_prior, na.rm = TRUE) - 0.3), 0.08)
  # The prior is at least as far from 0 as pure ML on average (it debiases the
  # attenuation, never adds it) -- a directional, not a point, check.
  expect_gte(mean(sig_prior, na.rm = TRUE) + 1e-8, mean(sig_ml, na.rm = TRUE))
})
