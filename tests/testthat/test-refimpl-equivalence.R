# Reference-implementation equivalence + a CI recovery gate for the abundance
# families the audit flagged as unevenly validated. The N-mixture family already
# gates head-to-head against unmarked::pcount (test-nmix-laplace.R); this adds the
# analogous removal check against unmarked::multinomPois and a CI-runnable recovery
# gate for the multispecies ms_abun family (whose only existing cross-check is an
# offline spAbundance probe).

# --------------------------------------------------------------------------- #
# removal: head-to-head against unmarked::multinomPois (removalPiFun).          #
# --------------------------------------------------------------------------- #

test_that("removal_laplace matches unmarked::multinomPois", {
  skip_if_fast()
  skip_if_not_installed("unmarked")

  sim <- simulate_removal(N = 200, K = 4, n_abund_covs = 1, n_det_covs = 1,
                          seed = 7)
  N <- nrow(sim$y); K <- ncol(sim$y)
  # removal_laplace consumes a long (site x occasion) response: abundance design
  # is per-site, detection design per-observation.
  fit <- removal_laplace(
    y = as.vector(sim$y), site_idx = rep(seq_len(N), times = K),
    X_lambda = cbind(1, sim$data$abund_cov1),
    X_p      = cbind(1, rep(sim$data$det_cov1, times = K)),
    mixture  = "P", max_iter = 100L, tol = 1e-8)

  umf <- unmarked::unmarkedFrameMPois(
    y = sim$y,
    siteCovs = data.frame(abund_cov1 = sim$data$abund_cov1,
                          det_cov1   = sim$data$det_cov1),
    type = "removal")
  um <- unmarked::multinomPois(~ det_cov1 ~ abund_cov1, data = umf)
  uc <- unmarked::coef(um)

  # tulpa matches unmarked's MLE essentially exactly here.
  expect_gte(fit$log_lik, as.numeric(unmarked::logLik(um)) - 1e-3)
  expect_lt(abs(fit$beta_lambda[1] - uc["lambda(Int)"]),        5e-3)
  expect_lt(abs(fit$beta_lambda[2] - uc["lambda(abund_cov1)"]), 5e-3)
  expect_lt(abs(fit$beta_p[1]      - uc["p(Int)"]),             5e-3)
  expect_lt(abs(fit$beta_p[2]      - uc["p(det_cov1)"]),        5e-3)
})

# --------------------------------------------------------------------------- #
# distance: head-to-head against unmarked::distsamp (halfnorm key).            #
# --------------------------------------------------------------------------- #

test_that("distance_laplace matches unmarked::distsamp", {
  skip_if_fast()
  skip_if_not_installed("unmarked")

  sim <- simulate_distance(N = 200, seed = 11)
  fit <- distance_laplace(
    y = sim$y, X_lambda = cbind(1, sim$data$abund_cov1),
    X_sigma = cbind(1, sim$data$sigma_cov1), cutpoints = sim$cutpoints,
    key = "halfnorm", transect = sim$truth$transect, mixture = "P",
    max_iter = 100L, tol = 1e-8)

  umf <- unmarked::unmarkedFrameDS(
    y = sim$y, survey = sim$truth$transect, dist.breaks = sim$cutpoints,
    unitsIn = "m",
    siteCovs = data.frame(abund_cov1 = sim$data$abund_cov1,
                          sigma_cov1 = sim$data$sigma_cov1),
    tlength = if (sim$truth$transect == "line") rep(1, nrow(sim$y)) else NULL)
  um <- unmarked::distsamp(~ sigma_cov1 ~ abund_cov1, data = umf,
                           keyfun = "halfnorm", output = "abund")
  uc <- unmarked::coef(um)

  expect_gte(fit$log_lik, as.numeric(unmarked::logLik(um)) - 1e-3)
  expect_lt(abs(fit$beta_lambda[1] - uc[1]), 5e-3)
  expect_lt(abs(fit$beta_lambda[2] - uc[2]), 5e-3)
  expect_lt(abs(fit$beta_sigma[1]  - uc[3]), 5e-3)
  expect_lt(abs(fit$beta_sigma[2]  - uc[4]), 5e-3)
})

# --------------------------------------------------------------------------- #
# ms_abun: CI-runnable recovery gate on a tiny fixture (the multispecies        #
# claim was only checked by an offline spAbundance probe; this catches a        #
# regression in CI).                                                            #
# --------------------------------------------------------------------------- #

test_that("ms_abun recovers the community-mean (lambda, p) intercepts (tiny fixture)", {
  skip_on_cran()
  skip_if_fast()
  # ms_abun is a community EM (multi-minute per fit), so this is a single small
  # fixture rather than a multi-seed sweep -- a CI smoke gate that catches a gross
  # regression in the multispecies recovery (the offline spAbundance probe in
  # dev_notes/ remains the precise cross-check).
  sim <- simulate_ms_abun(n_species = 8L, N = 70L, J = 4L, seed = 21L)
  fit <- tobs(~ 1, data = data.frame(s = seq_len(70L)), family = ms_abun(),
              detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(8L)),
              method = "laplace", control = list(verbose = FALSE))
  # Community means (truth: mu_lambda[1] = 1.1 on log scale, mu_p[1] = 0.3 logit).
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - 1.1), 0.4)
  expect_lt(abs(fit$means[["p_(Intercept)"]]      - 0.3), 0.4)
})
