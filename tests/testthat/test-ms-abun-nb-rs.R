# Community / multispecies N-mixture with a PER-SPECIES negative-binomial
# dispersion random effect:
#
#   log_r_s ~ N(mu_log_r, sigma_log_r)   (r_s ~ LogNormal)
#
# so (mu_log_r, sigma_log_r) are community hyperparameters alongside
# (mu_lambda, mu_p, Sigma_lambda, Sigma_p). The oracle widens the per-species RE
# vector to (b_lambda_s, b_p_s, b_logr_s); the third covariance block is the
# scalar sigma_log_r^2. These are parameter-recovery / coverage gates, not smoke.

test_that("ms_abun(negbin) floors the log_r block at the converged AGHQ order", {
  skip_on_cran()
  skip_if_fast()
  # . Two Gauss-Hermite nodes leave no freedom to represent curvature, so where
  # the 1-D log_r posterior is not near-Gaussian the marginal can come out
  # arbitrarily sharp -- a 17x collapse of fit$sds[["log_r"]] with the fit
  # converged, mu_log_r ordinary and nothing warned. The rule is converged at
  # three nodes, so the scalar order is a floor and a request for two does not
  # get two. Asserted on the per-block grid the fit reports rather than on an SE,
  # because the collapse is data-dependent: it needs a fixture whose realised
  # dispersion spread is wide, and only some seeds are.
  sim <- simulate_ms_abun(n_species = 3, N = 12, J = 2,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(3), 0.2), mu_p = c(0.5, -0.2),
                          sd_lambda = 0.3, sd_p = 0.3,
                          mixture = "negbin", size = 5, sigma_logr = 0.4,
                          seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, n.quad = 3L, n.quad.scalar = 2L,
                             max.iter = 5L))

  # Blocks are ordered [lambda | p | log_r]. The coefficient blocks take the
  # requested n.quad; the scalar block is floored, not clamped down to 2.
  grid <- fit$ms_community$n_quad_grid
  expect_length(grid, 3L)
  expect_identical(as.integer(grid[1:2]), c(3L, 3L))
  expect_identical(as.integer(grid[[3L]]), .TOBS_MIN_SCALAR_NQUAD)
  expect_gte(as.integer(grid[[3L]]), 3L)
})

test_that("ms_abun(negbin) recovers mu_log_r, sigma_log_r, and per-species r_s", {
  skip_on_cran()
  skip_if_fast()
  set.seed(14)
  sim <- simulate_ms_abun(n_species = 22, N = 130, J = 5,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.4, sd_p = 0.35,
                          mixture = "negbin", size = 5, sigma_logr = 0.5,
                          seed = 14)

  # n.quad = 3 keeps the AGHQ grid tractable: the per-species RE dimension is
  # p_lambda + p_p + 1 = 5, so 3^5 = 243 nodes per species (5^5 = 3125 at the
  # default). The small-cluster debias the test cares about survives.
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, n.quad = 3L))

  expect_true(isTRUE(fit$convergence$converged))
  expect_equal(fit$mixture, "negbin")

  d <- fit$ms_dispersion
  expect_false(is.null(d))
  expect_true(all(c("mu_log_r", "sigma_log_r", "r_s", "r") %in% names(d)))
  expect_length(d$r_s, 22L)
  expect_true(all(is.finite(d$r_s)) && all(d$r_s > 0))

  # mu_log_r: within ~2.5 marginal SE of truth (the log_r SE is on fit$sds).
  z_mu <- abs(d$mu_log_r - sim$truth$mu_log_r) / fit$sds[["log_r"]]
  expect_lt(z_mu, 2.5)
  # And in absolute terms within ~0.4 on the log scale at this fixture size.
  expect_lt(abs(d$mu_log_r - sim$truth$mu_log_r), 0.4)

  # sigma_log_r: recovered in the right ballpark (per-species variance components
  # are noisy at S ~ 20, so allow a wide band but reject collapse-to-0 / blow-up).
  expect_gt(d$sigma_log_r, 0.5 * sim$truth$sigma_log_r)
  expect_lt(d$sigma_log_r, 1.8 * sim$truth$sigma_log_r)

  # Per-species r_s: track the simulated truth (the point of the per-species RE).
  # The dispersion RE is only weakly identified from replicated counts, so the
  # correlation gate is modest (probe at S = 20: cor ~ 0.6) but well above chance.
  expect_gt(cor(d$r_s, sim$truth$r_s), 0.45)
  # Most species' r_s within a factor ~2.5 of truth on the log scale.
  log_err <- abs(log(d$r_s) - log(sim$truth$r_s))
  expect_gt(mean(log_err < log(2.5)), 0.7)

  # Community coefficient means still recovered alongside the dispersion RE.
  truth_m <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  names(truth_m) <- setdiff(names(fit$means), "log_r")
  for (nm in names(truth_m))
    expect_lt(abs(fit$means[[nm]] - truth_m[[nm]]) / fit$sds[[nm]], 3.0)

  # ranef() carries the per-species log-dispersion deviation as a "logr" arm.
  re <- ranef(fit)
  expect_true("logr" %in% re$arm)
  expect_equal(sum(re$arm == "logr"), 22L)
})

test_that("ms_abun(negbin) mu_log_r 95% CI covers at the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  covered <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 18, N = 100, J = 5,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                            sd_lambda = 0.4, sd_p = 0.35,
                            mixture = "negbin", size = 5, sigma_logr = 0.5,
                            seed = 500 + s)
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
                family = ms_abun(mixture = "negbin"),
                detection = ~ det_cov1,
                species = sim$species, method = "laplace",
                control = list(verbose = FALSE, n.quad = 3L))
    lr   <- fit$means[["log_r"]]
    selr <- fit$sds[["log_r"]]
    lo <- lr - 1.96 * selr; hi <- lr + 1.96 * selr
    covered[s] <- sim$truth$mu_log_r >= lo && sim$truth$mu_log_r <= hi
  }
  # Nominal 95%; allow Monte-Carlo slack over 20 seeds.
  #
  # This measures 0.800 and so does not pass. The threshold carries no
  # measurement (79e5eb7 recorded none) and the shortfall is understood rather
  # than mysterious: coverage here is conditional on the dispersion variance
  # being recovered. Across 39 fits at 8 and 36 species, those returning
  # sigma_log_r >= 0.30 covered 33/34 while those below covered 2/5, with the
  # point estimate 2.2x further out and the interval 28% narrower together. The
  # misses are variance-component collapse, not a uniform SE miss -- dropping the
  # four largest |z| from the 59-fit pool takes sqrt(mean(z^2)) from 1.189 to
  # 1.003, and the body's robust scale is 0.885 (IQR) / 0.931 (MAD).
  #
  # What this assertion should be instead needs a decision about what the family
  # claims, and the conditional form needs sigma_log_r recorded at S = 18, which
  # this fixture's earlier runs did not keep. Left as-is deliberately rather than
  # retuned to pass. See and NOTES_measurements.md.
  expect_gt(mean(covered), 0.8)
})

test_that("ms_abun(negbin) matches spAbundance::msNMix (validation only)", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("spAbundance")
  set.seed(77)
  sim <- simulate_ms_abun(n_species = 16, N = 100, J = 5,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.4, sd_p = 0.35,
                          mixture = "negbin", size = 5, sigma_logr = 0.5,
                          seed = 77)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, n.quad = 3L))
  # Validation: the tulpaObs community-mean log-dispersion should sit in the
  # neighbourhood of the simulated truth (head-to-head is informational, not a
  # hard gate -- spAbundance's NB parameterization / MCMC default differs).
  expect_lt(abs(fit$ms_dispersion$mu_log_r - sim$truth$mu_log_r), 0.6)
})
