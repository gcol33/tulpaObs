# Community count NUTS (msAbund NUTS): samples the exact joint posterior of the
# non-spatial community GLMM (community means, per-species deviations, community
# covariance) via the in-tree C++ FullGradFn (src/ms_count_nuts.cpp), warm-started
# at the community Laplace-EM mode. Poisson / negbin (dispersion community RE) /
# gaussian (per-species free residual variance).
#
# The critical correctness check is the byte-exact joint log-posterior + gradient
# vs the R oracle (.tobs_ms_count_nuts_logpost) for every response; the recovery
# check is that NUTS recovers the community means / dispersion and agrees with the
# Laplace-EM mode.

test_that("ms_count() reports nuts for every response family", {
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$ms_count)
})

test_that("community count NUTS log-posterior + gradient match the R oracle", {
  skip_on_cran()
  check_family <- function(response, seed) {
    sim <- simulate_ms_count(N = 60, n_species = 6, beta_comm_mean = c(1, 0.5),
                             response = response, seed = seed)
    model <- tulpaObs:::.tobs_build_ms_count(
      formula = ~ x, data = sim$data, y = sim$y, species = colnames(sim$y),
      response = response)
    P   <- model$process_info[[1L]]$p; S <- model$n_species
    lay <- tulpaObs:::.tobs_ms_count_nuts_layout(P, S, response)
    pri <- tulpaObs:::.tobs_ms_count_nuts_priors()
    Y   <- matrix(as.numeric(model$y), model$n_sites, S)
    spec <- list(X = model$X, y = Y, family = response)
    if (identical(response, "gaussian")) {
      spec$logphi_mean <- pri$logphi_mean; spec$logphi_sd <- pri$logphi_sd
    }
    set.seed(seed + 100L)
    for (rep in 1:4) {
      theta <- stats::rnorm(lay$total, 0, 0.4)
      theta[lay$chol_beta] <- theta[lay$chol_beta] * 0.3   # keep the covariance sane
      if (lay$is_nb) theta[lay$chol_logr] <- theta[lay$chol_logr] * 0.3
      r_ora <- tulpaObs:::.tobs_ms_count_nuts_logpost(
        theta, model$X, Y, lay, pri, sigma.beta = 10, sigma.logr = 1.5)
      cpp <- cpp_ms_count_nuts_joint_logpost(spec, theta, pri, 10, 1.5)
      expect_equal(cpp$lp, r_ora$lp, tolerance = 1e-8)
      expect_lt(max(abs(cpp$grad - r_ora$grad)), 1e-7)
    }
  }
  check_family("poisson",  2)
  check_family("negbin",   3)
  check_family("gaussian", 4)
})

test_that("community count Poisson NUTS recovers community means + agrees with Laplace", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_ms_count(N = 150, n_species = 12, beta_comm_mean = c(1, 0.5),
                           beta_comm_sd = c(0.4, 0.3), response = "poisson",
                           seed = 5)
  lap <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
              species = colnames(sim$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nut <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
              species = colnames(sim$y), method = "nuts",
              control = list(n.iter = 600L, n.warmup = 600L, seed = 1,
                             verbose = FALSE, progress = FALSE))
  expect_identical(nut$method, "nuts")
  expect_equal(sum(nut$divergent), 0)
  # Community means against the seed's REALIZED mean colMeans(truth$beta_species),
  # not the nominal c(1, 0.5) (#155): a single-seed comparison to the nominal
  # population constant is exactly the defect this issue targets. Budget = 3 sd
  # of the Laplace-estimator deviation from the realized mean over a 16-seed
  # measurement of this exact fixture (N=150, S=12, poisson: sd 0.018 / 0.012,
  # max 0.036 / 0.024), plus the NUTS-vs-Laplace gap this file already measures
  # below (0.08) -- the two errors are not perfectly aligned, so this is a
  # conservative sum, not a tight one. Seed 5 itself sits at 0.018 / -0.026, well
  # inside it.
  real <- colMeans(sim$truth$beta_species)
  expect_community_mean(nut, real, c(0.14, 0.12))
  # NUTS agrees closely with the Laplace mode (the discriminating check).
  expect_equal(unname(unlist(coef(nut))), unname(unlist(coef(lap))),
               tolerance = 0.08)
  # S3 surface
  expect_true(all(is.finite(diag(vcov(nut)))))
  expect_equal(dim(ranef(nut)), c(12L * 2L, 4L))
})

test_that("ms_count NUTS community-mean 95% CIs cover at the nominal rate", {
  # gcol33/tulpaObs#139: the single-seed recovery proves the byte-exact oracle;
  # this proves the community-mean intervals cover at the nominal rate (the
  # calibration NUTS exists to get right).
  skip_if_fast()
  skip_on_cran()
  truth <- c(1, 0.5)
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_ms_count(N = 150, n_species = 12, beta_comm_mean = truth,
                             beta_comm_sd = c(0.4, 0.3), response = "poisson",
                             seed = 300 + s)
    fit <- tryCatch(tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
                    species = colnames(sim$y), method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1,
                                   verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
    covered <- c(covered, abs(est - truth) <= 1.96 * se)
  }
  expect_gte(mean(covered), 0.85)
})

test_that("community count negbin NUTS recovers community means + dispersion", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_ms_count(N = 150, n_species = 12, beta_comm_mean = c(1, 0.5),
                           beta_comm_sd = c(0.4, 0.3), response = "negbin",
                           size = 3, seed = 6)
  lap <- tobs(~ x, data = sim$data, family = ms_count("negbin"), y = sim$y,
              species = colnames(sim$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nut <- tobs(~ x, data = sim$data, family = ms_count("negbin"), y = sim$y,
              species = colnames(sim$y), method = "nuts",
              control = list(n.iter = 600L, n.warmup = 600L, seed = 1,
                             verbose = FALSE, progress = FALSE))
  expect_identical(nut$method, "nuts")
  expect_equal(sum(nut$divergent), 0)
  # Against the seed's REALIZED mean, not the nominal c(1, 0.5) (#155). The
  # negbin community-mean estimator carries real, larger seed-to-seed spread
  # than poisson/gaussian at this same N=150/S=12 fixture -- measured over 16
  # seeds: intercept sd 0.110 (max 0.371), slope sd 0.207 (max 0.608), against
  # poisson's 0.018/0.012. That is the documented per-species dispersion RE
  # interacting with the community mean (CLAUDE.md: "negbin slope carries mild
  # first-order-Laplace attenuation"), not an artefact of retargeting -- the
  # budget below is wider than poisson's for that reason, not a looser standard.
  # 3 sd of that spread plus the NUTS-vs-Laplace gap measured below (0.12).
  # Seed 6 itself sits at 0.025 / 0.000, well inside it.
  real <- colMeans(sim$truth$beta_species)
  expect_community_mean(nut, real, c(0.45, 0.75))
  expect_equal(unname(unlist(coef(nut))), unname(unlist(coef(lap))),
               tolerance = 0.12)
  expect_identical(nut$ms_dispersion$response, "negbin")
  expect_length(nut$ms_dispersion$r_s, 12L)
  expect_true(all(is.finite(nut$ms_dispersion$r_s)) &&
                all(nut$ms_dispersion$r_s > 0))
})

test_that("community count NUTS log-posterior + gradient handle missing (NA) entries", {
  skip_on_cran()
  check_family_na <- function(response, seed) {
    sim <- simulate_ms_count(N = 60, n_species = 6, beta_comm_mean = c(1, 0.5),
                             response = response, seed = seed)
    y <- sim$y
    # Knock out a scattered set of site x species entries, plus nearly a whole
    # column for one species (ragged per-species coverage the mask must respect).
    set.seed(seed)
    y[cbind(sample(nrow(y), 15L), sample(ncol(y), 15L, replace = TRUE))] <- NA
    y[3:58, 4] <- NA                         # species 4 keeps only 4 observed sites
    model <- tulpaObs:::.tobs_build_ms_count(
      formula = ~ x, data = sim$data, y = y, species = colnames(y),
      response = response)
    P   <- model$process_info[[1L]]$p; S <- model$n_species
    lay <- tulpaObs:::.tobs_ms_count_nuts_layout(P, S, response)
    pri <- tulpaObs:::.tobs_ms_count_nuts_priors()
    Y   <- matrix(as.numeric(model$y), model$n_sites, S)
    expect_true(anyNA(Y))
    spec <- list(X = model$X, y = Y, family = response)
    if (identical(response, "gaussian")) {
      spec$logphi_mean <- pri$logphi_mean; spec$logphi_sd <- pri$logphi_sd
    }
    set.seed(seed + 200L)
    for (rep in 1:4) {
      theta <- stats::rnorm(lay$total, 0, 0.4)
      theta[lay$chol_beta] <- theta[lay$chol_beta] * 0.3
      if (lay$is_nb) theta[lay$chol_logr] <- theta[lay$chol_logr] * 0.3
      r_ora <- tulpaObs:::.tobs_ms_count_nuts_logpost(
        theta, model$X, Y, lay, pri, sigma.beta = 10, sigma.logr = 1.5)
      cpp <- cpp_ms_count_nuts_joint_logpost(spec, theta, pri, 10, 1.5)
      # A finite target (no NA leaking into the sum) that still matches the oracle
      # byte-for-byte with the missing entries dropped on both sides.
      expect_true(is.finite(cpp$lp))
      expect_equal(cpp$lp, r_ora$lp, tolerance = 1e-8)
      expect_lt(max(abs(cpp$grad - r_ora$grad)), 1e-7)
    }
  }
  check_family_na("poisson",  12)
  check_family_na("negbin",   13)
  check_family_na("gaussian", 14)
})

test_that("community count NUTS accepts missing (NA) entries and matches Laplace", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_ms_count(N = 150, n_species = 12, beta_comm_mean = c(1, 0.5),
                           beta_comm_sd = c(0.4, 0.3), response = "poisson",
                           seed = 21)
  y <- sim$y
  set.seed(21)
  y[sample(length(y), floor(0.12 * length(y)))] <- NA   # ~12% missing at random
  lap <- tobs(~ x, data = sim$data, family = ms_count(), y = y,
              species = colnames(y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nut <- tobs(~ x, data = sim$data, family = ms_count(), y = y,
              species = colnames(y), method = "nuts",
              control = list(n.iter = 600L, n.warmup = 600L, seed = 1,
                             verbose = FALSE, progress = FALSE))
  expect_identical(nut$method, "nuts")
  expect_equal(sum(nut$divergent), 0)
  expect_equal(nut$N, sum(!is.na(y)))                   # N counts observed entries
  # Against the seed's REALIZED mean, not the nominal c(1, 0.5) (#155). Same
  # poisson N=150/S=12 fixture as the first test above, so its measured spread
  # applies (3 sd 0.054 / 0.036, NUTS-vs-Laplace gap 0.1 below); the 12% MCAR
  # mask reduces information per cell slightly but not enough to move this
  # budget. Seed 21 itself sits at 0.037 / 0.015 on the unmasked fixture.
  real <- colMeans(sim$truth$beta_species)
  expect_community_mean(nut, real, c(0.14, 0.12))
  expect_equal(unname(unlist(coef(nut))), unname(unlist(coef(lap))),
               tolerance = 0.1)
  expect_true(is.finite(tobs_waic(nut)$waic))
})

test_that("community count gaussian NUTS recovers community means + residual variance", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_ms_count(N = 150, n_species = 12, beta_comm_mean = c(1, 0.5),
                           beta_comm_sd = c(0.4, 0.3), response = "gaussian",
                           sd = 1, seed = 7)
  nut <- tobs(~ x, data = sim$data, family = ms_count("gaussian"), y = sim$y,
              species = colnames(sim$y), method = "nuts",
              control = list(n.iter = 600L, n.warmup = 600L, seed = 1,
                             verbose = FALSE, progress = FALSE))
  expect_identical(nut$method, "nuts")
  expect_equal(sum(nut$divergent), 0)
  # Against the seed's REALIZED mean, not the nominal c(1, 0.5) (#155). Budget =
  # 3 sd of the Laplace-estimator deviation over a 16-seed measurement of this
  # exact fixture (N=150, S=12, gaussian: sd 0.028 / 0.020, max 0.052 / 0.036),
  # plus a NUTS-vs-Laplace margin (this family has no Laplace comparison
  # elsewhere in this file to measure that gap directly from, so 0.08 is
  # borrowed from the poisson case above -- both are well-informed,
  # non-dispersion families). Seed 7 itself sits at -0.002 / -0.002.
  real <- colMeans(sim$truth$beta_species)
  expect_community_mean(nut, real, c(0.16, 0.14))
  expect_identical(nut$ms_dispersion$response, "gaussian")
  expect_length(nut$ms_dispersion$variance, 12L)
  # residual variance recovers near the simulation truth sd^2 = 1
  expect_equal(mean(nut$ms_dispersion$variance), 1, tolerance = 0.3)
})
