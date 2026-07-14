# Community count NUTS (msAbund NUTS): samples the exact joint posterior of the
# non-spatial community Poisson GLMM (community means, per-species deviations,
# community covariance) via the in-tree C++ FullGradFn (src/ms_count_nuts.cpp),
# warm-started at the community Laplace-EM mode. Poisson.
#
# The critical correctness check is the byte-exact joint log-posterior + gradient
# vs the R oracle (.tobs_ms_count_nuts_logpost); the recovery check is that NUTS
# recovers the community means and agrees with the Laplace-EM mode.

test_that("ms_count() reports nuts + the negbin/gaussian NUTS gate", {
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$ms_count)
  sim <- simulate_ms_count(N = 40, n_species = 5, response = "negbin", seed = 1)
  expect_error(
    tobs(~ x, data = sim$data, family = ms_count("negbin"), y = sim$y,
         species = colnames(sim$y), method = "nuts"),
    "Poisson-only")
})

test_that("community count NUTS log-posterior + gradient match the R oracle", {
  skip_on_cran()
  sim <- simulate_ms_count(N = 60, n_species = 6, beta_comm_mean = c(1, 0.5),
                           response = "poisson", seed = 2)
  model <- tulpaObs:::.tobs_build_ms_count(
    formula = ~ x, data = sim$data, y = sim$y, species = colnames(sim$y),
    response = "poisson")
  P <- model$process_info[[1L]]$p; S <- model$n_species
  lay <- tulpaObs:::.tobs_ms_count_nuts_layout(P, S)
  pri <- tulpaObs:::.tobs_ms_count_nuts_priors()
  Y   <- matrix(as.numeric(model$y), model$n_sites, S)
  spec <- list(X = model$X, y = Y)

  set.seed(7)
  for (rep in 1:4) {
    theta <- stats::rnorm(lay$total, 0, 0.4)
    theta[lay$chol] <- theta[lay$chol] * 0.3          # keep the covariance sane
    r_ora <- tulpaObs:::.tobs_ms_count_nuts_logpost(theta, model$X, Y, lay, pri,
                                                    sigma.beta = 10)
    cpp   <- cpp_ms_count_nuts_joint_logpost(spec, theta, pri, 10)
    expect_equal(cpp$lp, r_ora$lp, tolerance = 1e-8)
    expect_lt(max(abs(cpp$grad - r_ora$grad)), 1e-7)
  }
})

test_that("community count NUTS recovers community means + agrees with Laplace", {
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
  # community means recover, and NUTS agrees with the Laplace mode
  expect_equal(unname(unlist(coef(nut))), c(1, 0.5), tolerance = 0.1)
  expect_equal(unname(unlist(coef(nut))), unname(unlist(coef(lap))),
               tolerance = 0.08)
  # S3 surface
  expect_true(all(is.finite(diag(vcov(nut)))))
  expect_equal(dim(ranef(nut)), c(12L * 2L, 4L))
})
