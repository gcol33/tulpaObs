# Value-only kernel path for ms_distance().
#
# The community-latent mode-adaptation backtracking line search
# (R/community_latent.R) calls the oracle's ll_cell() / data_ll() to evaluate a
# trial point's log-likelihood only -- it never reads the gradient/Fisher
# fields working() needs. cpp_distance_site_sweep(value_only = TRUE) skips the
# per-bin quadrature second-derivative accumulation and the whole detection-arm
# Louis/Fisher block in distance_kernel.h, filling only log_lik (plus the
# mean_N / var_N / boundary_weight the log-likelihood sum already produces as a
# side effect). What is asserted here is that this is a pure performance path:
# log_lik must be bit-for-bit identical to the full sweep, at both detection
# keys and with the oracle's site-subset (idx) contract.

test_that("cpp_distance_site_sweep(value_only = TRUE) matches the full sweep exactly", {
  set.seed(101)
  sim <- simulate_ms_distance(n_species = 4, N = 30, key = "halfnorm",
                              n_abund_covs = 1, seed = 101)
  model <- tulpaObs:::.tobs_build_ms_distance(
    ~ abund_cov1, ~ 1, sim$data, sim$y, sim$species, sim$cutpoints,
    key = "halfnorm", transect = "line")
  eng <- tulpaObs:::.tobs_ms_distance_engine(model)

  n <- model$n_sites
  eta_lam <- rnorm(n, log(20), 0.3)
  eta_sig <- rnorm(n, log(40), 0.2)

  full  <- eng$sweep(1, eta_lam, eta_sig, 0, value_only = FALSE)
  value <- eng$sweep(1, eta_lam, eta_sig, 0, value_only = TRUE)

  expect_identical(as.numeric(full$log_lik), as.numeric(value$log_lik))
  expect_identical(full$n_inadmissible, value$n_inadmissible)
})

test_that("value_only path matches under the hazard-rate key (nd == 2 detection block)", {
  set.seed(102)
  sim <- simulate_ms_distance(n_species = 4, N = 25, key = "hazard", shape = 0.3,
                              n_abund_covs = 1, seed = 102)
  model <- tulpaObs:::.tobs_build_ms_distance(
    ~ abund_cov1, ~ 1, sim$data, sim$y, sim$species, sim$cutpoints,
    key = "hazard", transect = "line")
  eng <- tulpaObs:::.tobs_ms_distance_engine(model)

  n <- model$n_sites
  eta_lam <- rnorm(n, log(20), 0.3)
  eta_sig <- rnorm(n, log(40), 0.2)

  full  <- eng$sweep(1, eta_lam, eta_sig, 0.3, value_only = FALSE)
  value <- eng$sweep(1, eta_lam, eta_sig, 0.3, value_only = TRUE)

  expect_identical(as.numeric(full$log_lik), as.numeric(value$log_lik))
})

test_that("oracle ll_cell() uses the value-only path and matches the site-subset contract", {
  set.seed(103)
  sim <- simulate_ms_distance(n_species = 3, N = 40, key = "halfnorm",
                              n_abund_covs = 1, seed = 103)
  model <- tulpaObs:::.tobs_build_ms_distance(
    ~ abund_cov1, ~ 1, sim$data, sim$y, sim$species, sim$cutpoints,
    key = "halfnorm", transect = "line")
  eng <- tulpaObs:::.tobs_ms_distance_engine(model)

  n <- model$n_sites; S <- model$n_species
  eta_sig <- rnorm(n, log(40), 0.2)
  eta_sig_list <- lapply(seq_len(S), function(s) eta_sig)
  oracle <- tulpaObs:::.tobs_ms_distance_oracle(eng, eta_sig_list, 0, n, S)

  eta_mat <- matrix(rnorm(n * S, log(20), 0.3), n, S)
  idx <- c(2L, 5L, 9L, n)
  ll_sub  <- oracle$ll_cell(eta_mat, idx = idx)
  ll_full <- oracle$ll_cell(eta_mat)

  expect_identical(ll_sub, ll_full[idx, , drop = FALSE])
  expect_equal(oracle$data_ll(eta_mat), sum(ll_full))
})
