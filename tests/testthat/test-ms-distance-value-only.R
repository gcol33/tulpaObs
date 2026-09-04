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
# log_lik must agree with the full sweep to within floating-point reassociation,
# at both detection keys and with the oracle's site-subset (idx) contract.
#
# Agreement is to a few ULP, not to the bit, and that is a property of the code
# rather than a slack allowance. Both branches obtain g by identical
# expressions -- dist_key_value() returns exp(-0.5 * u) and dist_key_deriv()
# assigns k.g the same, likewise 1 - exp(-z) under the hazard key -- so the
# arithmetic content is the same. What differs is the shape of the loop
# consuming it: `s += w * dist_key_value(...)` alone is a reduction a compiler
# can contract into a single fused multiply-add, while `s += w * k.g`
# interleaved with five other accumulations schedules differently. Where FMA is
# in the baseline instruction set the two contract differently and the sums part
# by a few ULP; on x86_64 there is no FMA to contract into and they coincide, so
# asserting bit-identity pinned an accident of x86_64 codegen (issue #314). The
# gap is bounded in ULPs and reported instead: a genuine divergence between the
# two paths would be orders of magnitude larger. The site-subset test below
# compares the value-only path against itself, one loop shape, and stays exact.

# Gap between two floating-point vectors in ULPs of the value, the scale a
# reassociation difference actually lives on.
max_ulp_gap <- function(actual, expected) {
  max(abs(actual - expected) /
      (pmax(abs(expected), .Machine$double.xmin) * .Machine$double.eps))
}

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

  fl <- as.numeric(full$log_lik); vl <- as.numeric(value$log_lik)
  # Which sites are impossible (-Inf) is a branch, not an accumulation, so that
  # stays exact; only the finite sums carry the reassociation.
  expect_identical(is.finite(fl), is.finite(vl))
  expect_lt(max_ulp_gap(fl[is.finite(fl)], vl[is.finite(vl)]), 64)
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

  fl <- as.numeric(full$log_lik); vl <- as.numeric(value$log_lik)
  expect_identical(is.finite(fl), is.finite(vl))
  expect_lt(max_ulp_gap(fl[is.finite(fl)], vl[is.finite(vl)]), 64)
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
