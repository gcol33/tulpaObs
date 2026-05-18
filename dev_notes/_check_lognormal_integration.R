devtools::load_all(quiet = TRUE)

# Single-seed sanity check on the new outer phi-grid integration for the
# lognormal cover-hurdle joint path. Expect:
#   * fit$sigma_pos finite and within 25% of truth (0.4)
#   * fit$sigma_pos_sd finite and > 0
#   * fit$sigma_pos +/- 3 * sigma_pos_sd brackets the truth
source("tests/testthat/test-cover-hurdle-nested-joint-recovery.R",
       echo = FALSE)

# `chain_adj_for_test` and `simulate_joint_lognormal_for_recovery` are now
# defined in the global env via source().
truth_sigma <- 0.4
n_s   <- 30L
adj   <- chain_adj_for_test(n_s)
seeds <- 1001:1010
sigma_hats <- numeric(length(seeds))
sigma_sds  <- numeric(length(seeds))
for (i in seq_along(seeds)) {
  sim <- simulate_joint_lognormal_for_recovery(
    N = 400, n_s = n_s, sigma_pos_true = truth_sigma, seed = seeds[i]
  )
  spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")
  fit <- tobs(
    formula  = ~ x,
    data     = sim$data,
    family   = cover("lognormal"),
    y        = sim$y,
    spatial  = spatial,
    engine   = "nested_laplace",
    control  = list(
      sigma_grid     = c(0.3, 0.6, 0.9),
      rho_grid       = c(0.5, 0.7, 0.9),
      sigma_pos_grid = c(0.3, 0.6, 0.9)
    )
  )
  sigma_hats[i] <- fit$sigma_pos
  sigma_sds[i]  <- fit$sigma_pos_sd
  cat(sprintf("seed %d  sigma_pos = %.4f  +/- %.4f  rel_err = %+5.1f%%\n",
              seeds[i], fit$sigma_pos, fit$sigma_pos_sd,
              100 * (fit$sigma_pos - truth_sigma) / truth_sigma))
}
cat(sprintf("\nmean(sigma_pos) = %.4f (truth %.4f)  mean rel err = %.1f%%\n",
            mean(sigma_hats), truth_sigma,
            100 * abs(mean(sigma_hats) - truth_sigma) / truth_sigma))
cat(sprintf("max rel err = %.1f%%\n",
            100 * max(abs(sigma_hats - truth_sigma) / truth_sigma)))
cat(sprintf("mean(sigma_pos_sd) = %.4f\n", mean(sigma_sds)))
