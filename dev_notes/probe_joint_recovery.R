# Probe the joint nested_laplace recovery tolerance for sigma_pos / phi_pos
# before pinning the assertion in tests/. Prints per-seed estimates.
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa",    quiet = TRUE, export_all = FALSE)
  load_all(".",            quiet = TRUE, export_all = FALSE)
})

source("tests/testthat/test-cover-hurdle-nested-joint-recovery.R")

n_s <- 30L
adj <- chain_adj_for_test(n_s)

cat("--- lognormal sigma_pos recovery (10 seeds, N=400) ---\n")
truth_sigma <- 0.4
sigma_hats <- numeric(10L)
for (r in seq_len(10L)) {
  sim <- simulate_joint_lognormal_for_recovery(
    N = 400, n_s = n_s, sigma_pos_true = truth_sigma, seed = 1000L + r
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
      sigma_grid = c(0.3, 0.6, 0.9),
      rho_grid   = c(0.5, 0.7, 0.9),
      alpha_grid = c(0.5, 1.0, 1.5)
    )
  )
  sigma_hats[r] <- fit$sigma_pos
  cat(sprintf("  seed %d: sigma_hat = %.4f  (rel err %.3f)\n",
              1000L + r, fit$sigma_pos,
              (fit$sigma_pos - truth_sigma) / truth_sigma))
}
cat(sprintf("  mean = %.4f, rel mean err = %.3f, max abs rel err = %.3f\n\n",
            mean(sigma_hats),
            (mean(sigma_hats) - truth_sigma) / truth_sigma,
            max(abs(sigma_hats - truth_sigma) / truth_sigma)))

cat("--- beta phi_pos recovery (10 seeds, N=600) ---\n")
truth_phi <- 30
phi_hats <- numeric(10L)
for (r in seq_len(10L)) {
  sim <- simulate_joint_beta_for_recovery(
    N = 600, n_s = n_s, phi = truth_phi, seed = 2000L + r
  )
  spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")
  fit <- tobs(
    formula  = ~ x,
    data     = sim$data,
    family   = cover("beta"),
    y        = sim$y,
    spatial  = spatial,
    engine   = "nested_laplace",
    control  = list(
      sigma_grid = c(0.3, 0.5, 0.8),
      rho_grid   = c(0.5, 0.7, 0.9),
      alpha_grid = c(0.5, 1.0, 1.5)
    )
  )
  phi_hats[r] <- fit$phi_pos
  cat(sprintf("  seed %d: phi_hat = %.4f  (rel err %.3f)\n",
              2000L + r, fit$phi_pos,
              (fit$phi_pos - truth_phi) / truth_phi))
}
cat(sprintf("  mean = %.4f, rel mean err = %.3f, max abs rel err = %.3f\n",
            mean(phi_hats),
            (mean(phi_hats) - truth_phi) / truth_phi,
            max(abs(phi_hats - truth_phi) / truth_phi)))
