# Probe the joint nested_laplace recovery numbers for sigma_pos / phi_pos.
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa",    quiet = TRUE, export_all = FALSE)
  load_all(".",            quiet = TRUE, export_all = FALSE)
})

chain_adj_for_test <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  adj
}

simulate_joint_beta <- function(N = 600, n_s = 30,
                                sigma = 0.5, rho = 0.7,
                                alpha = 1.0, phi = 30,
                                beta_occ = c(0.2, 0.7),
                                beta_pos = c(0.4, -0.5), seed = 23) {
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  phi_f       <- rnorm(n_s, 0, 1)
  theta_f     <- rnorm(n_s, 0, 1)
  w_s         <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)
  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
  mu_pos  <- plogis(eta_pos)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos]  <- rbeta(sum(is_pos),
                     mu_pos[is_pos] * phi,
                     (1 - mu_pos[is_pos]) * phi)
  y[!is_pos] <- 0
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x, region = factor(spatial_idx)), y = y,
       n_pos = sum(is_pos))
}

n_s <- 30L
adj <- chain_adj_for_test(n_s)

cat("--- beta phi_pos recovery (10 seeds, N=600) ---\n")
truth_phi <- 30
phi_hats <- numeric(10L)
n_pos_v  <- integer(10L)
for (r in seq_len(10L)) {
  sim <- simulate_joint_beta(N = 600, n_s = n_s, phi = truth_phi, seed = 2000L + r)
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
  n_pos_v[r]  <- sim$n_pos
  cat(sprintf("  seed %d  n_pos %d  phi_hat = %.3f  (rel err %.3f)\n",
              2000L + r, sim$n_pos, fit$phi_pos,
              (fit$phi_pos - truth_phi) / truth_phi))
}
cat(sprintf("  mean = %.3f, rel mean err = %.3f, max abs rel err = %.3f\n",
            mean(phi_hats),
            (mean(phi_hats) - truth_phi) / truth_phi,
            max(abs(phi_hats - truth_phi) / truth_phi)))
