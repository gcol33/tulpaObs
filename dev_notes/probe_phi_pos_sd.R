# Quick smoke for the phi_pos_sd addition. Loads tulpaObs and runs the
# minimal fit from the new test, prints the new field.
suppressPackageStartupMessages({
  library(devtools)
  pkgbuild::clean_dll(".")
  load_all("../tulpa",    quiet = TRUE, export_all = FALSE)
  load_all(".",           quiet = TRUE, export_all = FALSE)
})

chain_adj <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  adj
}

simulate_joint_beta <- function(N = 600, n_s = 25, sigma = 0.5, rho = 0.7,
                                alpha = 1.0, phi = 30,
                                beta_occ = c(0.2, 0.7),
                                beta_pos = c(0.4, -0.5), seed = 3001) {
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  phi_f   <- rnorm(n_s); theta_f <- rnorm(n_s)
  w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)
  x       <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
  mu_pos  <- plogis(eta_pos)
  y <- numeric(N); is_pos <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos),
                     mu_pos[is_pos] * phi,
                     (1 - mu_pos[is_pos]) * phi)
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x, region = factor(spatial_idx)), y = y)
}

adj <- chain_adj(25L)
sim <- simulate_joint_beta(N = 600, n_s = 25, phi = 30, seed = 3001)
spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

t0 <- Sys.time()
fit <- tobs(
  formula  = ~ x,
  data     = sim$data,
  family   = cover("beta"),
  y        = sim$y,
  spatial  = spatial,
  engine   = "nested_laplace",
  control  = list(
    sigma_grid     = c(0.3, 0.5, 0.8),
    rho_grid       = c(0.5, 0.7, 0.9),
    sigma_pos_grid = c(0.25, 0.5, 0.75)
  )
)
cat(sprintf("Fit elapsed: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("Joint beta fit:\n")
cat(sprintf("  phi_pos     = %.4f  (truth 30)\n", fit$phi_pos))
cat(sprintf("  phi_pos_sd  = %.4f\n", fit$phi_pos_sd))
cat(sprintf("  95%% CI      = [%.2f, %.2f]\n",
            fit$phi_pos - 1.96 * fit$phi_pos_sd,
            fit$phi_pos + 1.96 * fit$phi_pos_sd))

cat("\nSeparate-hurdle beta fit (should be NA):\n")
fit_sep <- tobs(formula = ~ x, data = sim$data,
                family = cover(positive = "beta"), y = sim$y)
cat(sprintf("  phi_pos     = %.4f\n", fit_sep$phi_pos))
cat(sprintf("  phi_pos_sd  = %s\n",   format(fit_sep$phi_pos_sd)))
