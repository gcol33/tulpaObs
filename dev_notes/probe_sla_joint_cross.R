# Cross-check probe: replicate the failing test scenario and print
# per-grid gamma + mixture contributions for the joint beta SLA path,
# then the standalone SLA gamma for comparison.
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",         quiet = TRUE, export_all = FALSE)
})

# Replicate .make_cover_data(seed=105, alpha=0.01, sigma=0.05, N=200, n_s=25).
set.seed(105L)
N <- 200L; n_s <- 25L; sigma <- 0.05; rho <- 0.7; alpha_true <- 0.01
beta_occ <- c(-0.3, 0.7); beta_pos <- c(0.4, -0.5); phi_disp <- 30
spatial_idx <- sample.int(n_s, N, replace = TRUE)
phi_f   <- rnorm(n_s, 0, 1); phi_f   <- phi_f   - mean(phi_f)
theta_f <- rnorm(n_s, 0, 1); theta_f <- theta_f - mean(theta_f)
w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)
x <- rnorm(N)
eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
occur   <- rbinom(N, 1, plogis(eta_occ))
eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha_true * w_s[spatial_idx]
mu_pos  <- plogis(eta_pos)
y <- numeric(N); is_pos <- occur == 1L
if (any(is_pos)) {
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi_disp,
                     (1 - mu_pos[is_pos]) * phi_disp)
}
y[!is_pos] <- 0; y <- pmin(pmax(y, 0), 1 - 1e-6)
dat <- data.frame(x = x, region = factor(spatial_idx, levels = seq_len(n_s)))

adj <- matrix(0L, n_s, n_s)
for (s in seq_len(n_s - 1L)) { adj[s, s+1L] <- 1L; adj[s+1L, s] <- 1L }

spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

cat("--- standalone SLA on cover('beta') ---\n")
fit_sep <- suppressMessages(tobs(
  formula = ~ x, data = dat, family = cover("beta"),
  y = y, approx = "simplified_laplace"
))
cat("  skew_pos = "); print(fit_sep$skew_pos)
cat("  skew_occ = "); print(fit_sep$skew_occ)
cat("  beta_pos = "); print(fit_sep$beta_pos)
cat("  se_pos   = "); print(fit_sep$se_pos)
cat("  phi_pos  = "); print(fit_sep$phi_pos)

cat("\n--- joint SLA on cover('beta') ---\n")
fit_joint <- suppressMessages(tobs(
  formula = ~ x, data = dat, family = cover("beta"),
  y = y, spatial = spatial, engine = "nested_laplace",
  approx = "simplified_laplace",
  control = list(sigma_grid = c(0.01,0.02,0.03), rho_grid = 0.5,
                 sigma_pos_grid = c(0.0,0.01,0.02))
))
cat("  skew_pos = "); print(fit_joint$skew_pos)
cat("  skew_occ = "); print(fit_joint$skew_occ)
cat("  beta_pos = "); print(fit_joint$beta_pos)
cat("  se_pos   = "); print(fit_joint$se_pos)
cat("  phi_pos  = "); print(fit_joint$phi_pos)
cat("  sla_status = ", fit_joint$sla_status, "\n")

# Per-grid debug for the joint fit.
fit <- fit_joint$joint
layout <- fit$arm_layout
p_occ <- layout$p[1]; p_pos <- layout$p[2]
bpos_idx <- layout$beta_start[2] + seq_len(p_pos)

cat("\n--- theta_grid columns: ", paste(colnames(fit$theta_grid), collapse=", "), "\n")
cat("  weights:  ", paste(round(fit$weights, 4), collapse=", "), "\n")
cat("  n_grid:   ", length(fit$weights), "\n")

enc_sla <- fit_joint$encoding
# Re-derive spi_full / spi_pos as in fit_cover_hurdle_joint_nested:
prior <- tulpa::prior_from_spec(spatial, dat[fit_joint$encoding$obs_keep, , drop=FALSE])
enc_sla$..spi_full <- as.integer(prior$spatial_idx)
enc_sla$..spi_pos  <- as.integer(prior$spatial_idx[fit_joint$encoding$idx_pos])

inner_var <- tulpaObs:::.joint_inner_var(fit, bpos_idx)
cat("  inner_var[, pos] (per grid):\n"); print(round(inner_var, 5))

for (k in seq_along(fit$weights)) {
  g <- tulpaObs:::.sla_inner_gamma_joint(k, bpos_idx, fit, enc_sla, "beta")
  cat(sprintf("  k=%d phi=%6.2f s_pos=%.3f w=%.3f  gamma_pos = ",
              k, fit$theta_grid[k, "phi_pos"],
              fit$theta_grid[k, "sigma_pos"], fit$weights[k]))
  cat(paste(round(g, 4), collapse=", "), "\n")
}

cat("\nCompare: standalone gamma_pos =", round(fit_sep$skew_pos, 4), "\n")
