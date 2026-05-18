# Decompose the mixture skewness into its three terms to identify which
# dominates the gap with standalone SLA at vanishing sigma.
suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa", quiet = TRUE, export_all = FALSE)
  load_all(".",         quiet = TRUE, export_all = FALSE)
})
set.seed(105L)
N <- 200L; n_s <- 25L; sigma <- 0.05; rho <- 0.7; alpha_true <- 0.01
beta_occ <- c(-0.3, 0.7); beta_pos <- c(0.4, -0.5); phi_disp <- 30
spatial_idx <- sample.int(n_s, N, replace = TRUE)
phi_f <- rnorm(n_s); phi_f <- phi_f - mean(phi_f)
theta_f <- rnorm(n_s); theta_f <- theta_f - mean(theta_f)
w_s <- sigma * (sqrt(rho)*phi_f + sqrt(1-rho)*theta_f)
x <- rnorm(N)
eta_occ <- beta_occ[1] + beta_occ[2]*x + w_s[spatial_idx]
occur <- rbinom(N, 1, plogis(eta_occ))
eta_pos <- beta_pos[1] + beta_pos[2]*x + alpha_true*w_s[spatial_idx]
mu_pos <- plogis(eta_pos)
y <- numeric(N); is_pos <- occur==1L
if(any(is_pos)) y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos]*phi_disp, (1-mu_pos[is_pos])*phi_disp)
y <- pmin(pmax(y,0), 1-1e-6)
dat <- data.frame(x = x, region = factor(spatial_idx, levels=seq_len(n_s)))
adj <- matrix(0L, n_s, n_s)
for (s in seq_len(n_s-1L)) { adj[s,s+1L]<-1L; adj[s+1L,s]<-1L }
spatial <- tulpa::spatial_bym2(adj, level="group", group_var="region")

fit_joint <- suppressMessages(tobs(
  formula = ~x, data = dat, family = cover("beta"), y = y,
  spatial = spatial, engine = "nested_laplace",
  approx = "simplified_laplace",
  control = list(sigma_grid=c(0.01,0.02,0.03), rho_grid=0.5,
                 sigma_pos_grid=c(0.0,0.01,0.02))
))
fit <- fit_joint$joint
layout <- fit$arm_layout
bpos_idx <- layout$beta_start[2] + seq_len(layout$p[2])
inner_var <- tulpaObs:::.joint_inner_var(fit, bpos_idx)
modes <- fit$modes[, bpos_idx, drop=FALSE]
weights <- fit$weights
mu_j <- as.numeric(crossprod(weights, modes))
vom <- as.numeric(crossprod(weights, modes^2)) - mu_j^2
mov <- as.numeric(crossprod(weights, inner_var))
sigma2_j <- vom + mov
sigma_j <- sqrt(sigma2_j)

cat("Per-arm betapos summary (joint):\n")
cat(sprintf("  mu_j   = %s\n", paste(round(mu_j, 5), collapse=",")))
cat(sprintf("  vom    = %s\n", paste(round(vom, 6), collapse=",")))
cat(sprintf("  mov    = %s\n", paste(round(mov, 6), collapse=",")))
cat(sprintf("  sigma_j = %s\n", paste(round(sigma_j, 5), collapse=",")))

# Compute per-grid gamma for the j=intercept and j=x.
enc_sla <- fit_joint$encoding
prior <- tulpa::prior_from_spec(spatial, dat[fit_joint$encoding$obs_keep,,drop=FALSE])
enc_sla$..spi_full <- as.integer(prior$spatial_idx)
enc_sla$..spi_pos  <- as.integer(prior$spatial_idx[fit_joint$encoding$idx_pos])

n_grid <- length(weights)
gamma_grid <- matrix(NA_real_, n_grid, 2)
for (k in seq_len(n_grid)) {
  gamma_grid[k,] <- tulpaObs:::.sla_inner_gamma_joint(k, bpos_idx, fit, enc_sla, "beta")
}

cat("\nMixture decomposition for each coefficient:\n")
for (j in seq_len(2L)) {
  dmu <- modes[,j] - mu_j[j]
  T1 <- sum(weights * dmu^3)                       # var-of-means cubed
  T2 <- sum(weights * 3 * dmu * inner_var[,j])     # cross term
  T3 <- sum(weights * gamma_grid[,j] * sqrt(pmax(inner_var[,j],0))^3)  # skew term
  M3 <- T1 + T2 + T3
  cat(sprintf("\nj=%d (intercept=1, x=2):\n", j))
  cat(sprintf("  T1 (var-of-means^3)   = %+.5e\n", T1))
  cat(sprintf("  T2 (cross)            = %+.5e\n", T2))
  cat(sprintf("  T3 (skew)             = %+.5e\n", T3))
  cat(sprintf("  M3                    = %+.5e\n", M3))
  cat(sprintf("  sigma_j^3             = %+.5e\n", sigma_j[j]^3))
  cat(sprintf("  gamma_j_total         = %+.5f\n", M3 / sigma_j[j]^3))
  cat(sprintf("  gamma_j from T3 only  = %+.5f\n", T3 / sigma_j[j]^3))
}
