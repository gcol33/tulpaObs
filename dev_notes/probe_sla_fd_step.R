# Test different FD step sizes on the problematic grid point.
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
p_pos <- layout$p[2]
bpos_idx <- layout$beta_start[2] + seq_len(p_pos)
enc_sla <- fit_joint$encoding
prior <- tulpa::prior_from_spec(spatial, dat[fit_joint$encoding$obs_keep,,drop=FALSE])
enc_sla$..spi_full <- as.integer(prior$spatial_idx)
enc_sla$..spi_pos  <- as.integer(prior$spatial_idx[fit_joint$encoding$idx_pos])

# Manually compute d3 at k=68 with various step scales.
k <- 68
n_x <- layout$n_x
Qp <- fit$Q_csc_p_per_grid; Qi <- fit$Q_csc_i_per_grid; Qx <- fit$Q_csc_x_per_grid
Qk_lt <- Matrix::sparseMatrix(
  i = as.integer(Qi[[k]])+1L, p = as.integer(Qp[[k]]), x = as.numeric(Qx[[k]]),
  dims = c(n_x,n_x), symmetric=FALSE, index1=TRUE
)
Qk <- Matrix::forceSymmetric(Qk_lt, uplo="L")
V <- tulpaObs:::.joint_constrained_solve_columns(Qk, layout, bpos_idx)

cat("--- k=68, phi=", fit$theta_grid[k,"phi_pos"],
    ", s_pos=", fit$theta_grid[k,"sigma_pos"], "---\n")
for (j in seq_len(p_pos)) {
  v_j <- V[,j]
  sigma_j <- sqrt(max(v_j[bpos_idx[j]],0))
  norm_v <- sqrt(sum(v_j^2))
  beta_share <- sqrt(sum(v_j[layout$beta_start[2]+seq_len(p_pos)]^2)) / norm_v
  cat(sprintf("\nj=%d: sigma_j=%.5f, ||v_j||=%.5f, beta_share=%.3f\n",
              j, sigma_j, norm_v, beta_share))
  cat(sprintf("  v_j components: beta_pos= %s, field_max=%.5f\n",
              paste(round(v_j[layout$beta_start[2]+seq_len(p_pos)], 5), collapse=","),
              max(abs(v_j[layout$phi_start+seq_len(2*n_s)]))))
  beta_hat <- fit$modes[k,]
  for (h_scale in c(0.1, 0.5, 1.0, 2.0, 5.0)) {
    eps_h <- .Machine$double.eps^(1/5)
    h <- h_scale * eps_h * sigma_j / norm_v
    Lf <- function(s) tulpaObs:::.loglik_cover_joint_at_grid(beta_hat + s*v_j, k, fit, enc_sla, "beta")
    L_p2 <- Lf(2*h); L_p1 <- Lf(h); L_m1 <- Lf(-h); L_m2 <- Lf(-2*h)
    d3 <- (L_p2 - 2*L_p1 + 2*L_m1 - L_m2) / (2*h^3)
    g <- d3 / sigma_j^3
    cat(sprintf("  h_scale=%-4g h=%.3e  d3=%+.4e  gamma=%+.4f\n",
                h_scale, h, d3, g))
  }
  # Also: step purely along e_{beta_j} (no field):
  e_j <- numeric(n_x); e_j[bpos_idx[j]] <- 1
  # Use sigma_j (marginal SD) as the step length along e_j.
  h_e <- .Machine$double.eps^(1/5) * sigma_j
  Lf <- function(s) tulpaObs:::.loglik_cover_joint_at_grid(beta_hat + s*e_j, k, fit, enc_sla, "beta")
  L_p2 <- Lf(2*h_e); L_p1 <- Lf(h_e); L_m1 <- Lf(-h_e); L_m2 <- Lf(-2*h_e)
  d3 <- (L_p2 - 2*L_p1 + 2*L_m1 - L_m2) / (2*h_e^3)
  cat(sprintf("  PURE-e_j step: h=%.3e d3=%+.4e gamma_proxy=%+.4f\n",
              h_e, d3, d3/sigma_j^3))  # NB: this is NOT correct gamma_j -- different direction.
}
