# Minimal reproducer for the tulpa_nested_laplace_joint() copy-arm
# regression with gaussian likelihood. Even at the *true* noise scale
# (phi = sigma_pos_true), the joint integrand assigns the highest marginal
# log-likelihood to the alpha = 0 cell, despite alpha_true = 1.
#
# Same engine on a beta copy arm (cover('beta')) recovers alpha at
# 80-97% nominal coverage (INLAabun Demo 3), so the bug is gaussian-specific.

suppressPackageStartupMessages(library(tulpa))

chain_adj <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) for (j in setdiff(c(s - 1L, s + 1L),
                                            c(0L, n_s + 1L))) adj[s, j] <- 1L
  adj
}

set.seed(50101L)
N <- 300L; n_s <- 25L
sigma_field_true <- 0.6; rho_true <- 0.7
alpha_true <- 1.0; sigma_pos_true <- 0.3
adj <- chain_adj(n_s)

region <- sample.int(n_s, N, replace = TRUE)
pf <- rnorm(n_s); tf <- rnorm(n_s)
w <- sigma_field_true * (sqrt(rho_true) * pf + sqrt(1 - rho_true) * tf)
x <- rnorm(N)
betaO <- c(-0.3, 0.7); betaP <- c(-1.5, 0.3)

occur <- rbinom(N, 1L, plogis(betaO[1] + betaO[2] * x + w[region]))
log_y <- rnorm(N, betaP[1] + betaP[2] * x + alpha_true * w[region],
               sigma_pos_true)
is_pos <- occur == 1L
y_pos <- log_y[is_pos]
X_pos <- model.matrix(~ x[is_pos])
spi_pos <- as.integer(region[is_pos])

X_occ <- model.matrix(~ x)
spi_occ <- as.integer(region)

prior <- list(
  type            = "bym2",
  n_spatial_units = n_s,
  scale_factor    = 1.0,
  sigma_grid      = c(0.3, 0.6, 0.9),
  rho_grid        = c(0.5, 0.7, 0.9)
)
# Adjacency CSR.
adj_rp <- integer(n_s + 1L); adj_ci <- integer(0L); adj_nb <- integer(n_s)
for (s in seq_len(n_s)) {
  nbr <- which(adj[s, ] == 1L) - 1L
  adj_ci <- c(adj_ci, nbr); adj_nb[s] <- length(nbr)
  adj_rp[s + 1L] <- length(adj_ci)
}
prior$adj_row_ptr <- adj_rp
prior$adj_col_idx <- adj_ci
prior$n_neighbors <- adj_nb

arm_occ <- list(y = occur, n_trials = rep(1L, N), X = X_occ,
                spatial_idx = spi_occ, re_idx = rep(0, N),
                n_re_groups = 0L, sigma_re = 1.0,
                family = "binomial", phi = 1.0)

cat(sprintf("\nTruth: alpha = %.1f, sigma_field = %.1f, rho = %.1f, sigma_pos = %.1f\n",
            alpha_true, sigma_field_true, rho_true, sigma_pos_true))
cat("Max log_marginal per (alpha) across (sigma, rho) at three phi settings:\n\n")
cat(sprintf("%-6s | %s\n", "phi",
            paste(sprintf("alpha=%.1f", c(0, 0.5, 1.0, 1.5)), collapse = " ")))
cat(strrep("-", 60), "\n", sep = "")

for (test_phi in c(1.0, 0.62, sigma_pos_true)) {
  arm_pos <- list(y = y_pos, n_trials = rep(1L, length(y_pos)),
                  X = X_pos, spatial_idx = spi_pos,
                  re_idx = rep(0, length(y_pos)),
                  n_re_groups = 0L, sigma_re = 1.0,
                  family = "gaussian", phi = test_phi)
  fit <- tulpa_nested_laplace_joint(
    responses = list(occ = arm_occ, pos = arm_pos),
    prior     = prior,
    copy      = list(arm = "pos", alpha_grid = c(0, 0.5, 1.0, 1.5)),
    max_iter  = 50L, tol = 1e-6, n_threads = 1L
  )
  df <- as.data.frame(cbind(fit$theta_grid, log_marginal = fit$log_marginal))
  by_alpha <- aggregate(log_marginal ~ alpha, df, max)
  cat(sprintf("%-6.2f | %s\n", test_phi,
              paste(sprintf("%.1f", by_alpha$log_marginal), collapse = "  ")))
}
cat("\nThe engine prefers alpha = 0 at every phi tested, including phi = sigma_pos_true.\n")
cat("Expected: highest log_marginal cell at (sigma=0.6, rho=0.7, alpha=1.0).\n")
