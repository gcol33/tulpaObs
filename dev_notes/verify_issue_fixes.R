## Quick verification that issues #1 and #2 are fixed.
## Run with: Rscript dev_notes/verify_issue_fixes.R

repo <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
suppressPackageStartupMessages(devtools::load_all(repo, quiet = TRUE))

cat("\n=== Issue #1: occu Laplace SEs vary with N ===\n")
for (NN in c(100, 600, 2000)) {
  sim <- simulate_occu(N = NN, J = 6, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.5, 1.2), beta_det = c(0, 0.8), seed = 1)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, engine = "laplace",
              control = list(verbose = FALSE))
  cat(sprintf("N=%4d: sds = %s\n", NN,
              paste(round(fit$sds, 4), collapse = ", ")))
}

cat("\n=== Issue #2: cover joint nested_laplace SE > var-of-means alone ===\n")
set.seed(1)
n_s <- 25; N_per <- 12
adj <- matrix(0L, n_s, n_s)
for (i in seq_len(n_s - 1)) { adj[i, i + 1] <- 1L; adj[i + 1, i] <- 1L }

site_id <- rep(seq_len(n_s), each = N_per)
N <- length(site_id)
phi_true <- rnorm(n_s, sd = 0.7)
x <- rnorm(N)
beta_occ_true <- c(0.4, 0.7); beta_pos_true <- c(0.4, -0.5)
alpha_true <- 1.0
eta_occ <- beta_occ_true[1] + beta_occ_true[2] * x + phi_true[site_id]
psi <- plogis(eta_occ)
y_occ <- rbinom(N, 1, psi)
eta_pos <- beta_pos_true[1] + beta_pos_true[2] * x + alpha_true * phi_true[site_id]
y_pos <- ifelse(y_occ == 1, exp(rnorm(N, eta_pos, 0.3)), 0)
y_pos <- pmin(y_pos, 1 - 1e-6)

dat <- data.frame(region = factor(site_id), x = x)
spec <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")
fit <- tobs(formula = ~ x, data = dat, family = cover("lognormal"),
            y = y_pos, spatial = spec, engine = "nested_laplace",
            control = list(alpha_grid = c(0, 0.5, 1, 1.5, 2)))
cat("se_occ:", round(fit$se_occ, 4), "\n")
cat("se_pos:", round(fit$se_pos, 4), "\n")

## Compare to var-of-means alone
joint <- fit$joint
layout <- joint$arm_layout
p_occ <- layout$p[1]; p_pos <- layout$p[2]
bocc_idx <- layout$beta_start[1] + seq_len(p_occ)
bpos_idx <- layout$beta_start[2] + seq_len(p_pos)
bo <- as.numeric(crossprod(joint$weights, joint$modes[, bocc_idx, drop = FALSE]))
bp <- as.numeric(crossprod(joint$weights, joint$modes[, bpos_idx, drop = FALSE]))
vom_occ <- sqrt(pmax(0, as.numeric(crossprod(joint$weights,
  joint$modes[, bocc_idx, drop = FALSE]^2)) - bo^2))
vom_pos <- sqrt(pmax(0, as.numeric(crossprod(joint$weights,
  joint$modes[, bpos_idx, drop = FALSE]^2)) - bp^2))
cat("var-of-means-only se_occ:", round(vom_occ, 4), "\n")
cat("var-of-means-only se_pos:", round(vom_pos, 4), "\n")
cat("ratio total/var-of-means (occ):", round(fit$se_occ / pmax(vom_occ, 1e-8), 2), "\n")
cat("ratio total/var-of-means (pos):", round(fit$se_pos / pmax(vom_pos, 1e-8), 2), "\n")
cat("(ratios > 1 => Mean-of-Var component contributes; old code had ratio = 1)\n")

cat("\nQ_csc_n =", joint$Q_csc_n,
    "n_grid =", length(joint$Q_csc_p_per_grid), "\n")
