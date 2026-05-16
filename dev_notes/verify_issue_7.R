#!/usr/bin/env Rscript
# Verify tulpaObs#7 fix: psi-arm Laplace SE should match MC sd of beta_psi_hat
# across seeds. Mirror the D1 sweep settings from the issue (smaller seed
# count for fast iteration during development).
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

N_SEEDS <- 30L
N <- 600L
J <- 6L
BETA_OCC <- c(0.5, 1.2)
BETA_DET <- c(0.0, 0.8)

beta_hat <- matrix(NA_real_, N_SEEDS, 4,
                   dimnames = list(NULL,
                                   c("psi_(Intercept)", "psi_occ_cov1",
                                     "p_(Intercept)",   "p_det_cov1")))
se_hat <- beta_hat

for (s in seq_len(N_SEEDS)) {
  sim <- simulate_occu(N = N, J = J,
                       n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = BETA_OCC,
                       beta_det = BETA_DET,
                       seed = 1000L + s)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, engine = "laplace",
              control = list(verbose = FALSE))
  beta_hat[s, ] <- fit$means[colnames(beta_hat)]
  se_hat[s, ]   <- fit$sds[colnames(se_hat)]
}

# Sanity: estimates should be unbiased, SE should match MC sd, coverage ~ 0.95.
truth <- c(BETA_OCC, BETA_DET)
out <- data.frame(
  param         = colnames(beta_hat),
  truth         = truth,
  mean_est      = colMeans(beta_hat),
  median_se     = apply(se_hat, 2, median),
  mc_sd         = apply(beta_hat, 2, sd),
  se_over_sd    = apply(se_hat, 2, median) / apply(beta_hat, 2, sd),
  cov_95        = colMeans(abs(beta_hat - matrix(truth, N_SEEDS, 4,
                                                  byrow = TRUE)) <
                            1.96 * se_hat),
  stringsAsFactors = FALSE
)
print(out, row.names = FALSE, digits = 3)
