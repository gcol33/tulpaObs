#!/usr/bin/env Rscript
suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
N_SEEDS <- 20L
beta_hat <- matrix(NA_real_, N_SEEDS, 4)
se_hat <- beta_hat
nms <- c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)", "p_det_cov1")
for (s in seq_len(N_SEEDS)) {
  sim <- simulate_occu(N = 600L, J = 6L, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.5, 1.2), beta_det = c(0, 0.8),
                       seed = 2000L + s)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, engine = "laplace",
              priors = FALSE, control = list(verbose = FALSE))
  beta_hat[s, ] <- fit$means[nms]
  se_hat[s, ]   <- fit$sds[nms]
}
res <- data.frame(
  param      = nms,
  truth      = c(0.5, 1.2, 0, 0.8),
  mean_est   = round(colMeans(beta_hat), 3),
  median_se  = round(apply(se_hat, 2, median), 3),
  mc_sd      = round(apply(beta_hat, 2, sd), 3),
  se_over_sd = round(apply(se_hat, 2, median) / apply(beta_hat, 2, sd), 3)
)
print(res, row.names = FALSE)
