# Verify issue #3 fix via the public tobs() API after detection-weighting fix.
# Expect: psi(0) close to 0.5, psi(1) close to 1.2, p(0) close to 0.0, p(1) close to 0.8.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

sim_occu <- function(N = 600, J = 6, beta_occ = c(0.5, 1.2),
                     beta_det = c(0, 0.8), seed = 42) {
  set.seed(seed)
  occ_cov1 <- rnorm(N)
  det_cov1 <- rnorm(N)
  X_occ <- cbind(1, occ_cov1)
  X_det <- cbind(1, det_cov1)
  psi <- plogis(as.vector(X_occ %*% beta_occ))
  p   <- plogis(as.vector(X_det %*% beta_det))
  z   <- rbinom(N, 1, psi)
  y   <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p[i])
  list(y = y, occ_cov1 = occ_cov1, det_cov1 = det_cov1)
}

cat("Truth: beta_occ = (0.500, 1.200), beta_det = (0.000, 0.800)\n\n")
cat(sprintf("%5s | %-32s\n", "seed", "tobs(occu(), engine='laplace')"))
cat(sprintf("%5s-+-%-32s\n", "-----", "--------------------------------"))
for (seed in c(42, 99, 123, 7, 2026)) {
  sim <- sim_occu(seed = seed)
  site_data <- data.frame(occ_cov1 = sim$occ_cov1, det_cov1 = sim$det_cov1)
  fit <- tobs(
    formula   = ~ occ_cov1,
    detection = ~ det_cov1,
    data      = site_data,
    family    = occu(),
    y         = sim$y,
    engine    = "laplace",
    control   = list(verbose = FALSE)
  )
  bo <- fit$means[grep("^psi_", names(fit$means))]
  bd <- fit$means[grep("^p_",   names(fit$means))]
  cat(sprintf("%5d | psi=(%.2f, %.2f) p=(%.2f, %.2f)\n",
              seed, bo[1], bo[2], bd[1], bd[2]))
}
