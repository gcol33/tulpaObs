# Smoke: confirm the n_quad wiring runs end-to-end through tobs() and that a
# higher n_quad lifts the (Laplace-attenuated) community SDs.
suppressMessages(devtools::load_all(".", quiet = TRUE))

sim <- simulate_ms_abun(n_species = 12, N = 80, J = 4,
                        n_abund_covs = 1, n_det_covs = 1,
                        mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.4),
                        sd_lambda = 0.5, sd_p = 0.4, seed = 2001)

fit_one <- function(nq) {
  t0 <- Sys.time()
  f <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
            detection = ~ det_cov1, species = sim$species, method = "laplace",
            control = list(re.aghq = nq > 1, n.quad = nq, verbose = FALSE))
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cm <- f$ms_community
  list(nq = nq, secs = dt, means = f$means,
       sd = c(cm$sd_lambda, cm$sd_p), stored_nq = cm$n_quad)
}

nm <- c("lambda_(Int)", "lambda_cov1", "p_(Int)", "p_cov1")
for (nq in c(1L, 3L)) {
  r <- fit_one(nq)
  cat(sprintf("\n=== n.quad = %d (stored %d, %.1fs) ===\n", nq, r$stored_nq, r$secs))
  cat("means:", paste(sprintf("%.3f", r$means), collapse = " "), "\n")
  cat("SDs  :", paste(sprintf("%s=%.3f", nm, r$sd), collapse = "  "), "\n")
}
cat("\npopulation SD: lambda 0.5/0.5, p 0.4/0.4\n")
