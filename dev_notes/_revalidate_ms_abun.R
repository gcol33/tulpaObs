# Re-validate ms_abun after the native NMixCommunityOracle: timing vs the 482s
# R-closure baseline, and the full test-ms-abun.R (incl. the 20-seed coverage,
# previously ~2h on the R-closure path).
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))

cat("installed tulpa:", as.character(utils::packageVersion("tulpa")), "\n")

sim <- simulate_ms_abun(n_species = 12, N = 80, J = 4,
                        n_abund_covs = 1, n_det_covs = 1,
                        mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.4),
                        sd_lambda = 0.5, sd_p = 0.4, seed = 2001)
for (nq in c(1L, 3L, 5L)) {
  t0 <- Sys.time()
  f <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
            detection = ~ det_cov1, species = sim$species, method = "laplace",
            control = list(re.aghq = nq > 1, n.quad = nq, verbose = FALSE))
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cm <- f$ms_community
  cat(sprintf("n.quad=%d  %.1fs (was ~482s/2951s on R-closure)  SDs %.3f/%.3f/%.3f/%.3f\n",
              nq, dt, cm$sd_lambda[1], cm$sd_lambda[2], cm$sd_p[1], cm$sd_p[2]))
}

cat("\n================ test-ms-abun.R (full, NOT_CRAN) ================\n")
testthat::test_file("tests/testthat/test-ms-abun.R", reporter = "summary")
