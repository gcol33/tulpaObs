# Re-validate ms_abun after the EM fast-path at n_quad=1 (the default). Timing
# vs the 482s (R-closure) and 201s (native-oracle joint-BFGS) baselines, plus the
# full test-ms-abun.R (incl. the 20-seed coverage).
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
cat("installed tulpa:", as.character(utils::packageVersion("tulpa")), "\n")

sim <- simulate_ms_abun(n_species = 12, N = 80, J = 4,
                        n_abund_covs = 1, n_det_covs = 1,
                        mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.4),
                        sd_lambda = 0.5, sd_p = 0.4, seed = 2001)
t0 <- Sys.time()
f <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
          detection = ~ det_cov1, species = sim$species, method = "laplace",
          control = list(verbose = FALSE))   # default: EM, n_quad=1
dt <- as.numeric(Sys.time() - t0, units = "secs")
cm <- f$ms_community
cat(sprintf("DEFAULT (EM, n_quad=1): %.1fs  (was 482s R-closure / 201s joint-BFGS)\n", dt))
cat(sprintf("  means: %s\n", paste(sprintf("%.3f", f$means), collapse = " ")))
cat(sprintf("  SDs  : %.3f/%.3f/%.3f/%.3f\n", cm$sd_lambda[1], cm$sd_lambda[2],
            cm$sd_p[1], cm$sd_p[2]))

cat("\n================ test-ms-abun.R (full, NOT_CRAN) ================\n")
tt <- Sys.time()
testthat::test_file("tests/testthat/test-ms-abun.R", reporter = "summary")
cat(sprintf("\ntest-ms-abun.R wall time: %.0fs\n",
            as.numeric(Sys.time() - tt, units = "secs")))
