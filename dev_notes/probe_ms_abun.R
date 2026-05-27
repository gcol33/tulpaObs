# Probe: community / multispecies N-mixture (ms_abun) through the public tobs()
# API -> tulpa's C++ Laplace-EM. Recovery of the community means (mu_lambda,
# mu_p), community SDs, and per-species coefficients vs simulated truth.
suppressMessages(devtools::load_all("."))

set.seed(2026)
sim <- simulate_ms_abun(n_species = 15, N = 120, J = 5,
                        n_abund_covs = 1, n_det_covs = 1,
                        mu_lambda = c(log(4), 0.5),
                        mu_p      = c(0.3, -0.4),
                        sd_lambda = 0.5, sd_p = 0.4,
                        mixture = "poisson", seed = 2026)

t0  <- Sys.time()
fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
            family = ms_abun(), detection = ~ det_cov1,
            species = sim$species, method = "laplace",
            control = list(verbose = FALSE))
cat(sprintf("fit time: %.1fs   converged: %s   method: %s\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            isTRUE(fit$convergence$converged), fit$method))

cm <- fit$ms_community
cat("\n=== community means (estimate vs truth) ===\n")
est_mu <- fit$means
truth_mu <- c(sim$truth$mu_lambda, sim$truth$mu_p)
names(truth_mu) <- names(est_mu)
print(round(rbind(estimate = est_mu, truth = truth_mu, se = fit$sds,
                  z = (est_mu - truth_mu) / fit$sds), 3))

cat("\n=== community SDs (estimate vs realized sample vs population) ===\n")
emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))
est_sd   <- c(cm$sd_lambda, cm$sd_p)
truth_sd <- c(sim$truth$sd_lambda, sim$truth$sd_p)
names(truth_sd) <- names(emp_sd) <- names(est_sd)
print(round(rbind(estimate = est_sd, `realized (sample)` = emp_sd,
                  population = truth_sd), 3))

cat("\n=== per-species coefficient recovery (cor with truth) ===\n")
cat("lambda arm:", round(diag(cor(cm$coef_lambda, sim$truth$beta_lambda)), 3), "\n")
cat("p arm:     ", round(diag(cor(cm$coef_p,      sim$truth$beta_p)), 3), "\n")

cat("\n=== 95% CI coverage of the community means ===\n")
lo <- est_mu - 1.96 * fit$sds; hi <- est_mu + 1.96 * fit$sds
print(data.frame(truth = round(truth_mu, 3),
                 covered = truth_mu >= lo & truth_mu <= hi))

# S3 smoke: coef/ranef/print should not error on the new family.
cat("\n=== S3 smoke ===\n")
cat("class:", paste(class(fit), collapse = "/"), "\n")
cat("coef() names:", paste(names(coef(fit)), collapse = ", "), "\n")
