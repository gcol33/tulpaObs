# Focused converged-reference check for the community-SD comparison.
# The 8-seed benchmark left spAbundance's variance-component Rhat high (up to
# 1.77 at 15k iters), so the SD comparison was confounded. Here: 2 seeds, a much
# longer chain, Rhat reported separately for the community MEANS vs the community
# VARIANCES, and the community SD compared tulpa vs converged-sp vs realized.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(spAbundance)
})

n_species <- 12L; N_sites <- 80L; J_visits <- 4L
mu_lambda <- c(log(4), 0.5); mu_p <- c(0.3, -0.4)
sd_lambda <- 0.5; sd_p <- 0.4

# ~3x the budget of the 8-seed run, 3 chains.
n.batch <- 1500L; batch.length <- 25L; n.chains <- 3L
n.burn <- 18750L; n.thin <- 10L
priors <- list(beta.comm.normal = list(mean = 0, var = 100),
               alpha.comm.normal = list(mean = 0, var = 100),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

for (seed in c(2001L, 2003L)) {
  sim <- simulate_ms_abun(n_species = n_species, N = N_sites, J = J_visits,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = mu_lambda, mu_p = mu_p,
                          sd_lambda = sd_lambda, sd_p = sd_p, seed = seed)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
              detection = ~ det_cov1, species = sim$species, method = "laplace",
              control = list(verbose = FALSE))
  cm <- fit$ms_community
  tulpa_sd <- c(cm$sd_lambda, cm$sd_p)

  sp_data <- list(y = aperm(sim$y, c(3, 1, 2)),
                  abund.covs = data.frame(abund_cov1 = sim$data$abund_cov1),
                  det.covs = list(det_cov1 = sim$data$det_cov1))
  t0 <- Sys.time()
  sp <- msNMix(~ abund_cov1, ~ det_cov1, data = sp_data, priors = priors,
               family = "Poisson", n.batch = n.batch, batch.length = batch.length,
               n.burn = n.burn, n.thin = n.thin, n.chains = n.chains, verbose = FALSE)
  t_sp <- as.numeric(Sys.time() - t0, units = "secs")

  sp_sd <- c(apply(sqrt(as.matrix(sp$tau.sq.beta.samples)),  2, median),
             apply(sqrt(as.matrix(sp$tau.sq.alpha.samples)), 2, median))
  emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))
  nm <- names(fit$means)

  rhat_mean <- max(c(sp$rhat$beta.comm, sp$rhat$alpha.comm), na.rm = TRUE)
  rhat_var  <- max(c(sp$rhat$tau.sq.beta, sp$rhat$tau.sq.alpha), na.rm = TRUE)

  cat(sprintf("\n=== seed %d (sp %.0fs, %d x %d x %d chains) ===\n",
              seed, t_sp, n.batch, batch.length, n.chains))
  cat(sprintf("Rhat: community MEANS max=%.3f | community VARIANCES max=%.3f\n",
              rhat_mean, rhat_var))
  cat("# Community SDs: population | realized(n=12) | tulpa | sp(converged)\n")
  print(round(data.frame(population = c(rep(sd_lambda, 2), rep(sd_p, 2)),
                         realized = emp_sd, tulpa = tulpa_sd, sp = sp_sd,
                         tulpa_vs_sp = tulpa_sd - sp_sd,
                         row.names = nm), 3))
}
