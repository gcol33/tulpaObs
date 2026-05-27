# Discriminating probe: does the AGHQ engine move the community SD AT ALL?
#
# The N=80/J=4 smoke showed n.quad=3 == n.quad=1 (no debias). Hypothesis: at
# 320 obs/species the per-species marginal is already sharply peaked, so the
# per-group Laplace integral is near-exact and AGHQ correctly does nothing --
# the SD gap vs spAbundance is the few-groups finite-sample / IG-prior spread,
# not a per-group Laplace bias. This probe creates the regime AGHQ IS meant to
# fix: tiny per-species data (N=15 sites, J=2 visits) and a low-dim (intercept-
# only, d=2) RE so n.quad=5 (25 nodes) is cheap. If AGHQ lifts the SD here but
# not at N=80, the engine works and the N=80 case legitimately needs no debias.
suppressMessages(devtools::load_all(".", quiet = TRUE))

sim <- simulate_ms_abun(n_species = 12, N = 15, J = 2,
                        n_abund_covs = 0, n_det_covs = 0,
                        mu_lambda = log(6), mu_p = 0.0,
                        sd_lambda = 0.8, sd_p = 0.8, seed = 707)
emp_sd <- c(sd(sim$truth$beta_lambda[, 1]), sd(sim$truth$beta_p[, 1]))

for (nq in c(1L, 5L)) {
  t0 <- Sys.time()
  f <- tobs(~ 1, data = sim$data, y = sim$y, family = ms_abun(),
            detection = ~ 1, species = sim$species, method = "laplace",
            control = list(re.aghq = nq > 1, n.quad = nq, verbose = FALSE))
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cm <- f$ms_community
  cat(sprintf("n.quad=%d (%.1fs): sd_lambda=%.3f  sd_p=%.3f\n",
              nq, dt, cm$sd_lambda[1], cm$sd_p[1]))
}
cat(sprintf("population SD 0.8/0.8 ; realized sd_lambda=%.3f sd_p=%.3f\n",
            emp_sd[1], emp_sd[2]))
