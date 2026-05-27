# Validate the AGHQ debias of the community N-mixture variance components.
#
# The 8-seed head-to-head confirmed tulpaObs ms_abun community SDs ran ~10-19%
# below a properly-converged spAbundance::msNMix reference (the small-cluster
# Laplace nAGQ=1 attenuation). This re-runs tulpa at n_quad in {1, 3, 5} against
# the same converged reference and measures (a) how much each n_quad lifts the
# SDs toward the reference and (b) the per-fit runtime, to pick the ms_abun
# default. Community MEANS and per-fit timing are reported alongside.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(spAbundance)
})

n_species <- 12L; N_sites <- 80L; J_visits <- 4L
mu_lambda <- c(log(4), 0.5); mu_p <- c(0.3, -0.4)
sd_lambda <- 0.5; sd_p <- 0.4
seeds   <- c(2001L, 2003L)
n_quads <- c(1L, 3L, 5L)

# Converged spAbundance reference (Rhat <= 1.006 in the earlier ref run).
n.batch <- 1500L; batch.length <- 25L; n.chains <- 3L
n.burn <- 18750L; n.thin <- 10L
priors <- list(beta.comm.normal = list(mean = 0, var = 100),
               alpha.comm.normal = list(mean = 0, var = 100),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

rows <- list()
for (seed in seeds) {
  sim <- simulate_ms_abun(n_species = n_species, N = N_sites, J = J_visits,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = mu_lambda, mu_p = mu_p,
                          sd_lambda = sd_lambda, sd_p = sd_p, seed = seed)
  nm  <- NULL
  emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))

  tulpa_sd <- list()
  for (nq in n_quads) {
    t0 <- Sys.time()
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
                detection = ~ det_cov1, species = sim$species, method = "laplace",
                control = list(re.aghq = nq > 1, n.quad = nq, verbose = FALSE))
    dt  <- as.numeric(Sys.time() - t0, units = "secs")
    cm  <- fit$ms_community
    nm  <- names(fit$means)
    tulpa_sd[[as.character(nq)]] <- list(sd = c(cm$sd_lambda, cm$sd_p),
                                         means = fit$means, secs = dt)
    cat(sprintf("seed %d  n.quad=%d  %.1fs  (stored n_quad=%d)\n",
                seed, nq, dt, cm$n_quad))
  }

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
  rhat_var <- max(c(sp$rhat$tau.sq.beta, sp$rhat$tau.sq.alpha), na.rm = TRUE)
  cat(sprintf("seed %d  spAbundance %.0fs  (var Rhat max=%.3f)\n", seed, t_sp, rhat_var))

  for (k in seq_along(nm)) {
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, param = nm[k], population = c(rep(sd_lambda, 2), rep(sd_p, 2))[k],
      realized = emp_sd[k],
      tulpa_nq1 = tulpa_sd[["1"]]$sd[k],
      tulpa_nq3 = tulpa_sd[["3"]]$sd[k],
      tulpa_nq5 = tulpa_sd[["5"]]$sd[k],
      sp_ref = sp_sd[k], sp_rhat_var = rhat_var,
      stringsAsFactors = FALSE)
  }
  attr(rows, "timing") <- rbind(attr(rows, "timing"),
    data.frame(seed = seed,
               nq1 = tulpa_sd[["1"]]$secs, nq3 = tulpa_sd[["3"]]$secs,
               nq5 = tulpa_sd[["5"]]$secs, sp = t_sp))
}

df <- do.call(rbind, rows)
out_tsv <- "dev_notes/_validate_ms_abun_aghq.tsv"
write.table(df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# Bias of each n_quad vs the converged reference, averaged over seeds/params.
rel <- function(col) mean(abs(df[[col]] - df$sp_ref) / df$sp_ref)
cat("\n================ SD comparison vs converged spAbundance ================\n")
print(round(df[, c("seed","param","realized","tulpa_nq1","tulpa_nq3","tulpa_nq5","sp_ref")], 3))
cat("\nmean |tulpa - sp_ref| / sp_ref :\n")
cat(sprintf("  n.quad=1 : %.1f%%\n", 100 * rel("tulpa_nq1")))
cat(sprintf("  n.quad=3 : %.1f%%\n", 100 * rel("tulpa_nq3")))
cat(sprintf("  n.quad=5 : %.1f%%\n", 100 * rel("tulpa_nq5")))
cat("\ntiming (s):\n"); print(attr(rows, "timing"))
cat("\nwrote", out_tsv, "\n")
