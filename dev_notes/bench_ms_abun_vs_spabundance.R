# =============================================================================
# Head-to-head: tulpaObs ms_abun() vs spAbundance::msNMix()
#
# The credibility check named in gcol33/tulpa#31's acceptance criteria:
# fit BOTH the C++ Laplace-EM (tulpaObs) and the reference Adaptive-MCMC
# (spAbundance) community N-mixture on a *shared* simulated panel, and compare
# their estimates of the community means (mu_lambda, mu_p), the community SDs
# (sqrt tau.sq), and the per-species coefficients -- both against each other and
# against the simulated truth.
#
# Community structure matches: simulate_ms_abun() draws per-species coefficients
# from independent (diagonal) Gaussian community priors, and msNMix estimates a
# diagonal community variance (tau.sq.beta / tau.sq.alpha). tulpaObs estimates a
# full Sigma; we compare its sqrt-diagonal.
#
# Data-format bridge: tulpaObs y is [sites x visits x species]; spAbundance y is
# [species x sites x visits]. det.covs is a named list (site-level vectors).
#
#   Rscript dev_notes/bench_ms_abun_vs_spabundance.R
# =============================================================================

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(spAbundance)
})

# ---- benchmark configuration ------------------------------------------------
N_SEEDS   <- 8L
n_species <- 12L
N_sites   <- 80L
J_visits  <- 4L
mu_lambda <- c(log(4), 0.5)     # (intercept, abund slope) on log scale
mu_p      <- c(0.3, -0.4)       # (intercept, det slope)   on logit scale
sd_lambda <- 0.5
sd_p      <- 0.4

# spAbundance Adaptive-MCMC budget (reference posterior).
n.batch      <- 600L            # 600 x 25 = 15000 iterations / chain
batch.length <- 25L
n.chains     <- 2L
n.burn       <- 7500L
n.thin       <- 5L              # -> (15000-7500)/5 = 1500 post draws / chain

# Vague community priors (spAbundance defaults are already weakly informative;
# state them explicitly so the reference is reproducible).
priors <- list(
  beta.comm.normal  = list(mean = 0, var = 100),
  alpha.comm.normal = list(mean = 0, var = 100),
  tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig   = list(a = 0.1, b = 0.1)
)

post_mean <- function(s) colMeans(as.matrix(s))
post_med  <- function(s) apply(as.matrix(s), 2, stats::median)

# ---- one shared-panel comparison -------------------------------------------
run_one <- function(seed, detailed = FALSE) {
  sim <- simulate_ms_abun(n_species = n_species, N = N_sites, J = J_visits,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = mu_lambda, mu_p = mu_p,
                          sd_lambda = sd_lambda, sd_p = sd_p,
                          mixture = "poisson", seed = seed)

  truth_mu <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  p_lam <- length(sim$truth$mu_lambda); p_p <- length(sim$truth$mu_p)

  # --- tulpaObs Laplace-EM ---
  t0 <- Sys.time()
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(), detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE))
  t_tulpa <- as.numeric(Sys.time() - t0, units = "secs")
  cm <- fit$ms_community
  tulpa_mu <- fit$means
  tulpa_sd <- c(cm$sd_lambda, cm$sd_p)
  tulpa_coef_lam <- cm$coef_lambda      # [S x p_lam]
  tulpa_coef_p   <- cm$coef_p           # [S x p_p]

  # --- spAbundance msNMix (Adaptive MCMC) ---
  sp_data <- list(
    y          = aperm(sim$y, c(3, 1, 2)),                  # [species x sites x visits]
    abund.covs = data.frame(abund_cov1 = sim$data$abund_cov1),
    det.covs   = list(det_cov1 = sim$data$det_cov1)
  )
  t0 <- Sys.time()
  sp <- msNMix(abund.formula = ~ abund_cov1, det.formula = ~ det_cov1,
               data = sp_data, priors = priors, family = "Poisson",
               n.batch = n.batch, batch.length = batch.length,
               n.burn = n.burn, n.thin = n.thin, n.chains = n.chains,
               verbose = FALSE)
  t_sp <- as.numeric(Sys.time() - t0, units = "secs")

  sp_mu <- c(post_mean(sp$beta.comm.samples), post_mean(sp$alpha.comm.samples))
  # community SD = posterior median of sqrt(tau.sq) per coefficient.
  sp_sd <- c(post_med(sqrt(as.matrix(sp$tau.sq.beta.samples))),
             post_med(sqrt(as.matrix(sp$tau.sq.alpha.samples))))
  # per-species coefs: colnames are coef-major (<coef>-sp1, <coef>-sp2, ...).
  sp_coef_lam <- matrix(post_mean(sp$beta.samples),  nrow = n_species, ncol = p_lam)
  sp_coef_p   <- matrix(post_mean(sp$alpha.samples), nrow = n_species, ncol = p_p)

  rhat_max <- suppressWarnings(max(unlist(sp$rhat), na.rm = TRUE))

  # realized (sample) community SD = SD of the actually-drawn per-species coefs.
  emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))

  if (detailed) {
    nm <- names(tulpa_mu)
    cat(sprintf("\n--- DETAILED PANEL (seed %d) ---\n", seed))
    cat(sprintf("converged: tulpa=%s   spAbundance max Rhat=%.3f\n",
                isTRUE(fit$convergence$converged), rhat_max))
    cat(sprintf("wall time: tulpa=%.2fs   spAbundance=%.1fs   (%.0fx)\n",
                t_tulpa, t_sp, t_sp / t_tulpa))
    cat("\n# Community means (mu): truth | tulpa | spAbundance\n")
    print(round(data.frame(truth = truth_mu, tulpa = tulpa_mu,
                           spAbundance = sp_mu,
                           tulpa_err = tulpa_mu - truth_mu,
                           sp_err = sp_mu - truth_mu,
                           row.names = nm), 3))
    cat("\n# Community SDs: population | realized | tulpa | spAbundance\n")
    print(round(data.frame(
      population = c(rep(sd_lambda, p_lam), rep(sd_p, p_p)),
      realized = emp_sd, tulpa = tulpa_sd, spAbundance = sp_sd,
      row.names = nm), 3))
    cat("\n# Per-species coef agreement (correlation across species)\n")
    for (c in seq_len(p_lam))
      cat(sprintf("  lambda[%d]: cor(tulpa,sp)=%.3f  cor(tulpa,truth)=%.3f  cor(sp,truth)=%.3f\n",
                  c, cor(tulpa_coef_lam[, c], sp_coef_lam[, c]),
                  cor(tulpa_coef_lam[, c], sim$truth$beta_lambda[, c]),
                  cor(sp_coef_lam[, c], sim$truth$beta_lambda[, c])))
    for (c in seq_len(p_p))
      cat(sprintf("  p[%d]:      cor(tulpa,sp)=%.3f  cor(tulpa,truth)=%.3f  cor(sp,truth)=%.3f\n",
                  c, cor(tulpa_coef_p[, c], sp_coef_p[, c]),
                  cor(tulpa_coef_p[, c], sim$truth$beta_p[, c]),
                  cor(sp_coef_p[, c], sim$truth$beta_p[, c])))
  }

  # per-species cross-method agreement summary (pooled over coefs/arms)
  coef_cor_tulpa_sp <- c(
    sapply(seq_len(p_lam), function(c) cor(tulpa_coef_lam[, c], sp_coef_lam[, c])),
    sapply(seq_len(p_p),   function(c) cor(tulpa_coef_p[, c],   sp_coef_p[, c])))
  coef_mad_tulpa_sp <- c(
    sapply(seq_len(p_lam), function(c) mean(abs(tulpa_coef_lam[, c] - sp_coef_lam[, c]))),
    sapply(seq_len(p_p),   function(c) mean(abs(tulpa_coef_p[, c]   - sp_coef_p[, c]))))

  data.frame(
    seed = seed,
    mu_max_abs_diff = max(abs(tulpa_mu - sp_mu)),
    mu_tulpa_rmse   = sqrt(mean((tulpa_mu - truth_mu)^2)),
    mu_sp_rmse      = sqrt(mean((sp_mu - truth_mu)^2)),
    sd_max_abs_diff = max(abs(tulpa_sd - sp_sd)),
    coef_cor_min    = min(coef_cor_tulpa_sp),
    coef_mad_max    = max(coef_mad_tulpa_sp),
    rhat_max        = rhat_max,
    t_tulpa = t_tulpa, t_sp = t_sp,
    stringsAsFactors = FALSE)
}

# ---- run ---------------------------------------------------------------------
cat("=============================================================\n")
cat(" ms_abun() [tulpaObs Laplace-EM]  vs  msNMix() [spAbundance MCMC]\n")
cat(sprintf(" S=%d species, N=%d sites, J=%d visits | %d seeds\n",
            n_species, N_sites, J_visits, N_SEEDS))
cat(sprintf(" MCMC: %d x %d iters x %d chains, burn %d, thin %d\n",
            n.batch, batch.length, n.chains, n.burn, n.thin))
cat("=============================================================\n")

res <- vector("list", N_SEEDS)
for (k in seq_len(N_SEEDS)) {
  res[[k]] <- run_one(seed = 2000L + k, detailed = (k == 1L))
  cat(sprintf("seed %d done (tulpa %.2fs / sp %.1fs)\n",
              2000L + k, res[[k]]$t_tulpa, res[[k]]$t_sp))
}
res <- do.call(rbind, res)

cat("\n=============================================================\n")
cat(" AGGREGATE OVER", N_SEEDS, "SEEDS (median [min, max])\n")
cat("=============================================================\n")
summ <- function(x) sprintf("%.3f  [%.3f, %.3f]", median(x), min(x), max(x))
cat(sprintf("Community-mean max |tulpa - sp|     : %s\n", summ(res$mu_max_abs_diff)))
cat(sprintf("Community-mean RMSE vs truth, tulpa : %s\n", summ(res$mu_tulpa_rmse)))
cat(sprintf("Community-mean RMSE vs truth, sp    : %s\n", summ(res$mu_sp_rmse)))
cat(sprintf("Community-SD  max |tulpa - sp|      : %s\n", summ(res$sd_max_abs_diff)))
cat(sprintf("Per-species coef cor(tulpa, sp) min : %s\n", summ(res$coef_cor_min)))
cat(sprintf("Per-species coef max |tulpa - sp|   : %s\n", summ(res$coef_mad_max)))
cat(sprintf("spAbundance max Rhat                : %s\n", summ(res$rhat_max)))
cat(sprintf("Wall time  tulpa (s)                : %s\n", summ(res$t_tulpa)))
cat(sprintf("Wall time  spAbundance (s)          : %s\n", summ(res$t_sp)))
cat(sprintf("Speedup  (sp / tulpa)               : %s\n", summ(res$t_sp / res$t_tulpa)))

out_tsv <- "dev_notes/_bench_ms_abun_vs_spabundance.tsv"
utils::write.table(res, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nper-seed table written to", out_tsv, "\n")
