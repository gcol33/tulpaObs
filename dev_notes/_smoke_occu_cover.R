# Quick smoke test for occu_cover() v1: simulate, fit, eyeball recovery.
# Run from the tulpaObs repo root:
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/_smoke_occu_cover.R

devtools::load_all(quiet = TRUE)

set.seed(1)
sim <- simulate_occu_cover(
  N = 200L, J = 5L,
  beta_occ = c(qlogis(0.4), 0.9),
  beta_p   = c(0.0,         0.6),
  beta_pos = c(log(0.10),  -0.4),
  sigma_pos = 0.4, positive = "lognormal", seed = 1L
)

cat(sprintf("Sites occupied (truth): %d / %d\n", sum(sim$truth$z), length(sim$truth$z)))
cat(sprintf("Plots detected: %d / %d\n", sum(sim$y), length(sim$y)))

long <- data.frame(
  site_id = rep(seq_len(nrow(sim$y)), each = ncol(sim$y)),
  visit   = rep(seq_len(ncol(sim$y)), times = nrow(sim$y)),
  y       = as.vector(t(sim$y)),
  det_cov1 = sim$visit_data$det_cov1,
  pos_cov1 = sim$visit_data$pos_cov1
)
od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                det.covs = c("det_cov1", "pos_cov1"))
cell_dat <- cbind(data.frame(site_id = seq_len(nrow(sim$y))), sim$data)

y_pos <- sim$y_pos
y_pos[is.na(y_pos)] <- 0

t0 <- Sys.time()
fit <- tobs(
  formula   = ~ occ_cov1, data = cell_dat,
  family    = occu_cover("lognormal"),
  detection = ~ det_cov1, positive = ~ pos_cov1,
  y         = od$y, y_pos = y_pos, visits = od$det.covs,
  method    = "laplace",
  control   = list(verbose = FALSE, max.iter = 500L)
)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("\nFit runtime: %.2f s\n", dt))
cat("Converged: ", fit$convergence$converged, "\n")
cat("\nEstimates vs truth:\n")
truth <- c(
  "psi_(Intercept)" = sim$truth$beta_occ[1],
  "psi_occ_cov1"    = sim$truth$beta_occ[2],
  "p_(Intercept)"   = sim$truth$beta_p[1],
  "p_det_cov1"      = sim$truth$beta_p[2],
  "pos_(Intercept)" = sim$truth$beta_pos[1],
  "pos_pos_cov1"    = sim$truth$beta_pos[2],
  "log_sigma_pos"   = log(sim$truth$sigma_pos)
)
print(data.frame(
  truth = truth,
  est   = fit$means[names(truth)],
  se    = fit$sds[names(truth)]
))
