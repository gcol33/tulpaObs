# Sanity: simulate with alpha = 0 (cover arm doesn't see the field).
# A correctly-specified fitter should recover alpha ~ 0 and sigma from psi arm alone.

devtools::load_all(quiet = TRUE)

set.seed(11)
N <- 50L; J <- 6L
adj <- matrix(0L, N, N)
for (s in seq_len(N)) {
  if (s > 1L)  adj[s, s - 1L] <- 1L
  if (s < N)   adj[s, s + 1L] <- 1L
}

sim <- simulate_occu_cover(
  N = N, J = J,
  beta_occ = c(qlogis(0.5),  0.7),
  beta_p   = c(0.0,           0.8),
  beta_pos = c(log(0.20),    -0.4),
  positive = "lognormal", sigma_pos = 0.35,
  adj = adj, sigma = 1.0, alpha = 0.0,   # <- alpha = 0
  seed = 11L
)

long <- data.frame(
  site_id = rep(seq_len(N), each = J),
  visit   = rep(seq_len(J), times = N),
  y       = as.vector(t(sim$y)),
  det_cov1 = sim$visit_data$det_cov1,
  pos_cov1 = sim$visit_data$pos_cov1
)
od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                det.covs = c("det_cov1", "pos_cov1"))
cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

t0 <- Sys.time()
fit <- tobs(
  formula   = ~ occ_cov1 + bym2(graph = adj),
  data      = cell_dat,
  family    = occu_cover("lognormal"),
  detection = ~ det_cov1, positive = ~ pos_cov1,
  y         = od$y, y_pos = y_pos, visits = od$det.covs,
  method    = "nested_laplace",
  control   = list(verbose = FALSE, max.iter = 500L)
)
cat(sprintf("Fit runtime: %.2f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\nalpha = 0 truth case:\n")
cat(sprintf("  alpha:      truth=0.00  est=%.3f  se=%.3f\n",
            fit$means["alpha"], fit$sds["alpha"]))
cat(sprintf("  sigma:      truth=1.00  est=%.3f\n",
            exp(fit$means["log_sigma"])))
cat(sprintf("  field cor:  %.3f\n", cor(fit$spatial_field, sim$truth$f)))
