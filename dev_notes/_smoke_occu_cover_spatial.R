# Spatial v2 smoke test for occu_cover(): generate on a chain adjacency,
# fit with nested_laplace, eyeball recovery of alpha + sigma + slopes + field.
# Run from the tulpaObs repo root:
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/_smoke_occu_cover_spatial.R

devtools::load_all(quiet = TRUE)

set.seed(11)
N <- 100L
J <- 6L

# Chain adjacency on N cells.
adj <- matrix(0L, N, N)
for (s in seq_len(N)) {
  if (s > 1L)  adj[s, s - 1L] <- 1L
  if (s < N)   adj[s, s + 1L] <- 1L
}

sim <- simulate_occu_cover(
  N = N, J = J,
  beta_occ = c(qlogis(0.5),  0.7),
  beta_p   = c(0.0,           0.8),     # stronger detection
  beta_pos = c(log(0.20),    -0.4),
  positive = "lognormal", sigma_pos = 0.35,
  adj = adj, sigma = 1.0, alpha = 1.0,  # stronger field signal
  seed = 11L
)

cat(sprintf("Cells occupied (truth): %d / %d\n", sum(sim$truth$z), N))
cat(sprintf("Plots detected: %d / %d\n", sum(sim$y), length(sim$y)))
cat(sprintf("Field sd (true): %.3f\n", sd(sim$truth$f)))

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
  detection = ~ det_cov1,
  positive  = ~ pos_cov1,
  y         = od$y, y_pos = y_pos, visits = od$det.covs,
  method    = "nested_laplace",
  control   = list(verbose = FALSE, max.iter = 500L)
)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("\nFit runtime: %.2f s\n", dt))
cat("Converged: ", fit$convergence$converged, "\n")

cat("\nHyperparameter estimates:\n")
hyper <- data.frame(
  parameter = c("alpha", "sigma (= exp(log_sigma))"),
  truth     = c(1.0, 1.0),
  est       = c(fit$means["alpha"], exp(fit$means["log_sigma"])),
  se        = c(fit$sds["alpha"],
                exp(fit$means["log_sigma"]) * fit$sds["log_sigma"])
)
print(hyper, row.names = FALSE)

cat("\nFixed-effect slope estimates:\n")
fe <- data.frame(
  parameter = c("psi_occ_cov1", "p_det_cov1", "pos_pos_cov1"),
  truth     = c(0.7, 0.5, -0.4),
  est       = c(fit$means["psi_occ_cov1"], fit$means["p_det_cov1"],
                fit$means["pos_pos_cov1"]),
  se        = c(fit$sds["psi_occ_cov1"], fit$sds["p_det_cov1"],
                fit$sds["pos_pos_cov1"])
)
print(fe, row.names = FALSE)

cat(sprintf("\nField recovery: cor(z_hat, f_true) = %.3f\n",
            cor(fit$spatial_field, sim$truth$f)))
