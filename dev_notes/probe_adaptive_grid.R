devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
set.seed(1)
N <- 30L; J <- 4L
adj <- matrix(0L, N, N)
for (s in seq_len(N)) {
  if (s > 1L) adj[s, s - 1L] <- 1L
  if (s < N)  adj[s, s + 1L] <- 1L
}
sim <- simulate_occu_cover(
  N = N, J = J, positive = "lognormal",
  adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
)
long <- data.frame(
  site_id  = rep(seq_len(N), each = J),
  visit    = rep(seq_len(J), times = N),
  y        = as.vector(t(sim$y)),
  det_cov1 = sim$visit_data$det_cov1,
  pos_cov1 = sim$visit_data$pos_cov1
)
od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                 det.covs = c("det_cov1", "pos_cov1"))
cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
t0 <- proc.time()
fit <- suppressWarnings(tobs(
  formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
  family = occu_cover("lognormal"),
  detection = ~ det_cov1, positive = ~ pos_cov1,
  y = od$y, y_pos = y_pos, visits = od$det.covs,
  method = "nested_laplace",
  control = list(verbose = FALSE, max.iter = 500L,
                 engine = "joint_coupled",
                 adaptive.grid = TRUE)
))
cat("elapsed=", (proc.time() - t0)[3], "s\n")
cat("n_grid=", nrow(fit$joint_fit$theta_grid), "\n")
cat("finite log_marginal=", sum(is.finite(fit$joint_fit$log_marginal)),
    "/", length(fit$joint_fit$log_marginal), "\n")
cat("sigma_mean=", fit$means[["sigma"]], "\n")
cat("alpha_mean=", fit$means[["alpha"]], "\n")
cat("intercept_sds=",
    round(fit$sds[c("psi_(Intercept)", "p_(Intercept)", "pos_(Intercept)")], 3),
    "\n")
