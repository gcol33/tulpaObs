# Probe: abun(mixture = "negbin") wiring through tulpaObs.
# Non-spatial NB recovery + the log_r coordinate surfacing in coef()/vcov().
suppressMessages(devtools::load_all("."))

set.seed(101)
sim <- simulate_abun(N = 300, J = 5, n_abund_covs = 2, n_det_covs = 1,
                     beta_lambda = c(log(4), 0.6, -0.4),
                     beta_p = c(0.3, 0.5),
                     mixture = "negbin", size = 1.5)

dat <- sim$data
dat$site <- factor(seq_len(nrow(dat)))

fit <- tobs(~ abund_cov1 + abund_cov2, data = dat, y = sim$y,
            family = abun(mixture = "negbin"),
            detection = ~ det_cov1, method = "laplace", verbose = FALSE)

cat("=== coef ===\n"); print(lapply(coef(fit), round, 3))
cat("\ntruth beta_lambda:", round(sim$truth$beta_lambda, 3), "\n")
cat("truth beta_p:     ", round(sim$truth$beta_p, 3), "\n")
cat("truth size r:     ", sim$truth$size, "\n\n")

cat("mixture:", fit$mixture, "\n")
cat("dispersion:\n"); print(fit$nmix_dispersion)
cat("\nlog_r in vcov dimnames:", "log_r" %in% rownames(fit$vcov), "\n")
cat("vcov dim:", paste(dim(fit$vcov), collapse = " x "), "\n")
cat("SE(log_r):", round(sqrt(fit$vcov["log_r","log_r"]), 4), "\n\n")

# Poisson sanity: same call, no log_r coordinate.
fitP <- tobs(~ abund_cov1 + abund_cov2, data = dat, y = sim$y,
             family = abun(mixture = "poisson"),
             detection = ~ det_cov1, method = "laplace", verbose = FALSE)
cat("Poisson mixture:", fitP$mixture, " log_r present:",
    "log_r" %in% rownames(fitP$vcov), " vcov dim:",
    paste(dim(fitP$vcov), collapse = "x"), "\n")

# simulate() under NB returns counts (smoke).
ys <- simulate(fit, nsim = 1)
cat("simulate() NB: dim", paste(dim(ys), collapse = "x"),
    " mean", round(mean(ys, na.rm = TRUE), 2), "\n")
