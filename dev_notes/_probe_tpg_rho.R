# Probe: (1) isolate the AR1 sigma/rho grid update on a long known eta series;
#        (2) full t_occu fit at higher T to see rho recovery scaling.
suppressMessages({
  library(tulpa); devtools::load_all(".", quiet = TRUE)
})
sink("dev_notes/_probe_tpg_rho.out")

## (1) Isolated AR1 hyperparameter update -----------------------------------
# Given a true AR1 eta series, does the grid update recover rho / sigma?
ar1_update_check <- function(T_s, rho_true, sigma_true, nsweep = 4000, seed = 1) {
  set.seed(seed)
  eta <- numeric(T_s)
  eta[1] <- rnorm(1, 0, sigma_true / sqrt(1 - rho_true^2))
  for (t in 2:T_s) eta[t] <- rho_true * eta[t - 1] + rnorm(1, 0, sigma_true)
  rho_grid <- seq(-0.95, 0.95, by = 0.05)
  ig_a <- 0.1; ig_b <- 0.1
  sigma2 <- 1; rho <- 0.0
  keep <- matrix(NA_real_, nsweep, 2)
  for (s in seq_len(nsweep)) {
    ss <- (1 - rho^2) * eta[1]^2 + sum((eta[-1] - rho * eta[-T_s])^2)
    sigma2 <- 1 / rgamma(1, ig_a + T_s / 2, ig_b + 0.5 * ss)
    lp <- vapply(rho_grid, function(rg) {
      s2 <- (1 - rg^2) * eta[1]^2 + sum((eta[-1] - rg * eta[-T_s])^2)
      0.5 * log(1 - rg^2) - 0.5 * s2 / sigma2
    }, numeric(1))
    w <- exp(lp - max(lp)); rho <- sample(rho_grid, 1, prob = w)
    keep[s, ] <- c(sqrt(sigma2), rho)
  }
  b <- keep[(nsweep / 2 + 1):nsweep, , drop = FALSE]
  c(sigma_hat = mean(b[, 1]), rho_hat = mean(b[, 2]))
}
cat("== isolated AR1 update (truth rho=0.6 sigma=0.7) ==\n")
for (T_s in c(8, 20, 50, 200)) {
  r <- ar1_update_check(T_s, 0.6, 0.7)
  cat(sprintf("  T=%3d  sigma_hat=%.3f  rho_hat=%.3f\n", T_s, r[1], r[2]))
}

## (2) Full t_occu fit, scaling T -------------------------------------------
cat("\n== full t_occu fit (truth beta0=0.2 p=0.4 rho=0.6 sigma=0.7) ==\n")
for (T_s in c(8, 15, 25)) {
  sim <- simulate_t_occu(N = 250, T_seasons = T_s, J = 4, beta_occ = c(0.2),
                         p = 0.4, rho = 0.6, sigma = 0.7, seed = 3)
  fit <- tobs(~ 1, family = t_occu(), detection = ~ 1, y = sim$y,
              data = sim$data, method = "pg_gibbs",
              control = list(n.iter = 2500L, n.warmup = 1200L,
                             n.chains = 2L, seed = 7L))
  m <- fit$means
  eta_cor <- suppressWarnings(cor(fit$temporal_field, sim$truth$eta))
  cat(sprintf(
    "  T=%2d  psi0=%+.3f(0.2)  p=%.3f(0.4)  sigma=%.3f(0.7)  rho=%.3f(0.6)  eta_cor=%.3f  maxRhat=%.3f\n",
    T_s, m[["psi_(Intercept)"]], plogis(m[["p_(Intercept)"]]),
    exp(m[["log_sigma_ar1"]]), m[["rho_ar1"]], eta_cor,
    max(fit$rhat, na.rm = TRUE)))
}
cat("\nDONE\n")
sink()
