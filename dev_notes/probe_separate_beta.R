suppressPackageStartupMessages({
  library(devtools)
  load_all("../tulpa",    quiet = TRUE, export_all = FALSE)
  load_all(".",            quiet = TRUE, export_all = FALSE)
})

simulate_beta_separate <- function(N = 400, beta_occ = c(-0.4, 0.8),
                                   beta_pos = c(0.4, -0.5), phi = 30, seed = 5000) {
  set.seed(seed)
  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  occur   <- rbinom(N, 1L, plogis(eta_occ))
  mu_pos  <- plogis(eta_pos)
  y       <- numeric(N); is_pos <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi,
                     (1 - mu_pos[is_pos]) * phi)
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x), y = y, n_pos = sum(is_pos))
}

cat("--- separate-hurdle beta phi_pos recovery (10 seeds, N=400) ---\n")
truth_phi <- 30
phi_hats <- numeric(10L)
for (r in seq_len(10L)) {
  sim <- simulate_beta_separate(N = 400, phi = truth_phi, seed = 5000L + r)
  fit <- tobs(formula = ~ x, data = sim$data,
              family = cover(positive = "beta"), y = sim$y)
  phi_hats[r] <- fit$phi_pos
  cat(sprintf("  seed %d  n_pos %d  phi_hat = %.3f  (rel err %.3f)\n",
              5000L + r, sim$n_pos, fit$phi_pos,
              (fit$phi_pos - truth_phi) / truth_phi))
}
cat(sprintf("  mean = %.3f, rel mean err = %.3f, max abs rel err = %.3f\n",
            mean(phi_hats),
            (mean(phi_hats) - truth_phi) / truth_phi,
            max(abs(phi_hats - truth_phi) / truth_phi)))
