# Reproduce test-sla-cover-hurdle.R:221 (deterministic) and print the actual
# magnitudes, to tell a genuine numerical regression from a too-tight heuristic.
suppressMessages(devtools::load_all("."))

simulate_beta_cover_local <- function(N = 400, beta_occ = c(-0.3, 0.7),
                                      beta_pos = c(0.4, -1.0), phi = 25,
                                      seed = 1) {
  set.seed(seed)
  x <- runif(N, -2, 2)
  occur <- rbinom(N, 1, plogis(beta_occ[1] + beta_occ[2] * x))
  mu <- plogis(beta_pos[1] + beta_pos[2] * x)
  y <- numeric(N); is_pos <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos), mu[is_pos] * phi, (1 - mu[is_pos]) * phi)
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos, phi = phi))
}

sim <- simulate_beta_cover_local(N = 400, seed = 25)
fit <- tobs(~ x, data = sim$data, family = cover("beta"), y = sim$y)

enc <- fit$encoding
beta_occ_sc <- tulpaObs:::.scale_beta_vec(fit$beta_occ, enc$scale_occ)
beta_pos_sc <- tulpaObs:::.scale_beta_vec(fit$beta_pos, enc$scale_pos)
ll_occ <- tulpaObs:::.loglik_cover_occ(beta_occ_sc, enc)
ll_pos <- tulpaObs:::.loglik_cover_pos_beta(beta_pos_sc, fit$phi_pos, enc)
ll_total <- ll_occ + ll_pos
lm_total <- fit$log_marginal["occ"] + fit$log_marginal["pos"]
ratio <- abs(ll_total - lm_total) / max(abs(ll_total), 1)

cat("=== components ===\n")
cat(sprintf("ll_occ = %.3f   ll_pos = %.3f   ll_total = %.3f\n", ll_occ, ll_pos, ll_total))
cat(sprintf("log_marginal: occ = %.3f  pos = %.3f  total = %.3f\n",
            fit$log_marginal["occ"], fit$log_marginal["pos"], lm_total))
cat(sprintf("relative diff = %.4f   (test threshold 0.5)\n", ratio))

cat("\n=== recovery (the substantive check) ===\n")
cat("beta_occ est:", round(fit$beta_occ, 3), " truth:", sim$truth$beta_occ, "\n")
cat("beta_pos est:", round(fit$beta_pos, 3), " truth:", sim$truth$beta_pos, "\n")
cat("phi_pos est:", round(fit$phi_pos, 2), " truth:", sim$truth$phi, "\n")
