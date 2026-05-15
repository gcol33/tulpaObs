# Verify issue #4 fix: cover('lognormal') joint nested_laplace should now
# recover alpha rather than collapse to ~0, and report sigma_pos close to
# sigma_pos_true rather than sqrt(sigma_pos^2 + alpha^2 * field_sigma^2).
#
# Reduced grid for speed (3 cells, 4 seeds) — full sweep in the issue body.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

chain_adj <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) for (j in setdiff(c(s - 1L, s + 1L),
                                            c(0L, n_s + 1L))) adj[s, j] <- 1L
  adj
}

cfg <- list(N = 300L, n_s = 25L, sigma = 0.6, rho = 0.7,
            betaO = c(-0.3, 0.7), betaP = c(-1.5, 0.3),
            n_seeds = 4L)
adj <- chain_adj(cfg$n_s)

simulate <- function(alpha_true, sigma_pos_true, seed) {
  set.seed(seed)
  region <- sample.int(cfg$n_s, cfg$N, replace = TRUE)
  pf <- rnorm(cfg$n_s, 0, 1); tf <- rnorm(cfg$n_s, 0, 1)
  w  <- cfg$sigma * (sqrt(cfg$rho) * pf + sqrt(1 - cfg$rho) * tf)
  x  <- rnorm(cfg$N)
  occur <- rbinom(cfg$N, 1L, plogis(cfg$betaO[1] + cfg$betaO[2] * x + w[region]))
  eta_pos <- cfg$betaP[1] + cfg$betaP[2] * x + alpha_true * w[region]
  log_y <- rnorm(cfg$N, eta_pos, sigma_pos_true)
  y <- numeric(cfg$N); is_pos <- occur == 1L
  y[is_pos] <- exp(log_y[is_pos])
  y <- pmin(y, 1 - 1e-6)
  list(data = data.frame(x = x, region = factor(region)), y = y)
}

extract_alpha <- function(fit) {
  k <- which(colnames(fit$joint$theta_grid) == "alpha")
  m <- sum(fit$joint$weights * fit$joint$theta_grid[, k])
  s <- sqrt(max(0, sum(fit$joint$weights *
                       fit$joint$theta_grid[, k]^2) - m^2))
  c(m, s)
}

cat(sprintf("%-8s %-9s | %-22s | %-12s\n",
            "alpha_t", "sigma_t", "alpha_hat (sd seeds)", "sigma_hat"))
cat(strrep("-", 60), "\n", sep = "")

grid <- expand.grid(alpha_true = c(0.0, 1.0, 1.5),
                    sigma_pos_true = c(0.3, 1.0))
for (g in seq_len(nrow(grid))) {
  a_t <- grid$alpha_true[g]; s_t <- grid$sigma_pos_true[g]
  alphas <- sigmas <- numeric(cfg$n_seeds)
  for (k in seq_len(cfg$n_seeds)) {
    sim <- simulate(a_t, s_t, 50000L + 100L * g + k)
    fit <- tryCatch(
      tobs(formula = ~ x, data = sim$data, family = cover("lognormal"),
           y = sim$y,
           spatial = tulpa::spatial_bym2(adj, level = "group",
                                         group_var = "region"),
           engine = "nested_laplace",
           control = list(sigma_grid = c(0.3, 0.6, 0.9),
                          rho_grid   = c(0.5, 0.7, 0.9),
                          alpha_grid = c(0.0, 0.5, 1.0, 1.5))),
      error = function(e) { message("seed ", k, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) { alphas[k] <- NA; sigmas[k] <- NA; next }
    alphas[k] <- extract_alpha(fit)[1]
    sigmas[k] <- fit$sigma_pos
  }
  cat(sprintf("%-8.2f %-9.2f | %.3f (sd=%.3f)        | %.3f\n",
              a_t, s_t, mean(alphas, na.rm = TRUE),
              sd(alphas, na.rm = TRUE),
              mean(sigmas, na.rm = TRUE)))
}
