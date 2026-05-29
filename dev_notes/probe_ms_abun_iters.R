# probe_ms_abun_iters.R - inspect EM iter counts and per-iter timing.
suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

probe <- function(S, N, J, max_iter = 100L, label = "") {
  sim <- simulate_ms_abun(n_species = S, N = N, J = J,
                          n_abund_covs = 1L, n_det_covs = 1L,
                          mu_lambda = c(log(3), 0.4),
                          mu_p      = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          mixture = "poisson", seed = 100L + S)
  t0 <- proc.time()[["elapsed"]]
  fit <- tobs(formula = ~ abund_cov1, detection = ~ det_cov1,
              data = sim$data,
              family = ms_abun(mixture = "poisson"),
              y = sim$y, species = sim$species,
              control = list(verbose = FALSE, max.iter = max_iter))
  dt <- proc.time()[["elapsed"]] - t0
  ni <- fit$convergence$n_iter %||% NA_integer_
  conv <- isTRUE(fit$convergence$converged)
  cat(sprintf("  %-12s  S=%2d N=%3d J=%d  max_iter=%3d  -> %.2fs  n_iter=%s converged=%s\n",
              label, S, N, J, max_iter, dt, as.character(ni), conv))
}

cat("== ms_abun() EM iter probe ==\n")
probe(S = 5L,  N = 20L, J = 3L, max_iter = 100L, label = "tiny default")
probe(S = 5L,  N = 20L, J = 3L, max_iter = 20L,  label = "tiny iter=20")
probe(S = 12L, N = 60L, J = 4L, max_iter = 100L, label = "coverage def")
probe(S = 12L, N = 60L, J = 4L, max_iter = 20L,  label = "coverage 20")
probe(S = 14L, N = 90L, J = 4L, max_iter = 100L, label = "recovery def")
