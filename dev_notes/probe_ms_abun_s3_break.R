suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
set.seed(3)
sim <- simulate_ms_abun(n_species = 8, N = 40, J = 3, seed = 3)
fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
            family = ms_abun(), detection = ~ det_cov1,
            species = sim$species, method = "laplace",
            control = list(verbose = FALSE))

cat("== fit summary ==\n")
cat("mixture:           ", fit$mixture, "\n")
cat("means:\n"); print(fit$means)
cat("sds:\n"); print(fit$sds)
cat("converged:         ", fit$convergence$converged, "\n")
cat("nrow(vcov):        ", nrow(fit$vcov), "\n")
cat("length(means):     ", length(fit$means), "\n")
cat("ms_dispersion:     ", capture.output(str(fit$ms_dispersion)), "\n")
cat("\n== coef(fit) ==\n")
cf <- coef(fit)
cat("class:", class(cf), "\n")
print(cf)
