suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})
Sys.setenv(NOT_CRAN = "true")
set.seed(3)
sim <- simulate_ms_abun(n_species = 8, N = 40, J = 3, seed = 3)
fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
            family = ms_abun(), detection = ~ det_cov1,
            species = sim$species, method = "laplace",
            control = list(verbose = FALSE))
cf <- coef(fit)
cat("class(cf):", class(cf), "\n")
cat("is.list(cf):", is.list(cf), "\n")
cat("names(cf):", names(cf), "\n")
cat("str(cf):\n"); str(cf)
