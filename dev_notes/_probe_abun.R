# Probe: non-spatial Poisson N-mixture wiring + recovery.
suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
})

set.seed(1)
beta_lambda <- c(log(4), 0.6, -0.4)   # intercept + 2 abundance slopes
beta_p      <- c(0.2, 0.5)            # intercept + 1 detection slope
sim <- simulate_abun(N = 300, J = 5,
                     n_abund_covs = 2, n_det_covs = 1,
                     beta_lambda = beta_lambda, beta_p = beta_p, seed = 7)

cat("mean count:", mean(sim$y), " max count:", max(sim$y), "\n")

fit <- tobs(
  formula   = ~ abund_cov1 + abund_cov2,
  data      = sim$data,
  family    = abun(),
  detection = ~ det_cov1,
  y         = sim$y,
  method    = "laplace"
)

print(fit)
cat("\n--- coefficients (truth in parens) ---\n")
truth <- c(beta_lambda, beta_p)
names(truth) <- names(fit$means)
est <- fit$means
se  <- fit$sds
for (nm in names(est)) {
  cat(sprintf("  %-22s est=% .3f  se=%.3f  truth=% .3f  z=% .2f\n",
              nm, est[nm], se[nm], truth[nm], (est[nm] - truth[nm]) / se[nm]))
}

cat("\n--- summary() ---\n")
print(summary(fit))

cat("\n--- fitted(): lambda/p/N heads ---\n")
fv <- fitted(fit)
cat("lambda:", round(head(fv$lambda), 2), "\n")
cat("p     :", round(head(fv$p), 2), "\n")
cat("E[N]  :", round(head(fv$N), 2), "\n")

cat("\n--- predict(terms='abund_cov1') head ---\n")
print(head(predict(fit, terms = "abund_cov1")))

cat("\n--- vcov dim:", paste(dim(vcov(fit)), collapse = "x"), "\n")
cat("--- nobs:", nobs(fit), "\n")
cat("--- logLik:", as.numeric(logLik(fit)), "\n")
cat("DONE\n")
