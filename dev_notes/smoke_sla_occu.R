# smoke_sla_occu.R
#
# End-to-end smoke test: fit a single-season occu model with
# approx = "simplified_laplace" and verify gamma is computed and stored.

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
devtools::load_all(quiet = TRUE)

set.seed(2026)
n_sites <- 60
n_visits <- 3
elev <- rnorm(n_sites)
effort <- rnorm(n_sites * n_visits)
beta_psi <- c(-0.5, 1.0)
beta_p   <- c(-1.2, 0.4)

eta_psi <- cbind(1, elev) %*% beta_psi
psi <- plogis(eta_psi)
z <- rbinom(n_sites, 1, psi)

eta_p <- cbind(1, effort) %*% beta_p
p <- plogis(eta_p)
p_mat <- matrix(p, n_sites, n_visits, byrow = TRUE)

y <- matrix(0L, n_sites, n_visits)
for (i in seq_len(n_sites)) {
  for (j in seq_len(n_visits)) {
    y[i, j] <- if (z[i] == 1) rbinom(1, 1, p_mat[i, j]) else 0L
  }
}

cat("prevalence (z):", round(mean(z), 2), "\n")
cat("any-detected rate:", round(mean(rowSums(y) > 0), 2), "\n\n")

# --- Gaussian Laplace (default) ----------------------------------------
fit_g <- tobs(
  formula   = ~ elev,
  data      = data.frame(elev = elev),
  family    = occu(),
  detection = ~ 1,
  y         = y,
  approx    = "gaussian_laplace",
  control   = list(verbose = FALSE)
)
cat("--- Gaussian Laplace ---\n")
cat("  means:", round(fit_g$means, 3), "\n")
cat("  sds  :", round(fit_g$sds, 3),   "\n")
cat("  skew :", if (is.null(fit_g$skew)) "NULL" else paste(round(fit_g$skew, 3), collapse = " "), "\n")
cat("  sla_status:", fit_g$sla_status, "\n\n")

# --- Simplified Laplace ------------------------------------------------
fit_s <- tobs(
  formula   = ~ elev,
  data      = data.frame(elev = elev),
  family    = occu(),
  detection = ~ 1,
  y         = y,
  approx    = "simplified_laplace",
  control   = list(verbose = FALSE)
)
cat("--- Simplified Laplace ---\n")
cat("  means:", round(fit_s$means, 3), "\n")
cat("  sds  :", round(fit_s$sds, 3),   "\n")
cat("  skew :", round(fit_s$skew, 3), "\n")
cat("  sla_status:", fit_s$sla_status, "\n\n")

# Marginal moments from the (SLA-resampled) draws — confirm they match
# the requested skewness
emp_skew <- apply(fit_s$draws, 2, function(d) {
  m <- mean(d); s <- sd(d)
  mean(((d - m) / s)^3)
})
cat("--- Draws cross-check (SLA fit) ---\n")
cat("  requested skew (fit$skew):", round(fit_s$skew, 3), "\n")
cat("  empirical skew of draws  :", round(emp_skew, 3), "\n")

# Compare CIs
ci_g <- t(apply(fit_g$draws, 2, quantile, c(0.025, 0.975)))
ci_s <- t(apply(fit_s$draws, 2, quantile, c(0.025, 0.975)))
cat("\n--- 95% CIs ---\n")
cat("  Gaussian Laplace:\n"); print(round(ci_g, 3))
cat("  Simplified Laplace:\n"); print(round(ci_s, 3))
