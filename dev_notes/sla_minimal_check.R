# sla_minimal_check.R
#
# Validates the upstream-tulpa SLA spec (dev_notes/upstream_tulpa_sla_spec.md §4)
# by running the SLA assembly end-to-end using only what tulpa_laplace()
# already exposes. If this works, no engine-level tulpa changes are needed
# for fixed-effects SLA marginals.

library(tulpa)

set.seed(1)
n <- 100
X <- cbind(intercept = 1, x = rnorm(n))
beta_true <- c(-1, 0.5)
eta <- X %*% beta_true
p <- plogis(eta)
y <- rbinom(n, 1, p)

cat("--- Inputs ---\n")
cat("  n =", n, ", p_fixed =", ncol(X), "\n")
cat("  beta_true =", beta_true, "\n")
cat("  prevalence =", mean(y), "\n\n")

fit <- tulpa::tulpa_laplace(y = y, n_trials = rep(1L, n), X = X,
                            family = "binomial", return_hessian = TRUE)

cat("--- tulpa_laplace() output ---\n")
str(fit, max.level = 1, give.attr = FALSE)
cat("\n")

# Required exposures for SLA:
stopifnot("mode missing"   = !is.null(fit$mode))
stopifnot("H_beta missing" = !is.null(fit$H_beta))

beta_hat <- fit$mode[seq_len(ncol(X))]
Sigma <- solve(fit$H_beta)
sigma_j <- sqrt(diag(Sigma))

cat("--- SLA assembly ---\n")
cat("  beta_hat   =", round(beta_hat, 4), "\n")
cat("  sigma_j    =", round(sigma_j, 4), "\n\n")

eta_hat <- as.numeric(X %*% beta_hat)
p_hat   <- plogis(eta_hat)
# n_trials = 1 here, so l3 per site:
l3 <- -p_hat * (1 - p_hat) * (1 - 2 * p_hat)

# v_{i,j} = (X Sigma)_{i,j}, then gamma_j = sigma_j^{-3} * sum_i l3_i * v_ij^3
XSig <- X %*% Sigma
gamma_j <- vapply(seq_len(ncol(X)), function(j) {
  sum(l3 * XSig[, j]^3) / sigma_j[j]^3
}, numeric(1))

names(gamma_j) <- colnames(X)
cat("  gamma_j (SLA) =", round(gamma_j, 4), "\n\n")

# Cross-check: independent posterior moments via brute-force MCMC ----------
# Use a simple Metropolis on the joint posterior for n_iter draws, compute
# marginal skewness empirically.

log_post <- function(beta) {
  eta <- as.numeric(X %*% beta)
  sum(y * eta - log1p(exp(eta))) + sum(dnorm(beta, 0, 5, log = TRUE))
}

n_iter <- 50000
burnin <- 5000
draws <- matrix(NA_real_, n_iter, length(beta_hat))
cur <- beta_hat
cur_lp <- log_post(cur)
prop_sd <- 0.25 * sigma_j  # ad-hoc proposal scale
acc <- 0
for (it in seq_len(n_iter)) {
  prop <- cur + rnorm(length(cur), 0, prop_sd)
  prop_lp <- log_post(prop)
  if (log(runif(1)) < prop_lp - cur_lp) {
    cur <- prop; cur_lp <- prop_lp; acc <- acc + 1
  }
  draws[it, ] <- cur
}
draws <- draws[(burnin + 1):n_iter, ]

emp_skew <- apply(draws, 2, function(d) {
  m <- mean(d); s <- sd(d)
  mean(((d - m) / s)^3)
})
names(emp_skew) <- colnames(X)

cat("--- Cross-check vs Metropolis posterior ---\n")
cat("  acceptance rate =", round(acc / n_iter, 3), "\n")
cat("  gamma_j (MCMC)  =", round(emp_skew, 4), "\n")
cat("  gamma_j (SLA)   =", round(gamma_j, 4), "\n")
cat("  abs diff        =", round(abs(gamma_j - emp_skew), 4), "\n")
cat("  sign match      =", all(sign(gamma_j) == sign(emp_skew)), "\n")

# Posterior sigma (cross-check sigma_j against MCMC SD)
mcmc_sd <- apply(draws, 2, sd)
cat("\n  sigma_j (Laplace) =", round(sigma_j, 4), "\n")
cat("  sigma_j (MCMC)    =", round(mcmc_sd, 4), "\n")
