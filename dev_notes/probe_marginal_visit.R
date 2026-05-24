# Visit-level detection: confirm the closed-form marginal MLE is calibrated for
# the issue8 config (per-visit detection covariate). Marginal likelihood:
#   no detection: L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)
#   any detection: L_i = psi_i * prod_j Bern(y_ij; p_ij)   (z_i = 1 forced)
nll_visit <- function(par, Y, V, psiX) {
  # par = (psi_int, p_int, p_slope); psiX intercept-only here.
  N <- nrow(Y); J <- ncol(Y)
  psi <- plogis(par[1])
  p <- plogis(pmin(pmax(par[2] + par[3] * V, -30), 30))  # N x J
  anydet <- rowSums(Y) > 0
  ll <- numeric(N)
  for (i in seq_len(N)) {
    if (anydet[i]) {
      ll[i] <- log(psi) + sum(Y[i, ] * log(p[i, ]) + (1 - Y[i, ]) * log(1 - p[i, ]))
    } else {
      ll[i] <- log(psi * prod(1 - p[i, ]) + (1 - psi))
    }
  }
  -sum(ll)
}
N <- 300L; J <- 6L; p0 <- 0.2; p1 <- 1.2; psi0 <- plogis(0.4)
E <- S <- matrix(NA_real_, 30, 3)
for (k in 1:30) {
  set.seed(5000L + k)
  z <- rbinom(N, 1, psi0); V <- matrix(rnorm(N * J), N, J)
  Y <- matrix(0L, N, J)
  for (i in seq_len(N)) Y[i, ] <- rbinom(J, 1, z[i] * plogis(p0 + p1 * V[i, ]))
  o <- optim(c(0, 0, 0), nll_visit, Y = Y, V = V, method = "BFGS", hessian = TRUE)
  E[k, ] <- o$par; S[k, ] <- sqrt(diag(solve(o$hessian)))
}
lab <- c("psi_int", "p_int", "p_visit_slope"); truth <- c(0.4, p0, p1)
cat(sprintf("== visit-level N=%d J=%d (analytic marginal MLE) ==\n", N, J))
for (j in 1:3)
  cat(sprintf("  %-14s est=%+.3f (truth %+.2f)  medSE=%.3f  mcSD=%.3f  ratio=%.2f\n",
              lab[j], mean(E[, j]), truth[j], median(S[, j]), sd(E[, j]),
              median(S[, j]) / sd(E[, j])))
