# Env-immune (base R only) check of the single-season occupancy SE calibration.
# The marginal likelihood is closed-form (z integrated out):
#   site i, any detection: L_i = psi_i * prod_j Bern(y_ij; p_ij)
#   site i, no detection:  L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)
# The exact observed information is the Hessian of this marginal log-lik at the
# MLE. If the analytic marginal-MLE SE matches the Monte-Carlo sd of the MLE
# across seeds (ratio ~ 1), then the right fix for tobs is to take SEs from the
# closed-form marginal Hessian -- and the deterministic-EM Louis/M-step SE is
# the approximation that under-disperses the detection slope.

nll <- function(par, y, Xo, Xd, nvalid, ndet, anydet) {
  po <- ncol(Xo); pd <- ncol(Xd)
  psi <- plogis(pmin(pmax(Xo %*% par[1:po], -30), 30))
  p   <- plogis(pmin(pmax(Xd %*% par[po + 1:pd], -30), 30))
  ll <- numeric(length(psi))
  det <- anydet
  ll[det]  <- log(psi[det]) + ndet[det] * log(p[det]) +
              (nvalid[det] - ndet[det]) * log(1 - p[det])
  ll[!det] <- log(psi[!det] * (1 - p[!det])^nvalid[!det] + (1 - psi[!det]))
  -sum(ll)
}

fit_marginal <- function(N, J, bo = c(0.5, 1.2), bd = c(0, 0.8), seed) {
  set.seed(seed)
  xo <- rnorm(N); xd <- rnorm(N)
  Xo <- cbind(1, xo); Xd <- cbind(1, xd)
  psi <- plogis(Xo %*% bo); p <- plogis(Xd %*% bd)
  z <- rbinom(N, 1, psi)
  Y <- matrix(0L, N, J)
  for (i in seq_len(N)) Y[i, ] <- rbinom(J, 1, z[i] * p[i])
  nvalid <- rep(J, N); ndet <- rowSums(Y); anydet <- ndet > 0
  o <- optim(c(bo, bd), nll, y = Y, Xo = Xo, Xd = Xd, nvalid = nvalid,
             ndet = ndet, anydet = anydet, method = "BFGS", hessian = TRUE)
  se <- sqrt(diag(solve(o$hessian)))
  list(est = o$par, se = se)  # par = (psi_int, psi_slope, p_int, p_slope)
}

for (cfg in list(c(600, 6), c(400, 16), c(200, 4))) {
  N <- cfg[1]; J <- cfg[2]
  E <- S <- matrix(NA_real_, 30, 4)
  for (k in 1:30) { r <- fit_marginal(N, J, seed = 7000 + k); E[k, ] <- r$est; S[k, ] <- r$se }
  lab <- c("psi_int", "psi_slope", "p_int", "p_slope")
  cat(sprintf("\n== N=%d J=%d (analytic marginal MLE) ==\n", N, J))
  for (j in 1:4)
    cat(sprintf("  %-10s est=%+.3f  medSE=%.3f  mcSD=%.3f  ratio=%.2f\n",
                lab[j], mean(E[, j]), median(S[, j]), sd(E[, j]),
                median(S[, j]) / sd(E[, j])))
}
