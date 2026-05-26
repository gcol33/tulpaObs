# Probe: AGHQ variance-component debias for occupancy RE.
# Compares raw EM (Laplace, nAGQ=1) vs AGHQ vs NUTS on sigma recovery, and
# checks the bias-vs-cluster-size direction (smaller per-group n -> more
# attenuation, more for AGHQ to recover).
devtools::load_all(".", quiet = TRUE)

sim_re_int <- function(seed, ng, per, J = 6L, b0 = 0.3, b1 = -0.6,
                       sigma = 0.9, p = 0.45) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  b <- rnorm(ng, 0, sigma)
  z <- rbinom(N, 1, plogis(b0 + b1 * x + b[g]))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), b = b)
}

sig_of <- function(fit) fit$means[[grep("^sigma_", names(fit$means), value = TRUE)[1]]]

cat("== sigma recovery (truth = 0.9), EM vs AGHQ ==\n")
for (per in c(8L, 15L, 30L)) {
  s <- sim_re_int(seed = 11, ng = 30L, per = per)
  fit_em <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                 family = occu(), method = "laplace",
                 control = list(re.aghq = FALSE, verbose = FALSE))
  fit_aq <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                 family = occu(), method = "laplace",
                 control = list(re.aghq = TRUE, n.quad = 9L, verbose = FALSE))
  cat(sprintf("  per=%2d  EM sigma=%.3f   AGHQ sigma=%.3f   (applied=%s)\n",
              per, sig_of(fit_em), sig_of(fit_aq),
              isTRUE(fit_aq$aghq$applied)))
}

cat("\n== correlated (1 + x | g), true sigma=(.77,.63), rho=+0.61 ==\n")
sim_corr <- function(seed = 402, ng = 40L, per = 12L, J = 6L) {
  set.seed(seed)
  N <- ng * per; g <- rep(seq_len(ng), each = per); x <- rnorm(N)
  Sig <- matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)
  U <- matrix(rnorm(ng * 2), ng, 2) %*% chol(Sig)
  z <- rbinom(N, 1, plogis(0.2 + U[g, 1] + (-0.4 + U[g, 2]) * x))
  y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * 0.5)
  list(y = y, d = data.frame(g = factor(g), x = x))
}
sc <- sim_corr(per = 10L)
for (aq in c(FALSE, TRUE)) {
  f <- tobs(~ x + (1 + x | g), data = sc$d, y = sc$y, detection = ~ 1,
            family = occu(), method = "laplace",
            control = list(re.aghq = aq, n.quad = 7L, verbose = FALSE))
  sg <- f$means[grep("^sigma_", names(f$means))]
  rho <- f$means[[grep("^cor_", names(f$means), value = TRUE)]]
  cat(sprintf("  aghq=%-5s  sigma=(%.3f, %.3f)  rho=%.3f\n",
              aq, sg[1], sg[2], rho))
}
