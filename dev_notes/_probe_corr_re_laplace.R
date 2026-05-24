# Multi-seed recovery probe for correlated random slopes (1 + x | g) on the
# deterministic Laplace path (gcol33/tulpaObs#11, consuming tulpa's cov_blocks).
# Reports sigma / correlation / BLUP recovery vs simulated truth so the test
# tolerances are set from observed behavior, not guessed.
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages({
  devtools::install("../tulpa", quick = TRUE, upgrade = FALSE, quiet = TRUE)
  devtools::load_all(quiet = TRUE)
})

sim_occu_re_corr <- function(seed = 404, ng = 40L, per = 25L, J = 6L,
                             b0 = 0.2, b1 = -0.4, p = 0.5,
                             Sigma = matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  U <- matrix(rnorm(ng * 2L), ng, 2L) %*% chol(Sigma)
  eta <- b0 + U[g, 1] + (b1 + U[g, 2]) * x
  z <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), U = U, Sigma = Sigma)
}

Sig <- matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)
rho_true <- Sig[1, 2] / sqrt(Sig[1, 1] * Sig[2, 2])
cat(sprintf("truth: sigma0=%.3f sigma1=%.3f rho=%.3f b0=0.2 b1=-0.4\n",
            sqrt(Sig[1, 1]), sqrt(Sig[2, 2]), rho_true))

res <- lapply(seq_len(6), function(sd) {
  s <- sim_occu_re_corr(seed = 400 + sd)
  fit <- tobs(~ x + (1 + x | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace", control = list(verbose = FALSE))
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  cor_nm <- grep("^cor_",   names(fit$means), value = TRUE)
  cf <- coef(fit)$psi
  re <- ranef(fit)
  b0h <- re$estimate[re$term == "(Intercept)"]
  b1h <- re$estimate[re$term == "x"]
  data.frame(
    seed = 400 + sd,
    sig0 = unname(fit$means[sig_nm[1]]),
    sig1 = unname(fit$means[sig_nm[2]]),
    rho  = unname(fit$means[cor_nm[1]]),
    b0   = cf[["(Intercept)"]],
    b1   = cf[["x"]],
    cor0 = cor(b0h, s$U[, 1]),
    cor1 = cor(b1h, s$U[, 2]))
})
res <- do.call(rbind, res)
print(round(res, 3))
cat("\nmeans:\n"); print(round(colMeans(res[-1]), 3))
