suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs"))
set.seed(101)
ng <- 30; per <- 20; N <- ng * per; J <- 8
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(N); z <- rnorm(N)
sig <- c(0.9, 0.7, 0.5)
R <- matrix(c(1, .5, .3,  .5, 1, .2,  .3, .2, 1), 3, 3)
Sig <- diag(sig) %*% R %*% diag(sig)
B <- MASS::mvrnorm(ng, mu = c(0, 0, 0), Sigma = Sig)
gi <- as.integer(g)
eta <- 0.3 + B[gi, 1] + B[gi, 2] * x + B[gi, 3] * z
psi <- plogis(eta); zocc <- rbinom(N, 1, psi); pdet <- 0.6
y <- matrix(0L, N, J); for (i in 1:N) y[i, ] <- rbinom(J, 1, zocc[i] * pdet)
d <- data.frame(g = g, x = x, z = z)

fit <- tobs(~ (1 + x + z | g), data = d, y = y, detection = ~ 1,
            family = occu(), engine = "nuts",
            control = list(iter = 500, warmup = 250, seed = 1, verbose = FALSE))

m <- fit$means
cat("total params:", length(m), "\n")
# Layout: [psi_int, p_int, log_sig(3), chol(3), z-effects(ng*3)]
sig_hat <- exp(m[3:5])
cat("sigma truth:", round(sig, 3), "\n")
cat("sigma hat  :", round(sig_hat, 3), "\n")

# Reconstruct correlated group effects: L = tanh-cholesky of chol_raw, then
# effect_gc = sigma_c * (L %*% z_g)_c.
chol_raw <- m[6:8]
L <- matrix(0, 3, 3)
idx <- 1
for (rr in 1:3) {
  s2 <- 0
  for (cc in seq_len(rr - 1)) { L[rr, cc] <- tanh(chol_raw[idx]); s2 <- s2 + L[rr, cc]^2; idx <- idx + 1 }
  L[rr, rr] <- sqrt(max(1 - s2, 1e-10))
}
zoff <- 8
Bhat <- matrix(0, ng, 3)
for (gg in seq_len(ng)) {
  zg <- m[zoff + (gg - 1) * 3 + 1:3]
  Bhat[gg, ] <- sig_hat * as.numeric(L %*% zg)
}
cat("cor(Bhat, Btrue) per coef:\n")
for (c in 1:3) cat("  coef", c, ":", round(cor(Bhat[, c], B[, c]), 3), "\n")

cat("\n-- tulpa::ranef --\n")
print(tryCatch(tulpa::ranef(fit), error = function(e) conditionMessage(e)))
cat("\n=== done ===\n")
