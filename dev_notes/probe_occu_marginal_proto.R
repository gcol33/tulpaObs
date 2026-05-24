# Pure-R prototype of the occupancy marginal-likelihood state posterior, to
# validate the math BEFORE building the engine `occupancy` family.
#
# Marginalising the latent occupancy state z gives, per site, a single binary
# response D_i = 1{>=1 detection} with mean mu_i = q_i * sigma(eta_i), where
# q_i = 1 - (1-p)^{J_i} is the per-site detection probability. The state field
# is then an ICAR-GMRF + this scaled-logistic Bernoulli likelihood -- NO EM,
# NO M-inflation. Held-out sites have J_i = 0 -> q_i = 0 -> no information.
#
# Direct GMRF Laplace over a tau grid, marginalise psi over the grid, check
# held-out coverage. If ~0.95, the math is right and the engine build follows.
#
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/probe_occu_marginal_proto.R

sig <- function(x) 1 / (1 + exp(-x))

make_grid <- function(gx = 10, gy = 10, J = 8, p_det = 0.5, seed = 1) {
  set.seed(seed)
  n <- gx * gy
  coord <- expand.grid(cx = seq_len(gx), cy = seq_len(gy))
  W <- matrix(0, n, n)
  idx_of <- function(i, j) (j - 1) * gx + i
  for (i in seq_len(gx)) for (j in seq_len(gy)) {
    a <- idx_of(i, j)
    if (i < gx) { b <- idx_of(i + 1, j); W[a, b] <- 1; W[b, a] <- 1 }
    if (j < gy) { b <- idx_of(i, j + 1); W[a, b] <- 1; W[b, a] <- 1 }
  }
  u <- 0.9 * scale(coord$cx)[, 1] + 0.7 * scale(coord$cy)[, 1] +
       1.2 * exp(-((coord$cx - 5)^2 + (coord$cy - 5)^2) / 6)
  u <- u - mean(u)
  psi <- sig(u)
  z <- rbinom(n, 1, psi)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i]) y[i, ] <- rbinom(J, 1, p_det)
  list(W = W, y = y, psi = psi, n = n, J = J, p = p_det)
}

# ICAR structure matrix R = D - W (intrinsic, rank n-1). Diffuse global level
# via a tiny ridge eps*I (stands in for the improper intercept).
fit_field <- function(D, q, R, tau, eps = 1e-4) {
  n <- length(D)
  Q <- tau * R + eps * diag(n)
  b <- rep(0, n)
  for (it in seq_len(100)) {
    s  <- sig(b)
    mu <- q * s
    mu <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
    grad_lik <- (D - q * s) * (1 - s) / (1 - mu)      # dl/db per site
    info     <- q * s * (1 - s)^2 / (1 - mu)          # Fisher info per site
    info[q == 0] <- 0; grad_lik[q == 0] <- 0
    H  <- diag(info) + Q
    rhs <- grad_lik - as.vector(Q %*% b)
    step <- solve(H, rhs)
    b <- b + step
    if (max(abs(step)) < 1e-8) break
  }
  s  <- sig(b); mu <- pmin(pmax(q * s, 1e-10), 1 - 1e-10)
  info <- q * s * (1 - s)^2 / (1 - mu); info[q == 0] <- 0
  H <- diag(info) + Q
  Hinv <- solve(H)
  # Laplace log-marginal: data ll + prior + 0.5 log|Q*| - 0.5 log|H|.
  ll  <- sum(ifelse(q == 0, 0, D * log(mu) + (1 - D) * log(1 - mu)))
  # |Q| via its eps-regularised determinant (proper here).
  logdetQ <- determinant(Q, logarithm = TRUE)$modulus
  logdetH <- determinant(H, logarithm = TRUE)$modulus
  lp  <- -0.5 * as.numeric(t(b) %*% Q %*% b)
  logmarg <- ll + lp + 0.5 * as.numeric(logdetQ) - 0.5 * as.numeric(logdetH)
  list(b = b, var = diag(Hinv), logmarg = logmarg)
}

run_one <- function(seed) {
  d <- make_grid(seed = seed)
  heldout <- seq(2, d$n, by = 4)
  Jvec <- rep(d$J, d$n); Jvec[heldout] <- 0L          # held-out: no visits
  q <- 1 - (1 - d$p)^Jvec                              # per-site detection prob
  Dbin <- as.integer(rowSums(d$y) > 0); Dbin[heldout] <- 0L
  deg <- rowSums(d$W); R <- diag(deg) - d$W

  tau_grid <- exp(seq(log(0.2), log(20), length.out = 15))
  fits <- lapply(tau_grid, function(tau) fit_field(Dbin, q, R, tau))
  lm   <- vapply(fits, function(f) f$logmarg, numeric(1))
  w    <- exp(lm - max(lm)); w <- w / sum(w)

  # Marginalise psi over the tau grid: eta mixture sum_k w_k N(b_k[i], var_k[i]).
  B   <- sapply(fits, function(f) f$b)                 # n x K
  V   <- sapply(fits, function(f) f$var)               # n x K
  gh  <- { n <- 15; k <- seq_len(n - 1); J <- matrix(0, n, n)
           J[cbind(k, k + 1)] <- sqrt(k / 2); J[cbind(k + 1, k)] <- sqrt(k / 2)
           e <- eigen(J, symmetric = TRUE); o <- order(e$values)
           list(x = e$values[o], wt = sqrt(pi) * e$vectors[1, o]^2) }
  psi_mean <- numeric(d$n); psi_lo <- numeric(d$n); psi_hi <- numeric(d$n)
  for (i in seq_len(d$n)) {
    m <- B[i, ]; s <- sqrt(pmax(V[i, ], 1e-12))
    # mean via GH over each cell's Gaussian, weighted by w
    em <- 0
    for (g in seq_along(gh$x))
      em <- em + (gh$wt[g] / sqrt(pi)) * sum(w * sig(m + sqrt(2) * s * gh$x[g]))
    psi_mean[i] <- em
    cdf <- function(t) sum(w * pnorm(t, m, s))
    lo <- min(m - 8 * s); hi <- max(m + 8 * s)
    q025 <- uniroot(function(t) cdf(t) - 0.025, c(lo, hi))$root
    q975 <- uniroot(function(t) cdf(t) - 0.975, c(lo, hi))$root
    psi_lo[i] <- sig(q025); psi_hi[i] <- sig(q975)
  }
  cov_i <- d$psi >= psi_lo & d$psi <= psi_hi
  list(ho_cov = mean(cov_i[heldout]), ob_cov = mean(cov_i[-heldout]),
       ho_w = median((psi_hi - psi_lo)[heldout]),
       ho_mae = mean(abs(psi_mean[heldout] - d$psi[heldout])),
       ess = 1 / sum(w^2))
}

res <- lapply(1:8, run_one)
cat(sprintf("held-out coverage : %.3f\n", mean(sapply(res, `[[`, "ho_cov"))))
cat(sprintf("observed coverage : %.3f\n", mean(sapply(res, `[[`, "ob_cov"))))
cat(sprintf("held-out CI width : %.3f\n", mean(sapply(res, `[[`, "ho_w"))))
cat(sprintf("held-out MAE      : %.3f\n", mean(sapply(res, `[[`, "ho_mae"))))
cat(sprintf("grid ESS (of 15)  : %.2f\n", mean(sapply(res, `[[`, "ess"))))
