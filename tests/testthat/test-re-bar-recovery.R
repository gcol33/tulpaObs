# Parameter-recovery tests for the lme4 bar-syntax random-effect blocks added
# in gcol33/tulpaObs#10: a correlated random intercept + 2 slopes (a 3x3
# block), and a slope-only block (`(0 + x | g)`, no group intercept).
#
# Random slopes are fit by the NUTS engine (the EM-Laplace path does not carry
# formula random effects; nested-Laplace rejects slopes). The NUTS parameter
# vector for a single RE term on the occupancy predictor with detection ~ 1 is
# laid out as
#   [psi_(Intercept), p_(Intercept), log_sigma(1..q), chol_raw(1..k(k-1)/2),
#    z_effects(group-major: g * q + c)]
# with q = re_n_coefs and the non-centered effects recovered as
#   b_{g,c} = sigma_c * (L %*% z_g)_c,  L = tanh-Cholesky(chol_raw).
# The tests assert the total parameter count (guards against layout drift, and
# pins the k(k-1)/2 Cholesky size from populate_helpers.h) and that the
# reconstructed group effects correlate with the simulated truth.

# Reconstruct the lower-triangular tanh-Cholesky factor from the off-diagonal
# raw parameters (strictly-lower, row-major), mirroring tulpa_priors_re.h.
tanh_chol <- function(raw, k) {
  L <- matrix(0, k, k)
  idx <- 1L
  for (r in seq_len(k)) {
    s2 <- 0
    for (cc in seq_len(r - 1L)) {
      L[r, cc] <- tanh(raw[idx]); s2 <- s2 + L[r, cc]^2; idx <- idx + 1L
    }
    L[r, r] <- sqrt(max(1 - s2, 1e-10))
  }
  L
}

# Draw ng x k correlated effects with covariance Sigma (no MASS dependency).
rmvn_rows <- function(ng, Sigma) {
  k <- nrow(Sigma)
  matrix(rnorm(ng * k), ng, k) %*% chol(Sigma)  # rows ~ N(0, Sigma)
}

test_that("(1 + x + z | g) recovers a 3x3 correlated RE block", {
  skip_on_cran()

  set.seed(101)
  ng <- 25L; per <- 20L; N <- ng * per; J <- 8L
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N); z <- rnorm(N)
  sig <- c(0.9, 0.7, 0.5)
  R <- matrix(c(1, .5, .3,  .5, 1, .2,  .3, .2, 1), 3, 3)
  Sigma <- diag(sig) %*% R %*% diag(sig)
  B <- rmvn_rows(ng, Sigma)                      # true (b0, b_x, b_z) per group
  eta <- 0.3 + B[g, 1] + B[g, 2] * x + B[g, 3] * z
  zocc <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, zocc[i] * 0.6)
  d <- data.frame(g = factor(g), x = x, z = z)

  fit <- tobs(~ (1 + x + z | g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.iter = 400, n.warmup = 200, seed = 1, verbose = FALSE))

  # Layout / structure: q = 3, k(k-1)/2 = 3 Cholesky params (NOT k(k+1)/2 = 6).
  expect_length(fit$re, 1L)
  expect_identical(fit$re[[1]]$type, "slope")
  expect_true(fit$re[[1]]$intercept)
  expect_equal(length(fit$means), 2L + 3L + 3L + ng * 3L)

  m <- fit$means
  sig_hat <- exp(m[3:5])
  expect_true(all(is.finite(sig_hat)) && all(sig_hat > 0))
  expect_true(all(sig_hat < 2))                  # not blown up

  L <- tanh_chol(m[6:8], 3L)
  zoff <- 8L
  Bhat <- t(vapply(seq_len(ng), function(gg) {
    as.numeric(sig_hat * (L %*% m[zoff + (gg - 1L) * 3L + 1:3]))
  }, numeric(3)))

  # The block is genuinely fit: reconstructed group effects track the truth
  # on all three coefficients (probe: 0.86 / 0.70 / 0.80 at N=600/iter=500).
  for (cc in 1:3) {
    expect_gt(cor(Bhat[, cc], B[, cc]), 0.45)
  }
})

test_that("(1 + x + z || g) is an uncorrelated multi-slope block (no Cholesky)", {
  skip_on_cran()

  set.seed(303)
  ng <- 20L; per <- 20L; N <- ng * per; J <- 8L
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N); z <- rnorm(N)
  sig <- c(0.8, 0.6, 0.5)
  B <- sapply(sig, function(s) rnorm(ng, 0, s))   # independent columns
  eta <- 0.3 + B[g, 1] + B[g, 2] * x + B[g, 3] * z
  zocc <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, zocc[i] * 0.6)
  d <- data.frame(g = factor(g), x = x, z = z)

  fit <- tobs(~ (1 + x + z || g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.iter = 300, n.warmup = 150, seed = 1, verbose = FALSE))

  # q = 3 but uncorrelated -> NO Cholesky block: [psi_int, p_int, sigma(3),
  # z_effects(ng*3)]. (Before per-term correlation handling, re_correlated was
  # left unsized for a fully-`||` block and the param layout read past its end.)
  expect_false(fit$re[[1]]$correlated)
  expect_equal(length(fit$means), 2L + 3L + ng * 3L)

  m <- fit$means
  sig_hat <- exp(m[3:5])
  expect_true(all(is.finite(sig_hat)) && all(sig_hat > 0) && all(sig_hat < 2))
  # Diagonal block: b_{g,c} = sigma_c * z_{g,c}.
  zoff <- 5L
  Bhat <- t(vapply(seq_len(ng), function(gg) {
    sig_hat * m[zoff + (gg - 1L) * 3L + 1:3]
  }, numeric(3)))
  for (cc in 1:3) expect_gt(cor(Bhat[, cc], B[, cc]), 0.45)
})

test_that("(0 + x | g) is a slope-only block with no group intercept", {
  skip_on_cran()

  set.seed(202)
  ng <- 25L; per <- 20L; N <- ng * per; J <- 8L
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  sig_slope <- 0.8
  b <- rnorm(ng, 0, sig_slope)                   # random slope, no intercept
  eta <- 0.3 + b[g] * x
  zocc <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, zocc[i] * 0.6)
  d <- data.frame(g = factor(g), x = x)

  fit <- tobs(~ (0 + x | g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.iter = 400, n.warmup = 200, seed = 1, verbose = FALSE))

  # Structural proof of "no intercept": q = 1 -> one effect per group, so the
  # param vector is [psi_int, p_int, log_sigma, z_effects(ng)]. An implicit
  # intercept would double the effect block to 2 * ng and add a 2nd sigma.
  expect_length(fit$re, 1L)
  expect_identical(fit$re[[1]]$type, "slope")
  expect_false(fit$re[[1]]$intercept)
  expect_equal(length(fit$means), 2L + 1L + ng)

  m <- fit$means
  sig_hat <- exp(m[3])
  expect_true(is.finite(sig_hat) && sig_hat > 0.3 && sig_hat < 1.8)

  # Uncorrelated single slope: b_g = sigma * z_g.
  bhat <- sig_hat * m[3 + seq_len(ng)]
  expect_gt(cor(bhat, b), 0.45)
})
