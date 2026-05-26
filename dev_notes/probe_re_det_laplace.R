# De-risk milestone 1 (RE on the detection predictor, Laplace path).
#
# Before wiring a detection-arm RE into .tobs_em_laplace_re, confirm the
# UPSTREAM capability it relies on: tulpa_laplace() must recover a known
# random intercept from a *weighted* binomial fit (the occupancy detection
# M-step is a weighted binomial -- rows weighted by P(z=1|y)), with
# return_re_cov = TRUE giving sane per-group posterior covariance blocks.
#
# This is a pure tulpa probe (no occupancy EM): G observer groups, each with
# n_per site-visits, logit(p) = b0 + b_g, b_g ~ N(0, sigma^2). Fit with an RE
# on the observer grouping and per-row weights, check b_g / sigma recovery.

suppressMessages(devtools::load_all("."))

set.seed(1)
G       <- 30L     # observer groups
n_per   <- 12L     # site-visits per observer
J       <- 4L      # binomial trials per row (visits)
b0_true <- -0.3
sig_true <- 0.8

run_one <- function(seed) {
  set.seed(seed)
  obs   <- rep(seq_len(G), each = n_per)
  b_g   <- rnorm(G, 0, sig_true)
  eta   <- b0_true + b_g[obs]
  p     <- plogis(eta)
  N     <- G * n_per
  y     <- rbinom(N, J, p)
  X     <- matrix(1, N, 1L)
  # Mild row weights in (0.2, 1] to mimic the occupancy posterior weighting.
  w     <- runif(N, 0.2, 1)

  re_list <- list(list(idx = obs, n_groups = G, n_coefs = 1L,
                       sigma = sig_true))   # fixed-sigma fit (one M-step)
  fit <- tulpa::tulpa_laplace(
    y = y, n_trials = rep(J, N), X = X, weights = w,
    re_list = re_list, family = "binomial",
    return_hessian = TRUE, return_re_cov = TRUE)

  b0_hat <- fit$mode[1]
  b_hat  <- fit$mode[-1]
  list(b0 = b0_hat,
       cor_b = cor(b_hat, b_g),
       rmse_b = sqrt(mean((b_hat - b_g)^2)),
       n_cov = length(fit$cov_blocks),
       cov_ok = !is.null(fit$cov_blocks) && length(fit$cov_blocks) == G)
}

res <- lapply(1:8, run_one)
cat(sprintf("b0 truth %.2f -> mean est %.3f\n",
            b0_true, mean(vapply(res, `[[`, numeric(1), "b0"))))
cat(sprintf("cor(b_hat, b_true): mean %.3f  (min %.3f)\n",
            mean(vapply(res, `[[`, numeric(1), "cor_b")),
            min(vapply(res, `[[`, numeric(1), "cor_b"))))
cat(sprintf("rmse(b_hat): mean %.3f\n",
            mean(vapply(res, `[[`, numeric(1), "rmse_b"))))
cat(sprintf("cov_blocks returned & length G: %s\n",
            all(vapply(res, `[[`, logical(1), "cov_ok"))))
