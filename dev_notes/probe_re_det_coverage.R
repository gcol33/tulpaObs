# Cleaner CI-coverage read for the detection-RE Laplace path (AGHQ default),
# 24 seeds. Targets the detection intercept (the parameter partially confounded
# with the RE mean) and the occupancy slope. Reports sigma_det recovery too.

suppressMessages(devtools::load_all("."))

sim_det_re <- function(seed, N = 400L, J = 6L, ng = 40L,
                       b0 = 0.4, b1 = -0.7, d0 = 0.2, sigma = 0.8) {
  set.seed(seed)
  occ_cov  <- rnorm(N)
  observer <- sample.int(ng, N, replace = TRUE)
  b_obs    <- rnorm(ng, 0, sigma)
  psi <- plogis(b0 + b1 * occ_cov); z <- rbinom(N, 1, psi)
  p   <- plogis(d0 + b_obs[observer])
  y   <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p[i])
  list(data = data.frame(occ_cov = occ_cov, observer = factor(observer)),
       y = y, sigma = sigma, d0 = d0, b0 = b0, b1 = b1)
}

seeds <- 1:24
sig <- numeric(length(seeds))
cov_d0 <- cov_b1 <- cov_b0 <- logical(length(seeds))
for (k in seq_along(seeds)) {
  s <- sim_det_re(seeds[k])
  fa <- tobs(~ occ_cov, detection = ~ (1 | observer), family = occu(),
             data = s$data, y = s$y, method = "laplace",
             control = list(verbose = FALSE))
  sig[k] <- fa$means[[grep("^sigma_", names(fa$means), value = TRUE)]]
  ci <- confint(fa)
  rr <- function(pat) grep(pat, rownames(ci))
  cov_d0[k] <- ci[rr("^p_\\(Intercept\\)$"), 1] <= s$d0 && s$d0 <= ci[rr("^p_\\(Intercept\\)$"), 2]
  cov_b0[k] <- ci[rr("psi_\\(Intercept\\)|^\\(Intercept\\)"), 1][1] <= s$b0 &&
               s$b0 <= ci[rr("psi_\\(Intercept\\)|^\\(Intercept\\)"), 2][1]
  cov_b1[k] <- ci[rr("occ_cov"), 1] <= s$b1 && s$b1 <= ci[rr("occ_cov"), 2]
}

cat(sprintf("sigma_det: mean %.3f (truth 0.80), bias %+.3f\n", mean(sig), mean(sig) - 0.8))
cat(sprintf("CI coverage over %d seeds:\n", length(seeds)))
cat(sprintf("  p_(Intercept)   %.0f%%\n", 100 * mean(cov_d0)))
cat(sprintf("  psi_(Intercept) %.0f%%\n", 100 * mean(cov_b0)))
cat(sprintf("  occ slope       %.0f%%\n", 100 * mean(cov_b1)))
