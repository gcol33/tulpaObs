# Validate the detection-RE Laplace path: (1) AGHQ reduces the small-cluster
# attenuation of sigma_det vs the raw nAGQ=1 EM, (2) the detection-intercept and
# occupancy fixed effects keep near-nominal CI coverage. Larger ng so per-group
# detection info is reasonable; truth sigma = 0.8.

suppressMessages(devtools::load_all("."))

sim_det_re <- function(seed, N = 400L, J = 6L, ng = 40L,
                       b0 = 0.4, b1 = -0.7, d0 = 0.2, sigma = 0.8) {
  set.seed(seed)
  occ_cov  <- rnorm(N)
  observer <- sample.int(ng, N, replace = TRUE)
  b_obs    <- rnorm(ng, 0, sigma)
  psi <- plogis(b0 + b1 * occ_cov)
  z   <- rbinom(N, 1, psi)
  p   <- plogis(d0 + b_obs[observer])
  y   <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p[i])
  list(data = data.frame(occ_cov = occ_cov, observer = factor(observer)),
       y = y, sigma = sigma, d0 = d0, b1 = b1)
}

fit_one <- function(s, aghq) {
  tobs(~ occ_cov, detection = ~ (1 | observer), family = occu(),
       data = s$data, y = s$y, method = "laplace",
       control = list(verbose = FALSE, re.aghq = aghq))
}

seeds <- 1:12
sig_em <- sig_aq <- numeric(length(seeds))
cov_d0 <- cov_b1 <- logical(length(seeds))
for (k in seq_along(seeds)) {
  s <- sim_det_re(seeds[k])
  fe <- fit_one(s, FALSE); fa <- fit_one(s, TRUE)
  sig_em[k] <- fe$means[[grep("^sigma_", names(fe$means), value = TRUE)]]
  sig_aq[k] <- fa$means[[grep("^sigma_", names(fa$means), value = TRUE)]]
  ci <- confint(fa)
  d0_row <- grep("^p_\\(Intercept\\)$", rownames(ci))
  b1_row <- grep("occ_cov", rownames(ci))
  cov_d0[k] <- ci[d0_row, 1] <= s$d0 && s$d0 <= ci[d0_row, 2]
  cov_b1[k] <- ci[b1_row, 1] <= s$b1 && s$b1 <= ci[b1_row, 2]
}

cat(sprintf("truth sigma_det = 0.80\n"))
cat(sprintf("EM   (nAGQ=1): mean %.3f  bias %+.3f\n", mean(sig_em), mean(sig_em) - 0.8))
cat(sprintf("AGHQ default : mean %.3f  bias %+.3f\n", mean(sig_aq), mean(sig_aq) - 0.8))
cat(sprintf("AGHQ |bias| < EM |bias|: %s\n",
            abs(mean(sig_aq) - 0.8) < abs(mean(sig_em) - 0.8)))
cat(sprintf("CI coverage (12 seeds): p_(Intercept) %.0f%%  occ slope %.0f%%\n",
            100 * mean(cov_d0), 100 * mean(cov_b1)))
