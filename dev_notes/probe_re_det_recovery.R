# End-to-end recovery for a DETECTION random intercept on the Laplace path
# (milestone 1: RE on the detection predictor). Occupancy psi = sigmoid(b0 +
# b1 occ_cov); detection p = sigmoid(d0 + b_obs[observer]), b_obs ~ N(0, sigma).
# Check: the fit succeeds (gate -> dual-arm EM -> AGHQ det branch), recovers the
# detection sigma and the fixed effects, and AGHQ reports arm = "det".

suppressMessages(devtools::load_all("."))

sim_det_re <- function(seed, N = 300L, J = 5L, ng = 30L,
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
       y = y, sigma = sigma, d0 = d0, b0 = b0, b1 = b1)
}

run <- function(seed) {
  s <- sim_det_re(seed)
  fit <- tobs(~ occ_cov, detection = ~ (1 | observer), family = occu(),
              data = s$data, y = s$y, method = "laplace",
              control = list(verbose = FALSE))
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  list(
    aghq_applied = isTRUE(fit$aghq$applied),
    arm = fit$aghq$arm %||% NA_character_,
    sig_name = sig_nm,
    sigma = unname(fit$means[sig_nm]),
    p0 = plogis(unname(fit$means[["p_(Intercept)"]])),
    b1 = unname(coef(fit)$psi[["occ_cov"]]),
    truth = s
  )
}

r1 <- run(1)
cat("AGHQ applied:", r1$aghq_applied, " arm:", r1$arm, "\n")
cat("sigma param name(s):", r1$sig_name, "\n")
cat(sprintf("sigma_det: est %.3f  truth %.3f\n", r1$sigma, r1$truth$sigma))
cat(sprintf("p intercept: est %.3f  truth %.3f\n", r1$p0, plogis(r1$truth$d0)))
cat(sprintf("occ slope b1: est %.3f  truth %.3f\n", r1$b1, r1$truth$b1))

cat("\n-- multi-seed sigma recovery (per-group n ~ 10) --\n")
sigs <- vapply(1:8, function(s) run(s)$sigma, numeric(1))
cat(sprintf("mean sigma_det over 8 seeds: %.3f (truth 0.8); range [%.3f, %.3f]\n",
            mean(sigs), min(sigs), max(sigs)))
