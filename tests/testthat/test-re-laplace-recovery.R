# Parameter-recovery tests for formula random effects fit under the
# DETERMINISTIC engine (gcol33/tulpaObs#11). The default Laplace engine now
# fits iid intercept RE and uncorrelated random slopes via a variance-component
# EM (R/em_laplace_re.R) instead of silently dropping them. Deterministic
# Laplace variance estimates for binary occupancy carry the usual small-cluster
# (PQL) bias, so the sigma tolerances are generous and the calibrated check is
# against the NUTS fit on the same data.

sim_occu_re_intercept <- function(seed = 101, ng = 30L, per = 25L, J = 6L,
                                  b0 = 0.3, b1 = -0.6, sigma = 0.9, p = 0.45) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  b_true <- rnorm(ng, 0, sigma)
  z <- rbinom(N, 1, plogis(b0 + b1 * x + b_true[g]))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), b_true = b_true)
}

# Single-season occupancy with a CORRELATED random intercept + slope:
# (b0_g, b1_g) ~ N(0, Sigma) per group, eta = b0 + b0_g + (b1 + b1_g) * x.
sim_occu_re_corr <- function(seed = 404, ng = 40L, per = 25L, J = 6L,
                             b0 = 0.2, b1 = -0.4, p = 0.5,
                             Sigma = matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  U <- matrix(rnorm(ng * 2L), ng, 2L) %*% chol(Sigma)  # rows ~ N(0, Sigma)
  eta <- b0 + U[g, 1] + (b1 + U[g, 2]) * x
  z <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), U = U, Sigma = Sigma)
}

test_that("iid intercept RE is fit (not dropped) by the default Laplace engine", {
  s <- sim_occu_re_intercept()
  fit <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace",
              control = list(verbose = FALSE))

  expect_identical(fit$method, "laplace")
  # The random effect is present, not silently dropped: a sigma hyperparameter
  # and per-group BLUPs are reported.
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 1L)
  expect_true(is.finite(fit$means[[sig_nm]]) && fit$means[[sig_nm]] > 0.4 &&
              fit$means[[sig_nm]] < 1.4)

  cf <- coef(fit)$psi
  expect_lt(abs(cf[["(Intercept)"]] - 0.3), 0.3)
  expect_lt(abs(cf[["x"]] - (-0.6)), 0.3)
  expect_lt(abs(plogis(fit$means[["p_(Intercept)"]]) - 0.45), 0.08)

  # Occupancy fixed-effect SE comes from the natural-scale observed info, not
  # the M-inflated M-step Hessian (which would be ~sqrt(M)=~31x too small).
  se_int <- fit$sds[["psi_(Intercept)"]]
  expect_true(is.finite(se_int) && se_int > 0.03 && se_int < 1)

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 30L)
  expect_gt(cor(re$estimate, s$b_true), 0.7)
})

test_that("uncorrelated random slopes (1 + x || g) recover under Laplace", {
  set.seed(202)
  ng <- 30L; per <- 30L; N <- ng * per; J <- 6L
  g <- rep(seq_len(ng), each = per); x <- rnorm(N)
  st <- c(0.7, 0.5)
  b0 <- rnorm(ng, 0, st[1]); b1 <- rnorm(ng, 0, st[2])
  z <- rbinom(N, 1, plogis(0.2 + b0[g] + (-0.4 + b1[g]) * x))
  y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * 0.5)
  d <- data.frame(g = factor(g), x = x)

  fit <- tobs(~ x + (1 + x || g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "laplace", control = list(verbose = FALSE))

  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 2L)   # intercept + slope sigma
  sig_hat <- fit$means[sig_nm]
  expect_true(all(is.finite(sig_hat) & sig_hat > 0.2 & sig_hat < 1.4))

  re <- ranef(fit)
  expect_equal(sort(unique(re$term)), c("(Intercept)", "x"))
  b0_hat <- re$estimate[re$term == "(Intercept)"]
  b1_hat <- re$estimate[re$term == "x"]
  expect_gt(cor(b0_hat, b0), 0.6)
  expect_gt(cor(b1_hat, b1), 0.4)
})

test_that("deterministic sigma and BLUPs track the NUTS fit on the same data", {
  skip_on_cran()
  s <- sim_occu_re_intercept(seed = 7, ng = 30L, per = 25L)

  fit_l <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                family = occu(), method = "laplace", control = list(verbose = FALSE))
  fit_n <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                family = occu(), method = "nuts",
                control = list(n.iter = 400, n.warmup = 200, seed = 1, verbose = FALSE))

  sig_l <- fit_l$means[[grep("^sigma_", names(fit_l$means), value = TRUE)]]
  sig_n <- exp(fit_n$means[[grep("^log_sigma_", names(fit_n$means), value = TRUE)]])
  # Deterministic Laplace sigma is in the NUTS ballpark (PQL bias is modest at
  # this cluster size; both should land within ~60% of each other).
  expect_lt(abs(sig_l - sig_n) / sig_n, 0.6)

  re_l <- ranef(fit_l); re_n <- ranef(fit_n)
  expect_gt(cor(re_l$estimate, re_n$estimate), 0.85)
})

test_that("correlated random slopes (1 + x | g) recover under Laplace", {
  skip_on_cran()
  s <- sim_occu_re_corr(seed = 402)

  fit <- tobs(~ x + (1 + x | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace", control = list(verbose = FALSE))

  expect_identical(fit$method, "laplace")

  # intercept + slope sigma, plus the off-diagonal reported as a correlation
  # (the block is genuinely full, not forced diagonal like `(1 + x || g)`).
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 2L)
  sig_hat <- fit$means[sig_nm]
  expect_true(all(is.finite(sig_hat) & sig_hat > 0.2 & sig_hat < 1.6))

  cor_nm <- grep("^cor_", names(fit$means), value = TRUE)
  expect_length(cor_nm, 1L)
  rho_hat <- fit$means[[cor_nm]]
  # Simulated correlation is +0.61; recover the sign / a positive association
  # (PQL attenuates the magnitude, so the bound is generous).
  expect_gt(rho_hat, 0.2)
  expect_lt(rho_hat, 0.99)

  cf <- coef(fit)$psi
  expect_lt(abs(cf[["(Intercept)"]] - 0.2), 0.3)
  expect_lt(abs(cf[["x"]] - (-0.4)), 0.3)

  re <- ranef(fit)
  expect_equal(sort(unique(re$term)), c("(Intercept)", "x"))
  b0_hat <- re$estimate[re$term == "(Intercept)"]
  b1_hat <- re$estimate[re$term == "x"]
  expect_gt(cor(b0_hat, s$U[, 1]), 0.6)
  expect_gt(cor(b1_hat, s$U[, 2]), 0.4)
})

test_that("RE forms the deterministic engine cannot fit error toward NUTS", {
  s <- sim_occu_re_intercept(seed = 3, ng = 12L, per = 12L)

  # Non-single model types route RE to NUTS (validator guard, exercised
  # directly to avoid building a full community dataset here).
  re_int <- tulpaObs:::.tobs_term_re(group = s$d$g, type = "intercept")
  stub <- list(model_type = "community", data = s$d, X_det_visit = NULL)
  expect_error(
    tulpaObs:::.validate_re_laplace(re_int, stub, NULL, "gaussian_laplace"),
    "single|nuts")

  # RE on the detection predictor is still NUTS-only.
  re_det <- tulpaObs:::.tobs_term_re(group = s$d$g, type = "intercept")
  re_det$shared <- c(FALSE, TRUE)
  stub2 <- list(model_type = "single", data = s$d, X_det_visit = NULL)
  expect_error(
    tulpaObs:::.validate_re_laplace(re_det, stub2, NULL, "gaussian_laplace"),
    "detection|nuts")
})
