# Tier A regression tests for the simplified Laplace path.
# See dev_notes/simplified_laplace_derivation.md §2 and §3.

test_that("approx = 'gaussian_laplace' is the default and produces NULL skew", {
  set.seed(101)
  n_sites <- 50; n_visits <- 3
  y <- matrix(rbinom(n_sites * n_visits, 1, 0.3), n_sites, n_visits)
  fit <- tobs(
    formula   = ~ 1,
    data      = data.frame(x = rep(0, n_sites)),
    family    = occu(),
    detection = ~ 1,
    y         = y,
    control   = list(verbose = FALSE)
  )
  expect_null(fit$skew)
  expect_identical(fit$sla_status, "off")
})

test_that("approx = 'simplified_laplace' attaches a numeric skew vector", {
  set.seed(102)
  n_sites <- 80; n_visits <- 4
  elev <- rnorm(n_sites)
  y <- matrix(rbinom(n_sites * n_visits, 1, plogis(0.2 + 0.5 * elev) * 0.6),
              n_sites, n_visits)
  fit <- tobs(
    formula   = ~ elev,
    data      = data.frame(elev = elev),
    family    = occu(),
    detection = ~ 1,
    y         = y,
    method    = "laplace_sla",
    control   = list(verbose = FALSE)
  )
  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew))
  expect_equal(length(fit$skew), length(fit$means))
  expect_named(fit$skew, names(fit$means))
  expect_true(all(is.finite(fit$skew)))
})

test_that("SLA draws have approximately the requested skew (small-gamma regime)", {
  set.seed(103)
  n_sites <- 200; n_visits <- 5
  elev <- rnorm(n_sites)
  y <- matrix(0L, n_sites, n_visits)
  psi <- plogis(0.3 + 0.8 * elev)
  z <- rbinom(n_sites, 1, psi)
  p <- 0.5
  for (i in seq_len(n_sites)) {
    y[i, ] <- if (z[i]) rbinom(n_visits, 1, p) else 0L
  }

  fit <- tobs(
    formula   = ~ elev,
    data      = data.frame(elev = elev),
    family    = occu(),
    detection = ~ 1,
    y         = y,
    method    = "laplace_sla",
    control   = list(verbose = FALSE)
  )

  emp <- apply(fit$draws, 2, function(d) {
    m <- mean(d); s <- sd(d)
    mean(((d - m) / s)^3)
  })

  # Only check parameters with |requested gamma| <= 0.5 — outside that
  # range the cumulant expansion saturates (see derivation §2.6).
  small <- which(abs(fit$skew) <= 0.5 & abs(fit$skew) > 0.02)
  if (length(small) > 0) {
    # 1000 SN draws have MC sigma_skew ~ sqrt(15 / N) ~ 0.12, so 0.2 is a
    # generous tolerance to account for MC noise + cumulant-expansion error.
    expect_true(all(abs(emp[small] - fit$skew[small]) < 0.25))
  }
})

test_that("SLA gracefully no-ops when EM has no detections at all", {
  set.seed(104)
  n_sites <- 30; n_visits <- 2
  y <- matrix(0L, n_sites, n_visits)
  fit <- tobs(
    formula   = ~ 1,
    data      = data.frame(x = rep(0, n_sites)),
    family    = occu(),
    detection = ~ 1,
    y         = y,
    method    = "laplace_sla",
    control   = list(verbose = FALSE)
  )
  expect_true(fit$sla_status %in% c("simplified_laplace",
                                     paste0("fallback_gaussian (", "EM fits missing for occ or det block)"),
                                     "fallback_gaussian (Louis I_obs (occ) not invertible)"))
})

test_that("Gaussian likelihood regression: l3_gaussian_identity == 0", {
  eta <- rnorm(10, 0, 2)
  expect_identical(tulpaObs:::.l3_gaussian_identity(eta), rep(0, 10))
})

test_that("Binomial l3 has expected sign structure", {
  # Vanishes at p = 1/2:
  expect_equal(tulpaObs:::.l3_binomial_logit(0, 1), 0)
  # Negative for p < 1/2, positive for p > 1/2 (since 1-2p flips sign):
  expect_true(tulpaObs:::.l3_binomial_logit(-2, 1) < 0)  # p < 1/2
  expect_true(tulpaObs:::.l3_binomial_logit( 2, 1) > 0)  # p > 1/2
  # Linear in n_trials:
  expect_equal(tulpaObs:::.l3_binomial_logit(1, 10),
               10 * tulpaObs:::.l3_binomial_logit(1, 1))
})

test_that("tulpa::sn_match round-trips for a moderately skewed distribution", {
  sn <- tulpa::sn_match(mu = 0, sigma = 1, gamma = 0.3)
  expect_false(is.null(sn))
  # Validate SN moments numerically via the local Azzalini sampler
  set.seed(1); s <- tulpaObs:::.sn_sample(50000, sn)
  expect_equal(mean(s), 0, tolerance = 0.02)
  expect_equal(sd(s), 1, tolerance = 0.02)
  m <- mean(s); ss <- sd(s)
  emp_skew <- mean(((s - m) / ss)^3)
  expect_equal(emp_skew, 0.3, tolerance = 0.05)
})

test_that("tulpa::sn_match returns NULL (with warning) above the SN ceiling", {
  # tulpa::sn_match warns + returns NULL above the SN representability ceiling
  expect_warning(res1 <- tulpa::sn_match(0, 1, 1.5), "exceeds skew-normal ceiling")
  expect_null(res1)
  expect_warning(res2 <- tulpa::sn_match(0, 1, -0.999), "exceeds skew-normal ceiling")
  expect_null(res2)
  # NaN gamma is a programmer bug upstream — sn_match errors rather than no-ops.
  expect_error(tulpa::sn_match(0, 1, NaN), "finite numeric")
})

test_that("FD gamma matches analytical gamma on a controlled single-season occu", {
  # Head-to-head numerical check: hold beta_hat and Sigma fixed by hand,
  # compute gamma both via the closed-form .sla_gamma_diag and via the
  # generic .sla_gamma_fd against .loglik_occu_single. If they disagree
  # by more than a few percent, the FD path is wrong and the non-diagonal
  # extensions can't be trusted.
  set.seed(2027)
  n_sites <- 200; n_visits <- 4
  elev <- rnorm(n_sites)
  X_occ <- cbind(`(Intercept)` = 1, elev = elev)
  X_det <- cbind(`(Intercept)` = rep(1, n_sites))
  beta_occ <- c(0.4, 0.7)
  beta_det <- -0.2
  psi_true <- plogis(X_occ %*% beta_occ)
  z <- rbinom(n_sites, 1, psi_true)
  p <- plogis(beta_det)
  y <- matrix(0L, n_sites, n_visits)
  for (i in seq_len(n_sites)) {
    y[i, ] <- if (z[i]) rbinom(n_visits, 1, p) else 0L
  }
  model <- structure(list(
    y           = y,
    X_processes = list(X_occ, X_det),
    n_sites     = n_sites,
    max_visits  = n_visits,
    process_info = list(
      list(name = "psi", p = 2L, coef_names = c("(Intercept)", "elev")),
      list(name = "p",   p = 1L, coef_names = "(Intercept)")
    ),
    model_type = "single"
  ), class = "tobs_model")

  # Hand-construct Sigma matching what the analytical SLA path would use:
  # full Louis I_obs for occ; for det, the original-likelihood weighted
  # binomial information at the mode.
  beta_hat <- c(beta_occ, beta_det)
  psi_hat  <- as.numeric(plogis(X_occ %*% beta_occ))
  p_hat    <- as.numeric(plogis(X_det %*% beta_det))

  # E-step weights for any-det and no-det sites
  any_det <- rowSums(y > 0) > 0
  n_valid <- rowSums(y >= 0)
  q_i     <- (1 - p_hat)^n_valid
  w_z1    <- ifelse(any_det, 1, psi_hat * q_i / (psi_hat * q_i + 1 - psi_hat))

  # Louis I_obs (occ): X' diag(psi(1-psi) - w(1-w)) X
  d_occ <- psi_hat * (1 - psi_hat) - w_z1 * (1 - w_z1)
  I_obs_occ <- t(X_occ) %*% diag(d_occ) %*% X_occ
  Sigma_occ <- solve(I_obs_occ)

  # Det H at original likelihood: only any-det sites contribute
  d_det <- ifelse(any_det, n_valid * p_hat * (1 - p_hat), 0)
  H_det <- t(X_det) %*% diag(d_det) %*% X_det
  Sigma_det <- solve(H_det)

  # Block-diagonal joint Sigma
  p_occ <- 2L; p_det <- 1L
  Sigma <- matrix(0, p_occ + p_det, p_occ + p_det)
  Sigma[seq_len(p_occ), seq_len(p_occ)] <- Sigma_occ
  Sigma[p_occ + seq_len(p_det), p_occ + seq_len(p_det)] <- Sigma_det

  # Analytical gamma: closed-form l3_occ (any-det vs no-det branches) +
  # .sla_gamma_diag on each block separately.
  pp     <- psi_hat * (1 - psi_hat)
  l3_occ <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    if (any_det[i]) {
      l3_occ[i] <- -psi_hat[i] * (1 - psi_hat[i]) * (1 - 2 * psi_hat[i])
    } else {
      one_m_q <- 1 - q_i[i]
      up      <- pp[i] * one_m_q
      upp     <- pp[i] * (1 - 2 * psi_hat[i]) * one_m_q
      uppp    <- pp[i] * (1 - 6 * pp[i]) * one_m_q
      v       <- 1 - psi_hat[i] * one_m_q
      vp      <- -up; vpp <- -upp; vppp <- -uppp
      l3_occ[i] <- vppp / v - 3 * vpp * vp / v^2 + 2 * (vp / v)^3
    }
  }
  gamma_occ_ana <- tulpaObs:::.sla_gamma_diag(l3_occ, X_occ, Sigma_occ)

  l3_det <- ifelse(any_det, -n_valid * p_hat * (1 - p_hat) * (1 - 2 * p_hat), 0)
  gamma_det_ana <- tulpaObs:::.sla_gamma_diag(l3_det, X_det, Sigma_det)
  gamma_ana <- c(gamma_occ_ana, gamma_det_ana)

  # FD gamma against the same Sigma + same likelihood
  log_lik_fn <- function(beta) tulpaObs:::.loglik_occu_single(beta, model)
  gamma_fd   <- tulpaObs:::.sla_gamma_fd(beta_hat, Sigma, log_lik_fn)

  # On the same Sigma + likelihood the two paths should agree to FD
  # truncation precision. Compare absolute |fd - ana| against a tolerance
  # that handles both regimes: for |ana| > 0.05, demand <= 5% relative;
  # for |ana| <= 0.05, demand <= 0.01 absolute (FD truncation floor).
  gamma_ana_u <- unname(gamma_ana)
  for (j in seq_along(gamma_ana_u)) {
    abs_err <- abs(gamma_fd[j] - gamma_ana_u[j])
    tol <- max(0.01, 0.05 * abs(gamma_ana_u[j]))
    expect_lt(abs_err, tol,
              label = sprintf("|gamma_fd - gamma_ana| at coef %d (fd = %.5f, ana = %.5f)",
                              j, gamma_fd[j], gamma_ana_u[j]))
  }
})

test_that("sla_replace_draws caps |gamma| > 0.95 and tracks it in attrs", {
  set.seed(105)
  draws <- matrix(rnorm(2000), 1000, 2)
  colnames(draws) <- c("a", "b")
  means <- c(0, 0); sds <- c(1, 1)
  gamma <- c(0.2, 1.5)        # 'b' over the cap (0.95)
  res <- tulpaObs:::.sla_replace_draws(draws, means, sds, gamma)
  expect_true("b" %in% attr(res, "sla_clipped"))
  # Empirical skew of 'a' should be in the right ballpark — 1000 SN draws
  # have MC sigma_skew ~ sqrt(15/N) ~ 0.12, so 0.25 is a loose tolerance
  # that still catches sign-flip or order-of-magnitude bugs.
  m <- mean(res[, "a"]); s <- sd(res[, "a"])
  expect_equal(mean(((res[, "a"] - m) / s)^3), 0.2, tolerance = 0.25)
})
