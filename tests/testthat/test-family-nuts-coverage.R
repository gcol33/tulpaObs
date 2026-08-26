# =============================================================================
# Family-level NUTS calibration: 20-seed 95% CI coverage on the coefficients the
# family NUTS path exists to calibrate. The single-seed recovery tests in each
# family's own test file prove the byte-exact C++<->R oracle and point recovery;
# these prove the intervals cover at the nominal rate. The estimand is the
# calibrated coefficient interval: NUTS samples the EXACT closed-form marginal, so
# an unbiased posterior should contain truth ~95% of the time even for
# weakly-identified coefficients (a wide but honest CI still covers).
#
# Pooled over all coefficients x 20 seeds, assert >= the 0.85 rubric floor (a
# calibrated sampler measures ~0.93-0.95; the floor absorbs Monte-Carlo slack).
# All heavy, so skip_if_fast()-gated -- zero cost in the fast suite.
#
# Each block also asserts the FIT count before scoring the rate (#276 item 1):
# `covered` pools several coefficients per fit, so a regression that makes
# most seeds error would otherwise leave the whole estimand measured on
# whichever seed or two survived and still read as a passing rate.
# =============================================================================

.nuts_ci_cover <- function(fit, truth) {
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  abs(est - truth) <= 1.96 * se        # per-coefficient 95% Wald CI containment
}

# The single-season occupancy NUTS path (cpp_occu_fit: occu / dyn_occu / int_occu)
# reports the pooled draws but not a per-parameter posterior SD, so read the 95%
# credible interval straight off the draws. `cols` selects the named coefficients
# to score (truth aligns to them positionally).
.nuts_ci_cover_draws <- function(fit, truth, cols = seq_along(truth)) {
  d <- fit$draws[, cols, drop = FALSE]
  q <- apply(d, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
  truth >= q[1, ] & truth <= q[2, ]
}

cuts5 <- seq(0, 1, length.out = 6)

test_that("distance NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(45), 0.4); beta_sigma <- c(log(0.45), 0.2)
  truth <- c(beta_lambda, beta_sigma)
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_distance(N = 120, cutpoints = cuts5, key = "halfnorm",
                             transect = "line", n_abund_covs = 1, n_sigma_covs = 1,
                             beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                             seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data,
                    family = distance(key = "halfnorm", transect = "line",
                                      cutpoints = sim$cutpoints),
                    detection = ~ sigma_cov1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

test_that("removal NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(7), 0.5); beta_p <- c(0.3, -0.3)
  truth <- c(beta_lambda, beta_p)
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_removal(N = 80, K = 5, n_abund_covs = 1, n_det_covs = 1,
                            beta_lambda = beta_lambda, beta_p = beta_p, seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data, family = removal(),
                    detection = ~ det_cov1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

test_that("fp_occu NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.6)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_fp_occu(N = 300, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                            p11 = 0.6, p10 = 0.05, b = 0.5, seed = 600 + s)
    fit <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                    detection = ~ 1, y = sim$y, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

test_that("dyn_abun NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(6), 0.4)
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2))
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_dyn_abun(N = 70, T = 3, J = 3, n_abund_covs = 1,
                             beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                             gamma = 1.2, seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 26),
                    detection = ~ 1, y = sim$y, method = "nuts",
                    control = list(n.iter = 250L, n.warmup = 250L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

test_that("abun NUTS coefficient 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  beta_lambda <- c(log(4), 0.5, -0.3); beta_p <- c(0.2, -0.4)
  truth <- c(beta_lambda, beta_p)
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_abun(N = 60, J = 4, n_abund_covs = 2, n_det_covs = 1,
                         beta_lambda = beta_lambda, beta_p = beta_p,
                         mixture = "poisson", seed = 600 + s)
    fit <- tryCatch(tobs(~ abund_cov1 + abund_cov2, data = sim$data, y = sim$y,
                    family = abun(), detection = ~ det_cov1, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   adapt.delta = 0.9, verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

# --- Single-season occupancy NUTS (cpp_occu_fit path) -----------------------
# The capability table lists dyn_occu / int_occu NUTS as "Yes", but NUTS
# appeared there only as a gated-error smoke check with no recovery test.
# These score the calibrated occupancy/transition intervals from the
# pooled draws (this path reports no per-parameter SD).

test_that("dyn_occu NUTS transition-parameter 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  # psi1 / p / gamma / eps intercepts on the logit scale.
  truth <- c(qlogis(0.5), qlogis(0.5), qlogis(0.3), qlogis(0.2))
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    set.seed(200 + s)
    n <- 160L; Tn <- 4L; J <- 4L; psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p <- 0.5
    z <- matrix(NA_integer_, n, Tn); z[, 1] <- rbinom(n, 1, psi1)
    for (t in 2:Tn) z[, t] <- ifelse(z[, t - 1] == 1, rbinom(n, 1, 1 - eps),
                                     rbinom(n, 1, gam))
    y <- array(0L, c(n, J, Tn))
    for (i in 1:n) for (t in 1:Tn) y[i, , t] <- if (z[i, t]) rbinom(J, 1, p) else 0L
    fit <- tryCatch(tobs(~ 1, data = data.frame(idx = 1:n), family = dyn_occu(),
                    y = y, detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
                    method = "nuts",
                    control = list(n.iter = 700L, n.warmup = 350L, seed = 1L,
                                   verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      covered <- c(covered, .nuts_ci_cover_draws(fit, truth))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})

test_that("int_occu NUTS occupancy 95% CIs cover at the nominal rate", {
  skip_on_cran(); skip_if_fast()
  truth <- c(0, 0.4)                      # psi_(Intercept), psi_x
  covered <- logical(0); n_fit <- 0L
  for (s in seq_len(20L)) {
    sim <- simulate_int_occu(N_total = 300, n_data = 1L, J = 5L,
                             beta_occ = c(0, 0.4), beta_det = list(c(0, -0.3)),
                             seed = 30L + s)
    fit <- tryCatch(tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1,
                    y = sim$y, method = "nuts",
                    control = list(n.iter = 700L, n.warmup = 350L, seed = 1L,
                                   verbose = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      n_fit <- n_fit + 1L
      cols <- match(c("psi_(Intercept)", "psi_x"), colnames(fit$draws))
      covered <- c(covered, .nuts_ci_cover_draws(fit, truth, cols = cols))
    }
  }
  expect_true(n_fit >= floor(0.8 * 20L))
  expect_gte(mean(covered), 0.85)
})
