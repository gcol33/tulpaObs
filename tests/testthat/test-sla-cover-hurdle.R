# Phase 3.5 SLA regression tests for cover(positive = "beta") and
# cover(positive = "lognormal").
#
# See R/sla_cover_hurdle.R and R/simplified_laplace.R.

# ---- shared simulators ----------------------------------------------------

simulate_beta_cover_local <- function(N = 400, beta_occ = c(-0.3, 0.7),
                                      beta_pos = c(0.4, -1.0), phi = 25,
                                      seed = 1) {
  set.seed(seed)
  x <- runif(N, -2, 2)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  mu      <- plogis(eta_pos)
  y <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos]  <- rbeta(sum(is_pos), mu[is_pos] * phi, (1 - mu[is_pos]) * phi)
  y[!is_pos] <- 0
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos, phi = phi))
}

simulate_lognormal_cover_local <- function(N = 400, beta_occ = c(-0.3, 0.6),
                                           beta_pos = c(-1.0, 0.4),
                                           sigma_pos = 0.5, seed = 2) {
  set.seed(seed)
  x <- runif(N, -2, 2)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  occur <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  is_pos <- occur == 1L
  y <- numeric(N)
  if (any(is_pos)) {
    lc <- rnorm(sum(is_pos), eta_pos[is_pos], sigma_pos)
    y[is_pos] <- pmin(exp(lc), 1 - 1e-6)
  }
  list(data = data.frame(x = x), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                    sigma_pos = sigma_pos))
}


# ---- smoke tests (beta) ----------------------------------------------------

test_that("cover(positive='beta') with simplified_laplace attaches skew + status", {
  sim <- simulate_beta_cover_local(N = 400, seed = 11)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y,
    approx  = "simplified_laplace"
  )
  expect_s3_class(fit, "cover_fit")
  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew_occ))
  expect_true(is.numeric(fit$skew_pos))
  expect_equal(length(fit$skew_occ), length(fit$beta_occ))
  expect_equal(length(fit$skew_pos), length(fit$beta_pos))
  expect_true(all(is.finite(fit$skew_occ)))
  expect_true(all(is.finite(fit$skew_pos)))
  expect_named(fit$skew_occ, colnames(fit$encoding$occ_data$X))
  expect_named(fit$skew_pos, colnames(fit$encoding$pos_data$X))
  # Draws should be filled in.
  expect_true(is.matrix(fit$draws_occ))
  expect_true(is.matrix(fit$draws_pos))
  expect_equal(nrow(fit$draws_occ), 1000L)
  expect_equal(nrow(fit$draws_pos), 1000L)
})


# ---- smoke tests (lognormal) ----------------------------------------------

test_that("cover(positive='lognormal') with simplified_laplace attaches skew + status", {
  sim <- simulate_lognormal_cover_local(N = 400, seed = 12)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("lognormal"),
    y       = sim$y,
    approx  = "simplified_laplace"
  )
  expect_s3_class(fit, "cover_fit")
  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew_occ))
  expect_true(is.numeric(fit$skew_pos))
  expect_true(all(is.finite(fit$skew_occ)))
  # The Lognormal arm's log-likelihood is Gaussian in beta_pos at fixed
  # sigma_pos, so its third derivative is identically zero. SLA should
  # therefore return gamma_pos == 0 to FD truncation precision.
  expect_true(all(abs(fit$skew_pos) < 1e-3))
})


# ---- default approx is gaussian_laplace -----------------------------------

test_that("cover() default leaves SLA off and produces no skew fields", {
  sim <- simulate_beta_cover_local(N = 300, seed = 13)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )
  expect_identical(fit$sla_status %||% "off", "off")
  expect_null(fit$skew_occ)
  expect_null(fit$skew_pos)
  expect_null(fit$draws_occ)
  expect_null(fit$draws_pos)
})


# ---- log-lik correctness (occurrence arm) ---------------------------------

test_that(".loglik_cover_occ is locally maximised at truth (Bernoulli arm)", {
  sim <- simulate_beta_cover_local(N = 800, seed = 21)
  enc <- tulpaObs:::encode_cover_hurdle(~ x, sim$data, sim$y,
                                        positive = "beta")
  ll_truth <- tulpaObs:::.loglik_cover_occ(sim$truth$beta_occ, enc)
  expect_true(is.finite(ll_truth))
  # Perturb each coefficient by +/- 0.3 and check the log-lik decreases.
  for (j in seq_along(sim$truth$beta_occ)) {
    b_up <- sim$truth$beta_occ; b_up[j] <- b_up[j] + 0.3
    b_dn <- sim$truth$beta_occ; b_dn[j] <- b_dn[j] - 0.3
    ll_up <- tulpaObs:::.loglik_cover_occ(b_up, enc)
    ll_dn <- tulpaObs:::.loglik_cover_occ(b_dn, enc)
    expect_lt(ll_up, ll_truth + 1e-6)
    expect_lt(ll_dn, ll_truth + 1e-6)
  }
})


# ---- log-lik correctness (beta positive arm) ------------------------------

test_that(".loglik_cover_pos_beta is locally maximised at truth", {
  sim <- simulate_beta_cover_local(N = 800, seed = 22)
  enc <- tulpaObs:::encode_cover_hurdle(~ x, sim$data, sim$y,
                                        positive = "beta")
  ll_truth <- tulpaObs:::.loglik_cover_pos_beta(
    sim$truth$beta_pos, sim$truth$phi, enc
  )
  expect_true(is.finite(ll_truth))
  for (j in seq_along(sim$truth$beta_pos)) {
    b_up <- sim$truth$beta_pos; b_up[j] <- b_up[j] + 0.3
    b_dn <- sim$truth$beta_pos; b_dn[j] <- b_dn[j] - 0.3
    ll_up <- tulpaObs:::.loglik_cover_pos_beta(b_up, sim$truth$phi, enc)
    ll_dn <- tulpaObs:::.loglik_cover_pos_beta(b_dn, sim$truth$phi, enc)
    expect_lt(ll_up, ll_truth + 1e-6)
    expect_lt(ll_dn, ll_truth + 1e-6)
  }
})


# ---- log-lik correctness (lognormal positive arm) -------------------------

test_that(".loglik_cover_pos_lognormal is locally maximised at truth", {
  sim <- simulate_lognormal_cover_local(N = 800, seed = 23)
  enc <- tulpaObs:::encode_cover_hurdle(~ x, sim$data, sim$y,
                                        positive = "lognormal")
  ll_truth <- tulpaObs:::.loglik_cover_pos_lognormal(
    sim$truth$beta_pos, sim$truth$sigma_pos, enc
  )
  expect_true(is.finite(ll_truth))
  for (j in seq_along(sim$truth$beta_pos)) {
    b_up <- sim$truth$beta_pos; b_up[j] <- b_up[j] + 0.3
    b_dn <- sim$truth$beta_pos; b_dn[j] <- b_dn[j] - 0.3
    ll_up <- tulpaObs:::.loglik_cover_pos_lognormal(b_up, sim$truth$sigma_pos, enc)
    ll_dn <- tulpaObs:::.loglik_cover_pos_lognormal(b_dn, sim$truth$sigma_pos, enc)
    expect_lt(ll_up, ll_truth + 1e-6)
    expect_lt(ll_dn, ll_truth + 1e-6)
  }
})


# ---- Sigma / gamma sanity --------------------------------------------------

test_that(".sla_compute_cover_hurdle returns finite gamma_occ / gamma_pos on a recoverable beta sim", {
  sim <- simulate_beta_cover_local(N = 500, seed = 24)
  enc <- tulpaObs:::encode_cover_hurdle(~ x, sim$data, sim$y,
                                        positive = "beta")
  fits <- tulpaObs:::fit_cover_hurdle(enc, positive = "beta")
  sla_res <- tulpaObs:::.sla_compute_cover_hurdle(fits, enc, "beta")
  expect_true(isTRUE(sla_res$valid))
  expect_true(all(is.finite(sla_res$gamma_occ)))
  expect_true(all(is.finite(sla_res$gamma_pos)))
})


# ---- Log-lik consistent with fit at the mode ------------------------------

test_that("log-likelihood evaluated at the fitted mode is finite and order-of-magnitude consistent with log_marginal", {
  sim <- simulate_beta_cover_local(N = 400, seed = 25)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )
  enc <- fit$encoding
  ll_occ <- tulpaObs:::.loglik_cover_occ(fit$beta_occ, enc)
  ll_pos <- tulpaObs:::.loglik_cover_pos_beta(fit$beta_pos, fit$phi_pos, enc)
  expect_true(is.finite(ll_occ))
  expect_true(is.finite(ll_pos))
  # ll_total and the log_marginal differ by the Laplace correction (log|H|
  # terms) but should remain in the same order of magnitude.
  ll_total <- ll_occ + ll_pos
  lm_total <- fit$log_marginal["occ"] + fit$log_marginal["pos"]
  expect_lt(abs(ll_total - lm_total) / max(abs(ll_total), 1), 0.5)
})
