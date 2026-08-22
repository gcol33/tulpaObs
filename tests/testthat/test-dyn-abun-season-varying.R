# Season-varying survival (omega) and recruitment (gamma) for the Dail-Madsen
# open N-mixture (dyn_abun()). The transition from season t-1 to t uses
# interval-(t-1) vital rates, so a season covariate on omega / gamma drives
# the dynamics. The forward HMM kernel carries length-(T-1) interval-indexed
# eta with per-interval forward-mode gradients; the constant-rate (scalar)
# path is the broadcast special case.
#
# Tests: an FD anchor on the per-interval analytic gradient, a no-regression
# guard that the scalar broadcast is bit-identical to a per-interval vector of
# the same constant, the long-form binder shape, and slope recovery for a
# season-varying omega and gamma fit against simulated truth.


test_that("per-interval omega/gamma analytic gradient matches finite differences", {
  set.seed(7)
  T <- 4L; J <- 3L; K <- 18L; nIv <- T - 1L
  lambda <- 6; p <- 0.5
  omega_iv <- c(0.45, 0.65, 0.55); gamma_iv <- c(0.8, 1.4, 0.6)
  Ni <- rpois(1, lambda); ymat <- matrix(0L, T, J)
  for (t in seq_len(T)) {
    if (t > 1L) Ni <- rbinom(1, Ni, omega_iv[t - 1L]) + rpois(1, gamma_iv[t - 1L])
    ymat[t, ] <- rbinom(J, Ni, p)
  }
  yflat <- as.integer(t(ymat))
  eta_l <- log(lambda); eta_p <- qlogis(p)
  eo <- qlogis(omega_iv); eg <- log(gamma_iv)
  ll <- function(o, g) tulpaObs:::cpp_dyn_abun_total_log_lik(
    yflat, 1L, T, J, K, eta_l, eta_p, as.numeric(o), as.numeric(g))$log_lik
  an <- tulpaObs:::cpp_dyn_abun_total_log_lik(yflat, 1L, T, J, K, eta_l, eta_p, eo, eg)
  h <- 1e-6
  g_om_fd <- numeric(nIv); g_gm_fd <- numeric(nIv)
  for (iv in seq_len(nIv)) {
    ep <- eo; em <- eo; ep[iv] <- ep[iv] + h; em[iv] <- em[iv] - h
    g_om_fd[iv] <- (ll(ep, eg) - ll(em, eg)) / (2 * h)
    ep <- eg; em <- eg; ep[iv] <- ep[iv] + h; em[iv] <- em[iv] - h
    g_gm_fd[iv] <- (ll(eo, ep) - ll(eo, em)) / (2 * h)
  }
  expect_equal(as.numeric(an$grad_eta_omega), g_om_fd, tolerance = 1e-5)
  expect_equal(as.numeric(an$grad_eta_gamma), g_gm_fd, tolerance = 1e-5)
})


test_that("scalar (constant-rate) kernel is bit-identical to the broadcast path", {
  # The per-site scalar eta and a per-interval vector of that same constant must
  # give the identical log-likelihood and (summed) score -- the no-regression
  # guard that the season-varying kernel did not perturb the constant-rate fit.
  set.seed(3)
  T <- 4L; J <- 3L; K <- 25L; nIv <- T - 1L
  lambda <- 7; p <- 0.5; om_c <- 0.55; gm_c <- 1.1
  Ni <- rpois(1, lambda); ymat <- matrix(0L, T, J)
  for (t in seq_len(T)) {
    if (t > 1L) Ni <- rbinom(1, Ni, om_c) + rpois(1, gm_c)
    ymat[t, ] <- rbinom(J, Ni, p)
  }
  yflat <- as.integer(t(ymat))
  sc <- tulpaObs:::cpp_dyn_abun_total_log_lik(
    yflat, 1L, T, J, K, log(lambda), qlogis(p), qlogis(om_c), log(gm_c))
  vc <- tulpaObs:::cpp_dyn_abun_total_log_lik(
    yflat, 1L, T, J, K, log(lambda), qlogis(p),
    rep(qlogis(om_c), nIv), rep(log(gm_c), nIv))
  expect_identical(sc$log_lik, vc$log_lik)
  expect_identical(as.numeric(sc$grad_eta_omega), sum(as.numeric(vc$grad_eta_omega)))
  expect_identical(as.numeric(sc$grad_eta_gamma), sum(as.numeric(vc$grad_eta_gamma)))
})


test_that("binder builds long-form omega design only when season-varying", {
  sim_sv <- simulate_dyn_abun(N = 20, T = 4, J = 2, n_abund_covs = 1,
                              beta_lambda = c(log(5), 0.3), p = 0.5, gamma = 1,
                              beta_omega = c(qlogis(0.6), 0.5), seed = 9)
  m_sv <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim_sv$data, sim_sv$y,
                                          omega_formula = ~ season_cov, K_max = 20)
  expect_true(m_sv$omega_season_varying)
  expect_false(m_sv$gamma_season_varying)
  expect_equal(nrow(m_sv$X_processes[[3]]), 20L * 3L)   # (site x interval) rows
  expect_equal(nrow(m_sv$X_processes[[4]]), 20L)        # per-site gamma

  # Constant-rate formula collapses to the per-site design (broadcast path).
  m_c <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim_sv$data, sim_sv$y,
                                         omega_formula = ~ 1, K_max = 20)
  expect_false(m_c$omega_season_varying)
  expect_equal(nrow(m_c$X_processes[[3]]), 20L)
})


test_that("dyn_abun recovers a season-varying survival slope", {
  skip_on_cran()
  skip_if_fast()
  beta_omega <- c(qlogis(0.6), 0.8)         # logit intercept + season slope
  est <- numeric(3); se <- numeric(3); int <- numeric(3)
  for (k in seq_len(3)) {
    sim <- simulate_dyn_abun(N = 300, T = 5, J = 3, n_abund_covs = 1,
                             beta_lambda = c(log(7), 0.3), p = 0.55, gamma = 1.2,
                             beta_omega = beta_omega, seed = 200 + k)
    fit <- tobs(formula = ~ abund_cov1, data = sim$data,
                family = dyn_abun(K_max = 35), detection = ~ 1, y = sim$y,
                omega = ~ season_cov, method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    if (k == 1L) {
      expect_s3_class(fit, "tobs_fit")
      expect_true("omega_season_cov" %in% names(fit$means))
      expect_true(all(is.finite(vcov(fit))))
    }
    est[k] <- fit$means[["omega_season_cov"]]
    se[k]  <- fit$sds[["omega_season_cov"]]
    int[k] <- fit$means[["omega_(Intercept)"]]
  }
  # Slope within 3 SE per seed; mean slope + intercept near truth.
  expect_true(all(abs(est - beta_omega[2]) / se < 3))
  expect_lt(abs(mean(est) - beta_omega[2]), 0.25)
  expect_lt(abs(mean(int) - beta_omega[1]), 0.25)
})


test_that("dyn_abun recovers a season-varying recruitment slope", {
  skip_on_cran()
  skip_if_fast()
  beta_gamma <- c(log(1.0), 0.6)            # log intercept + season slope
  est <- numeric(3); se <- numeric(3); int <- numeric(3)
  for (k in seq_len(3)) {
    sim <- simulate_dyn_abun(N = 300, T = 5, J = 3, n_abund_covs = 1,
                             beta_lambda = c(log(7), 0.3), p = 0.55, omega = 0.6,
                             beta_gamma = beta_gamma, seed = 300 + k)
    fit <- tobs(formula = ~ abund_cov1, data = sim$data,
                family = dyn_abun(K_max = 35), detection = ~ 1, y = sim$y,
                gamma = ~ season_cov, method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    if (k == 1L) {
      expect_true("gamma_season_cov" %in% names(fit$means))
      expect_equal(nrow(fit$model$X_processes[[4]]), 300L * 4L)  # long-form gamma
    }
    est[k] <- fit$means[["gamma_season_cov"]]
    se[k]  <- fit$sds[["gamma_season_cov"]]
    int[k] <- fit$means[["gamma_(Intercept)"]]
  }
  expect_true(all(abs(est - beta_gamma[2]) / se < 3))
  expect_lt(abs(mean(est) - beta_gamma[2]), 0.2)
  expect_lt(abs(mean(int) - beta_gamma[1]), 0.25)
})


test_that("constant-rate dyn_abun fit is unchanged (no-regression at fit level)", {
  skip_on_cran()
  skip_if_fast()
  # A constant-rate (intercept-only omega/gamma) fit must recover truth exactly as
  # the pre-season-varying kernel did; the simulator's constant-rate RNG stream is
  # preserved, so this anchors the whole constant-rate path end to end.
  sim <- simulate_dyn_abun(N = 250, T = 4, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(6), 0.4), p = 0.5, omega = 0.6,
                           gamma = 1.2, seed = 11)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 35),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  truth <- c(log(6), 0.4, qlogis(0.5), qlogis(0.6), log(1.2))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                            "p_(Intercept)", "omega_(Intercept)", "gamma_(Intercept)"))
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(abs(est[4] - qlogis(0.6)), 0.2)
})
