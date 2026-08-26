# Dail-Madsen open-population N-mixture (Poisson initial abundance, binomial
# survival, Poisson recruitment), non-spatial Laplace (analytic-gradient BFGS over
# the exact HMM forward marginal) + NUTS.
#
# Recovery-grade tests: point recovery against simulated truth + 95% CI coverage
# across seeds, plus a correctness anchor (the C++ forward log-lik against an
# independent R forward recursion), an FD check of the analytic gradient (which is
# forward-mode differentiation through the scaled forward algorithm), and a
# byte-level C++ <-> R oracle cross-check.

# Independent R forward recursion for one site (the correctness anchor).
dm_fwd_R <- function(y_mat, lambda, p, omega, gamma, K) {
  S <- K + 1; n <- 0:K; T <- nrow(y_mat)
  obs <- function(t) {
    yt <- y_mat[t, ]; yt <- yt[!is.na(yt)]
    vapply(n, function(nn) if (!length(yt)) 1 else prod(dbinom(yt, nn, p)), numeric(1))
  }
  tr <- matrix(0, S, S)
  for (nn in 0:K) for (s in 0:nn) {
    g <- 0:(K - s)
    tr[nn + 1, s + g + 1] <- tr[nn + 1, s + g + 1] + dbinom(s, nn, omega) * dpois(g, gamma)
  }
  a <- dpois(n, lambda) * obs(1); ll <- log(sum(a)); a <- a / sum(a)
  for (t in 2:T) { pre <- as.vector(t(tr) %*% a); a <- pre * obs(t); ll <- ll + log(sum(a)); a <- a / sum(a) }
  ll
}

test_that("dyn_abun() family is wired and reports its supported methods", {
  f <- dyn_abun()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "dyn_abun")
  expect_identical(f$status, "working")
  expect_true(all(c("laplace", "nuts") %in% tulpaObs:::.tobs_family_methods$dyn_abun))
})

test_that("the C++ forward marginal matches an independent R forward recursion", {
  set.seed(1)
  K <- 40
  for (rep in 1:6) {
    T <- sample(2:4, 1); J <- sample(2:3, 1)
    lambda <- runif(1, 2, 8); p <- runif(1, 0.3, 0.7)
    omega <- runif(1, 0.4, 0.8); gamma <- runif(1, 0.3, 2)
    ymat <- matrix(rpois(T * J, lambda * p), T, J)
    yflat <- as.integer(t(ymat))                  # j + J*t layout
    out <- tulpaObs:::cpp_dyn_abun_total_log_lik(
      yflat, 1L, T, J, K, log(lambda), qlogis(p), qlogis(omega), log(gamma))
    expect_equal(out$log_lik, dm_fwd_R(ymat, lambda, p, omega, gamma, K),
                 tolerance = 1e-9, info = paste("rep", rep))
  }
})

test_that("analytic gradient matches finite differences (forward-mode diff)", {
  sim <- simulate_dyn_abun(N = 50, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(5), 0.3), p = 0.5, omega = 0.6,
                           gamma = 1, seed = 2)
  model <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim$data, sim$y, K_max = 30)
  lay  <- tulpaObs:::.tobs_dyn_abun_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_dyn_abun_nuts_marginal(model)
  theta <- c(log(5), 0.2, 0.1, qlogis(0.6), log(1))
  f <- function(th) tulpaObs:::.tobs_dyn_abun_nuts_logpost(th, marg, lay)$lp
  g_an <- tulpaObs:::.tobs_dyn_abun_nuts_logpost(theta, marg, lay)$grad
  g_fd <- sapply(seq_along(theta), function(j) {
    h <- 1e-5; tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  expect_equal(g_an, g_fd, tolerance = 1e-4)
})

test_that("C++ dyn_abun NUTS log-posterior matches the R oracle byte-for-byte", {
  sim <- simulate_dyn_abun(N = 40, T = 3, J = 2, n_abund_covs = 1, seed = 12)
  model <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim$data, sim$y, K_max = 25)
  lay  <- tulpaObs:::.tobs_dyn_abun_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_dyn_abun_nuts_marginal(model)
  theta <- c(log(5), 0.2, 0.1, qlogis(0.6), log(0.8))
  spec <- list(y = as.integer(model$y_flat), n_sites = model$n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
               X_omega = model$X_processes[[3]], X_gamma = model$X_processes[[4]])
  r_out <- tulpaObs:::.tobs_dyn_abun_nuts_logpost(theta, marg, lay, sigma.beta = 10)
  c_out <- tulpaObs:::cpp_dyn_abun_nuts_joint_logpost(spec, theta, 10)
  expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9)
  expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9)
})

test_that("dyn_abun Laplace recovers truth", {
  skip_if_fast()
  beta_lambda <- c(log(6), 0.4)
  sim <- simulate_dyn_abun(N = 250, T = 4, J = 3, n_abund_covs = 1,
                           beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                           gamma = 1.2, seed = 11)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 35),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(abs(est[2] - beta_lambda[2]), 0.2)
  expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                            "p_(Intercept)", "omega_(Intercept)", "gamma_(Intercept)"))
  expect_equal(unname(fit$intercepts$lambda), exp(est[1]), tolerance = 1e-8)
  expect_equal(unname(fit$intercepts$omega), plogis(est[4]), tolerance = 1e-8)
  expect_equal(unname(fit$intercepts$gamma), exp(est[5]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(6), 0.4)
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2))
  n_seed <- 20L
  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_abun(N = 150, T = 3, J = 3, n_abund_covs = 1,
                             beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                             gamma = 1.2, seed = 500 + s)
    fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 30),
                detection = ~ 1, y = sim$y, method = "laplace",
                control = list(verbose = FALSE))
    lo <- fit$means - 1.96 * fit$sds; hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.8),
              info = paste(round(cover_rate, 2), collapse = " | "))
})

test_that("S3 surface works for dyn_abun fits", {
  skip_if_fast()
  sim <- simulate_dyn_abun(N = 120, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(6), 0.4), seed = 3)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 28),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), sum(!is.na(sim$y)))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("lambda", "p", "omega", "gamma", "EN"))
  expect_length(fv$lambda, 120L)
  expect_true(all(fv$lambda > 0))
  expect_equal(dim(fv$EN), c(120L, 3L))

  pr <- predict(fit, X.0 = cbind(1, c(-1, 0, 1)), type = "lambda")
  expect_true(all(c("mean", "sd", "q2.5", "q50", "q97.5") %in% names(pr)))
  expect_true(all(pr$mean > 0) && all(diff(pr$mean) > 0))
  expect_true(all(pr$q2.5 <= pr$mean & pr$mean <= pr$q97.5))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim >= 0))

  rr <- residuals(fit, type = "pearson")
  expect_null(rr$occ)
  expect_equal(dim(rr$det), dim(sim$y))
})

test_that("dyn_abun rejects single-season data", {
  sim <- simulate_dyn_abun(N = 20, T = 4, J = 3, seed = 5)
  y1 <- sim$y[, , 1, drop = FALSE]
  expect_error(
    tobs(formula = ~ 1, data = sim$data, family = dyn_abun(), detection = ~ 1,
         y = y1, method = "laplace"),
    ">= 2 primary seasons")
})

test_that("dyn_abun NB analytic gradient (incl log_r) matches finite differences", {
  sim <- simulate_dyn_abun(N = 50, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = c(log(5), 0.3), p = 0.5, omega = 0.6,
                           gamma = 1, mixture = "negbin", r = 2, seed = 2)
  model <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim$data, sim$y,
                                           mixture = "negbin", K_max = 35)
  lay  <- tulpaObs:::.tobs_dyn_abun_nuts_layout(2L, 1L, 1L, 1L, use_nb = TRUE)
  expect_equal(lay$total, 6L)
  marg <- tulpaObs:::.tobs_dyn_abun_nuts_marginal(model)
  theta <- c(log(5), 0.2, 0.1, qlogis(0.6), log(1), log(2))   # last coord = log r
  f <- function(th) tulpaObs:::.tobs_dyn_abun_nuts_logpost(th, marg, lay)$lp
  g_an <- tulpaObs:::.tobs_dyn_abun_nuts_logpost(theta, marg, lay)$grad
  g_fd <- sapply(seq_along(theta), function(j) {
    h <- 1e-5; tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  expect_equal(g_an, g_fd, tolerance = 1e-4)
})

test_that("C++ dyn_abun NB NUTS log-posterior matches the R oracle byte-for-byte", {
  sim <- simulate_dyn_abun(N = 40, T = 3, J = 2, n_abund_covs = 1,
                           mixture = "negbin", r = 2.5, seed = 12)
  model <- tulpaObs:::.tobs_build_dyn_abun(~ abund_cov1, ~ 1, sim$data, sim$y,
                                           mixture = "negbin", K_max = 25)
  lay  <- tulpaObs:::.tobs_dyn_abun_nuts_layout(2L, 1L, 1L, 1L, use_nb = TRUE)
  marg <- tulpaObs:::.tobs_dyn_abun_nuts_marginal(model)
  theta <- c(log(5), 0.2, 0.1, qlogis(0.6), log(0.8), log(2.5))
  spec <- list(y = as.integer(model$y_flat), n_sites = model$n_sites,
               T = model$n_seasons, J = model$max_visits, K_max = model$K_max,
               X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
               X_omega = model$X_processes[[3]], X_gamma = model$X_processes[[4]],
               use_nb = TRUE)
  r_out <- tulpaObs:::.tobs_dyn_abun_nuts_logpost(theta, marg, lay, sigma.beta = 10)
  c_out <- tulpaObs:::cpp_dyn_abun_nuts_joint_logpost(spec, theta, 10)
  expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9)
  expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9)
})

test_that("dyn_abun negbin Laplace recovers truth (incl dispersion)", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(7), 0.4)
  sim <- simulate_dyn_abun(N = 300, T = 4, J = 3, n_abund_covs = 1,
                           beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                           gamma = 1.2, mixture = "negbin", r = 3, seed = 21)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data,
              family = dyn_abun(K_max = 45, mixture = "negbin"),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true("log_r" %in% names(fit$means))
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2), log(3))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(abs(est[2] - beta_lambda[2]), 0.25)
  # dispersion summary present and on the right order of magnitude
  expect_true(is.finite(fit$dispersion$r))
  expect_lt(abs(log(fit$dispersion$r) - log(3)), 0.8)
})

test_that("dyn_abun NUTS recovers truth and scores WAIC", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(6), 0.4)
  sim <- simulate_dyn_abun(N = 70, T = 3, J = 3, n_abund_covs = 1,
                           beta_lambda = beta_lambda, p = 0.5, omega = 0.6,
                           gamma = 1.2, seed = 31)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 26),
              detection = ~ 1, y = sim$y, method = "nuts",
              control = list(n.iter = 250L, n.warmup = 250L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))
  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 250L)
  truth <- c(beta_lambda, qlogis(0.5), qlogis(0.6), log(1.2))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 4))
  expect_lt(mean(fit$nuts$divergent), 0.2)
  w <- waic(fit)
  expect_true(is.finite(w$waic))
})


# --- NUTS + random effect ------------------------------------

# Dail-Madsen data with a per-site intercept RE on the initial-abundance arm.
sim_dyn_abun_lambda_re <- function(N, T, J, ngrp, beta_lambda, p, omega, gamma,
                                   sigma_b, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  lambda <- exp(as.numeric(model.matrix(~ x1, data) %*% beta_lambda) + b[grp])
  y <- array(0L, dim = c(N, J, T))
  for (i in seq_len(N)) {
    Ni <- stats::rpois(1L, lambda[i])
    for (t in seq_len(T)) {
      if (t > 1L) Ni <- stats::rbinom(1L, Ni, omega) + stats::rpois(1L, gamma)
      y[i, , t] <- stats::rbinom(J, Ni, p)
    }
  }
  list(y = y, data = data, sigma_b = sigma_b, beta_lambda = beta_lambda)
}

test_that("dyn_abun() NUTS samples a single initial-abundance RE and recovers it", {
  skip_on_cran()
  skip_if_fast()
  s <- sim_dyn_abun_lambda_re(N = 60, T = 3, J = 3, ngrp = 6,
                              beta_lambda = c(log(6), 0.3), p = 0.5, omega = 0.6,
                              gamma = 1, sigma_b = 0.5, seed = 9)
  fit <- tobs(formula = ~ x1 + (1 | g), data = s$data,
              family = dyn_abun(K_max = 28), detection = ~ 1, y = s$y,
              method = "nuts", verbose = FALSE,
              control = list(n.iter = 350L, n.warmup = 250L, seed = 1L))
  expect_identical(fit$method, "nuts")
  expect_identical(fit$re$arm, "lambda")
  expect_equal(fit$re$n_groups, 6L)
  expect_lt(abs(fit$re$sigma - 0.5), 0.45)
  expect_gt(fit$re$sigma_sd, 0)
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(6)), 0.4)
  expect_lt(mean(fit$nuts$divergent), 0.2)
})

# --- Laplace + random effect --------------------------------

test_that("dyn_abun() Laplace AGHQ recovers a site-grouped initial-abundance RE", {
  skip_on_cran()
  skip_if_fast()
  # The initial-abundance RE shifts only eta_lambda, which enters the season-1
  # initial distribution; the per-site conditional likelihood c(N_1) is integrated
  # against the RE by AGHQ on the exact HMM-forward marginal. The open model's
  # later seasons dilute the RE signal, so some residual sigma shrinkage remains
  # at small per-group n; the tolerance reflects that.
  est_sigma <- numeric(3); est_b0 <- numeric(3); est_b1 <- numeric(3)
  for (k in seq_len(3)) {
    s <- sim_dyn_abun_lambda_re(N = 120, T = 3, J = 3, ngrp = 15,
                                beta_lambda = c(log(4), 0.4), p = 0.5, omega = 0.5,
                                gamma = 0.5, sigma_b = 0.6, seed = 20 + k)
    fit <- tobs(formula = ~ x1 + (1 | g), data = s$data,
                family = dyn_abun(K_max = 30), detection = ~ 1, y = s$y,
                method = "laplace", verbose = FALSE,
                control = list(progress = FALSE, n.quad = 5L))
    if (k == 1L) {
      expect_identical(fit$method, "laplace")
      expect_identical(fit$dyn_abun_re$arm, "lambda")
      expect_true("sigma_g1_(Intercept)" %in% names(fit$means))
      re <- ranef(fit)
      expect_true(is.data.frame(re))
      expect_equal(nrow(re), 15L)
    }
    est_b0[k] <- fit$means[["lambda_(Intercept)"]]
    est_b1[k] <- fit$means[["lambda_x1"]]
    est_sigma[k] <- fit$means[["sigma_g1_(Intercept)"]]
  }
  expect_lt(abs(mean(est_b0) - log(4)), 0.3)
  expect_lt(abs(mean(est_b1) - 0.4), 0.2)
  expect_lt(abs(mean(est_sigma) - 0.6), 0.2)
  expect_gt(mean(est_sigma), 0.3)   # not collapsed to zero
})

# --- detection-arm random effect ----------------------------

# Dail-Madsen data with a per-site intercept RE on the detection (p) arm.
sim_dyn_abun_p_re <- function(N, T, J, ngrp, beta_lambda, beta_p, omega, gamma,
                              sigma_b, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  lambda <- exp(as.numeric(model.matrix(~ x1, data) %*% beta_lambda))
  pdet   <- stats::plogis(beta_p + b[grp])
  y <- array(0L, dim = c(N, J, T))
  for (i in seq_len(N)) {
    Ni <- stats::rpois(1L, lambda[i])
    for (t in seq_len(T)) {
      if (t > 1L) Ni <- stats::rbinom(1L, Ni, omega) + stats::rpois(1L, gamma)
      y[i, , t] <- stats::rbinom(J, Ni, pdet[i])
    }
  }
  list(y = y, data = data, sigma_b = sigma_b, beta_p = beta_p)
}

test_that("dyn_abun() Laplace AGHQ recovers a site-grouped detection RE", {
  skip_on_cran()
  skip_if_fast()
  # The detection RE shifts only eta_p, which enters every season's observation
  # pmf, so (unlike the initial-abundance arm) the per-site marginal is
  # re-evaluated through the exact HMM forward recursion per quadrature node with a
  # closed-form second-order eta_p forward-mode pass. Each site contributes T*J
  # detection observations, so the detection RE is well identified and AGHQ debias
  # recovers sigma_p with little attenuation.
  est_sigma <- numeric(3); est_p0 <- numeric(3)
  for (k in seq_len(3)) {
    s <- sim_dyn_abun_p_re(N = 100, T = 3, J = 4, ngrp = 12,
                           beta_lambda = c(log(8), 0.3), beta_p = stats::qlogis(0.5),
                           omega = 0.6, gamma = 1, sigma_b = 0.6, seed = 40 + k)
    fit <- tobs(formula = ~ x1, data = s$data, family = dyn_abun(K_max = 45),
                detection = ~ (1 | g), y = s$y, method = "laplace", verbose = FALSE,
                control = list(progress = FALSE, n.quad = 5L))
    if (k == 1L) {
      expect_identical(fit$method, "laplace")
      expect_identical(fit$dyn_abun_re$arm, "p")
      expect_true("sigma_p1_(Intercept)" %in% names(fit$means))
      re <- ranef(fit)
      expect_true(is.data.frame(re))
      expect_equal(nrow(re), 12L)
    }
    est_p0[k]    <- fit$means[["p_(Intercept)"]]
    est_sigma[k] <- fit$means[["sigma_p1_(Intercept)"]]
  }
  expect_lt(abs(mean(est_p0) - stats::qlogis(0.5)), 0.25)
  expect_lt(abs(mean(est_sigma) - 0.6), 0.2)
  expect_gt(mean(est_sigma), 0.35)   # not collapsed to zero
})

test_that("dyn_abun() NUTS samples a single detection RE and recovers it", {
  skip_on_cran()
  skip_if_fast()
  s <- sim_dyn_abun_p_re(N = 60, T = 3, J = 3, ngrp = 6,
                         beta_lambda = c(log(6), 0.3), beta_p = stats::qlogis(0.5),
                         omega = 0.6, gamma = 1, sigma_b = 0.5, seed = 14)
  fit <- tobs(formula = ~ x1, data = s$data, family = dyn_abun(K_max = 28),
              detection = ~ (1 | g), y = s$y, method = "nuts", verbose = FALSE,
              control = list(n.iter = 350L, n.warmup = 250L, n.chains = 1L,
                             seed = 1L))
  expect_identical(fit$method, "nuts")
  expect_identical(fit$re$arm, "p")
  expect_equal(fit$re$n_groups, 6L)
  expect_lt(abs(fit$re$sigma - 0.5), 0.45)
  expect_gt(fit$re$sigma_sd, 0)
  expect_true(any(grepl("^log_sigma_p_", colnames(fit$draws))))  # det-arm sigma
  expect_lt(mean(fit$nuts$divergent), 0.2)
})

test_that("dyn_abun() RE is one arm at a time (lambda OR p, not both)", {
  s <- sim_dyn_abun_p_re(N = 25, T = 3, J = 2, ngrp = 4,
                         beta_lambda = c(log(5), 0.2), beta_p = stats::qlogis(0.5),
                         omega = 0.6, gamma = 1, sigma_b = 0.4, seed = 2)
  # A random effect on BOTH the initial-abundance and detection arms in one fit is
  # rejected (the AGHQ path integrates one arm at a time).
  expect_error(
    tobs(formula = ~ x1 + (1 | g), data = s$data, family = dyn_abun(K_max = 25),
         detection = ~ (1 | g), y = s$y, method = "laplace",
         control = list(progress = FALSE)),
    "BOTH|one arm")
  # NUTS: a shared term across the two arms also routes to the one-arm error.
  expect_error(
    tobs(formula = ~ x1 + (1 | g), data = s$data, family = dyn_abun(K_max = 25),
         detection = ~ (1 | g), y = s$y, method = "nuts",
         control = list(n.iter = 20L, n.warmup = 10L)),
    "ONE arm|one arm|BOTH")
})


# --- areal spatial (ICAR / proper-CAR) on the initial-abundance arm (#51) ------

.sim_da_spatial <- function(adj, Tn = 3L, J = 2L, sd_phi = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ng <- nrow(adj)
  phi <- as.numeric(scale(rnorm(ng)))
  for (rep in 1:3) { pn <- phi; for (i in seq_len(ng)) { nbi <- which(adj[i,]==1L); pn[i] <- 0.5*phi[i]+0.5*mean(phi[nbi]) }; phi <- pn }
  phi <- sd_phi * as.numeric(scale(phi)); phi <- phi - mean(phi)
  x <- rnorm(ng); lam <- exp(log(6) + 0.5*x + phi)
  p <- 0.5; omega <- 0.5; gamma <- 0.6
  y <- array(0L, c(ng, J, Tn))
  for (i in seq_len(ng)) {
    Ni <- rpois(1, lam[i])
    for (t in seq_len(Tn)) {
      if (t > 1L) Ni <- rbinom(1, Ni, omega) + rpois(1, gamma)
      y[i, , t] <- rbinom(J, Ni, p)
    }
  }
  list(y = y, data = data.frame(abund_cov1 = x), phi = phi)
}

test_that("dyn_abun() areal ICAR recovers the initial-abundance slope + field", {
  skip_on_cran()
  skip_if_fast()
  # The field on log lambda_1 is fit by BFGS over the exact forward-HMM marginal
  # (cpp_dyn_abun_total_log_lik gradient) + the ICAR prior, with an FD-Hessian
  # observed-info Laplace marginal integrated over tau.
  adj <- rook_adj(5L)   # 25 sites
  slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
  for (s in 1:2) {
    sim <- .sim_da_spatial(adj, Tn = 3L, J = 2L, seed = 500 + s)
    fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
                family = dyn_abun(K_max = 25), detection = ~ 1, y = sim$y,
                method = "nested_laplace", control = list(progress = FALSE, verbose = FALSE))
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_true("lambda_abund_cov1" %in% names(fit$means))
      V <- vcov(fit)
      expect_true(all(is.finite(V)))
      expect_true(all(eigen(V, only.values = TRUE)$values > 0))
      expect_false(is.null(fit$spatial_field))
    }
    est <- fit$means[["lambda_abund_cov1"]]; se <- fit$sds[["lambda_abund_cov1"]]
    slopes <- c(slopes, est)
    slope_ok <- c(slope_ok, abs(est - 0.5) / se < 3.5)
    field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
  }
  expect_true(all(slope_ok))
  expect_lt(abs(mean(slopes) - 0.5), 0.2)
  expect_gt(mean(field_cor), 0.6)
})

test_that("dyn_abun() areal spatial: bym2 fits; nuts+icar samples (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(4L)
  s <- .sim_da_spatial(adj, Tn = 3L, J = 2L, seed = 3)
  fit <- tobs(formula = ~ abund_cov1 + bym2(graph = adj), data = s$data,
              family = dyn_abun(K_max = 20), detection = ~ 1, y = s$y,
              method = "nested_laplace", control = list(progress = FALSE, verbose = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_true(all(is.finite(vcov(fit))))
  expect_false(is.null(fit$spatial_field))
  # NUTS + areal now samples an intrinsic icar() field via the #71 sum-to-zero
  # reparameterisation (full recovery lives in test-count-spatial-nuts.R).
  fit_icar <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = s$data,
    family = dyn_abun(K_max = 20), detection = ~ 1, y = s$y, method = "nuts",
    control = list(n.iter = 40L, n.warmup = 40L, max.treedepth = 8L,
                   verbose = FALSE, progress = FALSE))
  expect_identical(fit_icar$method, "nuts")
  expect_lt(abs(mean(fit_icar$spatial_field)), 1e-6)
})

test_that("dyn_abun() bym2 + proper-CAR recover the initial-abundance field + slope (#131)", {
  skip_on_cran()
  skip_if_fast()
  # bym2 fits the rho-mixed unit field z = sqrt(rho) * phi + sqrt(1 - rho) * theta
  # over the exact forward-HMM marginal; proper-CAR (absent from the dyn_abun
  # suite before #131) a full-rank precision. Same 25-site fixture as the icar
  # recovery above; the forward-HMM marginal keeps this to two seeds.
  adj <- rook_adj(5L)
  for (term in c("bym2", "car_proper")) {
    tf <- if (term == "bym2") (~ abund_cov1 + bym2(graph = adj)) else
                              (~ abund_cov1 + car_proper(graph = adj))
    slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
    for (s in 1:2) {
      sim <- .sim_da_spatial(adj, Tn = 3L, J = 2L, seed = 500 + s)
      fit <- tobs(formula = tf, data = sim$data, family = dyn_abun(K_max = 25),
                  detection = ~ 1, y = sim$y, method = "nested_laplace",
                  control = list(progress = FALSE, verbose = FALSE))
      est <- fit$means[["lambda_abund_cov1"]]; se <- fit$sds[["lambda_abund_cov1"]]
      slopes    <- c(slopes, est)
      slope_ok  <- c(slope_ok, abs(est - 0.5) / se < 3.5)
      field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
    }
    expect_true(all(slope_ok), info = term)
    expect_lt(abs(mean(slopes) - 0.5), 0.2)
    expect_gt(mean(field_cor), 0.6)
  }
})

test_that("dyn_abun() temporal()-only field recovers the AR1 field + slope (#114)", {
  skip_on_cran()
  skip_if_fast()
  # A temporal() term on its own (no areal field) runs the shared areal-BFGS
  # driver with a single temporal block on the initial-abundance arm (#114).
  Tt <- 6L; per_t <- 20L; N <- Tt * per_t; Tocc <- 3L; J <- 2L
  fcor <- slope <- rep(NA_real_, 6L)
  for (s in seq_len(6L)) {
    set.seed(800L + s)
    period <- rep(seq_len(Tt), each = per_t)
    rho <- 0.7; sig <- 0.6; u <- numeric(Tt)
    u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
    for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
    u <- u - mean(u)
    x <- stats::rnorm(N)
    lam1 <- exp(log(6) + 0.5 * x + u[period])
    y <- array(0L, c(N, J, Tocc)); Ncur <- stats::rpois(N, lam1)
    for (t in 1:Tocc) {
      for (j in 1:J) y[, j, t] <- stats::rbinom(N, Ncur, 0.5)
      if (t < Tocc) Ncur <- stats::rbinom(N, Ncur, 0.6) + stats::rpois(N, 2)
    }
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period),
                         family = dyn_abun(K_max = 30L), detection = ~ 1, y = y,
                         method = "nested_laplace",
                         control = list(verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_null(fit$spatial_field)                 # temporal-only: no areal field
      expect_length(fit$temporal_field, Tt)
    }
    slope[s] <- fit$means[["lambda_x"]]
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.66)
  expect_lt(abs(mean(slope[ok]) - 0.5), 0.12)        # initial-abundance slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.8)       # AR1 temporal field recovered
})

test_that("dyn_abun() NUTS + temporal field samples the AR1 field, 0 divergences (#114)", {
  skip_on_cran()
  skip_if_fast()
  # A fixed-hyper non-centered temporal field rides the SAME dyn_abun NUTS field
  # block as the areal field (nuts_field_block.h), with field_map = period index
  # and a temporal whitened loading fixed at the nested-Laplace temporal-only
  # posterior mean. Structurally identical to the recovery-tested NUTS + areal
  # path (test-count-spatial-nuts.R); the forward-HMM makes long chains costly, so
  # this asserts it samples cleanly (0 divergences) and recovers the AR1 field --
  # a fuller recovery lives in the probe (dev_notes; cor 0.999). ar1 is full rank,
  # so the field is NOT sum-to-zero constrained (unlike the intrinsic areal case).
  Tt <- 6L; per_t <- 20L; N <- Tt * per_t; Tocc <- 3L; J <- 2L
  set.seed(801L)
  period <- rep(seq_len(Tt), each = per_t)
  rho <- 0.7; sig <- 0.6; u <- numeric(Tt)
  u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
  for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
  u <- u - mean(u)
  x <- stats::rnorm(N)
  lam1 <- exp(log(6) + 0.5 * x + u[period])
  y <- array(0L, c(N, J, Tocc)); Ncur <- stats::rpois(N, lam1)
  for (t in 1:Tocc) {
    for (j in 1:J) y[, j, t] <- stats::rbinom(N, Ncur, 0.5)
    if (t < Tocc) Ncur <- stats::rbinom(N, Ncur, 0.6) + stats::rpois(N, 2)
  }
  fit <- tobs(~ x + temporal(period, type = "ar1"),
              data = data.frame(x = x, period = period),
              family = dyn_abun(K_max = 30L), detection = ~ 1, y = y, method = "nuts",
              control = list(n.iter = 40L, n.warmup = 40L, max.treedepth = 6L,
                             verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nuts")
  expect_equal(fit$nuts$divergent_total, 0)
  expect_length(fit$temporal_field, Tt)
  expect_null(fit$spatial_field)
  expect_gt(abs(stats::cor(fit$temporal_field, u)), 0.6)
})

test_that("dyn_abun() NUTS + temporal + areal errors (temporal-only under NUTS, #114)", {
  # NUTS + temporal runs the fixed-hyper temporal field ON ITS OWN; a simultaneous
  # areal field composes with temporal only under nested_laplace.
  adj <- rook_adj(4L); ng <- nrow(adj)
  set.seed(9); period <- rep(1:4, length.out = ng)
  y <- array(rpois(ng * 2L * 3L, 3), c(ng, 2L, 3L))
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj) + temporal(period, type = "ar1"),
         data = data.frame(abund_cov1 = rnorm(ng), period = period),
         family = dyn_abun(K_max = 20L), detection = ~ 1, y = y, method = "nuts",
         control = list(verbose = FALSE, progress = FALSE)),
    "temporal"
  )
})
