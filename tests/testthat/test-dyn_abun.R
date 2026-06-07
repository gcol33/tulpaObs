# Dail-Madsen open-population N-mixture (Poisson initial abundance, binomial
# survival, Poisson recruitment), non-spatial Laplace (analytic-gradient BFGS over
# the exact HMM forward marginal) + NUTS (gcol33/tulpaObs#37).
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
  expect_true(all(pr > 0) && all(diff(pr) > 0))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim >= 0))

  rr <- residuals(fit, type = "pearson")
  expect_equal(dim(rr), dim(sim$y))
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
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
})
