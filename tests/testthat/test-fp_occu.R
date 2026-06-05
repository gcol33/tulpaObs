# Multistate false-positive occupancy (Miller et al. 2011 confirmed-detection
# design), non-spatial Laplace (analytic-gradient BFGS) + NUTS (gcol33/tulpaObs#40).
#
# Recovery-grade tests: point recovery against simulated truth + 95% CI coverage
# across seeds, plus a closed-form correctness anchor (the two-state marginal
# against a direct R computation, and the certain-detection identity), an FD check
# of the analytic gradient, and a byte-level C++ <-> R oracle cross-check.

fp_R_marginal <- function(y, psi, p11, p10, b) {
  # Direct two-state marginal for one site (y in {0,1,2}).
  g1 <- ifelse(y == 0, 1 - p11, ifelse(y == 1, p11 * (1 - b), p11 * b))
  g0 <- ifelse(y == 0, 1 - p10, ifelse(y == 1, p10, 0))
  log(psi * prod(g1) + (1 - psi) * prod(g0))
}

test_that("fp_occu() family is wired and reports its supported methods", {
  f <- fp_occu()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "fp_occu")
  expect_identical(f$status, "working")
  expect_true(all(c("laplace", "nuts") %in% tulpaObs:::.tobs_family_methods$fp_occu))
})

test_that("the marginal matches a direct two-state computation", {
  set.seed(1)
  for (rep in 1:8) {
    J <- sample(3:6, 1)
    psi <- runif(1, 0.2, 0.8); p11 <- runif(1, 0.3, 0.8)
    p10 <- runif(1, 0.01, 0.15); b <- runif(1, 0.2, 0.8)
    y <- sample(0:2, J, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    out <- tulpaObs:::cpp_fp_occu_total_log_lik(
      as.integer(y), rep(1L, J), qlogis(psi), qlogis(p11), qlogis(p10), qlogis(b))
    expect_equal(out$log_lik, fp_R_marginal(y, psi, p11, p10, b), tolerance = 1e-9)
  }
  # Certain-detection identity: any y == 2 makes the unoccupied term zero, so the
  # posterior occupancy is exactly 1 and L = psi * prod(g1).
  y <- c(0L, 2L, 1L)
  out <- tulpaObs:::cpp_fp_occu_total_log_lik(
    y, rep(1L, 3), qlogis(0.4), qlogis(0.6), qlogis(0.05), qlogis(0.5))
  expect_equal(out$w1, 1, tolerance = 1e-12)
  g1 <- c(1 - 0.6, 0.6 * 0.5, 0.6 * 0.5)
  expect_equal(out$log_lik, log(0.4) + sum(log(g1)), tolerance = 1e-9)
})

test_that("analytic gradient matches finite differences", {
  sim <- simulate_fp_occu(N = 100, J = 5, n_occ_covs = 1, seed = 7)
  model <- tulpaObs:::.tobs_build_fp_occu(~ occ_cov1, ~ 1, sim$data, sim$y)
  lay  <- tulpaObs:::.tobs_fp_occu_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_fp_occu_nuts_marginal(model)
  theta <- c(0.2, 0.4, 0.3, qlogis(0.06), 0.1)
  f <- function(th) tulpaObs:::.tobs_fp_occu_nuts_logpost(th, marg, lay)$lp
  g_an <- tulpaObs:::.tobs_fp_occu_nuts_logpost(theta, marg, lay)$grad
  g_fd <- sapply(seq_along(theta), function(j) {
    h <- 1e-5; tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  expect_equal(g_an, g_fd, tolerance = 1e-4)
})

test_that("C++ fp_occu NUTS log-posterior matches the R oracle byte-for-byte", {
  sim <- simulate_fp_occu(N = 60, J = 5, n_occ_covs = 1, seed = 12)
  model <- tulpaObs:::.tobs_build_fp_occu(~ occ_cov1, ~ 1, sim$data, sim$y)
  lay  <- tulpaObs:::.tobs_fp_occu_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_fp_occu_nuts_marginal(model)
  theta <- c(0.1, 0.3, 0.2, qlogis(0.05), -0.1)
  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = model$X_processes[[1]], X_p11 = model$X_processes[[2]],
               X_p10 = model$X_processes[[3]], X_b = model$X_processes[[4]],
               n_sites = model$n_sites)
  r_out <- tulpaObs:::.tobs_fp_occu_nuts_logpost(theta, marg, lay, sigma.beta = 10)
  c_out <- tulpaObs:::cpp_fp_occu_nuts_joint_logpost(spec, theta, 10)
  expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9)
  expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9)
})

test_that("fp_occu Laplace recovers truth", {
  skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.7)
  sim <- simulate_fp_occu(N = 600, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                          p11 = 0.6, p10 = 0.05, b = 0.5, seed = 11)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(abs(est[2] - beta_psi[2]), 0.2)
  expect_named(fit$means, c("psi_(Intercept)", "psi_occ_cov1",
                            "p11_(Intercept)", "p10_(Intercept)", "b_(Intercept)"))
  expect_equal(unname(fit$intercepts$psi), plogis(est[1]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_psi <- c(qlogis(0.55), 0.6)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  n_seed <- 30L
  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_fp_occu(N = 400, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                            p11 = 0.6, p10 = 0.05, b = 0.5, seed = 400 + s)
    fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
                detection = ~ 1, y = sim$y, method = "laplace",
                control = list(verbose = FALSE))
    lo <- fit$means - 1.96 * fit$sds; hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.85),
              info = paste(round(cover_rate, 2), collapse = " | "))
})

test_that("fp_occu supports a covariate on the false-positive arm", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  N <- 600; J <- 6
  x <- rnorm(N)
  psi <- plogis(0.3 + 0.6 * x); p11 <- 0.6; b <- 0.5
  p10 <- plogis(qlogis(0.05) + 0.5 * x)        # covariate-varying false positives
  z <- rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- rbinom(J, 1L, p11); cert <- rbinom(J, 1L, b)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else y[i, ] <- rbinom(J, 1L, p10[i])
  }
  dat <- data.frame(x = x)
  fit <- tobs(formula = ~ x, data = dat, family = fp_occu(), detection = ~ 1,
              y = y, fp_formula = ~ x, method = "laplace",
              control = list(verbose = FALSE))
  est <- fit$means
  expect_true("p10_x" %in% names(est))
  expect_lt(abs(unname(est["p10_x"]) - 0.5), 0.3)
  expect_lt(abs(unname(est["psi_x"]) - 0.6), 0.2)
})

test_that("S3 surface works for fp_occu fits", {
  skip_if_fast()
  sim <- simulate_fp_occu(N = 300, J = 5, n_occ_covs = 1,
                          beta_psi = c(qlogis(0.5), 0.6), seed = 3)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), sum(!is.na(sim$y)))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "p11", "p10", "b", "z"))
  expect_length(fv$psi, 300L)
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$z >= 0 & fv$z <= 1))

  pr <- predict(fit, X.0 = cbind(1, c(-1, 0, 1)), type = "psi")
  expect_length(pr, 3L)
  expect_true(all(diff(pr) > 0))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim %in% 0:2))

  rr <- residuals(fit, type = "pearson")
  expect_equal(dim(rr), dim(sim$y))
})

test_that("fp_occu rejects invalid detection states", {
  sim <- simulate_fp_occu(N = 30, J = 4, seed = 5)
  y_bad <- sim$y; y_bad[1, 1] <- 3L
  expect_error(
    tobs(formula = ~ 1, data = sim$data, family = fp_occu(), detection = ~ 1,
         y = y_bad, method = "laplace"),
    "detection states")
})

test_that("fp_occu NUTS recovers truth and scores WAIC", {
  skip_on_cran()
  skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.6)
  sim <- simulate_fp_occu(N = 300, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                          p11 = 0.6, p10 = 0.05, b = 0.5, seed = 31)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))
  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 500L)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(mean(fit$nuts$divergent), 0.2)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
})
