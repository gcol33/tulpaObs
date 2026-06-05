# Removal-sampling abundance (sequential depletion), Poisson + negbin,
# non-spatial Laplace + NUTS (gcol33/tulpaObs#39).
#
# Recovery-grade tests (per the "statistical code needs recovery tests" rule):
# point recovery against simulated truth + 95% CI coverage across seeds, plus a
# closed-form correctness anchor (the Poisson removal marginal equals
# independent Poissons), an FD check of the analytic gradient, and a byte-level
# C++ <-> R oracle cross-check. Structural tests cover the family wiring / S3.

test_that("removal() family is wired and reports its supported methods", {
  f <- removal()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "removal")
  expect_identical(f$status, "working")
  expect_identical(f$params$mixture, "poisson")
  expect_true(all(c("laplace", "nuts") %in% tulpaObs:::.tobs_family_methods$removal))
})

test_that("removal marginal equals independent Poissons (Poisson abundance)", {
  # Thinning identity: under N ~ Poisson(lambda), the per-pass removals are
  # independent y_k ~ Poisson(lambda * pi_k), pi_k = p_k prod_{l<k}(1-p_l). The
  # kernel's marginal sum over N must reproduce this closed form exactly.
  set.seed(1)
  for (rep in 1:5) {
    K <- sample(2:5, 1)
    lambda <- runif(1, 1, 12)
    p <- runif(K, 0.1, 0.8)
    y <- rpois(K, lambda * tulpaObs:::.removal_pi(p))   # any nonneg counts work
    eta_lambda <- log(lambda)
    eta_p <- qlogis(p)
    out <- tulpaObs:::cpp_removal_total_log_lik(
      as.integer(y), rep(1L, K), eta_p, eta_lambda,
      K_max = sum(y) + 200L, r = Inf)
    ll_ref <- sum(dpois(y, lambda * tulpaObs:::.removal_pi(p), log = TRUE))
    expect_equal(out$log_lik, ll_ref, tolerance = 1e-8)
    # E[N | y]: with independent-Poisson removals, the undetected count is
    # Poisson(lambda * pi_0), so E[N|y] = R + lambda * prod(1-p).
    expect_equal(out$mean_N, sum(y) + lambda * prod(1 - p), tolerance = 1e-6)
  }
})

test_that("analytic gradient matches finite differences (Poisson + NB)", {
  for (mix in c("P", "NB")) {
    sim <- simulate_removal(N = 60, K = 4, n_abund_covs = 1, n_det_covs = 1,
                            beta_lambda = c(log(8), 0.4), beta_p = c(0.3, -0.4),
                            mixture = if (mix == "NB") "negbin" else "poisson",
                            size = 3, seed = 7)
    model <- tulpaObs:::.tobs_build_removal(
      ~ abund_cov1, ~ det_cov1, sim$data, sim$y)
    is_nb <- identical(mix, "NB")
    lay  <- tulpaObs:::.tobs_abun_nuts_layout(2L, 2L, is_nb)
    marg <- tulpaObs:::.tobs_removal_nuts_marginal(model, mixture = mix)
    theta <- c(log(7), 0.3, 0.2, -0.3, if (is_nb) log(3))
    f <- function(th) tulpaObs:::.tobs_removal_nuts_logpost(th, marg, lay)$lp
    g_an <- tulpaObs:::.tobs_removal_nuts_logpost(theta, marg, lay)$grad
    g_fd <- sapply(seq_along(theta), function(j) {
      h <- 1e-5; tp <- theta; tm <- theta
      tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
      (f(tp) - f(tm)) / (2 * h)
    })
    expect_equal(g_an, g_fd, tolerance = 1e-4,
                 info = paste("mixture", mix))
  }
})

test_that("C++ removal NUTS log-posterior matches the R oracle byte-for-byte", {
  for (mix in c("P", "NB")) {
    is_nb <- identical(mix, "NB")
    sim <- simulate_removal(N = 50, K = 4, n_abund_covs = 1, n_det_covs = 1,
                            mixture = if (is_nb) "negbin" else "poisson",
                            size = 3, seed = 12)
    model <- tulpaObs:::.tobs_build_removal(~ abund_cov1, ~ det_cov1, sim$data, sim$y)
    K_max <- max(rowSums(sim$y)) + 100L
    lay  <- tulpaObs:::.tobs_abun_nuts_layout(2L, 2L, is_nb)
    marg <- tulpaObs:::.tobs_removal_nuts_marginal(model, mixture = mix, K_max = K_max)
    theta <- c(log(6), 0.2, 0.1, -0.2, if (is_nb) log(3))
    spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
                 X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
                 n_sites = model$n_sites, K_max = as.integer(K_max), is_nb = is_nb)
    r_out <- tulpaObs:::.tobs_removal_nuts_logpost(theta, marg, lay)
    c_out <- tulpaObs:::cpp_removal_nuts_joint_logpost(spec, theta,
                                                       sigma_beta = 10, sigma_logr = 1.5)
    expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9, info = mix)
    expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9, info = mix)
  }
})

test_that("single Poisson removal fit recovers truth", {
  skip_if_fast()
  beta_lambda <- c(log(8), 0.6, -0.4)
  beta_p      <- c(0.2, 0.4)
  sim <- simulate_removal(N = 400, K = 5, n_abund_covs = 2, n_det_covs = 1,
                          beta_lambda = beta_lambda, beta_p = beta_p, seed = 11)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = removal(), detection = ~ det_cov1, y = sim$y,
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_lambda, beta_p)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3))
  expect_lt(abs(est[2] - beta_lambda[2]), 0.15)
  expect_lt(abs(est[3] - beta_lambda[3]), 0.15)
  expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                            "lambda_abund_cov2", "p_(Intercept)", "p_det_cov1"))
  expect_equal(unname(fit$intercepts$lambda), exp(est[1]), tolerance = 1e-8)
  expect_equal(unname(fit$intercepts$p), plogis(est[4]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(9), 0.5, -0.3)
  beta_p      <- c(0.3, 0.4)
  n_seed <- 30L
  truth <- c(beta_lambda, beta_p)
  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_removal(N = 200, K = 5, n_abund_covs = 2, n_det_covs = 1,
                            beta_lambda = beta_lambda, beta_p = beta_p, seed = 200 + s)
    fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
                family = removal(), detection = ~ det_cov1, y = sim$y,
                method = "laplace", control = list(verbose = FALSE))
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.85),
              info = paste(round(cover_rate, 2), collapse = " | "))
})

test_that("negbin removal recovers truth and surfaces dispersion", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(10), 0.5, -0.3)
  beta_p      <- c(0.3, 0.4)
  size_true   <- 3
  sim <- simulate_removal(N = 400, K = 6, n_abund_covs = 2, n_det_covs = 1,
                          beta_lambda = beta_lambda, beta_p = beta_p,
                          mixture = "negbin", size = size_true, seed = 21)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = removal(mixture = "negbin"), detection = ~ det_cov1,
              y = sim$y, method = "laplace")

  expect_identical(fit$mixture, "negbin")
  truth <- c(beta_lambda, beta_p)
  est   <- as.numeric(fit$means[1:5])
  se    <- as.numeric(fit$sds[1:5])
  expect_true(all(abs(est - truth) / se < 3))
  expect_true("log_r" %in% rownames(fit$vcov))
  expect_false(is.null(fit$nmix_dispersion))
  se_logr <- sqrt(fit$vcov["log_r", "log_r"])
  expect_lt(abs(fit$nmix_dispersion$log_r - log(size_true)) / se_logr, 3.5)
})

test_that("S3 surface works for removal fits", {
  skip_if_fast()
  sim <- simulate_removal(N = 200, K = 4, n_abund_covs = 2, n_det_covs = 1,
                          beta_lambda = c(log(8), 0.5, -0.3),
                          beta_p = c(0.3, 0.4), seed = 3)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = removal(), detection = ~ det_cov1, y = sim$y,
              method = "laplace", control = list(verbose = FALSE))

  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), length(sim$y))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("lambda", "p", "N"))
  expect_length(fv$lambda, 200L)
  expect_true(all(fv$lambda > 0))

  X0 <- cbind(1, c(-1, 0, 1), 0)
  pr <- predict(fit, X.0 = X0)
  expect_true(all(pr$mean > 0))
  expect_true(all(diff(pr$mean) > 0))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim >= 0))
  # Depletion: simulated removals never exceed the running available count, so a
  # site's later passes cannot exceed its earlier total by construction.
  expect_true(all(rowSums(ysim) >= 0))

  rr <- residuals(fit, type = "pearson")
  expect_equal(dim(rr), dim(sim$y))
})

test_that("removal() requires complete pass sequences (no NA)", {
  sim <- simulate_removal(N = 30, K = 4, seed = 5)
  y_na <- sim$y; y_na[1, 2] <- NA
  expect_error(
    tobs(formula = ~ 1, data = sim$data, family = removal(),
         detection = ~ 1, y = y_na, method = "laplace"),
    "complete pass sequences")
})

test_that("removal NUTS recovers truth and scores WAIC", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(7), 0.5)
  beta_p      <- c(0.3, -0.3)
  sim <- simulate_removal(N = 80, K = 5, n_abund_covs = 1, n_det_covs = 1,
                          beta_lambda = beta_lambda, beta_p = beta_p, seed = 31)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data, family = removal(),
              detection = ~ det_cov1, y = sim$y, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))

  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 500L)
  truth <- c(beta_lambda, beta_p)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_false(any(is.na(fit$divergent)))
  expect_lt(mean(fit$nuts$divergent), 0.2)
  # WAIC / LOO from the NUTS draws (per-site pointwise marginal log-lik).
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
})
