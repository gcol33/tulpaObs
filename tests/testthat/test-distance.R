# Binned distance sampling (half-normal / hazard-rate key, line / point
# transect), Poisson + negbin, non-spatial Laplace + NUTS (gcol33/tulpaObs#38).
#
# Recovery-grade tests (per the "statistical code needs recovery tests" rule):
# point recovery against simulated truth + 95% CI coverage across seeds, plus a
# closed-form correctness anchor (the Poisson distance marginal equals
# independent per-bin Poissons by Poisson thinning), an FD check of the analytic
# gradient, an analytic-observed-information vs FD-Hessian check (the Louis
# curvature, including the second-derivative bin quadrature), and a byte-level
# C++ <-> R oracle cross-check. Structural tests cover the family wiring / S3.

cuts5 <- seq(0, 1, length.out = 6)   # five distance bins out to 1

test_that("distance() family is wired and reports its supported methods", {
  f <- distance(cutpoints = cuts5)
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "distance")
  expect_identical(f$status, "working")
  expect_identical(f$params$key, "halfnorm")
  expect_identical(f$params$transect, "line")
  expect_true(all(c("laplace", "nuts") %in% tulpaObs:::.tobs_family_methods$distance))
})

test_that("distance marginal equals independent per-bin Poissons (Poisson abundance)", {
  # Poisson thinning: under N ~ Poisson(lambda), the multinomial split into the
  # detection cells (pi_1, ..., pi_B, 1 - p) gives INDEPENDENT bin counts
  # y_b ~ Poisson(lambda * pi_b). The kernel's marginal sum over latent N must
  # reproduce this closed form exactly, and E[N|y] = R + lambda * (1 - p).
  set.seed(1)
  for (rep in 1:5) {
    for (trans in c("line", "point")) {
      lambda <- runif(1, 5, 40)
      sigma  <- runif(1, 0.25, 0.7)
      pi_b   <- tulpaObs:::.distance_pi(sigma, cuts5, "halfnorm", trans)
      y      <- rpois(length(pi_b), lambda * pi_b)
      out <- tulpaObs:::cpp_distance_total_log_lik(
        matrix(as.integer(y), 1L), log(lambda), log(sigma), 0,
        cuts5, if (trans == "point") 1L else 0L, 0L,
        K_max = sum(y) + 400L, r = Inf, quad_order = 64L)
      ll_ref <- sum(dpois(y, lambda * pi_b, log = TRUE))
      expect_equal(out$log_lik, ll_ref, tolerance = 1e-7,
                   info = paste(trans, rep))
      expect_equal(out$mean_N, sum(y) + lambda * (1 - sum(pi_b)),
                   tolerance = 1e-5, info = paste(trans, rep))
      expect_equal(out$p_det, sum(pi_b), tolerance = 1e-7)
    }
  }
})

test_that("analytic gradient matches finite differences (halfnorm, hazard, NB)", {
  cases <- list(
    list(key = "halfnorm", mix = "P"),
    list(key = "hazard",   mix = "P"),
    list(key = "halfnorm", mix = "NB"))
  for (cs in cases) {
    is_nb <- identical(cs$mix, "NB"); hazard <- identical(cs$key, "hazard")
    sim <- simulate_distance(N = 80, key = cs$key, transect = "line",
                             beta_lambda = c(log(30), 0.3),
                             beta_sigma  = c(log(0.45), 0.2), shape = 3,
                             mixture = if (is_nb) "negbin" else "poisson",
                             size = 5, seed = 7)
    model <- tulpaObs:::.tobs_build_distance(~ abund_cov1, ~ sigma_cov1, sim$data,
                                             sim$y, cutpoints = sim$cutpoints,
                                             key = cs$key, mixture = cs$mix)
    lay  <- tulpaObs:::.tobs_distance_nuts_layout(2L, 2L, hazard, is_nb)
    marg <- tulpaObs:::.tobs_distance_nuts_marginal(model, mixture = cs$mix)
    theta <- c(log(28), 0.25, log(0.45), 0.18,
               if (hazard) log(3), if (is_nb) log(5))
    f <- function(th) tulpaObs:::.tobs_distance_nuts_logpost(th, marg, lay)$lp
    g_an <- tulpaObs:::.tobs_distance_nuts_logpost(theta, marg, lay)$grad
    g_fd <- sapply(seq_along(theta), function(j) {
      h <- 1e-5; tp <- theta; tm <- theta
      tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
      (f(tp) - f(tm)) / (2 * h)
    })
    expect_equal(g_an, g_fd, tolerance = 1e-4, info = paste(cs$key, cs$mix))
  }
})

test_that("analytic observed information matches the FD Hessian (Louis curvature)", {
  # The Laplace vcov is the inverse of the analytic marginal observed information
  # (E[I_c|y] minus the rank-1 Var[N|y] coupling, with the second-derivative bin
  # quadrature). It must equal minus the finite-difference Jacobian of the exact
  # marginal gradient at the mode.
  for (key in c("halfnorm", "hazard")) {
    hazard <- identical(key, "hazard")
    sim <- simulate_distance(N = 120, key = key, transect = "line",
                             beta_lambda = c(log(35), 0.3),
                             beta_sigma  = c(log(0.45), 0.2), shape = 3, seed = 5)
    Xl <- model.matrix(~ abund_cov1, sim$data)
    Xs <- model.matrix(~ sigma_cov1, sim$data)
    raw <- tulpaObs:::distance_laplace(sim$y, Xl, Xs, sim$cutpoints, key = key,
                                       transect = "line", mixture = "P",
                                       verbose = FALSE)
    theta <- c(raw$beta_lambda, raw$beta_sigma, if (hazard) raw$eta_b)
    pl <- length(raw$beta_lambda); ps <- length(raw$beta_sigma)
    grad_at <- function(th) {
      bl <- th[seq_len(pl)]; bs <- th[pl + seq_len(ps)]
      eb <- if (hazard) th[pl + ps + 1L] else 0
      o <- tulpaObs:::cpp_distance_total_log_lik(
        sim$y, as.numeric(Xl %*% bl), as.numeric(Xs %*% bs), eb,
        sim$cutpoints, 0L, if (hazard) 1L else 0L, raw$K_max, Inf, 64L)
      g <- c(as.numeric(crossprod(Xl, o$grad_eta_lambda)),
             as.numeric(crossprod(Xs, o$grad_eta_sigma)))
      if (hazard) g <- c(g, o$grad_eta_b)
      g
    }
    p <- length(theta); Hfd <- matrix(0, p, p)
    for (i in seq_len(p)) {
      th1 <- theta; th2 <- theta; th1[i] <- th1[i] + 1e-5; th2[i] <- th2[i] - 1e-5
      Hfd[, i] <- -(grad_at(th1) - grad_at(th2)) / 2e-5
    }
    Hobs <- as.matrix(raw$H_obs)
    expect_lt(max(abs(Hobs - Hfd)) / max(abs(Hfd)), 1e-5, label = key)
  }
})

test_that("C++ distance NUTS log-posterior matches the R oracle byte-for-byte", {
  for (cs in list(list(key = "halfnorm", mix = "P"),
                  list(key = "hazard",   mix = "P"),
                  list(key = "halfnorm", mix = "NB"))) {
    is_nb <- identical(cs$mix, "NB"); hazard <- identical(cs$key, "hazard")
    sim <- simulate_distance(N = 60, key = cs$key, transect = "point",
                             mixture = if (is_nb) "negbin" else "poisson",
                             size = 4, seed = 12)
    model <- tulpaObs:::.tobs_build_distance(~ abund_cov1, ~ sigma_cov1, sim$data,
                                             sim$y, cutpoints = sim$cutpoints,
                                             key = cs$key, transect = "point",
                                             mixture = cs$mix)
    lay  <- tulpaObs:::.tobs_distance_nuts_layout(2L, 2L, hazard, is_nb)
    marg <- tulpaObs:::.tobs_distance_nuts_marginal(model, mixture = cs$mix)
    theta <- c(log(30), 0.2, log(0.45), 0.15,
               if (hazard) log(3), if (is_nb) log(4))
    spec <- list(y = model$y, X_lambda = model$X_processes[[1]],
                 X_sigma = model$X_processes[[2]], cutpoints = model$cutpoints,
                 transect = 1L, key = if (hazard) 1L else 0L,
                 K_max = marg$K_max, is_nb = is_nb,
                 quad_order = model$quad_order)
    r_out <- tulpaObs:::.tobs_distance_nuts_logpost(theta, marg, lay)
    c_out <- tulpaObs:::cpp_distance_nuts_joint_logpost(spec, theta, 10, 1.5, 1.5)
    expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9, info = paste(cs$key, cs$mix))
    expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9,
                 info = paste(cs$key, cs$mix))
  }
})

test_that("half-normal distance fit recovers truth", {
  skip_if_fast()
  beta_lambda <- c(log(50), 0.5, -0.3)
  beta_sigma  <- c(log(0.45), 0.25)
  sim <- simulate_distance(N = 400, cutpoints = cuts5, key = "halfnorm",
                           transect = "line", n_abund_covs = 2, n_sigma_covs = 1,
                           beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                           seed = 11)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = distance(key = "halfnorm", transect = "line",
                                cutpoints = sim$cutpoints),
              detection = ~ sigma_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_lambda, beta_sigma)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3))
  expect_lt(abs(est[2] - beta_lambda[2]), 0.12)
  expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                            "lambda_abund_cov2", "sigma_(Intercept)", "sigma_sigma_cov1"))
  expect_equal(unname(fit$intercepts$lambda), exp(est[1]), tolerance = 1e-8)
  expect_equal(unname(fit$intercepts$sigma), exp(est[4]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(45), 0.4, -0.3)
  beta_sigma  <- c(log(0.45), 0.2)
  n_seed <- 30L
  truth <- c(beta_lambda, beta_sigma)
  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_distance(N = 250, cutpoints = cuts5, key = "halfnorm",
                             transect = "line", n_abund_covs = 2, n_sigma_covs = 1,
                             beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                             seed = 300 + s)
    fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
                family = distance(key = "halfnorm", transect = "line",
                                  cutpoints = sim$cutpoints),
                detection = ~ sigma_cov1, y = sim$y, method = "laplace",
                control = list(verbose = FALSE))
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.85),
              info = paste(round(cover_rate, 2), collapse = " | "))
})

test_that("hazard-rate distance recovers truth and the shape", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(55), 0.4)
  beta_sigma  <- c(log(0.5), 0.2)
  shape_true  <- 3
  sim <- simulate_distance(N = 400, cutpoints = cuts5, key = "hazard",
                           transect = "line", n_abund_covs = 1, n_sigma_covs = 1,
                           beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                           shape = shape_true, seed = 17)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data,
              family = distance(key = "hazard", transect = "line",
                                cutpoints = sim$cutpoints),
              detection = ~ sigma_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  truth <- c(beta_lambda, beta_sigma)
  est   <- as.numeric(fit$means[1:4])
  se    <- as.numeric(fit$sds[1:4])
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_true("log_shape" %in% rownames(fit$vcov))
  se_shape <- sqrt(fit$vcov["log_shape", "log_shape"])
  expect_lt(abs(fit$distance_shape$eta_b - log(shape_true)) / se_shape, 3.5)
})

test_that("negbin distance recovers truth and surfaces dispersion", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(50), 0.5)
  beta_sigma  <- c(log(0.45), 0.2)
  size_true   <- 6
  sim <- simulate_distance(N = 400, cutpoints = cuts5, key = "halfnorm",
                           transect = "line", n_abund_covs = 1, n_sigma_covs = 1,
                           beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                           mixture = "negbin", size = size_true, seed = 21)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data,
              family = distance(key = "halfnorm", transect = "line",
                                cutpoints = sim$cutpoints, mixture = "negbin"),
              detection = ~ sigma_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_identical(fit$mixture, "negbin")
  truth <- c(beta_lambda, beta_sigma)
  est   <- as.numeric(fit$means[1:4])
  se    <- as.numeric(fit$sds[1:4])
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_true("log_r" %in% rownames(fit$vcov))
  expect_false(is.null(fit$nmix_dispersion))
  se_logr <- sqrt(fit$vcov["log_r", "log_r"])
  expect_lt(abs(fit$nmix_dispersion$log_r - log(size_true)) / se_logr, 3.5)
})

test_that("point-transect distance recovers truth", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(60), 0.4)
  beta_sigma  <- c(log(0.5), 0.2)
  sim <- simulate_distance(N = 400, cutpoints = cuts5, key = "halfnorm",
                           transect = "point", n_abund_covs = 1, n_sigma_covs = 1,
                           beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                           seed = 23)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data,
              family = distance(key = "halfnorm", transect = "point",
                                cutpoints = sim$cutpoints),
              detection = ~ sigma_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  truth <- c(beta_lambda, beta_sigma)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
})

test_that("S3 surface works for distance fits", {
  skip_if_fast()
  sim <- simulate_distance(N = 200, cutpoints = cuts5, key = "halfnorm",
                           transect = "line", n_abund_covs = 2, n_sigma_covs = 1,
                           beta_lambda = c(log(45), 0.5, -0.3),
                           beta_sigma = c(log(0.45), 0.2), seed = 3)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = distance(key = "halfnorm", transect = "line",
                                cutpoints = sim$cutpoints),
              detection = ~ sigma_cov1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))

  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), length(sim$y))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("lambda", "sigma", "p"))
  expect_length(fv$lambda, 200L)
  expect_true(all(fv$lambda > 0))
  expect_true(all(fv$p > 0 & fv$p <= 1))

  X0 <- cbind(1, c(-1, 0, 1), 0)
  pr <- predict(fit, X.0 = X0, type = "lambda")
  expect_length(pr, 3L)
  expect_true(all(pr > 0) && all(diff(pr) > 0))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim >= 0))

  rr <- residuals(fit, type = "pearson")
  expect_equal(dim(rr), dim(sim$y))
})

test_that("distance() validates cutpoints and rejects NA", {
  sim <- simulate_distance(N = 30, cutpoints = cuts5, seed = 5)
  expect_error(
    tobs(formula = ~ 1, data = sim$data, family = distance(),
         detection = ~ 1, y = sim$y, method = "laplace"),
    "cutpoints")
  y_na <- sim$y; y_na[1, 2] <- NA
  expect_error(
    tobs(formula = ~ 1, data = sim$data,
         family = distance(cutpoints = cuts5), detection = ~ 1, y = y_na,
         method = "laplace"),
    "must not contain NA")
})

test_that("distance NUTS recovers truth and scores WAIC", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(45), 0.4)
  beta_sigma  <- c(log(0.45), 0.2)
  sim <- simulate_distance(N = 120, cutpoints = cuts5, key = "halfnorm",
                           transect = "line", n_abund_covs = 1, n_sigma_covs = 1,
                           beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                           seed = 31)
  fit <- tobs(formula = ~ abund_cov1, data = sim$data,
              family = distance(key = "halfnorm", transect = "line",
                                cutpoints = sim$cutpoints),
              detection = ~ sigma_cov1, y = sim$y, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))

  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 500L)
  truth <- c(beta_lambda, beta_sigma)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_false(any(is.na(fit$divergent)))
  expect_lt(mean(fit$nuts$divergent), 0.2)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
})


# --- NUTS + random effect (tulpaObs#51) ------------------------------------

# Line-transect half-normal distance data with a per-site intercept RE on the
# abundance arm. lambda_i = exp(b0 + b1 x_i + b[group_i]); each individual is
# detected with prob exp(-d^2 / (2 sigma^2)) and binned by its distance.
sim_distance_lambda_re <- function(N, ngrp, cutpoints, beta_lambda, sigma,
                                   sigma_b, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- max(cutpoints); n_bins <- length(cutpoints) - 1L
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  x   <- stats::rnorm(N)
  data <- data.frame(x1 = x, g = factor(grp))
  lambda <- exp(beta_lambda[1] + beta_lambda[2] * x + b[grp])
  y <- matrix(0L, N, n_bins)
  for (i in seq_len(N)) {
    Ni <- stats::rpois(1L, lambda[i])
    if (Ni == 0L) next
    d <- stats::runif(Ni, 0, W)
    det <- stats::runif(Ni) < exp(-d^2 / (2 * sigma^2))
    if (any(det)) {
      bins <- findInterval(d[det], cutpoints, rightmost.closed = TRUE)
      bins <- bins[bins >= 1L & bins <= n_bins]
      for (bb in bins) y[i, bb] <- y[i, bb] + 1L
    }
  }
  list(y = y, data = data, sigma_b = sigma_b, beta_lambda = beta_lambda,
       sigma = sigma)
}

test_that("distance() NUTS samples a single abundance RE and recovers sigma + betas", {
  skip_on_cran()
  skip_if_fast()
  cuts <- seq(0, 1, length.out = 5)
  s <- sim_distance_lambda_re(N = 80, ngrp = 8, cutpoints = cuts,
                              beta_lambda = c(log(30), 0.3), sigma = 0.4,
                              sigma_b = 0.5, seed = 9)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1, data = s$data, y = s$y,
              family = distance(cutpoints = cuts, key = "halfnorm",
                                transect = "line", K_max = 120L),
              method = "nuts", verbose = FALSE,
              control = list(n.iter = 400L, n.warmup = 300L, seed = 1L))
  expect_identical(fit$method, "nuts")
  expect_identical(fit$re$arm, "lambda")
  expect_equal(fit$re$n_groups, 8L)
  expect_lt(abs(fit$re$sigma - 0.5), 0.4)
  expect_gt(fit$re$sigma_sd, 0)
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(30)), 0.4)
  expect_lt(abs(fit$means[["lambda_x1"]] - 0.3), 0.25)
  expect_lt(mean(fit$nuts$divergent), 0.15)
})

test_that("distance() Laplace AGHQ recovers a site-grouped abundance RE", {
  skip_on_cran()
  skip_if_fast()
  cuts <- seq(0, 1, length.out = 5)
  s <- sim_distance_lambda_re(N = 120, ngrp = 12, cutpoints = cuts,
                              beta_lambda = c(log(40), 0.3), sigma = 0.4,
                              sigma_b = 0.6, seed = 3)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1, data = s$data, y = s$y,
              family = distance(cutpoints = cuts, key = "halfnorm",
                                transect = "line", K_max = 250L),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = 5L))
  expect_identical(fit$method, "laplace")
  expect_identical(fit$nmix_re$arm, "lambda")
  expect_true("sigma_g1_(Intercept)" %in% names(fit$means))
  # AGHQ-debiased variance component near truth (some small-cluster attenuation
  # remains at this group size); fixed effects recover tightly.
  expect_lt(abs(fit$means[["sigma_g1_(Intercept)"]] - 0.6), 0.25)
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(40)), 0.3)
  expect_lt(abs(fit$means[["lambda_x1"]] - 0.3), 0.2)
  # Per-group BLUPs surface through ranef() (one row per group level).
  re <- ranef(fit)
  expect_true(is.data.frame(re))
  expect_equal(nrow(re), 12L)
})

test_that("distance() Laplace RE rejects the hazard key, detection arm, both arms", {
  cuts <- seq(0, 1, length.out = 5)
  s <- sim_distance_lambda_re(N = 30, ngrp = 5, cutpoints = cuts,
                              beta_lambda = c(log(20), 0.2), sigma = 0.4,
                              sigma_b = 0.4, seed = 2)
  # Hazard-rate key carries a global shape coordinate the count-family theta
  # layout cannot express, so the grouped-RE path is half-normal only.
  expect_error(
    tobs(formula = ~ x1 + (1 | g), detection = ~ 1, data = s$data, y = s$y,
         family = distance(cutpoints = cuts, key = "hazard"), method = "laplace"),
    "half-normal")
  # A detection-scale (sigma-arm) RE couples a site's bins through the latent N.
  expect_error(
    tobs(formula = ~ x1, detection = ~ (1 | g), data = s$data, y = s$y,
         family = distance(cutpoints = cuts), method = "laplace"),
    "abundance arm only")
  # Both arms at once is not yet supported.
  expect_error(
    tobs(formula = ~ x1 + (1 | g), detection = ~ (1 | g), data = s$data, y = s$y,
         family = distance(cutpoints = cuts), method = "laplace"),
    "abundance arm only|BOTH")
  # Detection-arm RE under NUTS is likewise abundance-arm only.
  expect_error(
    tobs(formula = ~ x1, detection = ~ (1 | g), data = s$data, y = s$y,
         family = distance(cutpoints = cuts), method = "nuts",
         control = list(n.iter = 20L, n.warmup = 10L)),
    "abundance arm only|single intercept")
})
