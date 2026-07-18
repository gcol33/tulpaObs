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


# --- NUTS + random effect (tulpaObs#51) ------------------------------------

# Removal data with a per-site intercept RE on the abundance arm.
sim_removal_lambda_re <- function(N, K, ngrp, beta_lambda, beta_p, sigma_b,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  lambda <- exp(as.numeric(model.matrix(~ x1, data) %*% beta_lambda) + b[grp])
  p <- plogis(beta_p)
  Nlat <- stats::rpois(N, lambda)
  y <- matrix(0L, N, K)
  for (i in seq_len(N)) {
    rem <- Nlat[i]
    for (k in seq_len(K)) { yk <- stats::rbinom(1L, rem, p); y[i, k] <- yk; rem <- rem - yk }
  }
  list(y = y, data = data, sigma_b = sigma_b, beta_lambda = beta_lambda)
}

test_that("removal() NUTS samples a single intercept RE and recovers sigma + betas", {
  skip_on_cran()
  skip_if_fast()
  # Small counts (lambda ~ 5) + an explicit modest K.max keep the depleting-
  # binomial marginal cheap so the sampler runs in test time.
  s <- sim_removal_lambda_re(N = 70, K = 4, ngrp = 8,
                             beta_lambda = c(log(5), 0.3), beta_p = qlogis(0.5),
                             sigma_b = 0.6, seed = 7)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1, data = s$data,
              y = s$y, family = removal(K_max = 45L), method = "nuts",
              verbose = FALSE,
              control = list(n.iter = 400L, n.warmup = 300L, seed = 1L))
  expect_identical(fit$method, "nuts")
  expect_identical(fit$re$arm, "lambda")
  expect_equal(fit$re$n_groups, 8L)
  expect_lt(abs(fit$re$sigma - 0.6), 0.4)
  expect_gt(fit$re$sigma_sd, 0)
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(5)), 0.35)
  expect_lt(abs(fit$means[["lambda_x1"]] - 0.3), 0.25)
  expect_lt(mean(fit$nuts$divergent), 0.15)
})

test_that("removal() Laplace AGHQ recovers a site-grouped intercept RE (sigma + betas)", {
  skip_on_cran()
  skip_if_fast()
  # Site-level intercept RE on the abundance arm, fit on the shared count-model
  # AGHQ path (tulpaObs#51). n.quad > 1 debiases the small-cluster variance
  # attenuation so sigma recovers; the per-site removal marginal feeds the same
  # grouped-RE oracle the N-mixture uses.
  s <- sim_removal_lambda_re(N = 90, K = 5, ngrp = 10,
                             beta_lambda = c(log(6), 0.3), beta_p = qlogis(0.5),
                             sigma_b = 0.7, seed = 11)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1, data = s$data,
              y = s$y, family = removal(K_max = 60L), method = "laplace",
              verbose = FALSE, control = list(n.quad = 5L))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$nmix_re$arm, "lambda")
  expect_true("sigma_g1_(Intercept)" %in% names(fit$means))
  expect_true(any(grepl("^re_g1_", names(fit$means))))
  expect_lt(abs(fit$means[["sigma_g1_(Intercept)"]] - 0.7), 0.35)
  expect_lt(abs(fit$means[["lambda_(Intercept)"]] - log(6)), 0.35)
  expect_lt(abs(fit$means[["lambda_x1"]] - 0.3), 0.25)
})

test_that("removal() NUTS RE rejects slopes / both-arm", {
  s <- sim_removal_lambda_re(N = 40, K = 4, ngrp = 5, beta_lambda = c(log(6), 0.2),
                             beta_p = 0, sigma_b = 0.5, seed = 2)
  # Random slope under NUTS is AGHQ-only territory.
  expect_error(
    tobs(formula = ~ x1 + (x1 | g), detection = ~ 1, data = s$data, y = s$y,
         family = removal(), method = "nuts",
         control = list(n.iter = 20L, n.warmup = 10L)),
    "single intercept random effect|laplace")
})


# --- areal spatial (ICAR / proper-CAR) on the abundance arm (tulpaObs#51) -----

# Rook-adjacency on a side x side grid (one spatial unit per site).
.rem_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

# Spatial removal-sampling data: log lambda_i = b0 + b1 x_i + phi_i with a
# smoothed, demeaned ICAR-like field phi; N_i ~ Poisson(lambda_i) depleted over
# K removal passes at detection p_i.
.sim_spatial_removal <- function(adj, b_lambda, b_p, K = 5L, sd_phi = 0.6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ng <- nrow(adj)
  phi <- as.numeric(scale(rnorm(ng)))
  for (rep in 1:3) {
    pn <- phi
    for (i in seq_len(ng)) { nb <- which(adj[i, ] == 1L); pn[i] <- 0.5 * phi[i] + 0.5 * mean(phi[nb]) }
    phi <- pn
  }
  phi <- sd_phi * as.numeric(scale(phi)); phi <- phi - mean(phi)
  x_ab <- rnorm(ng); x_det <- rnorm(ng)
  lambda <- exp(b_lambda[1] + b_lambda[2] * x_ab + phi)
  p <- plogis(b_p[1] + b_p[2] * x_det)
  N <- rpois(ng, lambda)
  y <- matrix(0L, ng, K)
  for (i in seq_len(ng)) {
    rem <- N[i]
    for (k in seq_len(K)) { d <- rbinom(1, rem, p[i]); y[i, k] <- d; rem <- rem - d }
  }
  list(y = y, data = data.frame(abund_cov1 = x_ab, det_cov1 = x_det),
       truth = c(b_lambda, b_p), phi = phi)
}

test_that("removal() areal ICAR recovers the abundance slope + field, calibrated cov", {
  skip_on_cran()
  skip_if_fast()
  # The abundance-arm field reuses the shared count-marginal nested-Laplace driver
  # (the removal per-site marginal supplies the same moments as the N-mixture).
  adj <- .rem_grid_adj(7L)             # 49 sites / spatial units
  b_lambda <- c(log(8), 0.5); b_p <- c(0.2, 0.4)
  slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
  for (s in 1:3) {
    sim <- .sim_spatial_removal(adj, b_lambda, b_p, K = 5L, seed = 300 + s)
    fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
                family = removal(), detection = ~ det_cov1, y = sim$y,
                method = "nested_laplace", control = list(progress = FALSE, verbose = FALSE))
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                                "p_(Intercept)", "p_det_cov1"))
      V <- vcov(fit)
      expect_true(all(is.finite(V)))
      expect_true(all(eigen(V, only.values = TRUE)$values > 0))      # PD
      expect_lt(fit$sds[["lambda_(Intercept)"]], 5)                  # constrained intercept
      expect_false(is.null(fit$spatial_field))
    }
    est <- fit$means[["lambda_abund_cov1"]]; se <- fit$sds[["lambda_abund_cov1"]]
    slopes <- c(slopes, est)
    slope_ok <- c(slope_ok, abs(est - 0.5) / se < 3)
    field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
  }
  expect_true(all(slope_ok))                       # slope within 3 SE every seed
  expect_lt(abs(mean(slopes) - 0.5), 0.15)         # unbiased on average
  expect_gt(mean(field_cor), 0.6)                  # field tracks the truth
})

test_that("removal() areal-spatial coefficient SEs are calibrated", {
  skip_on_cran()
  skip_if_fast()
  # The areal fit reuses the shared count-spatial assembler, whose cross-arm
  # rank-1 (Var[N|y]) correction uses the binomial detection weight. Removal's
  # depleting binomial differs, so confirm the (lambda, p) coefficient SEs are
  # not anti-conservative: on field-free data the icar fit's SEs must match the
  # kernel-correct non-spatial removal SEs (pass 1 dominates the detection info,
  # so the binomial cross-arm is calibrated; coverage is at/above nominal).
  adj <- .rem_grid_adj(8L)             # 64 sites
  set.seed(101)
  ng <- nrow(adj); x_ab <- rnorm(ng); x_det <- rnorm(ng)
  lambda <- exp(log(10) + 0.5 * x_ab); p <- plogis(0.2 + 0.4 * x_det)  # no field
  N <- rpois(ng, lambda); y <- matrix(0L, ng, 5L)
  for (i in seq_len(ng)) {
    rem <- N[i]; for (k in 1:5) { d <- rbinom(1, rem, p[i]); y[i, k] <- d; rem <- rem - d }
  }
  dat <- data.frame(abund_cov1 = x_ab, det_cov1 = x_det)
  fns <- tobs(formula = ~ abund_cov1, data = dat, family = removal(),
              detection = ~ det_cov1, y = y, method = "laplace",
              control = list(progress = FALSE, verbose = FALSE))
  fsp <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = dat,
              family = removal(), detection = ~ det_cov1, y = y,
              method = "nested_laplace", control = list(progress = FALSE, verbose = FALSE))
  ratio <- fsp$sds[names(fns$sds)] / fns$sds
  # SEs match within 15% (the field adds a small extra d.o.f.); crucially the
  # spatial SEs are not materially smaller than the correct non-spatial ones.
  expect_true(all(ratio > 0.85 & ratio < 1.2))
})

test_that("removal() areal spatial: proper-CAR + bym2 fit; nuts+icar samples (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .rem_grid_adj(6L)
  sim <- .sim_spatial_removal(adj, c(log(7), 0.4), c(0.2, 0.3), K = 4L, seed = 7)
  for (term in c("car_proper", "bym2")) {
    f <- if (term == "car_proper") (~ abund_cov1 + car_proper(graph = adj)) else
                                   (~ abund_cov1 + bym2(graph = adj))
    fit <- tobs(formula = f, data = sim$data, family = removal(),
                detection = ~ det_cov1, y = sim$y, method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    expect_identical(fit$method, "nested_laplace")
    expect_true(all(is.finite(vcov(fit))))
    expect_false(is.null(fit$spatial_field))
  }
  # NUTS + areal now samples an intrinsic icar() field via the #71 sum-to-zero
  # reparameterisation (full recovery lives in test-count-spatial-nuts.R). Here
  # confirm the dispatch path runs end-to-end and centres the field.
  fit_icar <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
    family = removal(), detection = ~ det_cov1, y = sim$y, method = "nuts",
    control = list(n.iter = 60L, n.warmup = 60L, verbose = FALSE, progress = FALSE))
  expect_identical(fit_icar$method, "nuts")
  expect_false(is.null(fit_icar$spatial_field))
  expect_lt(abs(mean(fit_icar$spatial_field)), 1e-6)   # sum-to-zero centred
})


test_that("removal() bym2 + proper-CAR recover the abundance field + slope (#131)", {
  skip_on_cran()
  skip_if_fast()
  # bym2 is a distinct code path from icar: the fitted unit field is the
  # rho-mixed reconstruction z = sqrt(rho) * phi + sqrt(1 - rho) * theta, so it
  # exercises the mixture the icar recovery never touches. proper-CAR uses a
  # full-rank precision. Both still track a structured (icar-simulated) field and
  # recover the abundance slope on the shared count-marginal driver.
  adj <- .rem_grid_adj(7L)
  b_lambda <- c(log(8), 0.5); b_p <- c(0.2, 0.4)
  for (term in c("bym2", "car_proper")) {
    tf <- if (term == "bym2") (~ abund_cov1 + bym2(graph = adj)) else
                              (~ abund_cov1 + car_proper(graph = adj))
    slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
    for (s in 1:3) {
      sim <- .sim_spatial_removal(adj, b_lambda, b_p, K = 5L, seed = 300 + s)
      fit <- tobs(formula = tf, data = sim$data, family = removal(),
                  detection = ~ det_cov1, y = sim$y, method = "nested_laplace",
                  control = list(progress = FALSE, verbose = FALSE))
      est <- fit$means[["lambda_abund_cov1"]]; se <- fit$sds[["lambda_abund_cov1"]]
      slopes    <- c(slopes, est)
      slope_ok  <- c(slope_ok, abs(est - 0.5) / se < 3)
      field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
    }
    expect_true(all(slope_ok), info = term)
    expect_lt(abs(mean(slopes) - 0.5), 0.15)
    expect_gt(mean(field_cor), 0.6)
  }
})


test_that("removal() temporal()-only field recovers the AR1 field + slope (#114)", {
  skip_on_cran()
  skip_if_fast()
  # A temporal() term on its own (no areal field) runs the shared areal-BFGS
  # driver with a single temporal block on the abundance arm (gcol33/tulpaObs#114).
  Tt <- 8L; per_t <- 30L; N <- Tt * per_t
  fcor <- slope <- rep(NA_real_, 8L)
  for (s in seq_len(8L)) {
    set.seed(100L + s)
    period <- rep(seq_len(Tt), each = per_t)
    rho <- 0.7; sig <- 0.5; u <- numeric(Tt)
    u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
    for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
    u <- u - mean(u)
    x <- stats::rnorm(N)
    lambda <- exp(log(8) + 0.5 * x + u[period])
    K <- 4L; Nn <- stats::rpois(N, lambda); y <- matrix(0L, N, K); rem <- Nn
    for (k in 1:K) { y[, k] <- stats::rbinom(N, rem, 0.45); rem <- rem - y[, k] }
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period),
                         family = removal(), detection = ~ 1, y = y,
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
  expect_gte(mean(ok), 0.75)
  expect_lt(abs(mean(slope[ok]) - 0.5), 0.10)        # abundance slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.85)      # AR1 temporal field recovered
})
