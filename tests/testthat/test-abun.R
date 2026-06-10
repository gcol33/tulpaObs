# N-mixture abundance (Royle 2004), Poisson, non-spatial Laplace.
#
# Recovery-grade tests (per the "statistical code needs recovery tests" rule):
# point recovery against simulated truth + 95% CI coverage across seeds. The
# structural tests below check the family wiring / S3 surface.

test_that("abun() family is wired and reports its supported methods", {
  f <- abun()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "abun")
  expect_identical(f$status, "working")
  expect_identical(f$params$mixture, "poisson")
})

test_that("single Poisson N-mixture fit recovers truth", {
  skip_if_fast()
  beta_lambda <- c(log(4), 0.6, -0.4)
  beta_p      <- c(0.2, 0.5)
  sim <- simulate_abun(N = 400, J = 5, n_abund_covs = 2, n_det_covs = 1,
                       beta_lambda = beta_lambda, beta_p = beta_p, seed = 11)

  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = abun(), detection = ~ det_cov1, y = sim$y,
              method = "laplace")

  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_lambda, beta_p)
  est   <- as.numeric(fit$means)
  se    <- as.numeric(fit$sds)

  # Each coefficient within ~3 SE of truth, and abundance slopes tight in
  # absolute terms (N=400 is well-identified).
  expect_true(all(abs(est - truth) / se < 3))
  expect_lt(abs(est[2] - beta_lambda[2]), 0.12)   # lambda_abund_cov1
  expect_lt(abs(est[3] - beta_lambda[3]), 0.12)   # lambda_abund_cov2

  # Coefficient layout + back-transformed intercepts.
  expect_named(fit$means, c("lambda_(Intercept)", "lambda_abund_cov1",
                            "lambda_abund_cov2", "p_(Intercept)", "p_det_cov1"))
  expect_equal(unname(fit$intercepts$lambda), exp(est[1]), tolerance = 1e-8)
  expect_equal(unname(fit$intercepts$p), plogis(est[4]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(5), 0.5, -0.3)
  beta_p      <- c(0.4, 0.4)
  n_seed <- 30L
  p_names <- c("lambda_(Intercept)", "lambda_abund_cov1", "lambda_abund_cov2",
               "p_(Intercept)", "p_det_cov1")
  truth <- c(beta_lambda, beta_p)

  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_abun(N = 150, J = 4, n_abund_covs = 2, n_det_covs = 1,
                         beta_lambda = beta_lambda, beta_p = beta_p, seed = 100 + s)
    fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
                family = abun(), detection = ~ det_cov1, y = sim$y,
                method = "laplace")
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  # Nominal 95%; allow Monte-Carlo slack at 30 seeds (>= 0.85 per coefficient).
  expect_true(all(cover_rate >= 0.85),
              info = paste(p_names, round(cover_rate, 2), collapse = " | "))
})

test_that("S3 surface works for N-mixture fits", {
  skip_if_fast()
  # Explicit betas so the abund_cov1 effect direction is deterministic.
  sim <- simulate_abun(N = 200, J = 4, n_abund_covs = 2, n_det_covs = 1,
                       beta_lambda = c(log(4), 0.5, -0.3),
                       beta_p = c(0.3, 0.4), seed = 3)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = abun(), detection = ~ det_cov1, y = sim$y,
              method = "laplace")

  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), sum(!is.na(sim$y)))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("lambda", "p", "N"))
  expect_length(fv$lambda, 200L)
  expect_true(all(fv$lambda > 0))

  # Design-matrix prediction on the abundance arm.
  X0 <- cbind(1, c(-1, 0, 1), 0)
  pr <- predict(fit, X.0 = X0)
  expect_true(all(pr$mean > 0))
  expect_true(all(diff(pr$mean) > 0))   # increasing in abund_cov1

  # simulate() returns a count matrix respecting the NA pattern.
  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim >= 0))

  rr <- residuals(fit, type = "pearson")
  expect_equal(dim(rr), dim(sim$y))
})

test_that("negbin N-mixture recovers truth, surfaces dispersion, covers CIs", {
  skip_on_cran()
  skip_if_fast()
  beta_lambda <- c(log(5), 0.6, -0.4)
  beta_p      <- c(0.4, 0.5)
  size_true   <- 2                       # Var(N) = lambda + lambda^2 / 2

  # --- point recovery on one well-identified fit ---
  sim <- simulate_abun(N = 400, J = 6, n_abund_covs = 2, n_det_covs = 1,
                       beta_lambda = beta_lambda, beta_p = beta_p,
                       mixture = "negbin", size = size_true, seed = 21)
  fit <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sim$data,
              family = abun(mixture = "negbin"), detection = ~ det_cov1,
              y = sim$y, method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$mixture, "negbin")

  truth <- c(beta_lambda, beta_p)
  est   <- as.numeric(fit$means[1:5])
  se    <- as.numeric(fit$sds[1:5])
  expect_true(all(abs(est - truth) / se < 3))

  # Dispersion is carried as the trailing log_r coordinate (with an SE) and
  # also summarized on the natural r scale; coef()/vcov() generics report the
  # two arms only (log_r is a nuisance, not an arm coefficient).
  expect_true("log_r" %in% rownames(fit$vcov))
  expect_equal(dim(fit$vcov), c(6L, 6L))
  expect_false(is.null(fit$nmix_dispersion))
  se_logr <- sqrt(fit$vcov["log_r", "log_r"])
  expect_lt(abs(fit$nmix_dispersion$log_r - log(size_true)) / se_logr, 3)
  expect_equal(fit$nmix_dispersion$r, exp(fit$nmix_dispersion$log_r),
               tolerance = 1e-6)
  cf <- coef(fit)
  expect_named(cf, c("lambda", "p"))
  expect_length(unlist(cf), 5L)

  # --- CI coverage across seeds (5 betas + log dispersion) ---
  n_seed <- 25L
  covered <- matrix(NA, n_seed, 6L)
  tr <- c(truth, log(size_true))
  for (s in seq_len(n_seed)) {
    sm <- simulate_abun(N = 200, J = 5, n_abund_covs = 2, n_det_covs = 1,
                        beta_lambda = beta_lambda, beta_p = beta_p,
                        mixture = "negbin", size = size_true, seed = 300 + s)
    ft <- tobs(formula = ~ abund_cov1 + abund_cov2, data = sm$data,
               family = abun(mixture = "negbin"), detection = ~ det_cov1,
               y = sm$y, method = "laplace", control = list(verbose = FALSE))
    es <- c(as.numeric(ft$means[1:5]), ft$nmix_dispersion$log_r)
    sd <- c(as.numeric(ft$sds[1:5]), sqrt(ft$vcov["log_r", "log_r"]))
    covered[s, ] <- abs(es - tr) <= 1.96 * sd
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.8),
              info = paste(c("lam0", "lam1", "lam2", "p0", "p1", "log_r"),
                           round(cover_rate, 2), collapse = " | "))
})

test_that("poisson and negbin paths differ only by the dispersion coordinate", {
  sim <- simulate_abun(N = 150, J = 4, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(4), 0.5), beta_p = c(0.3, 0.4),
                       mixture = "negbin", size = 2, seed = 7)
  base <- list(formula = ~ abund_cov1, data = sim$data, detection = ~ det_cov1,
               y = sim$y, method = "laplace", control = list(verbose = FALSE))

  fitP  <- do.call(tobs, c(base, list(family = abun(mixture = "poisson"))))
  fitNB <- do.call(tobs, c(base, list(family = abun(mixture = "negbin"))))

  expect_identical(fitP$mixture, "poisson")
  expect_null(fitP$nmix_dispersion)
  expect_false("log_r" %in% rownames(fitP$vcov))
  expect_equal(dim(fitP$vcov), c(4L, 4L))

  expect_identical(fitNB$mixture, "negbin")
  expect_equal(dim(fitNB$vcov), c(5L, 5L))   # 4 betas + log_r
  expect_true(is.finite(fitNB$nmix_dispersion$r))
})

test_that("unsupported methods are rejected with the supported set", {
  skip_if_fast()
  sim <- simulate_abun(N = 60, J = 3, seed = 9)
  # abun() supports laplace / nested_laplace / nuts; the EM-correction routes
  # (gibbs / mi) are not wired for the N-mixture marginal.
  expect_error(
    tobs(formula = ~ abund_cov1, data = sim$data, family = abun(),
         detection = ~ det_cov1, y = sim$y, method = "laplace_gibbs"),
    "not available for abun"
  )
})

# Build a rook-adjacency grid graph of `side` x `side` cells.
.grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng)) {
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L)
      adj[i, j] <- 1L
  }
  adj
}

# Simulate spatial N-mixture: log lambda_i = b0 + b1 x_i + phi_i with a
# smoothed, demeaned ICAR-like offset; counts y_ij ~ Binomial(N_i, p_i).
# `size = NULL` draws N ~ Poisson; a finite `size` draws N ~ NegBin(mu, size).
.sim_spatial_nmix <- function(adj, b_lambda, b_p, J = 5L, sd_phi = 0.6,
                              size = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ng <- nrow(adj)
  phi <- as.numeric(scale(rnorm(ng)))
  for (rep in 1:3) {
    pn <- phi
    for (i in seq_len(ng)) {
      nb <- which(adj[i, ] == 1L)
      pn[i] <- 0.5 * phi[i] + 0.5 * mean(phi[nb])
    }
    phi <- pn
  }
  phi <- sd_phi * as.numeric(scale(phi)); phi <- phi - mean(phi)
  x_ab <- rnorm(ng); x_det <- rnorm(ng)
  lambda <- exp(b_lambda[1] + b_lambda[2] * x_ab + phi)
  p <- plogis(b_p[1] + b_p[2] * x_det)
  N <- if (is.null(size)) rpois(ng, lambda) else rnbinom(ng, size = size, mu = lambda)
  y <- matrix(NA_integer_, ng, J)
  for (i in seq_len(ng)) y[i, ] <- rbinom(J, N[i], p[i])
  list(y = y, data = data.frame(abund_cov1 = x_ab, det_cov1 = x_det),
       truth = c(b_lambda, b_p))
}

test_that("areal-spatial negbin N-mixture fits and reports grid-integrated r", {
  skip_on_cran()
  skip_if_fast()
  adj <- .grid_adj(6L)
  b_lambda <- c(log(6), 0.5); b_p <- c(0.4, 0.4)
  sim <- .sim_spatial_nmix(adj, b_lambda, b_p, J = 6L, size = 3, seed = 77)
  fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
              family = abun(mixture = "negbin"), detection = ~ det_cov1,
              y = sim$y, method = "nested_laplace", control = list(verbose = FALSE))

  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$mixture, "negbin")

  # The NB size r is integrated over the outer grid (alongside tau), so it is
  # reported as a hyperparameter, not as a vcov coordinate.
  expect_false("log_r" %in% rownames(fit$vcov))
  expect_false(is.null(fit$nmix_dispersion))
  expect_true(is.finite(fit$nmix_dispersion$r))
  expect_true(is.finite(fit$nmix_hyper$r[["mean"]]))

  # Abundance slope recovers within ~3 SE under the spatial offset.
  expect_lt(abs(fit$means[["lambda_abund_cov1"]] - 0.5) /
              fit$sds[["lambda_abund_cov1"]], 3)
})

test_that("areal-spatial (icar/bym2) N-mixture fits with calibrated coef cov", {
  skip_on_cran()
  skip_if_fast()
  adj <- .grid_adj(6L)
  b_lambda <- c(log(5), 0.5); b_p <- c(0.3, 0.4)
  sim <- .sim_spatial_nmix(adj, b_lambda, b_p, J = 5L, seed = 42)
  truth <- c("lambda_(Intercept)", "lambda_abund_cov1",
             "p_(Intercept)", "p_det_cov1")

  for (mod in c("icar", "bym2")) {
    f <- if (mod == "icar") (~ abund_cov1 + icar(graph = adj)) else
                            (~ abund_cov1 + bym2(graph = adj))
    fit <- tobs(formula = f, data = sim$data, family = abun(),
                detection = ~ det_cov1, y = sim$y, method = "nested_laplace",
                control = list(verbose = FALSE))
    expect_identical(fit$method, "nested_laplace")

    V <- vcov(fit)
    expect_true(all(is.finite(V)))
    expect_true(all(eigen(V, only.values = TRUE)$values > 0))   # PD
    # Constrained intercept SE must be finite and sane (not the 1e3+ blow-up of
    # the unconstrained improper-prior Hessian).
    expect_lt(fit$sds[["lambda_(Intercept)"]], 5)
    # Covariate slopes recover within ~3 SE.
    expect_lt(abs(fit$means[["lambda_abund_cov1"]] - 0.5) /
                fit$sds[["lambda_abund_cov1"]], 3)
    expect_lt(abs(fit$means[["p_det_cov1"]] - 0.4) /
                fit$sds[["p_det_cov1"]], 3)
  }
})

test_that("spatial-slope 95% CI covers truth across seeds (calibration)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .grid_adj(6L)
  b_lambda <- c(log(5), 0.5); b_p <- c(0.3, 0.4)
  n_seed <- 20L
  cov_slope <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- .sim_spatial_nmix(adj, b_lambda, b_p, J = 5L, seed = 200 + s)
    fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
                family = abun(), detection = ~ det_cov1, y = sim$y,
                method = "nested_laplace", control = list(verbose = FALSE))
    est <- fit$means[["lambda_abund_cov1"]]; se <- fit$sds[["lambda_abund_cov1"]]
    cov_slope[s] <- abs(est - 0.5) <= 1.96 * se
  }
  # Calibrated SE -> ~95% coverage; allow Monte-Carlo slack at 20 seeds.
  expect_gte(mean(cov_slope), 0.8)
})

test_that("spatial N-mixture carries cross-arm (lambda,p) covariance", {
  skip_if_fast()
  # tulpaObs#19: under the spatial path the coefficient covariance must NOT be
  # block-diagonal across the abundance and detection arms -- the cross-arm
  # (lambda, p) block, folded through the shared field, has to be non-zero so
  # derived quantities combining the two arms propagate the correlation.
  adj <- .grid_adj(7L)
  b_lambda <- c(log(6), 0.6); b_p <- c(0.3, 0.5)
  sim <- .sim_spatial_nmix(adj, b_lambda, b_p, J = 6L, seed = 19)
  fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
              family = abun(), detection = ~ det_cov1, y = sim$y,
              method = "nested_laplace", control = list(verbose = FALSE))
  nm <- rownames(fit$vcov)
  cross <- fit$vcov[grep("^lambda_", nm), grep("^p_", nm), drop = FALSE]
  expect_gt(max(abs(cross)), 1e-6)
  # The draws used for derived quantities inherit that covariance.
  cor_draws <- cor(fit$draws[, grep("^lambda_", nm)[1]],
                   fit$draws[, grep("^p_", nm)[1]])
  expect_gt(abs(cor_draws), 0.02)
})

test_that("spatial N-mixture expected-count (lambda*p) CI is calibrated", {
  skip_on_cran()
  skip_if_fast()
  # The derived expected count mu = lambda * p combines BOTH arms, so its CI is
  # only calibrated if the posterior draws carry the cross-arm covariance
  # (tulpaObs#19). We evaluate mu per draw at a fixed design point and check
  # 95% quantile-CI coverage vs the simulated truth across seeds.
  adj <- .grid_adj(6L)
  b_lambda <- c(log(5), 0.5); b_p <- c(0.3, 0.4)
  x_lam <- c(1, 0.4); x_p <- c(1, -0.3)               # fixed evaluation point
  mu_true <- exp(sum(x_lam * b_lambda)) * plogis(sum(x_p * b_p))
  n_seed <- 20L
  covered <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- .sim_spatial_nmix(adj, b_lambda, b_p, J = 5L, seed = 500 + s)
    fit <- tobs(formula = ~ abund_cov1 + icar(graph = adj), data = sim$data,
                family = abun(), detection = ~ det_cov1, y = sim$y,
                method = "nested_laplace", control = list(verbose = FALSE))
    nm <- rownames(fit$vcov)
    bl <- fit$draws[, grep("^lambda_", nm), drop = FALSE]
    bp <- fit$draws[, grep("^p_", nm), drop = FALSE]
    mu_draws <- exp(as.numeric(bl %*% x_lam)) * plogis(as.numeric(bp %*% x_p))
    ci <- quantile(mu_draws, c(0.025, 0.975))
    covered[s] <- ci[1] <= mu_true && mu_true <= ci[2]
  }
  expect_gte(mean(covered), 0.8)
})


# --- NUTS (gcol33/tulpaObs#41) ----------------------------------------------
# The non-spatial N-mixture sampler: the in-tree C++ FullGradFn (src/abun_nuts.cpp)
# over the closed-form marginal, driving tulpa's NUTS engine. The R target
# .tobs_abun_nuts_logpost (R/abun_nuts.R) is the oracle; the FAST test
# cross-checks the C++ port against it (no fit), the SLOW tests sample.

test_that("cpp_abun_nuts_joint_logpost matches the R oracle (byte-exact)", {
  set.seed(7)
  for (is_nb in c(FALSE, TRUE)) {
    mix <- if (is_nb) "negbin" else "poisson"
    sim <- simulate_abun(N = 40, J = 4, n_abund_covs = 1, n_det_covs = 1,
                         mixture = mix, size = 3, seed = 19)
    y <- sim$y; n_sites <- nrow(y); mv <- ncol(y)
    valid  <- as.vector(t(!is.na(y)))
    smat   <- matrix(seq_len(n_sites), n_sites, mv)
    y_long <- as.integer(as.vector(t(y))[valid])
    site_idx <- as.vector(t(smat))[valid]
    X_lambda <- model.matrix(~ abund_cov1, sim$data)
    X_p <- model.matrix(~ det_cov1, sim$data)[site_idx, , drop = FALSE]
    p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
    K_max <- as.integer(max(y_long) + 80L)
    model <- list(y_long = y_long, site_idx = site_idx,
                  X_processes = list(X_lambda, X_p))
    marg <- tulpaObs:::.tobs_abun_nuts_marginal(
      model, mixture = if (is_nb) "NB" else "P", K_max = K_max)
    lay  <- tulpaObs:::.tobs_abun_nuts_layout(p_lam, p_p, is_nb)
    spec <- list(y = y_long, site_idx = as.integer(site_idx),
                 X_lambda = X_lambda, X_p = X_p, n_sites = n_sites,
                 K_max = K_max, is_nb = is_nb)
    theta <- c(sim$truth$beta_lambda, sim$truth$beta_p)
    if (is_nb) theta <- c(theta, log(sim$truth$size))
    theta <- theta + rnorm(length(theta), 0, 0.2)
    rr <- tulpaObs:::.tobs_abun_nuts_logpost(theta, marg, lay,
                                             sigma.beta = 10, sigma.logr = 1.5)
    cc <- tulpaObs:::cpp_abun_nuts_joint_logpost(spec, theta, 10, 1.5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-7)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-7)
  }
})

test_that("abun() declares nuts among its methods and samples an icar field", {
  skip_on_cran(); skip_if_fast()
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$abun)
  # An intrinsic (icar) field + nuts now samples via the sum-to-zero
  # reparameterisation (#71): its rank-deficient precision's constant direction is
  # dropped, so the field is centred and the geometry is well conditioned (the
  # recovery test above checks calibration + 0 divergences; here just confirm the
  # public path returns a nuts fit with a centred field).
  sim <- simulate_abun(N = 12, J = 3, n_abund_covs = 1, n_det_covs = 1, seed = 5)
  adj <- matrix(0L, 12, 12)
  for (i in 1:11) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  fit <- tobs(~ abund_cov1 + icar(graph = adj), data = sim$data, y = sim$y,
              family = abun(), detection = ~ det_cov1, method = "nuts",
              control = list(n.iter = 100L, n.warmup = 100L, seed = 1L,
                             verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nuts")
  expect_false(is.null(fit$spatial_field))
  expect_lt(abs(sum(fit$spatial_field)), 1e-6)          # sum-to-zero centred
})

test_that("tobs(abun(), method='nuts') recovers truth and scores WAIC", {
  skip_on_cran(); skip_if_fast()
  set.seed(101)
  sim <- simulate_abun(N = 60, J = 4, n_abund_covs = 2, n_det_covs = 1,
                       beta_lambda = c(log(4), 0.5, -0.3), beta_p = c(0.2, -0.4),
                       mixture = "poisson", seed = 101)
  fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data, y = sim$y,
              family = abun(), detection = ~ det_cov1, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))
  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 500L)
  truth <- c(sim$truth$beta_lambda, sim$truth$beta_p)
  expect_true(all(abs(as.numeric(fit$means) - truth) / as.numeric(fit$sds) < 3.5))
  # real sampler diagnostics replaced the NA placeholders
  expect_false(any(is.na(fit$divergent)))
  expect_false(any(is.na(fit$accept_prob)))
  expect_lt(mean(fit$divergent), 0.2)
  # WAIC / LOO from the NUTS draws (per-site marginal pointwise log-lik)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
  expect_lt(w$p_waic, nrow(sim$y))            # p_waic < n_sites
})

test_that("tobs(abun(mixture='negbin'), method='nuts') recovers dispersion", {
  skip_on_cran(); skip_if_fast()
  set.seed(101)
  sim <- simulate_abun(N = 60, J = 4, n_abund_covs = 2, n_det_covs = 1,
                       beta_lambda = c(log(4), 0.5, -0.3), beta_p = c(0.2, -0.4),
                       mixture = "negbin", size = 3, seed = 101)
  fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data, y = sim$y,
              family = abun(mixture = "negbin"), detection = ~ det_cov1,
              method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))
  expect_true("log_r" %in% names(fit$means))
  expect_false(is.null(fit$nmix_dispersion))
  expect_lt(abs(log(fit$nmix_dispersion$r) - log(3)), log(2.5))
  expect_lt(mean(fit$divergent), 0.2)
})


# --- NUTS + areal proper-CAR field on the abundance arm (tulpaObs#51) ----------

test_that("abun() NUTS + proper-CAR areal field recovers + calibrates to nested-Laplace", {
  skip_on_cran()
  skip_if_fast()
  # Fixed-hyper non-centered field (precision fixed at the nested-Laplace estimate):
  # the proper-CAR full-rank precision gives a clean non-centered geometry, so NUTS
  # samples without divergences and the coefficient SDs match the nested-Laplace SEs.
  adj <- .grid_adj(6L)
  sim <- .sim_spatial_nmix(adj, c(log(6), 0.5), c(0.2, 0.4), J = 5L, seed = 42)
  nl <- tobs(formula = ~ abund_cov1 + car_proper(graph = adj), data = sim$data,
             family = abun(), detection = ~ det_cov1, y = sim$y,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tobs(formula = ~ abund_cov1 + car_proper(graph = adj), data = sim$data,
             family = abun(), detection = ~ det_cov1, y = sim$y, method = "nuts",
             verbose = FALSE,
             control = list(n.iter = 400L, n.warmup = 400L, seed = 7L, progress = FALSE))
  expect_identical(nu$method, "nuts")
  expect_false(is.null(nu$spatial_field))
  expect_lt(mean(nu$nuts$divergent), 0.05)                          # clean geometry
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
              nu$sds[["lambda_abund_cov1"]], 3)                     # slope recovered
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.8)           # field agrees with NL
  # NUTS SDs calibrate to the nested-Laplace SEs (fixed-hyper, same field precision).
  rr <- nu$sds[["lambda_abund_cov1"]] / nl$sds[["lambda_abund_cov1"]]
  expect_gt(rr, 0.6); expect_lt(rr, 1.6)
})

test_that("abun() NUTS + icar/bym2 sum-to-zero field recovers, centred, 0 divergences (#71)", {
  skip_on_cran()
  skip_if_fast()
  # The intrinsic icar / bym2 fields are sampled via the sum-to-zero
  # reparameterisation: the whitened loading drops the precision null-space
  # (constant) direction, so the field-mean is removed and the NUTS geometry is
  # well conditioned (no flat direction maxing the tree depth). icar -> n_raw =
  # n - 1, automatically centred; bym2 -> structured (centred) + unstructured.
  for (typ in c("icar", "bym2")) {
    adj <- .grid_adj(6L)
    sim <- .sim_spatial_nmix(adj, c(log(6), 0.5), c(0.2, 0.4), J = 5L,
                             seed = if (typ == "icar") 21 else 22)
    f <- if (typ == "icar") (~ abund_cov1 + icar(graph = adj)) else
                            (~ abund_cov1 + bym2(graph = adj))
    nl <- tobs(formula = f, data = sim$data, family = abun(), detection = ~ det_cov1,
               y = sim$y, method = "nested_laplace",
               control = list(verbose = FALSE, progress = FALSE))
    nu <- tobs(formula = f, data = sim$data, family = abun(), detection = ~ det_cov1,
               y = sim$y, method = "nuts", verbose = FALSE,
               control = list(n.iter = 400L, n.warmup = 400L, seed = 3L, progress = FALSE))
    expect_identical(nu$method, "nuts")
    expect_equal(mean(nu$nuts$divergent), 0)                          # clean geometry
    expect_lt(max(nu$nuts$treedepth), 9)                              # no flat-mean blowup
    expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.95)          # reproduces NL field
    expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
                nu$sds[["lambda_abund_cov1"]], 3)                     # slope recovered
    if (typ == "icar") expect_lt(abs(sum(nu$spatial_field)), 1e-6)    # sum-to-zero
    rr <- nu$sds[["lambda_abund_cov1"]] / sqrt(diag(nl$vcov))[["lambda_abund_cov1"]]
    expect_gt(rr, 0.6); expect_lt(rr, 1.6)                            # SD calibrated to NL
  }
})
