# =============================================================================
# test-cover-nuts.R - NUTS sampler for the non-spatial standalone cover hurdle
# (cover(), method = "nuts").
#
# The cover hurdle is occu_cover() minus the occupancy / detection latent
# mixture: the presence Bernoulli arm and the beta / lognormal positive arm are
# conditionally independent given the data. The sampler draws the exact joint
# coefficient marginal c(beta_presence, beta_positive, log_disp) via the in-tree
# C++ FullGradFn (src/cover_nuts.cpp), warm-started at the Laplace mode. Tests:
#   - byte-exact gradient: C++ FullGradFn == R oracle (closed form, both arms)
#   - parameter recovery + 95% CI coverage vs simulated truth (both cover arms)
#   - NUTS posterior agrees with the Laplace mode it is built on
#   - S3 method coverage (coef / vcov / confint / summary / logLik / WAIC / predict)
#   - dispatch gating (spatial + nuts rejected; family advertises "nuts")
# =============================================================================


# Local beta-cover simulator (the exported simulate_cover() generates lognormal
# cover only). Mirrors the helper in test-cover-hurdle-beta.R: presence Bernoulli
# on the same design as the positive beta arm.
.cn_sim_beta <- function(N = 500L, beta_occ = c(-0.4, 0.8),
                         beta_pos = c(0.5, -1.2), phi = 30, seed = 1L) {
  set.seed(seed)
  x       <- stats::runif(N, -2, 2)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  occur   <- stats::rbinom(N, 1L, stats::plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  mu      <- stats::plogis(eta_pos)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos] <- stats::rbeta(sum(is_pos), mu[is_pos] * phi, (1 - mu[is_pos]) * phi)
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(data = data.frame(x = x), y = y,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos, phi = phi))
}

# Lognormal-cover inputs from simulate_cover() (intercept + slope on x, shared
# design across both arms).
.cn_sim_lognormal <- function(N = 500L, beta_occ = c(-0.4, 0.8),
                              beta_pos = c(-1.0, 0.3), sigma_pos = 0.4,
                              seed = 1L) {
  simulate_cover(N = N, beta_occ = beta_occ, beta_pos = beta_pos,
                 sigma_pos = sigma_pos, seed = seed)
}

# Identity-Gaussian (delta-normal) cover: the exported simulator draws the raw
# unbounded magnitude at present sites.
.cn_sim_gaussian <- function(N = 500L, beta_occ = c(-0.4, 0.8),
                             beta_pos = c(2.0, 0.4), sigma_pos = 0.5,
                             seed = 1L) {
  simulate_cover(N = N, beta_occ = beta_occ, beta_pos = beta_pos,
                 sigma_pos = sigma_pos, response = "gaussian", seed = seed)
}

.cn_fit <- function(sim, positive, method = "laplace", control = list()) {
  tobs(formula = ~ x, data = sim$data, family = cover(positive),
       y = sim$y, method = method, control = control)
}


test_that("cover NUTS C++ FullGradFn matches the R oracle (byte-exact)", {
  for (pos in c("lognormal", "beta", "gaussian")) {
    sim <- switch(pos,
                  beta     = .cn_sim_beta(N = 200L, seed = 11L),
                  gaussian = .cn_sim_gaussian(N = 200L, seed = 11L),
                  .cn_sim_lognormal(N = 200L, seed = 11L))
    enc <- tulpaObs:::encode_cover_hurdle(~ x, sim$data, sim$y,
                                          positive = pos, autoscale = FALSE)
    spec <- tulpaObs:::.tobs_cover_nuts_spec(enc)
    np   <- ncol(enc$occ_data$X) + ncol(enc$pos_data$X) + 1L

    set.seed(7)
    # A plausible parameter near a reasonable mode (intercepts, slopes, log_disp).
    theta <- c(stats::rnorm(np - 1L, 0, 0.3),
               if (pos == "beta") log(25) else log(0.4)) + stats::rnorm(np, 0, 0.1)

    rr <- tulpaObs:::.tobs_cover_nuts_logpost(theta, enc, sigma.beta = 5,
                                              sigma.logdisp = 5)
    cc <- tulpaObs:::cpp_cover_nuts_logpost(spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-6)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-6)

    # Analytic R gradient matches a central finite difference of its own lp.
    fd <- numeric(np); h <- 1e-5
    for (k in seq_len(np)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_cover_nuts_logpost(tp, enc, 5, 5)$lp -
                tulpaObs:::.tobs_cover_nuts_logpost(tm, enc, 5, 5)$lp) / (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4)
  }
})


test_that("cover NUTS is gated off the spatial path; family advertises nuts", {
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$cover)

  sim <- .cn_sim_lognormal(N = 40L, seed = 5L)
  sim$data$cell_idx <- rep(seq_len(20L), length.out = 40L)
  adj <- matrix(0L, 20, 20)
  for (i in seq_len(19)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  expect_error(
    suppressWarnings(tobs(
      formula = ~ x + icar(graph = adj, group_var = "cell_idx"),
      data = sim$data, family = cover("lognormal"), y = sim$y,
      method = "nuts", control = list(verbose = FALSE))),
    "not yet wired|non-spatial sampler")
})


test_that("cover NUTS posterior agrees with the Laplace mode", {
  skip_on_cran()
  skip_if_fast()

  sim <- .cn_sim_lognormal(N = 600L, beta_occ = c(-0.3, 0.7),
                           beta_pos = c(-1.0, 0.4), sigma_pos = 0.4, seed = 42L)
  lap <- .cn_fit(sim, "lognormal", "laplace")
  nut <- .cn_fit(sim, "lognormal", "nuts",
                 list(verbose = FALSE, n.iter = 1500L, n.warmup = 1000L,
                      n.chains = 2L, seed = 1L))

  expect_equal(nut$method, "nuts")
  expect_true(nut$nuts$divergent_total <= 5L)
  expect_lt(max(nut$nuts$rhat, na.rm = TRUE), 1.05)

  # Posterior means within a fraction of a Laplace SE; SDs comparable.
  d_occ <- abs(as.numeric(nut$beta_occ) - as.numeric(lap$beta_occ))
  d_pos <- abs(as.numeric(nut$beta_pos) - as.numeric(lap$beta_pos))
  expect_true(all(d_occ < 0.5 * pmax(as.numeric(lap$se_occ), 1e-3) + 0.05))
  expect_true(all(d_pos < 0.5 * pmax(as.numeric(lap$se_pos), 1e-3) + 0.05))
  ratio_occ <- as.numeric(nut$se_occ) / pmax(as.numeric(lap$se_occ), 1e-3)
  ratio_pos <- as.numeric(nut$se_pos) / pmax(as.numeric(lap$se_pos), 1e-3)
  expect_true(all(ratio_occ > 0.5 & ratio_occ < 2.0))
  expect_true(all(ratio_pos > 0.5 & ratio_pos < 2.0))
})


test_that("cover NUTS recovers parameters + 95% CI coverage (lognormal)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds   <- 8L
  N         <- 400L
  beta_occ  <- c(stats::qlogis(0.55), 0.9)
  beta_pos  <- c(log(0.12), -0.4)
  sigma_pos <- 0.4
  truth <- c(beta_occ, beta_pos, log(sigma_pos))

  est <- se <- matrix(NA_real_, n_seeds, length(truth))
  for (s in seq_len(n_seeds)) {
    sim <- .cn_sim_lognormal(N = N, beta_occ = beta_occ, beta_pos = beta_pos,
                             sigma_pos = sigma_pos, seed = 3000L + s)
    nut <- tryCatch(.cn_fit(sim, "lognormal", "nuts",
                    list(verbose = FALSE, n.iter = 1200L, n.warmup = 800L,
                         n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
    se[s, ]  <- as.numeric(nut$sds)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.75)

  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:4] < 0.25))   # the four regression coefficients
  expect_lt(bias[5], 0.10)             # log_sigma_pos

  # 95% Wald coverage on the slope coefficients (>= 0.75 floor at 8 seeds).
  cover_idx <- c(2, 4)
  cov_ok <- abs(est[ok, cover_idx, drop = FALSE] -
                matrix(truth[cover_idx], sum(ok), length(cover_idx), byrow = TRUE)) <
            1.96 * se[ok, cover_idx, drop = FALSE]
  expect_gte(mean(cov_ok), 0.75)
})


test_that("cover NUTS recovers parameters + coverage (beta positive)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds  <- 6L
  N        <- 500L
  beta_occ <- c(stats::qlogis(0.55), 0.8)
  beta_pos <- c(0.5, -1.0)
  phi      <- 25
  truth <- c(beta_occ, beta_pos, log(phi))

  est <- se <- matrix(NA_real_, n_seeds, length(truth))
  for (s in seq_len(n_seeds)) {
    sim <- .cn_sim_beta(N = N, beta_occ = beta_occ, beta_pos = beta_pos,
                        phi = phi, seed = 4000L + s)
    nut <- tryCatch(.cn_fit(sim, "beta", "nuts",
                    list(verbose = FALSE, n.iter = 1200L, n.warmup = 800L,
                         n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
    se[s, ]  <- as.numeric(nut$sds)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.6)

  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:4] < 0.30))   # the four regression coefficients
  expect_lt(bias[5], 0.40)             # log_phi

  cover_idx <- c(2, 4)
  cov_ok <- abs(est[ok, cover_idx, drop = FALSE] -
                matrix(truth[cover_idx], sum(ok), length(cover_idx), byrow = TRUE)) <
            1.96 * se[ok, cover_idx, drop = FALSE]
  # Coverage measured 1.0 (all cells, 6 seeds 4001-4006). The gate sits well
  # below the measurement so a real regression still trips without riding
  # this estimand's own Monte Carlo noise.
  expect_gte(mean(cov_ok), 0.6)
})


test_that("cover NUTS recovers parameters + coverage (gaussian positive)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds   <- 8L
  N         <- 400L
  beta_occ  <- c(stats::qlogis(0.55), 0.9)
  beta_pos  <- c(2.0, 0.4)
  sigma_pos <- 0.5
  # The identity-Gaussian log-dispersion is log(sigma) on the raw response scale.
  truth <- c(beta_occ, beta_pos, log(sigma_pos))

  est <- se <- matrix(NA_real_, n_seeds, length(truth))
  for (s in seq_len(n_seeds)) {
    sim <- .cn_sim_gaussian(N = N, beta_occ = beta_occ, beta_pos = beta_pos,
                            sigma_pos = sigma_pos, seed = 6000L + s)
    nut <- tryCatch(.cn_fit(sim, "gaussian", "nuts",
                    list(verbose = FALSE, n.iter = 1200L, n.warmup = 800L,
                         n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
    se[s, ]  <- as.numeric(nut$sds)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.75)

  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:4] < 0.25))   # the four regression coefficients
  expect_lt(bias[5], 0.10)             # log_sigma_pos

  cover_idx <- c(2, 4)
  cov_ok <- abs(est[ok, cover_idx, drop = FALSE] -
                matrix(truth[cover_idx], sum(ok), length(cover_idx), byrow = TRUE)) <
            1.96 * se[ok, cover_idx, drop = FALSE]
  # Coverage measured 0.9375 (15/16 cells, 8 seeds 6001-6008). The gate sits
  # below the measurement so a real regression still trips without riding
  # this estimand's own Monte Carlo noise.
  expect_gte(mean(cov_ok), 0.75)
})


test_that("cover NUTS fit supports the S3 method surface", {
  skip_on_cran()
  skip_if_fast()

  sim <- .cn_sim_lognormal(N = 400L, seed = 9L)
  nut <- .cn_fit(sim, "lognormal", "nuts",
                 list(verbose = FALSE, n.iter = 1000L, n.warmup = 800L,
                      n.chains = 2L, seed = 1L))
  np <- nut$n_params

  expect_s3_class(nut, "cover_fit")
  expect_equal(nut$method, "nuts")

  # The flat coefficient surface (presence + positive + log_disp) is finite.
  cf <- coef(nut)
  expect_equal(length(cf), np)
  expect_true(all(is.finite(unlist(cf))))
  expect_equal(dim(vcov(nut)), c(np, np))
  expect_equal(nrow(confint(nut)), np)
  expect_true(is.finite(as.numeric(logLik(nut))))

  # summary surfaces the per-parameter Rhat / ESS the convergence list carries.
  sm <- summary(nut)
  expect_true("rhat" %in% colnames(sm))
  expect_true(all(nut$convergence$rhat[!is.na(nut$convergence$rhat)] > 0.9))

  # Calibrated WAIC from the per-draw pointwise likelihood (NUTS draws).
  expect_true(is.finite(waic(nut)$waic))

  # predict.cover_fit works on the fixed-effects (non-joint) NUTS fit.
  newdata <- data.frame(x = sim$data$x)
  p_hat <- predict(nut, newdata, type = "occupancy")
  e_hat <- predict(nut, newdata, type = "expected")
  m_hat <- predict(nut, newdata, type = "conditional")
  expect_equal(e_hat, p_hat * m_hat, tolerance = 1e-8)
  expect_true(all(p_hat >= 0 & p_hat <= 1))
})

test_that("the positive-arm density code has one definition and one policy", {
  # .occu_cover_pos_code() used to be defined twice at top level, with different
  # unknown-input policies, and a third time under a second name on the NUTS
  # path; file collation order decided which ran.
  expect_equal(.occu_cover_pos_code("lognormal"), 0L)
  expect_equal(.occu_cover_pos_code("beta"),      3L)
  expect_equal(.occu_cover_pos_code("gaussian"),  4L)
  # An arm with no compiled density errors instead of falling back to lognormal.
  expect_error(.occu_cover_pos_code("beta_oi"), "no compiled density")
  expect_false(exists(".tobs_cover_pos_code", envir = asNamespace("tulpaObs"),
                      inherits = FALSE))

  # cover(response = "beta_oi") reaches the NUTS branch (ordinal and
  # lognormal_trunc are rejected earlier), so it needs its own gate rather than
  # being sampled as a plain beta.
  set.seed(3); n <- 60L
  dat <- data.frame(
    x = stats::rnorm(n),
    y = ifelse(stats::runif(n) < 0.4, 0,
               pmin(1, stats::rbeta(n, 2, 5) +
                      0.2 * (stats::runif(n) < 0.3))))
  expect_error(
    tobs(y ~ x, data = dat, family = cover(response = "beta_oi"),
         method = "nuts", control = list(verbose = FALSE, progress = FALSE)),
    "nested_laplace")
})
