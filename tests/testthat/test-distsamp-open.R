# distsamp_open() -- open-population distance sampling (unmarked distsampOpen).
# A Dail-Madsen open N-mixture (as dyn_abun()) with a distance-bin multinomial
# emission at each primary period. The band allocation is conditional on the
# period total, so it factors out of the abundance HMM: the marginal reuses the
# validated dyn_abun forward kernel (eta_p = logit(pdist)) plus the per-period band
# multinomials, maximised by optim BFGS with an analytic gradient (one kernel call
# per evaluation) and an observed-information vcov.
#
# Recovery-grade: the fixed effects recover with near-nominal 95% coverage across
# seeds. Structural tests cover the registry, the gates, the analytic gradient
# (FD-checked), and S3.

cutp <- c(0, 10, 20, 30, 40)

# --- (1) registry ----------------------------------------------------------

test_that("distsamp_open() reports laplace as its backend", {
  expect_identical(tulpaObs:::.tobs_family_methods$distsamp_open, "laplace")
})


# --- (2) analytic gradient matches finite differences ----------------------

test_that("the distsamp_open analytic gradient matches finite differences", {
  skip_on_cran()
  sim <- simulate_distsamp_open(N = 120, cutpoints = cutp, n_seasons = 3L,
           beta_lambda = c(log(10), 0.3), beta_sigma = c(log(18), 0.1),
           omega = 0.6, gamma = 2, seed = 4)
  m <- tulpaObs:::.tobs_build_distsamp_open(~ abund_cov1, ~ det_cov1, ~ 1, ~ 1,
         sim$data, sim$y, cutp, "line")
  th <- c(log(9), 0.2, log(17), 0.05, stats::qlogis(0.55), log(2.3))
  an <- tulpaObs:::.dso_grad(th, m)
  h  <- 1e-5; fd <- numeric(length(th))
  for (j in seq_along(th)) {
    tp <- th; tp[j] <- tp[j] + h; tm <- th; tm[j] <- tm[j] - h
    fd[j] <- -(tulpaObs:::.dso_negll(tp, m) - tulpaObs:::.dso_negll(tm, m)) / (2 * h)
  }
  expect_lt(max(abs(an - fd)), 1e-4)
  expect_gt(cor(an, fd), 0.9999)
})


# --- (3) dispatch + input gates --------------------------------------------

test_that("distsamp_open() gates its required inputs", {
  sim <- simulate_distsamp_open(N = 30, cutpoints = cutp, n_seasons = 3L, seed = 1)
  fam <- distsamp_open(cutpoints = cutp)

  expect_error(tobs(~ 1, data = sim$data, y = sim$y, family = fam),
               "detection")
  expect_error(tobs(~ 1, data = sim$data, family = fam, detection = ~ 1),
               "requires `y`")
  expect_error(tobs(~ 1, data = sim$data, y = sim$y, family = distsamp_open(),
                    detection = ~ 1), "cutpoints")
  expect_error(tobs(~ 1, data = sim$data, y = sim$y, family = fam,
                    detection = ~ 1, method = "nuts"), "laplace")
  # single primary period -> use distance()
  y1 <- sim$y[, , 1, drop = FALSE]
  expect_error(tobs(~ 1, data = sim$data, y = y1, family = fam, detection = ~ 1),
               ">= 2 primary periods|distance\\(\\)")
})


# --- (4) single fit + S3 ---------------------------------------------------

test_that("a distsamp_open fit recovers a single data set and wires S3", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_distsamp_open(N = 100, cutpoints = cutp, n_seasons = 4L,
           beta_lambda = c(log(8), 0.3), beta_sigma = c(log(18), 0.1),
           omega = 0.6, gamma = 1.8, seed = 301)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = distsamp_open(cutpoints = cutp), detection = ~ det_cov1,
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "laplace")
  expect_true(fit$convergence$converged)

  truth <- c(sim$truth$beta_lambda, sim$truth$beta_sigma,
             stats::qlogis(0.6), log(1.8))
  expect_true(all(abs(fit$means - truth) / fit$sds < 3.5))

  expect_length(unlist(coef(fit)), 6L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  fv <- fitted(fit)
  expect_length(fv$lambda, 100L); expect_length(fv$sigma, 100L)
  expect_length(fv$omega, 100L);  expect_length(fv$gamma, 100L)
  expect_length(predict(fit, type = "abundance"), 100L)
  expect_length(predict(fit, type = "survival"), 100L)
  expect_length(residuals(fit)$occ, 100L)
  expect_true(is.finite(waic(fit)$waic))
  expect_identical(dim(simulate(fit)), dim(sim$y))

  # nobs() counts the observed (site, band, period) counts.
  expect_identical(nobs(fit), sum(!is.na(sim$y)))
})


# --- (5) multi-seed recovery + coverage ------------------------------------

test_that("distsamp_open recovers lambda / sigma / omega / gamma with ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  tt <- c(log(8), 0.3, log(18), 0.1, stats::qlogis(0.6), log(1.8))
  np <- length(tt); n_seed <- 20L
  est <- matrix(NA_real_, n_seed, np); cov <- matrix(FALSE, n_seed, np)
  conv <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_distsamp_open(N = 100, cutpoints = cutp, n_seasons = 4L,
             beta_lambda = tt[1:2], beta_sigma = tt[3:4],
             omega = stats::plogis(tt[5]), gamma = exp(tt[6]), seed = 300 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, y = sim$y,
           family = distsamp_open(cutpoints = cutp), detection = ~ det_cov1,
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    b <- fit$means; se <- fit$sds
    est[s, ]  <- b; conv[s] <- isTRUE(fit$convergence$converged)
    cov[s, ]  <- (tt >= b - 1.96 * se) & (tt <= b + 1.96 * se)
  }
  # Measured (T = 4, 20 seeds): all converge, max coefficient bias ~0.068, every
  # per-parameter 95% coverage >= 0.90, pooled ~0.933. The survival / recruitment
  # ridge is resolved with 3 transitions (n_seasons = 4).
  expect_gte(sum(conv), n_seed - 1L)
  expect_equal(colMeans(est, na.rm = TRUE), tt, tolerance = 0.12)
  expect_true(all(colMeans(cov) >= 0.85))
})


# --- (6) negative-binomial initial abundance ------------------------------

test_that("the distsamp_open NB analytic gradient (incl. log_r) matches FD", {
  skip_on_cran()
  sim <- simulate_distsamp_open(N = 80, cutpoints = cutp, n_seasons = 3L,
           beta_lambda = c(log(6), 0.3), beta_sigma = c(log(18), 0.1),
           omega = 0.6, gamma = 1.2, mixture = "negbin", size = 5, seed = 4)
  m <- tulpaObs:::.tobs_build_distsamp_open(~ abund_cov1, ~ det_cov1, ~ 1, ~ 1,
         sim$data, sim$y, cutp, "line", mixture = "negbin")
  th <- c(log(6), 0.2, log(17), 0.05, stats::qlogis(0.55), log(1.1), log(4))
  an <- tulpaObs:::.dso_grad(th, m)
  expect_length(an, 7L)                       # the 4 arms + trailing log_r
  h  <- 1e-5; fd <- numeric(length(th))
  for (j in seq_along(th)) {
    tp <- th; tp[j] <- tp[j] + h; tm <- th; tm[j] <- tm[j] - h
    fd[j] <- -(tulpaObs:::.dso_negll(tp, m) - tulpaObs:::.dso_negll(tm, m)) / (2 * h)
  }
  expect_lt(max(abs(an - fd)), 1e-4)
  expect_gt(cor(an, fd), 0.9999)
})

test_that("a distsamp_open(negbin) fit recovers abundance / scale and surfaces r", {
  skip_if_fast()
  skip_on_cran()
  # Modest lambda keeps the cubic-in-K forward tractable; NB size / omega / gamma
  # sit on a weakly identified ridge at short series, so the recovery targets are
  # the well-identified abundance + distance-scale arms (and a finite r / log_r).
  sim <- simulate_distsamp_open(N = 100, cutpoints = cutp, n_seasons = 4L,
           beta_lambda = c(log(8), 0.3), beta_sigma = c(log(18), 0.1),
           omega = 0.6, gamma = 1.5, mixture = "negbin", size = 6, seed = 401)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = distsamp_open(cutpoints = cutp, mixture = "negbin"),
              detection = ~ det_cov1, method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$convergence$converged)
  expect_true("log_r" %in% names(fit$means))
  expect_true(is.finite(fit$r) && fit$r > 0)
  expect_length(unlist(coef(fit)), 6L)        # log_r is a trailing coord, not an arm
  b <- fit$means; se <- fit$sds
  tgt <- c("lambda_(Intercept)", "lambda_abund_cov1", "sigma_(Intercept)")
  tru <- c(log(8), 0.3, log(18))
  expect_true(all(abs(b[tgt] - tru) / se[tgt] < 3.5))
  expect_true(is.finite(waic(fit)$waic))
  expect_identical(dim(simulate(fit)), dim(sim$y))
})


# --- (7) zero-inflated initial abundance (ZIP / ZINB) ---------------------

test_that("a distsamp_open(zip) fit recovers abundance / scale and the ZI share", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_distsamp_open(N = 120, cutpoints = cutp, n_seasons = 4L,
           beta_lambda = c(log(8), 0.3), beta_sigma = c(log(18), 0.1),
           omega = 0.55, gamma = 1.2, mixture = "zip", zi = 0.35, seed = 411)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = distsamp_open(cutpoints = cutp, mixture = "zip"),
              detection = ~ det_cov1, omega = ~ 1, gamma = ~ 1,
              method = "laplace", control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$convergence$converged)
  expect_true(isTRUE(fit$zero_inflated))
  expect_true("zi_logit" %in% names(fit$means))
  # The structural-zero share is the distinctive ZIP quantity and recovers well.
  expect_equal(fit$zi_omega, 0.35, tolerance = 0.12)
  b <- fit$means; se <- fit$sds
  tgt <- c("lambda_(Intercept)", "lambda_abund_cov1", "sigma_(Intercept)")
  tru <- c(log(8), 0.3, log(18))
  expect_true(all(abs(b[tgt] - tru) / se[tgt] < 3.5))
  expect_true(is.finite(waic(fit)$waic))     # the per-site log_lik_site path
  expect_identical(dim(simulate(fit)), dim(sim$y))
})

test_that("distsamp_open(zinb) recovers the structural-zero share across seeds (#137)", {
  skip_if_fast()
  skip_on_cran()
  # Was smoke-only ("fits and names zi_logit / log_r"). The distinctive ZINB quantity
  # is the structural-zero share; over a few seeds it recovers to the same tolerance
  # the ZIP path meets (calibrated dev_notes/_calib_137_min.R: zi ~0.30, well within
  # 0.12 of truth). lambda / sigma / r sit on the usual ridge and are checked loosely.
  zi <- numeric(0)
  for (s in seq_len(3L)) {
    sim <- simulate_distsamp_open(N = 110, cutpoints = cutp, n_seasons = 4L,
             beta_lambda = c(log(9), 0.3), beta_sigma = c(log(16), 0.1),
             omega = 0.6, gamma = 1.2, mixture = "zinb", size = 8, zi = 0.3,
             seed = 80L + s)
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
                family = distsamp_open(cutpoints = cutp, mixture = "zinb"),
                detection = ~ det_cov1, omega = ~ 1, gamma = ~ 1,
                method = "laplace", control = list(verbose = FALSE))
    expect_s3_class(fit, "tobs_fit")
    expect_true(all(c("zi_logit", "log_r") %in% names(fit$means)))
    expect_true(is.finite(fit$r) && fit$r > 0)
    expect_true(is.finite(waic(fit)$waic))
    zi <- c(zi, fit$zi_omega)
  }
  # The structural-zero share recovers (the ZINB estimand); mean over seeds is tight.
  expect_lt(abs(mean(zi) - 0.3), 0.12)
})


# --- (8) alternative population dynamics (tp2..tp5) ------------------------

test_that("distsamp_open fits every alternative dynamics and recovers lambda/sigma", {
  skip_if_fast()
  skip_on_cran()
  # The density-regulated transitions (trend / autoreg / ricker / gompertz) use the
  # value-only forward-HMM kernel with a numeric gradient. Their transition rates
  # are per-capita (autoreg gamma), geometric (trend gamma), or bounded by a
  # carrying capacity (ricker / gompertz K), so the truth is set to keep abundance
  # low and the cubic-in-K forward cheap; an explicit small K_max bounds the cost.
  # The well-identified recovery target is the abundance / distance-scale arms --
  # survival / recruitment / K / r sit on the usual short-series ridge.
  cp3 <- c(0, 10, 20, 30)
  bl  <- c(log(4), 0.3); bs <- c(log(16), 0.1)
  spec <- list(
    notrend  = list(arm = "omega_(Intercept)",
                    sim = list(omega = 0.6, gamma = 1.5)),
    trend    = list(arm = "gamma_(Intercept)",
                    sim = list(gamma = 1.05)),
    autoreg  = list(arm = "gamma_(Intercept)",
                    sim = list(omega = 0.5, gamma = 0.25)),
    ricker   = list(arm = "K_(Intercept)",
                    sim = list(r = 0.3, K = 8)),
    gompertz = list(arm = "K_(Intercept)",
                    sim = list(r = 0.3, K = 8)))

  for (d in names(spec)) {
    args <- c(list(N = 40, cutpoints = cp3, n_seasons = 3L,
                   beta_lambda = bl, beta_sigma = bs, dynamics = d, seed = 51),
              spec[[d]]$sim)
    sim <- do.call(simulate_distsamp_open, args)
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
                family = distsamp_open(cutpoints = cp3, dynamics = d, K_max = 22L),
                detection = ~ det_cov1, method = "laplace",
                control = list(verbose = FALSE))
    info <- paste0("dynamics = ", d)
    expect_s3_class(fit, "tobs_fit")
    expect_true(isTRUE(fit$convergence$converged), info = info)
    # The dynamics-specific third / fourth arm is present and named.
    expect_true(spec[[d]]$arm %in% names(fit$means), info = info)
    # Abundance + distance-scale recover (rough tolerance at N = 40, bounded K).
    b <- fit$means
    expect_lt(abs(b[["lambda_(Intercept)"]] - bl[1]), 0.35)
    expect_lt(abs(b[["sigma_(Intercept)"]]  - bs[1]), 0.35)
    # S3 surface works on every dynamics.
    fv <- fitted(fit)
    expect_length(fv$lambda, 40L); expect_length(fv$sigma, 40L)
    expect_length(predict(fit, type = "abundance"), 40L)
    expect_true(is.finite(waic(fit)$waic), info = info)
    expect_identical(dim(simulate(fit)), dim(sim$y))
  }
})

test_that("distsamp_open density-dependent dynamics params have nominal CI coverage (#137)", {
  skip_if_fast()
  skip_on_cran()
  # Section 8 above only asserts lambda / sigma recovery per dynamics; the survival /
  # recruitment / carrying-capacity / growth params it introduces (omega, gamma, K, r)
  # were untested. They sit on the short-series ridge, so the POINT estimate can be
  # biased (e.g. ricker/gompertz r ~0.46 vs 0.3), but the observed-Fisher CI widens on
  # the ridge and COVERS the truth -- the honest, testable claim. Coverage is pooled
  # across every introduced param and seed (measured 38/40 = 0.95, per-arm 0.80-1.00;
  # dev_notes/_calib_137_min.R) so the floor is robust to the ridge point-bias.
  cp3 <- c(0, 10, 20, 30)
  bl  <- c(log(4), 0.3); bs <- c(log(16), 0.1)
  qlog <- stats::qlogis
  spec <- list(
    notrend  = list(sim = list(omega = 0.6, gamma = 1.5), ns = 3L, K = 22L,
                    arms = c("omega_(Intercept)"),
                    truth = function(s) c(qlog(s$truth$omega))),
    trend    = list(sim = list(gamma = 1.05), ns = 3L, K = 22L,
                    arms = c("gamma_(Intercept)"),
                    truth = function(s) c(log(s$truth$gamma))),
    autoreg  = list(sim = list(omega = 0.5, gamma = 0.25), ns = 3L, K = 22L,
                    arms = c("omega_(Intercept)", "gamma_(Intercept)"),
                    truth = function(s) c(qlog(s$truth$omega), log(s$truth$gamma))),
    ricker   = list(sim = list(K = 8, r = 0.3), ns = 4L, K = 24L,
                    arms = c("K_(Intercept)", "r_(Intercept)"),
                    truth = function(s) c(log(s$truth$K), s$truth$r)),
    gompertz = list(sim = list(K = 8, r = 0.3), ns = 4L, K = 24L,
                    arms = c("K_(Intercept)", "r_(Intercept)"),
                    truth = function(s) c(log(s$truth$K), s$truth$r)))

  covered <- logical(0)
  for (d in names(spec)) {
    sp <- spec[[d]]; info <- paste0("dynamics = ", d)
    lam_ok <- sig_ok <- logical(0)
    for (seed in seq_len(5L)) {
      args <- c(list(N = 40L, cutpoints = cp3, n_seasons = sp$ns,
                     beta_lambda = bl, beta_sigma = bs, dynamics = d, seed = 50L + seed),
                sp$sim)
      sim <- do.call(simulate_distsamp_open, args)
      fit <- tryCatch(tobs(~ abund_cov1, data = sim$data, y = sim$y,
                       family = distsamp_open(cutpoints = cp3, dynamics = d, K_max = sp$K),
                       detection = ~ det_cov1, method = "laplace",
                       control = list(verbose = FALSE)), error = function(e) NULL)
      if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
      est <- fit$means[sp$arms]; se <- fit$sds[sp$arms]; tv <- sp$truth(sim)
      covered <- c(covered, abs(est - tv) <= 1.96 * se)
      lam_ok <- c(lam_ok, abs(fit$means[["lambda_(Intercept)"]] - bl[1]) < 0.4)
      sig_ok <- c(sig_ok, abs(fit$means[["sigma_(Intercept)"]]  - bs[1]) < 0.4)
    }
    # Abundance / distance-scale recover on the majority of seeds per dynamics.
    expect_gte(mean(lam_ok), 0.6, label = paste0("lambda recovery ", info))
    expect_gte(mean(sig_ok), 0.6, label = paste0("sigma recovery ", info))
  }
  # The dynamics-introduced params (omega / gamma / K / r) cover at the nominal rate
  # when pooled -- wide ridge CIs cover even where the point estimate is biased.
  expect_gte(mean(covered), 0.75)
})
