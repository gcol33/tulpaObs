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
  expect_true(is.finite(tobs_waic(fit)$waic))
  expect_identical(dim(simulate(fit)), dim(sim$y))
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
