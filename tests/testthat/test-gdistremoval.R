# gdistremoval() -- joint distance + removal sampling (unmarked gdistremoval,
# Amundson et al. 2014). Single-season: the detected birds are cross-classified by
# a distance band AND a removal period. Binomial thinning of a Poisson N is closed
# under thinning, so the total detected is Poisson(lambda*pdist*prem) and the band
# / period allocations are two conditional multinomials -- a closed-form marginal,
# maximised by optim BFGS with an observed-information vcov (the double_observer()
# recipe, with a distance multinomial + a depleting-removal multinomial).
#
# Recovery-grade: the fixed effects recover to near-nominal 95% coverage across
# seeds. Structural tests cover the registry, the dispatch gates, the input
# cross-classification check, the closed-form == integrate cell probs, and S3.

cutp <- c(0, 10, 20, 30, 40)

# --- (1) registry ----------------------------------------------------------

test_that("gdistremoval() reports laplace as its backend", {
  expect_identical(tulpaObs:::.tobs_family_methods$gdistremoval, "laplace")
})


# --- (2) closed-form half-normal cell probs == the integrate path ----------

test_that("the closed-form half-normal band probabilities match numeric integration", {
  for (tr in c("line", "point")) {
    cf <- as.numeric(tulpaObs:::.gdr_dist_cp(18, cutp, tr))
    ig <- tulpaObs:::.distance_pi(18, cutp, "halfnorm", tr, NULL)
    expect_lt(max(abs(cf - ig)), 1e-10)
  }
})


# --- (3) dispatch + input gates --------------------------------------------

test_that("gdistremoval() gates its required inputs", {
  sim <- simulate_gdistremoval(N = 40, cutpoints = cutp, n_periods = 3L, seed = 1)
  ok  <- list(data = sim$data, y = sim$y, y_rem = sim$y_rem,
              family = gdistremoval(cutpoints = cutp))

  # missing detection / y / y_rem / cutpoints
  expect_error(tobs(~ 1, data = ok$data, y = ok$y, y_rem = ok$y_rem,
                    family = ok$family), "detection")
  expect_error(tobs(~ 1, data = ok$data, y_rem = ok$y_rem, family = ok$family,
                    detection = ~ 1), "requires `y`")
  expect_error(tobs(~ 1, data = ok$data, y = ok$y, family = ok$family,
                    detection = ~ 1), "y_rem")
  expect_error(tobs(~ 1, data = ok$data, y = ok$y, y_rem = ok$y_rem,
                    family = gdistremoval(), detection = ~ 1), "cutpoints")

  # engine gate: only laplace
  expect_error(tobs(~ 1, data = ok$data, y = ok$y, y_rem = ok$y_rem,
                    family = ok$family, detection = ~ 1, method = "nuts"),
               "laplace")

  # per-site totals of the two responses must match (same detected birds)
  yr_bad <- sim$y_rem; yr_bad[1, 1] <- yr_bad[1, 1] + 5L
  expect_error(tobs(~ 1, data = ok$data, y = ok$y, y_rem = yr_bad,
                    family = ok$family, detection = ~ 1),
               "totals of y .* must match|cross-classified")
})


# --- (4) single fit + S3 ---------------------------------------------------

test_that("a gdistremoval fit recovers a single data set and wires S3", {
  skip_on_cran()
  set.seed(11)
  sim <- simulate_gdistremoval(N = 300, cutpoints = cutp, n_periods = 4L,
           beta_lambda = c(log(30), 0.3), beta_sigma = c(log(18), 0.1),
           beta_r = c(stats::qlogis(0.4), -0.2), seed = 11)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, y_rem = sim$y_rem,
              family = gdistremoval(cutpoints = cutp), detection = ~ det_cov1,
              removal = ~ rem_cov1, method = "laplace",
              control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "laplace")
  expect_true(fit$convergence$converged)

  truth <- c(sim$truth$beta_lambda, sim$truth$beta_sigma, sim$truth$beta_r)
  expect_true(all(abs(fit$means - truth) / fit$sds < 3))

  # S3 surface runs
  expect_length(unlist(coef(fit)), 6L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  fv <- fitted(fit)
  expect_length(fv$lambda, 300L); expect_length(fv$sigma, 300L)
  expect_length(fv$r, 300L)
  expect_length(predict(fit, type = "abundance"), 300L)
  expect_length(predict(fit, type = "distance"), 300L)
  expect_length(predict(fit, type = "removal"), 300L)
  expect_length(residuals(fit)$occ, 300L)
  expect_true(is.finite(waic(fit)$waic))
  expect_true(is.matrix(simulate(fit)$yDist))

  # nobs() counts both response tables: the band allocation and the period
  # allocation are separate multinomial factors of the marginal.
  expect_identical(nobs(fit), sum(!is.na(sim$y)) + sum(!is.na(sim$y_rem)))
})


# --- (5) multi-seed recovery + coverage ------------------------------------

test_that("gdistremoval recovers lambda / sigma / r with ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  bl <- c(log(30), 0.3); bs <- c(log(18), 0.1); br <- c(stats::qlogis(0.4), -0.2)
  truth  <- c(bl, bs, br); np <- length(truth)
  n_seed <- 20L
  est <- matrix(NA_real_, n_seed, np); cov <- matrix(FALSE, n_seed, np)
  conv <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_gdistremoval(N = 300, cutpoints = cutp, n_periods = 4L,
             beta_lambda = bl, beta_sigma = bs, beta_r = br, seed = 100 + s)
    fit <- tryCatch(
      tobs(~ abund_cov1, data = sim$data, y = sim$y, y_rem = sim$y_rem,
           family = gdistremoval(cutpoints = cutp), detection = ~ det_cov1,
           removal = ~ rem_cov1, method = "laplace",
           control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    b <- fit$means; se <- fit$sds
    est[s, ]  <- b; conv[s] <- isTRUE(fit$convergence$converged)
    cov[s, ]  <- (truth >= b - 1.96 * se) & (truth <= b + 1.96 * se)
  }
  # every fit converges; the coefficient means are ~unbiased and the Wald CIs
  # cover at the nominal rate (measured: bias <= 0.011, pooled coverage ~0.975).
  expect_gte(sum(conv), n_seed - 1L)
  expect_equal(colMeans(est, na.rm = TRUE), truth, tolerance = 0.05)
  expect_true(all(colMeans(cov) >= 0.85))
})
