test_that("abun/removal/distance criteria integrate the fit's own K.max (#420)", {
  skip_on_cran()
  ctl <- list(verbose = FALSE, progress = FALSE)

  sa <- simulate_abun(N = 60, J = 3, n.abund.covs = 1, n.det.covs = 1,
                      beta.lambda = c(log(4), 0.3), beta.p = c(0.2, 0), seed = 3)
  ka <- as.integer(max(sa$y, na.rm = TRUE) + 250L)
  fa <- tobs(~ abund_cov1, data = sa$data, family = abun(K.max = ka),
             detection = ~ 1, y = sa$y, method = "laplace", control = ctl)

  sr <- simulate_removal(N = 60, K = 3, n.abund.covs = 1, n.det.covs = 1,
                         beta.lambda = c(log(4), 0.3), beta.p = c(0.2, 0), seed = 3)
  kr <- as.integer(max(rowSums(sr$y)) + 250L)
  fr <- tobs(~ abund_cov1, data = sr$data, family = removal(K.max = kr),
             detection = ~ 1, y = sr$y, method = "laplace", control = ctl)

  sd <- simulate_distance(N = 60, key = "halfnorm", transect = "line",
                          n.abund.covs = 1, n.sigma.covs = 1,
                          beta.lambda = c(log(10), 0.3),
                          beta.sigma = c(log(0.45), 0), seed = 3)
  kd <- as.integer(3L * max(rowSums(sd$y)) + 250L)
  fd <- tobs(~ abund_cov1, data = sd$data,
             family = distance(key = "halfnorm", transect = "line",
                               cutpoints = sd$cutpoints, K.max = kd),
             detection = ~ 1, y = sd$y, method = "laplace", control = ctl)

  for (cs in list(list(fa, ka), list(fr, kr), list(fd, kd))) {
    fit <- cs[[1]]; k <- cs[[2]]
    expect_identical(as.integer(fit$K_max), k)
    expect_identical(tulpaObs:::.tobs_model_for_scoring(fit)$K_max, k)
    expect_true(all(is.finite(as.numeric(waic(fit)$estimates))))
  }

  # Without K.max the scoring truncation is the one the fitter resolved.
  f0 <- tobs(~ abund_cov1, data = sa$data, family = abun(), detection = ~ 1,
             y = sa$y, method = "laplace", control = ctl)
  expect_identical(tulpaObs:::.tobs_model_for_scoring(f0)$K_max,
                   as.integer(f0$K_max))
})
