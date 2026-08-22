# gcol33/tulpaObs#223: two pre-existing gaps in fitted()/predict() found while
# fixing #218 (detection-arm spatial field). Neither is a spatial-field issue:
#
#   1. fitted(fit)$z on an int_occu() fit reported the prior marginal (psi),
#      never conditioning on any source's detection history.
#   2. predict(fit, X.0=, type="detection"/"both") never predicted detection
#      at all for single-season/int_occu() fits, silently returning only
#      occupancy regardless of `type`.

.iofp_sim <- function(seed, n_sites = 80, J1 = 4L, J2 = 3L,
                      beta_occ = c(0.2, 0.7), beta_det = c(0.4, 0.5)) {
  set.seed(seed)
  x_cov <- rnorm(n_sites); det_cov <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(beta_occ[1] + beta_occ[2] * x_cov))
  mk <- function(J, p0) {
    p <- plogis(p0 + beta_det[2] * det_cov)
    y <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
    y
  }
  list(data = data.frame(occ_cov = x_cov, det_cov = det_cov),
       y = list(src1 = mk(J1, beta_det[1]), src2 = mk(J2, beta_det[1] - 0.3)),
       z_true = z)
}

.iofp_fit <- function(sim) {
  suppressWarnings(tobs(~ occ_cov, data = sim$data, family = int_occu(),
                        detection = ~ det_cov, y = sim$y, method = "laplace",
                        control = list(verbose = FALSE, progress = FALSE)))
}

test_that("fitted() z posterior on an int_occu() fit pools all sources, not the prior", {
  sim <- .iofp_sim(seed = 1)
  fit <- .iofp_fit(sim)
  f <- fitted(fit)

  expect_false(isTRUE(all.equal(f$z, f$psi)))

  any_det <- (rowSums(sim$y$src1 > 0) > 0) | (rowSums(sim$y$src2 > 0) > 0)
  expect_true(all(f$z[any_det] == 1))
  # Sites with no detection anywhere pull posterior occupancy DOWN from the
  # prior marginal (more non-detection evidence than the marginal alone).
  expect_true(mean(f$z[!any_det]) < mean(f$psi[!any_det]))
  expect_true(all(f$z >= 0 & f$z <= 1))

  # nobs() counts every surveyed (site, visit) cell, pooled over the sources.
  expect_identical(nobs(fit),
                   sum(vapply(sim$y, function(m) sum(!is.na(m)), integer(1))))
})

test_that("fitted() z posterior matches a hand-rolled multi-source Bayes update", {
  sim <- .iofp_sim(seed = 2)
  fit <- .iofp_fit(sim)
  f <- fitted(fit)

  psi <- f$psi; p1 <- f$p$src1; p2 <- f$p$src2
  y1 <- sim$y$src1; y2 <- sim$y$src2
  z_manual <- numeric(nrow(y1))
  for (i in seq_along(z_manual)) {
    if (any(y1[i, ] == 1) || any(y2[i, ] == 1)) { z_manual[i] <- 1; next }
    prod_1mp <- (1 - p1[i])^ncol(y1) * (1 - p2[i])^ncol(y2)
    z_manual[i] <- psi[i] * prod_1mp / (psi[i] * prod_1mp + (1 - psi[i]))
  }
  expect_equal(f$z, z_manual, tolerance = 1e-8)
})

test_that("predict() detection design mode works for single-season occu()", {
  sim <- simulate_occu(N = 100, J = 5, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.3, 1.0), beta_det = c(0.7, 0.6), seed = 1)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
              y = sim$y, method = "laplace", control = list(verbose = FALSE))

  X_det <- fit$model$X_processes[[2L]]
  X_occ <- fit$model$X_processes[[1L]]

  pr <- predict(fit, type = "detection", X_det.0 = X_det[1:5, , drop = FALSE])
  expect_s3_class(pr, "data.frame")
  expect_identical(nrow(pr), 5L)
  expect_true(all(pr$mean >= 0 & pr$mean <= 1))

  both <- predict(fit, X.0 = X_occ[1:5, , drop = FALSE], type = "both",
                  X_det.0 = X_det[1:5, , drop = FALSE])
  expect_identical(names(both), c("occupancy", "detection"))
  expect_equal(both$detection, pr)

  # A design matrix for the type NOT requested is not required.
  expect_error(predict(fit, X.0 = X_occ[1:5, , drop = FALSE], type = "both"),
              "X_det.0")
  expect_error(predict(fit, X_det.0 = X_det[1:5, , drop = FALSE], type = "both"),
              "X.0")

  # No design matrix at all is still the in-sample fallback (matches
  # type = "occupancy" with no X.0), not an error.
  expect_identical(predict(fit, type = "detection"), fitted(fit))
})

test_that("predict() detection design mode is per-source for int_occu()", {
  sim <- .iofp_sim(seed = 3)
  fit <- .iofp_fit(sim)
  f <- fitted(fit)

  Xd1 <- fit$model$X_processes[[2L]][1:5, , drop = FALSE]
  Xd2 <- fit$model$X_processes[[3L]][1:5, , drop = FALSE]

  pr <- predict(fit, type = "detection", X_det.0 = list(src1 = Xd1, src2 = Xd2))
  expect_identical(names(pr), c("src1", "src2"))
  expect_identical(nrow(pr$src1), 5L)
  # Matches the in-sample posterior-mean predictor closely (draw-averaged vs
  # plugged-in mean, same relationship as the existing occupancy design mode).
  expect_lt(max(abs(pr$src1$mean - f$p$src1[1:5])), 0.02)
  expect_lt(max(abs(pr$src2$mean - f$p$src2[1:5])), 0.02)

  # Unnamed list in source order works too.
  pr2 <- predict(fit, type = "detection", X_det.0 = list(Xd1, Xd2))
  expect_equal(pr2, pr)

  expect_error(
    predict(fit, type = "detection", X_det.0 = Xd1),
    "list"
  )
  expect_error(
    predict(fit, type = "detection", X_det.0 = list(src1 = Xd1)),
    "name every source|unnamed entries"
  )
})
