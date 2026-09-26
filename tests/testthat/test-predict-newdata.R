# predict(newdata = ) on the families that predict from a design matrix: the
# frame is expanded through the predicted arm's fitted formula, so a call with
# `newdata` equals the same call with the design matrix it expands to.

.pnd_ctl <- list(verbose = FALSE, progress = FALSE)

.pnd_same <- function(a, b) {
  expect_equal(unclass(a), unclass(b), ignore_attr = TRUE)
}

test_that("occu(): newdata predicts per row, population-level with a random effect", {
  set.seed(3); ng <- 12; per <- 30; N <- ng * per; J <- 5
  lab <- sprintf("grp%02d", sample.int(ng)); gi <- rep(seq_len(ng), each = per)
  b <- rnorm(ng); x <- rnorm(N); w <- rnorm(N)
  z <- rbinom(N, 1, plogis(0.2 + 0.5 * x + b[gi]))
  y <- matrix(rbinom(N * J, 1, z * 0.7), N, J)
  d <- data.frame(x = x, w = w, g = factor(lab[gi]))
  fit <- suppressWarnings(tobs(~ x + (1 | g), data = d, y = y,
                               detection = ~ w, family = occu(),
                               method = "laplace", control = .pnd_ctl))
  nd <- data.frame(x = c(0, 1), w = c(0, 0.5),
                   g = factor(c(lab[1], "NEVER_SEEN"),
                              levels = c(lab, "NEVER_SEEN")))

  occ <- predict(fit, newdata = nd, type = "occupancy")
  expect_equal(nrow(occ), 2L)
  .pnd_same(occ, predict(fit, X.0 = cbind(1, nd$x), type = "occupancy"))
  # Positional data frame takes the same route.
  .pnd_same(predict(fit, nd), occ)

  both <- predict(fit, newdata = nd, type = "both")
  .pnd_same(both$occupancy, occ)
  .pnd_same(both$detection,
            predict(fit, X_det.0 = cbind(1, nd$w), type = "detection"))

  expect_error(predict(fit, newdata = nd, X.0 = cbind(1, 0)), "not several")
  expect_error(predict(fit, newdata = nd, terms = "x"), "not several")
  expect_error(predict(fit, newdata = nd, type = "state"), "in-sample")
  expect_error(predict(fit, newdata = data.frame(w = 1)), "covariate(s) x", fixed = TRUE)
})

test_that("occu(): newdata re-expands a factor on the fitted levels", {
  set.seed(11); N <- 200; J <- 4
  hab <- factor(sample(c("forest", "meadow", "scrub"), N, TRUE))
  z <- rbinom(N, 1, plogis(c(forest = 0.8, meadow = -0.4, scrub = 0.2)[hab]))
  y <- matrix(rbinom(N * J, 1, z * 0.6), N, J)
  fit <- suppressWarnings(tobs(~ hab, data = data.frame(hab = hab), y = y,
                               detection = ~ 1, family = occu(),
                               method = "laplace", control = .pnd_ctl))
  # A one-level frame still gets the fitted three-level contrast columns.
  nd <- data.frame(hab = "scrub")
  .pnd_same(predict(fit, newdata = nd),
            predict(fit, X.0 = cbind(1, 0, 1)))
})

test_that("int_occu(): newdata detection expands once per source", {
  set.seed(5); n <- 80
  oc <- rnorm(n); dc <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 + 0.7 * oc))
  mk <- function(J) matrix(rbinom(n * J, 1, z * plogis(0.3 + 0.5 * dc)), n, J)
  fit <- suppressWarnings(tobs(~ occ_cov, data = data.frame(occ_cov = oc,
                                                            det_cov = dc),
                               family = int_occu(), detection = ~ det_cov,
                               y = list(src1 = mk(4), src2 = mk(3)),
                               method = "laplace", control = .pnd_ctl))
  nd <- data.frame(occ_cov = c(-1, 1), det_cov = c(0, 2))
  det <- predict(fit, newdata = nd, type = "detection")
  expect_named(det, c("src1", "src2"))
  X <- cbind(1, nd$det_cov)
  .pnd_same(det, predict(fit, X_det.0 = list(src1 = X, src2 = X),
                         type = "detection"))
})

test_that("abun(): newdata predicts abundance and detection per row", {
  sim <- simulate_abun(seed = 1)
  fit <- suppressWarnings(tobs(~ abund_cov1 + abund_cov2, data = sim$data,
                               y = sim$y, detection = ~ det_cov1,
                               family = abun(), control = .pnd_ctl))
  nd <- sim$data[c(3, 7), , drop = FALSE]
  lam <- predict(fit, newdata = nd)
  expect_equal(nrow(lam), 2L)
  .pnd_same(lam, predict(fit, X.0 = cbind(1, nd$abund_cov1, nd$abund_cov2)))
  .pnd_same(predict(fit, newdata = nd, type = "detection"),
            predict(fit, X.0 = cbind(1, nd$det_cov1), type = "detection"))
})

test_that("distance() / fp_occu(): newdata reaches the design-matrix predictor", {
  sim <- simulate_distance(N = 80, seed = 18L)
  fd <- suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                              family = distance(cutpoints = sim$cutpoints,
                                                key = "halfnorm",
                                                transect = "line"),
                              detection = ~ sigma_cov1, y = sim$y,
                              method = "laplace", control = .pnd_ctl))
  nd <- sim$data[1:3, , drop = FALSE]
  .pnd_same(predict(fd, newdata = nd),
            predict(fd, X.0 = cbind(1, nd$abund_cov1)))
  .pnd_same(predict(fd, newdata = nd, type = "sigma"),
            predict(fd, X.0 = cbind(1, nd$sigma_cov1), type = "sigma"))

  sf <- simulate_fp_occu(N = 120, J = 4, seed = 7)
  ff <- suppressWarnings(tobs(~ occ_cov1, data = sf$data, family = fp_occu(),
                              detection = ~ 1, y = sf$y, method = "laplace",
                              control = .pnd_ctl))
  nf <- sf$data[1:2, , drop = FALSE]
  pf <- predict(ff, newdata = nf)
  expect_equal(nrow(pf), 2L)
  .pnd_same(pf, predict(ff, X.0 = cbind(1, nf$occ_cov1)))
})
