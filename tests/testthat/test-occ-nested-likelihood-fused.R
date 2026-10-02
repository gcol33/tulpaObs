# The marginalized single-season occupancy LikelihoodSpec (occ_make_nested_
# likelihood) fills tulpa's fused slot ll_eta_weights_fn. The fused callback has
# to return what ll_double returns and write what eta_weights_fn writes, and the
# weights have to be the density's own derivatives: the score is d ll / d eta
# and the working weight is the expected information q s (1 - s)^2 / (1 - q s).

.occ_fused_eval <- function(y, q, eta) {
  tulpaObs:::occ_nested_likelihood_eval(as.numeric(y), as.numeric(q),
                                        as.numeric(eta))
}

test_that("fused occupancy callback is the split pair, bit for bit", {
  set.seed(41)
  n <- 400L
  eta <- c(stats::rnorm(n - 6L, 0, 3), -40, 40, 0, 1e-8, -1e-8, 12)
  q   <- c(stats::runif(n - 6L), 0, 1, 1, 0.5, 0.999999, 0)
  y   <- stats::rbinom(n, 1L, 0.5)
  ev  <- .occ_fused_eval(y, q, eta)

  expect_identical(ev[, "ll_fused"],       ev[, "ll_split"])
  expect_identical(ev[, "grad_fused"],     ev[, "grad_split"])
  expect_identical(ev[, "neg_hess_fused"], ev[, "neg_hess_split"])
  expect_true(all(is.finite(ev)))
  # a site with no valid visit contributes nothing
  none <- q == 0
  expect_true(all(ev[none, ] == 0))
})

test_that("fused occupancy score matches finite differences of the value", {
  set.seed(42)
  n   <- 300L
  eta <- stats::rnorm(n, 0, 2)
  q   <- stats::runif(n, 0.05, 1)
  y   <- stats::rbinom(n, 1L, 0.5)
  h   <- 1e-6
  up  <- .occ_fused_eval(y, q, eta + h)[, "ll_fused"]
  dn  <- .occ_fused_eval(y, q, eta - h)[, "ll_fused"]
  ev  <- .occ_fused_eval(y, q, eta)
  fd  <- (up - dn) / (2 * h)
  expect_equal(ev[, "grad_fused"], fd, tolerance = 1e-7)

  # the working weight is the expected information of the scaled Bernoulli
  s  <- stats::plogis(eta)
  mu <- q * s
  expect_equal(ev[, "neg_hess_fused"], q * s * (1 - s)^2 / (1 - mu),
               tolerance = 1e-12)
  # q = 1 is a plain logit Bernoulli: score y - s, information s (1 - s)
  ev1 <- .occ_fused_eval(y, rep(1, n), eta)
  expect_equal(ev1[, "grad_fused"], y - s, tolerance = 1e-12)
  expect_equal(ev1[, "neg_hess_fused"], s * (1 - s), tolerance = 1e-12)
  expect_equal(ev1[, "ll_fused"], stats::dbinom(y, 1L, s, log = TRUE),
               tolerance = 1e-12)
})
