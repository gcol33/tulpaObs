# The quadrature behind the community factor magnitude (R/community_latent.R,
# gcol33/tulpaObs#153). The factor update fixes the loadings only up to one
# overall magnitude; that magnitude is set by the JOINT site marginal
#
#   L_i = integral prod_s f(y_is | eta_is + lambda_s' z) N(z; 0, I_Q) dz
#
# so these pin the pieces it is built from against independent references, rather
# than leaving them to be inferred from a downstream recovery run.

test_that("Gauss-Hermite nodes integrate the standard normal", {
  for (n in c(3L, 5L, 7L, 15L)) {
    gh <- tulpaObs:::.tobs_gh_nodes(n)
    expect_length(gh$x, n)
    expect_equal(sum(gh$w), 1, tolerance = 1e-12)
    # a rule with n nodes is exact to degree 2n - 1, so the moments it claims to
    # reproduce are the odd ones (0) and the even ones (1, 3, 15, ...)
    expect_equal(sum(gh$w * gh$x), 0, tolerance = 1e-10)
    expect_equal(sum(gh$w * gh$x^2), 1, tolerance = 1e-10)
    expect_equal(sum(gh$w * gh$x^3), 0, tolerance = 1e-10)
    if (n >= 3L) expect_equal(sum(gh$w * gh$x^4), 3, tolerance = 1e-10)
    if (n >= 5L) expect_equal(sum(gh$w * gh$x^6), 15, tolerance = 1e-9)
  }
})


test_that("the joint site marginal matches direct numerical integration", {
  # Q = 1 so the truth is a one-dimensional integral stats::integrate() can do.
  set.seed(4)
  Ns <- 6L; S <- 4L
  y   <- matrix(rbinom(Ns * S, 1, 0.5), Ns, S)
  eta <- matrix(rnorm(Ns * S, 0, 0.7), Ns, S)
  lam <- matrix(rnorm(S, 0, 0.9), S, 1L)
  or  <- tulpaObs:::.tobs_ms_count_oracle(y, link = "logit")

  # site i: integral prod_s Bern(y_is | plogis(eta_is + lam_s z)) phi(z) dz
  ref <- sum(vapply(seq_len(Ns), function(i) {
    f <- function(z) vapply(z, function(zz) {
      p <- stats::plogis(eta[i, ] + as.numeric(lam) * zz)
      prod(ifelse(y[i, ] > 0, p, 1 - p)) * stats::dnorm(zz)
    }, numeric(1))
    log(stats::integrate(f, -12, 12, rel.tol = 1e-10)$value)
  }, numeric(1)))

  gh <- tulpaObs:::.tobs_gh_nodes(25L)
  expect_equal(tulpaObs:::.tobs_latent_joint_marginal(or, eta, lam, gh), ref,
               tolerance = 1e-7)

  # Q = 2 against a dense tensor reference (no closed form to appeal to). The
  # tensor rule converges more slowly than the one-dimensional one, so assert a
  # relative tolerance and let the scale test below carry the question that
  # actually matters -- whether the node count moves the ARGMAX.
  lam2 <- matrix(rnorm(S * 2L, 0, 0.8), S, 2L)
  dense <- tulpaObs:::.tobs_latent_joint_marginal(or, eta, lam2,
                                                  tulpaObs:::.tobs_gh_nodes(35L))
  expect_equal(tulpaObs:::.tobs_latent_joint_marginal(
    or, eta, lam2, tulpaObs:::.tobs_gh_nodes(9L)), dense, tolerance = 1e-4)
})


test_that("the factor scale recovers a known inflation of the loadings", {
  skip_if_fast()
  skip_on_cran()
  # Simulate at a known loading matrix, hand the scale step lambda INFLATED by a
  # known factor, and check it discounts by roughly that factor. This is the
  # defect in #153 in miniature: the update lands on an over-large magnitude and
  # the joint marginal is what pulls it back.
  set.seed(11)
  Ns <- 400L; S <- 14L; Q <- 2L
  lam  <- matrix(rnorm(S * Q, 0, 0.8), S, Q)
  zeta <- matrix(rnorm(Ns * Q), Ns, Q)
  eta  <- matrix(rnorm(Ns * S, 0, 0.5), Ns, S)
  y    <- matrix(rbinom(Ns * S, 1, stats::plogis(eta + tcrossprod(zeta, lam))),
                 Ns, S)
  or <- tulpaObs:::.tobs_ms_count_oracle(y, link = "logit")
  gh <- tulpaObs:::.tobs_gh_nodes(5L)
  for (infl in c(1.5, 2)) {
    cc <- tulpaObs:::.tobs_latent_factor_scale(or, eta, lam * infl, gh)
    expect_equal(cc, 1 / infl, tolerance = 0.35)
  }
  # a degenerate (all-zero) loading matrix has no magnitude to set
  expect_identical(
    tulpaObs:::.tobs_latent_factor_scale(or, eta, matrix(0, S, Q), gh), 1)

  # The driver runs this on a 5-node rule. The raw tensor integral is not
  # converged there, but the scale is the ARGMAX over a grid, and the quadrature
  # error is smooth in c and largely cancels: the estimate has to be stable in the
  # node count, which is what licenses the cheap default.
  ref <- tulpaObs:::.tobs_latent_factor_scale(or, eta, lam * 1.5,
                                              tulpaObs:::.tobs_gh_nodes(21L))
  for (n in c(5L, 9L, 15L)) {
    expect_equal(tulpaObs:::.tobs_latent_factor_scale(
      or, eta, lam * 1.5, tulpaObs:::.tobs_gh_nodes(n)), ref, tolerance = 0.01)
  }
})
