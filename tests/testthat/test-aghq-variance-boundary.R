# Testing a scalar community variance component against its lower boundary
# (gcol33/tulpaObs#250 item 3).
#
# A variance component that has collapsed toward zero is invisible from the fit
# it arrives on: the optimizer converges, the point estimate is ordinary, and
# nothing warns. The gate is the component's own uncertainty rather than a cut
# on `sigma_hat`, which is a number with nothing behind it and does not transfer
# between fixtures.
#
# A 1x1 covariance block's integration coordinate is log(sigma), so
# tulpa_re_aghq()'s `re_par_se` is SE(log sigma) directly, and the Wald
# statistic for H0: sigma = 0 reduces to 1 / SE(log sigma) by the delta method.
# The critical value is the ordinary one-sided normal quantile because sigma = 0
# is a boundary point, where the one-sided null is 0.5 chi^2_0 + 0.5 chi^2_1.
#
# These drive the helper on constructed engine returns. The arithmetic is what
# is under test, so nothing here fits a model.

# One engine return carrying a single scalar covariance block.
vb_ref <- function(sigma, se_log, nc = 1L, full = FALSE) {
  list(Sigma_list = list(matrix(sigma^2, 1L, 1L)),
       re_par_se = se_log,
       re_par_layout = list(list(label = "logr", nc = nc, full = full,
                                 index = 1L, coord = "log_sd_1")))
}

test_that("the statistic is 1 / SE(log sigma), free of sigma's own scale", {
  a <- .tobs_aghq_variance_boundary(vb_ref(0.5, 0.25), 1L)
  b <- .tobs_aghq_variance_boundary(vb_ref(5.0, 0.25), 1L)

  expect_true(a$available)
  expect_equal(a$statistic, 1 / 0.25)
  # sigma_hat cancels, so two fits differing only in scale get the same verdict.
  # This is what makes the gate transferable between fixtures.
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$distinguishable, b$distinguishable)
})

test_that("the statistic equals sigma / SE(sigma) computed the long way", {
  # The delta method the helper's derivation relies on, done independently:
  # SE(sigma) = sigma * SE(log sigma), so sigma / SE(sigma) must reproduce it.
  sigma <- 0.37; se_log <- 0.42
  expect_equal(.tobs_aghq_variance_boundary(vb_ref(sigma, se_log), 1L)$statistic,
               sigma / (sigma * se_log))
})

test_that("the critical value is the boundary-aware one-sided quantile", {
  r <- .tobs_aghq_variance_boundary(vb_ref(0.5, 0.25), 1L)
  # For a one-sided W >= 0 under 0.5 chi^2_0 + 0.5 chi^2_1, P(W > c) = 1 - Phi(c),
  # so the boundary value and the ordinary one-sided normal quantile coincide.
  expect_equal(r$critical, stats::qnorm(1 - .TOBS_VC_BOUNDARY_ALPHA))
  expect_equal(r$alpha, .TOBS_VC_BOUNDARY_ALPHA)
  # Not the two-sided value, which is the easy thing to write by mistake.
  expect_false(isTRUE(all.equal(r$critical, stats::qnorm(1 - .TOBS_VC_BOUNDARY_ALPHA / 2))))
})

test_that("the verdict flips at the critical value and follows alpha", {
  crit <- stats::qnorm(1 - .TOBS_VC_BOUNDARY_ALPHA)
  expect_true(.tobs_aghq_variance_boundary(vb_ref(1, 1 / crit), 1L)$distinguishable)
  expect_false(.tobs_aghq_variance_boundary(vb_ref(1, 1.001 / crit), 1L)$distinguishable)

  # A stricter level demands more curvature, so a fit admitted at 5% can be
  # refused at 0.1%. The level is an argument, not a constant in the body.
  r <- vb_ref(1, 0.5)
  expect_true(.tobs_aghq_variance_boundary(r, 1L, alpha = 0.05)$distinguishable)
  expect_false(.tobs_aghq_variance_boundary(r, 1L, alpha = 0.001)$distinguishable)
})

test_that("every decline returns the same shape, never a silent verdict", {
  flds <- c("sigma", "se_log", "statistic", "critical", "alpha",
            "distinguishable", "available", "reason")
  cases <- list(
    no_re_par_layout = list(),
    block_not_scalar = vb_ref(0.5, 0.25, nc = 2L),
    no_re_par_se     = vb_ref(0.5, 0.25)[c("Sigma_list", "re_par_layout")],
    curvature_unavailable = vb_ref(0.5, NA_real_))
  for (nm in names(cases)) {
    r <- .tobs_aghq_variance_boundary(cases[[nm]], 1L)
    expect_identical(names(r), flds, info = nm)
    expect_false(r$available, info = nm)
    # NA, never FALSE: a block the curvature could not be read for has not been
    # shown to be at its boundary, and reporting FALSE would warn about it.
    expect_true(is.na(r$distinguishable), info = nm)
    expect_identical(r$reason, nm, info = nm)
  }
  # A zero or negative log-scale SE is curvature that was not read either.
  expect_false(.tobs_aghq_variance_boundary(vb_ref(0.5, 0), 1L)$available)
  # A block index past the layout declines rather than subscripting out of range.
  expect_false(.tobs_aghq_variance_boundary(vb_ref(0.5, 0.25), 2L)$available)
})

test_that("the warning fires once and names every component that failed", {
  recs <- list(sigma_log_r = .tobs_aghq_variance_boundary(vb_ref(0.05, 2.0), 1L),
               sigma_omega = .tobs_aghq_variance_boundary(vb_ref(0.02, 3.0), 1L))
  msg <- tryCatch({ .tobs_warn_variance_boundary(recs); NULL },
                  warning = function(w) conditionMessage(w))
  expect_true(is.character(msg))
  expect_match(msg, "sigma_log_r")
  expect_match(msg, "sigma_omega")
  expect_match(msg, "not distinguishable from zero")
  # One warning for both, so a reader does not correlate two.
  expect_warning(.tobs_warn_variance_boundary(recs), NULL)
})

test_that("a determined component, and an undecidable one, stay silent", {
  expect_silent(.tobs_warn_variance_boundary(
    list(sigma_log_r = .tobs_aghq_variance_boundary(vb_ref(0.5, 0.20), 1L))))
  # available = FALSE is not evidence of collapse and must not warn.
  expect_silent(.tobs_warn_variance_boundary(
    list(sigma_log_r = .tobs_aghq_variance_boundary(list(), 1L))))
  expect_silent(.tobs_warn_variance_boundary(list()))
})
