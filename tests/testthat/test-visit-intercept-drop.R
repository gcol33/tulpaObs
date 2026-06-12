# The visit-level design matrix is stacked onto a site-level detection arm that
# already carries an intercept, so the shared builder .tobs_build_visit_X() drops
# the visit (Intercept) to keep the combined design full rank. occu(), abun() and
# removal() previously inlined their own builders and kept the visit intercept;
# these tests lock the now-shared convention so all four families agree.

test_that(".tobs_build_visit_X drops the visit intercept by default", {
  vd <- data.frame(x = c(1, 2, 3, 4))            # 2 units x 2 slots

  X1 <- tulpaObs:::.tobs_build_visit_X(~ x, vd, 2, 2, "detection")
  expect_equal(colnames(X1), "x")
  expect_false("(Intercept)" %in% colnames(X1))

  # `~ x` (intercept dropped) and `~ 0 + x` describe the same visit design.
  X0 <- tulpaObs:::.tobs_build_visit_X(~ 0 + x, vd, 2, 2, "detection")
  expect_equal(as.vector(X1), as.vector(X0))

  # Opt out and the intercept stays.
  Xk <- tulpaObs:::.tobs_build_visit_X(~ x, vd, 2, 2, "detection",
                                       drop_intercept = FALSE)
  expect_true("(Intercept)" %in% colnames(Xk))

  # An intercept-only visit formula adds nothing beyond the site intercept.
  expect_null(tulpaObs:::.tobs_build_visit_X(~ 1, vd, 2, 2, "detection"))

  # Absent formula or data short-circuits to NULL.
  expect_null(tulpaObs:::.tobs_build_visit_X(NULL, vd, 2, 2, "detection"))
  expect_null(tulpaObs:::.tobs_build_visit_X(~ x, NULL, 2, 2, "detection"))
})

test_that(".tobs_build_visit_X zero-fills NA cells and checks the row count", {
  vdna <- data.frame(x = c(1, NA, 3, 4))
  Xna  <- tulpaObs:::.tobs_build_visit_X(~ x, vdna, 2, 2, "detection")
  expect_equal(Xna[2, 1], 0)

  expect_error(
    tulpaObs:::.tobs_build_visit_X(~ x, vdna, 3, 2, "detection"),
    "one row per unit-visit"
  )
})

test_that("occu / abun / removal drop the visit intercept (no double intercept)", {
  set.seed(1)
  n <- 6L; J <- 3L
  d  <- data.frame(z = rnorm(n))
  vd <- data.frame(x = rnorm(n * J))             # unit-major, n * J rows

  yb     <- matrix(rbinom(n * J, 1L, 0.5), n, J)
  m_occu <- tulpaObs:::.tobs_build_single(~ z, ~ 1, d, yb, ~ x, vd)
  expect_equal(m_occu$det_visit_names, "x")
  expect_false("(Intercept)" %in% m_occu$det_visit_names)

  yc     <- matrix(rpois(n * J, 2L), n, J)
  m_abun <- tulpaObs:::.tobs_build_abun(~ z, ~ 1, d, yc, ~ x, vd)
  dn_a   <- m_abun$process_info[[2L]]$coef_names
  expect_equal(sum(dn_a == "(Intercept)"), 1L)   # only the site arm's intercept
  expect_true("x" %in% dn_a)

  yr     <- matrix(rpois(n * J, 1L), n, J)        # complete: removal forbids NA
  m_rem  <- tulpaObs:::.tobs_build_removal(~ z, ~ 1, d, yr, ~ x, vd)
  dn_r   <- m_rem$process_info[[2L]]$coef_names
  expect_equal(sum(dn_r == "(Intercept)"), 1L)
  expect_true("x" %in% dn_r)
})
