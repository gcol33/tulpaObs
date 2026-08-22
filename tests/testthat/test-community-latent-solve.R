# =============================================================================
# test-community-latent-solve.R -- one ridge ladder behind every
# singular-Hessian retry in R/community_latent.R, and the field Newton holds a
# non-finite step instead of applying it.
#
# The file used to carry the ladder three times with different ridges
# (relative vs absolute), different tier counts and different failure values.
# The relative one guarded `is.finite(d) && d > 0` on its second tier and then
# reused the same `d` unguarded on its third, so a non-finite mean diagonal
# gave `max(NaN, 1e-6)` = NaN as the ridge -- in exactly the near-singular case
# the retry exists for.
#
# Separately, the field Newton applied its step unconditionally: a non-finite
# step poisoned `F`, made `max(abs(step)) < tol` NA so the loop never broke
# early, and reached the tau M-step with no error and no convergence flag.
# =============================================================================

test_that("the ridge ladder solves, retries and never builds a NaN ridge", {
  H <- Matrix::Diagonal(4, c(2, 3, 4, 5))
  g <- c(1, 1, 1, 1)
  expect_equal(.tobs_ridge_solve(H, g), as.numeric(solve(as.matrix(H), g)))
  # The inverse form (no `g`) is what the covariance solve uses.
  expect_equal(as.matrix(.tobs_ridge_solve(H)), solve(as.matrix(H)))

  # A singular H must still return a finite step rather than throw.
  Hs <- Matrix::Matrix(matrix(1, 2, 2), sparse = TRUE)
  v <- .tobs_ridge_solve(Hs, c(1, 1))
  expect_true(all(is.finite(v)))
  expect_true(all(is.finite(as.matrix(.tobs_ridge_solve(Hs)))))

  # The ridge itself is finite whatever the diagonal: a NaN mean diagonal falls
  # back to a unit scale rather than producing NaN * 1e-4.
  src <- paste(deparse(body(.tobs_ridge_solve)), collapse = " ")
  expect_true(grepl("is.finite(d) && d > 0", src, fixed = TRUE))
  expect_false(grepl("max(1e-4 * d,", src, fixed = TRUE))
})

test_that("the field Newton holds a non-finite step", {
  # The guard reads on the step, before it reaches F, and breaks the loop.
  src <- paste(deparse(.tobs_latent_field_solve), collapse = " ")
  i_guard <- regexpr("all(is.finite(step))", src, fixed = TRUE)
  i_apply <- regexpr("F + matrix(step", src, fixed = TRUE)
  expect_gt(i_guard, 0)
  expect_gt(i_apply, i_guard)
})

test_that("the three retries route through the one ladder", {
  for (fn in c(".tobs_latent_field_solve", ".tobs_latent_factor_update")) {
    src <- paste(deparse(get(fn, envir = asNamespace("tulpaObs"))),
                 collapse = " ")
    expect_true(grepl(".tobs_ridge_solve", src, fixed = TRUE), info = fn)
    # No second ladder left behind.
    expect_false(grepl("diag(1e-6, nrow(H))", src, fixed = TRUE), info = fn)
  }
})
