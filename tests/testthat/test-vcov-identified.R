# .tobs_vcov_identified(): the covariance of a fit whose observed information
# is singular or indefinite along some direction keeps the standard errors of
# the coefficients that do not load on it, and names the ones that do (#449).

test_that("a positive-definite information is inverted unchanged, silently", {
  H <- matrix(c(4, 1, 0, 1, 3, 0.5, 0, 0.5, 2), 3L,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  expect_silent(vc <- tulpaObs:::.tobs_vcov_identified(H))
  expect_equal(unname(vc$vcov), unname(solve(H)))
  expect_identical(vc$vcov, vc$vcov_draws)
  expect_length(vc$unidentified, 0L)
})

test_that("a ridge between two coefficients leaves the third's SE intact", {
  # b and c enter only through b + c: the information is singular along
  # (0, 1, -1) / sqrt(2), and a is resolved on its own.
  H <- matrix(0, 3L, 3L, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  H["a", "a"] <- 4
  H[c("b", "c"), c("b", "c")] <- 2
  expect_warning(vc <- tulpaObs:::.tobs_vcov_identified(H, what = "f()"),
                 "f\\(\\).*loading on b, c")
  expect_identical(vc$unidentified, c("b", "c"))
  expect_equal(vc$vcov["a", "a"], 1 / 4)
  expect_true(all(is.na(vc$vcov[c("b", "c"), ])))
  expect_true(all(is.na(vc$vcov[, c("b", "c")])))
  # The draws' covariance is finite: the resolved direction b + c at its
  # information, the unresolved one at the prior variance.
  vc2 <- suppressWarnings(
    tulpaObs:::.tobs_vcov_identified(H, prior_var = 25))
  Vd <- vc2$vcov_draws
  expect_true(all(is.finite(Vd)))
  u0 <- c(0, 1, -1) / sqrt(2)
  u1 <- c(0, 1, 1) / sqrt(2)
  expect_equal(drop(t(u0) %*% Vd %*% u0), 25)
  expect_equal(drop(t(u1) %*% Vd %*% u1), 1 / 4)
})

test_that("an indefinite information (a stop off the maximum) is flagged", {
  H <- diag(c(5, 2, -0.3))
  dimnames(H) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_warning(vc <- tulpaObs:::.tobs_vcov_identified(H), "loading on c")
  expect_equal(diag(vc$vcov)[c("a", "b")], c(a = 1 / 5, b = 1 / 2))
  expect_true(is.na(vc$vcov["c", "c"]))
})

test_that("a non-finite information reports every SE as NA, with a warning", {
  H <- diag(2); H[1, 1] <- NaN
  dimnames(H) <- list(c("a", "b"), c("a", "b"))
  expect_warning(vc <- tulpaObs:::.tobs_vcov_identified(H), "loading on a, b")
  expect_true(all(is.na(vc$vcov)))
})
