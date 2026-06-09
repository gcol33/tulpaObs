## Response on the top formula LHS (gcol33/tulpaObs#66).
##
## A single-vector-response family (the cover hurdle) may carry its response on
## the formula left-hand side and drop `y =`; the LHS form must produce the same
## fit as the equivalent `y =` call. Matrix / array / list response families
## keep `y =`, and a two-sided formula for those errors. The validation paths
## (response given twice; matrix family with a two-sided formula) run always; the
## equality-of-fits check is a real fit and is gated for fast iteration.

simulate_lhs_cover <- function(N = 600, beta_occ = c(-0.4, 0.8),
                               beta_pos = c(0.5, -1.2), phi = 30, seed = 7) {
  set.seed(seed)
  x       <- runif(N, -2, 2)
  occur   <- rbinom(N, 1, plogis(beta_occ[1] + beta_occ[2] * x))
  mu      <- plogis(beta_pos[1] + beta_pos[2] * x)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos]  <- rbeta(sum(is_pos), mu[is_pos] * phi, (1 - mu[is_pos]) * phi)
  y[!is_pos] <- 0
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  data.frame(cover = y, x = x)
}


test_that("cover(): LHS-response fit equals the y= fit (same coefficients)", {
  skip_if_fast()
  skip_on_cran()

  dat <- simulate_lhs_cover(N = 600, seed = 7)

  fit_lhs <- tobs(cover ~ x, data = dat, family = cover("beta"))
  fit_y   <- tobs(~ x, data = dat, family = cover("beta"), y = dat$cover)

  expect_s3_class(fit_lhs, "cover_fit")
  # The LHS form is exactly the y= form: identical coefficients and dispersion.
  expect_equal(fit_lhs$beta_occ, fit_y$beta_occ)
  expect_equal(fit_lhs$beta_pos, fit_y$beta_pos)
  expect_equal(fit_lhs$phi_pos,  fit_y$phi_pos)
})


test_that("cover(): an LHS expression is evaluated against data", {
  skip_if_fast()
  skip_on_cran()

  dat <- simulate_lhs_cover(N = 400, seed = 3)

  # A bare-expression LHS (here the identity, to keep cover in [0, 1]) resolves
  # against `data`, matching the equivalent y = I(cover) call.
  fit_expr <- tobs(I(cover) ~ x, data = dat, family = cover("beta"))
  fit_y    <- tobs(~ x, data = dat, family = cover("beta"), y = dat$cover)

  expect_equal(fit_expr$beta_occ, fit_y$beta_occ)
  expect_equal(fit_expr$beta_pos, fit_y$beta_pos)
})


test_that("cover(): two-sided formula plus y= errors (response given twice)", {
  dat <- simulate_lhs_cover(N = 50, seed = 1)
  expect_error(
    tobs(cover ~ x, data = dat, family = cover("beta"), y = dat$cover),
    "response is given twice"
  )
})


test_that("occu(): a two-sided formula errors and points to y=", {
  # A matrix-response family rejects an LHS response before any fit.
  y <- matrix(rbinom(30, 1, 0.4), nrow = 10, ncol = 3)
  dat <- data.frame(elev = rnorm(10))
  expect_error(
    tobs(elev_response ~ elev, data = dat, family = occu(),
         detection = ~ 1, y = y),
    "matrix / array supplied via `y =`"
  )
})


test_that("cover(): one-sided formula with y= still works (smoke)", {
  skip_if_fast()
  skip_on_cran()

  dat <- simulate_lhs_cover(N = 200, seed = 9)
  fit <- tobs(~ x, data = dat, family = cover("beta"), y = dat$cover)
  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)
})


test_that("family$response declares the response kind (single source of truth)", {
  expect_identical(cover()$response, "vector")
  expect_identical(occu()$response,  "matrix")
  expect_identical(abun()$response,  "matrix")
  expect_identical(occu_cover()$response, "matrix")
})
