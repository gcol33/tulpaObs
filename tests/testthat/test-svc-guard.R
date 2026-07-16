# svc() (the continuous NNGP spatially-varying coefficient) is consumed only by
# the single-season occupancy NUTS path. On every other family it used to be
# extracted and silently dropped (gcol33/tulpaObs#118) -- the model fit with no
# spatially-varying coefficient, no error. It must now error with a pointer to
# the recovery-tested areal-bar route. These are fast error paths (the guard
# fires before any fit), so they run in every tier.
#
# `prior_range` is supplied even though these calls are expected to fail: without
# it svc() errors in its own constructor, and that message also matches the regex
# below -- the tests would pass while never reaching the family guard.

test_that("svc() on abun() errors instead of silently dropping (tulpaObs#118)", {
  set.seed(1)
  y  <- matrix(stats::rpois(30 * 3, 3), 30, 3)
  df <- data.frame(lon = stats::runif(30), lat = stats::runif(30))
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)), data = df,
         family = abun(), detection = ~ 1, y = y),
    "svc|areal|spatially-varying")
})

test_that("svc() on occu() under laplace errors (would be silently dropped)", {
  set.seed(2)
  y  <- matrix(stats::rbinom(30 * 3, 1, 0.4), 30, 3)
  df <- data.frame(lon = stats::runif(30), lat = stats::runif(30))
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.1, 0.05)), data = df,
         family = occu(), detection = ~ 1, y = y, method = "laplace"),
    "svc|areal|spatially-varying")
})
