# The two ways of asking for a spatially-varying coefficient go through one verb
# (gcol33/tulpaObs#146): the areal flavour is the weighted bar
# `spatial(~ 1 + w || node, graph)`, the continuous flavour is
# `spatial(lon, lat, model = "svc", coefficients = )`. `svc()` stays as the
# direct constructor, the way `icar()` does alongside `spatial(model = "icar")`.
#
# The continuous form also selects its coefficients BY NAME, matched against the
# design of the arm the surfaces load on; `indices = ` (column positions) is the
# lower-level form and keeps working. These are constructor / resolver tests, so
# they run in every tier -- the recovery of the surfaces themselves lives in
# test-occu-svc-laplace-recovery.R and test-svc-families-recovery.R.

pr <- c(0.1, 0.05)

# ---------------------------------------------------------------------------
# Umbrella dispatch
# ---------------------------------------------------------------------------

test_that("spatial(model = \"svc\") builds the same term svc() does", {
  co <- cbind(stats::runif(20), stats::runif(20))
  direct   <- .tobs_term_svc(coords = co, indices = 2L, nn = 5, prior_range = pr)
  umbrella <- .tobs_term_spatial(coords = co, model = "svc", indices = 2L,
                                 nn = 5, prior_range = pr)
  expect_s3_class(umbrella, "tobs_svc")
  # The umbrella stamps `field_name` (NULL here); everything else is identical.
  umbrella$field_name <- NULL
  expect_equal(umbrella, direct)
})

test_that("svc is offered by the spatial umbrella and dispatches by name", {
  expect_true("svc" %in% .tobs_spatial_models)
  co <- cbind(stats::runif(12), stats::runif(12))
  tm <- .tobs_term_spatial(coords = co, model = "svc",
                           coefficients = "elev", prior_range = pr, nn = 4)
  expect_identical(tm$coefficients, "elev")
  expect_identical(tm$n_svc, 1L)
})

test_that("spatial(model = \"svc\") rejects an unknown argument", {
  co <- cbind(stats::runif(12), stats::runif(12))
  expect_error(
    .tobs_term_spatial(coords = co, model = "svc", coefficients = "elev",
                       prior_range = pr, graph = diag(3)),
    "unknown argument")
})

test_that("a bar with model = \"svc\" points at the continuous form", {
  expect_error(
    .tobs_term_spatial(~ 1 + w || cell, graph = matrix(0, 3, 3), model = "svc"),
    "coefficients")
})

# ---------------------------------------------------------------------------
# Coefficient selector
# ---------------------------------------------------------------------------

test_that("svc() requires exactly one coefficient selector", {
  co <- cbind(stats::runif(12), stats::runif(12))
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4),
               "coefficients")
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                              coefficients = "elev", indices = 2L),
               "not both")
})

test_that("svc() rejects a malformed or duplicated selector", {
  co <- cbind(stats::runif(12), stats::runif(12))
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                              coefficients = c("elev", "elev")),
               "duplicate")
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                              coefficients = 2), "names")
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                              indices = c(1L, 1L)), "duplicate")
  expect_error(.tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                              indices = 0L), "positive")
})

test_that("names and positions resolve to the same design columns", {
  co <- cbind(stats::runif(10), stats::runif(10))
  X  <- matrix(0, 10, 3,
               dimnames = list(NULL, c("(Intercept)", "elev", "w")))
  by_name <- .tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                            coefficients = c("(Intercept)", "w"))
  by_pos  <- .tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                            indices = c(1L, 3L))
  expect_identical(.tobs_svc_columns(by_name, X, "occu"), c(1L, 3L))
  expect_identical(.tobs_svc_columns(by_pos, X, "occu"), c(1L, 3L))
})

test_that("an unmatched coefficient name names the available columns", {
  co <- cbind(stats::runif(10), stats::runif(10))
  X  <- matrix(0, 10, 2, dimnames = list(NULL, c("(Intercept)", "elev")))
  tm <- .tobs_term_svc(coords = co, prior_range = pr, nn = 4,
                       coefficients = "elevation")
  err <- expect_error(.tobs_svc_columns(tm, X, "occu"), "elevation")
  expect_match(conditionMessage(err), "elev")
  expect_match(conditionMessage(err), "Available")
})

test_that("an out-of-range position still errors against the design width", {
  co <- cbind(stats::runif(10), stats::runif(10))
  X  <- matrix(0, 10, 2, dimnames = list(NULL, c("(Intercept)", "elev")))
  tm <- .tobs_term_svc(coords = co, prior_range = pr, nn = 4, indices = 5L)
  expect_error(.tobs_svc_columns(tm, X, "occu"), "out of range")
})

# ---------------------------------------------------------------------------
# End to end: a name-selected surface fits and reports its resolved column
# ---------------------------------------------------------------------------

test_that("a name-selected continuous SVC fits occu() and reports its column", {
  skip_if_fast()
  set.seed(11)
  n <- 60L; J <- 5L
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n),
                   w = stats::rnorm(n))
  psi <- stats::plogis(-0.2 + 0.5 * df$w)
  z   <- stats::rbinom(n, 1, psi)
  y   <- matrix(stats::rbinom(n * J, 1, 0.5) * z, n, J)

  fit <- tobs(~ w + spatial(lon, lat, model = "svc", coefficients = "w",
                            nn = 8, prior_range = c(0.2, 0.05)),
              data = df, family = occu(), detection = ~ 1, y = y,
              method = "laplace", verbose = FALSE)

  expect_s3_class(fit$svc, "tobs_svc")
  # "w" is the second occupancy design column (intercept first).
  expect_identical(fit$svc_indices, 2L)
  expect_identical(fit$svc_coefficients, "w")
  expect_length(fit$svc_field, n)
  expect_true(all(is.finite(fit$svc_field)))
})
