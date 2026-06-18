# Single source of truth for response / site / visit input handling (R/inputs.R):
# the site-count cross-check every family used to hand-roll, the canonical totals
# reader, and the tobs_data frame standing in for the (data, y, visits) triple.

test_that(".tobs_check_site_count is the single validation policy", {
  # Matrix families count "rows" (preserves the historical occu message).
  expect_error(tulpaObs:::.tobs_check_site_count(3, 5, "rows"),
               "y has 3 rows but data has 5 rows")
  # 3D / community families count "sites".
  expect_error(tulpaObs:::.tobs_check_site_count(3, 5, "sites"),
               "y has 3 sites but data has 5 rows")
  # Cover hurdle vector counts "values".
  expect_error(tulpaObs:::.tobs_check_site_count(3, 5, "values"),
               "y has 3 values but data has 5 rows")
  # Match -> silent.
  expect_silent(tulpaObs:::.tobs_check_site_count(5, 5))
})

test_that("the matrix family check still fires through tobs()", {
  expect_error(
    tobs(~ 1, data.frame(x = 1:5), family = occu(),
         detection = ~ 1, y = matrix(0, 3, 2)),
    "y has 3 rows but data has 5 rows"
  )
})

test_that(".tobs_input_dims reads totals from the response shape", {
  # Matrix: sites x visits, no sources.
  d <- tulpaObs:::.tobs_input_dims(matrix(0, 40, 3))
  expect_identical(d$n_sites, 40L)
  expect_identical(d$max_visits, 3L)
  expect_true(is.na(d$n_sources))

  # 3D array: leading two dims are sites x visits.
  d3 <- tulpaObs:::.tobs_input_dims(array(0, dim = c(12, 4, 6)))
  expect_identical(d3$n_sites, 12L)
  expect_identical(d3$max_visits, 4L)

  # List of per-source matrices: source count + first source's shape.
  dl <- tulpaObs:::.tobs_input_dims(list(matrix(0, 30, 4), matrix(0, 20, 3)))
  expect_identical(dl$n_sources, 2L)
  expect_identical(dl$n_sites, 30L)

  # Bare vector (cover hurdle): sites only.
  dv <- tulpaObs:::.tobs_input_dims(c(0.1, 0.2, 0.3))
  expect_identical(dv$n_sites, 3L)
  expect_true(is.na(dv$max_visits))

  expect_true(is.na(tulpaObs:::.tobs_input_dims(NULL)$n_sites))
})

test_that("a tobs_data frame stands in for (data, y, visits)", {
  sim <- simulate_occu(N = 30, J = 3, seed = 7)
  fr  <- tobs_format(y = sim$y, occ.covs = sim$data)

  # Frame carries the response; passing it again via y =/visits = is ambiguous.
  expect_error(
    tobs(~ occ_cov1, fr, family = occu(), detection = ~ 1, y = sim$y),
    "not `y =`"
  )
  expect_error(
    tobs(~ occ_cov1, fr, family = occu(), detection = ~ 1,
         visits = list(a = matrix(0, 30, 3))),
    "not `visits =`"
  )

  # A frame with no response errors with a clear pointer.
  empty <- structure(list(occ.covs = sim$data), class = "tobs_data")
  expect_error(
    tobs(~ occ_cov1, empty, family = occu(), detection = ~ 1),
    "no response"
  )
})

test_that("frame input fits identically to raw input", {
  skip_on_cran(); skip_if_fast()
  sim <- simulate_occu(N = 40, J = 3, seed = 11)

  # The laplace EM draws pseudo-samples for its vcov from the global RNG, so a
  # fixed seed before each fit isolates the comparison to the inputs: a frame and
  # the raw triple it unpacks to must reach the fitter byte-identical.
  set.seed(1)
  fit_raw <- tobs(~ occ_cov1 + occ_cov2, sim$data, family = occu(),
                  detection = ~ det_cov1, y = sim$y, method = "laplace")

  fr <- tobs_format(y = sim$y, occ.covs = sim$data)
  set.seed(1)
  fit_frame <- tobs(~ occ_cov1 + occ_cov2, fr, family = occu(),
                    detection = ~ det_cov1, method = "laplace")

  expect_equal(coef(fit_frame), coef(fit_raw))
})

test_that("control$verbose reports totals and the fit stores dims", {
  skip_on_cran(); skip_if_fast()
  sim <- simulate_occu(N = 40, J = 3, seed = 13)

  expect_message(
    fit <- tobs(~ occ_cov1, sim$data, family = occu(), detection = ~ 1,
                y = sim$y, method = "laplace", control = list(verbose = TRUE)),
    "fitting occu on 40 sites x 3 visits"
  )
  expect_identical(fit$dims$n_sites, 40L)
  expect_identical(fit$dims$max_visits, 3L)
})
