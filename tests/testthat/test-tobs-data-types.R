# =============================================================================
# test-tobs-data-types.R - tobs_data() typed response reshape.
#
# Guards the backbone fix: the response matrix type and validation follow the
# `type` argument. Regression for the bug where a continuous cover response was
# coerced with as.integer(), silently truncating every value < 1 to zero.
# =============================================================================

make_long <- function(yvals) {
  # 3 sites x 2 visits, first-appearance site order deliberately non-sorted.
  data.frame(
    site  = c(3, 3, 1, 1, 2, 2),
    visit = c(1, 2, 1, 2, 1, 2),
    y     = yvals,
    xcov  = seq_along(yvals) + 0.5
  )
}

test_that("type = 'cover' preserves continuous values (no integer truncation)", {
  cov <- c(0.02, 0.63, 0.001, 0.5, 0.0, 0.31)
  od <- tobs_data(make_long(cov), y = "y", site = "site", visit = "visit",
                  type = "cover")
  expect_true(is.double(od$y))
  # row 1 = site 3 (first appearance), values land by site/visit match
  expect_equal(unname(od$y[1, 1]), 0.02)
  expect_equal(unname(od$y[1, 2]), 0.63)
  expect_equal(sort(as.vector(od$y)), sort(cov))
  expect_false(any(as.vector(od$y) == 0 & sort(cov)[1] > 0))  # nothing truncated
})

test_that("type = 'occurrence' builds an integer 0/1 matrix", {
  od <- tobs_data(make_long(c(1, 0, 1, 1, 0, 0)), y = "y", site = "site",
                  visit = "visit", type = "occurrence")
  expect_true(is.integer(od$y))
  expect_setequal(stats::na.omit(as.vector(od$y)), c(0L, 1L))
})

test_that("type = 'abundance' keeps integer counts > 1", {
  od <- tobs_data(make_long(c(0, 3, 7, 1, 12, 0)), y = "y", site = "site",
                  visit = "visit", type = "abundance")
  expect_true(is.integer(od$y))
  expect_true(max(od$y, na.rm = TRUE) == 12L)
})

test_that("default type is occurrence (back-compatible)", {
  od_default <- tobs_data(make_long(c(1, 0, 1, 1, 0, 0)), y = "y",
                          site = "site", visit = "visit")
  od_occ <- tobs_data(make_long(c(1, 0, 1, 1, 0, 0)), y = "y", site = "site",
                      visit = "visit", type = "occurrence")
  expect_identical(od_default$y, od_occ$y)
})

test_that("type validation rejects out-of-domain responses", {
  expect_error(tobs_data(make_long(c(0, 2, 1, 0, 1, 0)), y = "y", site = "site",
                         visit = "visit", type = "occurrence"), "0/1")
  expect_error(tobs_data(make_long(c(0.1, 1.4, 0.2, 0, 0.3, 0)), y = "y",
                         site = "site", visit = "visit", type = "cover"),
               "\\[0, 1\\]")
  expect_error(tobs_data(make_long(c(0, -1, 2, 0, 1, 0)), y = "y", site = "site",
                         visit = "visit", type = "abundance"), "non-negative")
})

test_that("occurrence and cover reshapes align position-for-position", {
  # The occu_cover use case: detection (0/1) and cover (continuous) built from
  # the same long frame must share a layout so y_pos > 0 wherever detection==1.
  cov <- c(0.02, 0.0, 0.63, 0.0, 0.0, 0.31)
  occ <- as.integer(cov > 0)
  df  <- make_long(cov); df$occ <- occ
  od     <- tobs_data(df, y = "occ", site = "site", visit = "visit",
                      type = "occurrence")
  od_cov <- tobs_data(df, y = "y",   site = "site", visit = "visit",
                      type = "cover")
  detected <- od$y == 1 & !is.na(od$y)
  expect_true(all(od_cov$y[detected] > 0))
})
