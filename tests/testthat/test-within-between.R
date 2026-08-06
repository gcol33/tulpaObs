# Tests for within_between(), the Mundlak within/between decomposition helper.

test_that("within_between() arithmetic is exact and group-respecting", {
  set.seed(1)
  d <- data.frame(
    plot = rep(letters[1:4], each = 5),
    year = rep(2000:2004, times = 4) + rep(c(0L, 5L, 10L, 15L), each = 5),
    z    = stats::rnorm(20)
  )

  out <- within_between(d, group = "plot", vars = c("year", "z"))

  expect_true(all(c("year_btw", "year_wtn", "z_btw", "z_wtn") %in% names(out)))

  # Decomposition is exact.
  expect_equal(out$year_btw + out$year_wtn, out$year)
  expect_equal(out$z_btw    + out$z_wtn,    out$z)

  # Per-group between values are constant within each group.
  for (g in unique(out$plot)) {
    idx <- which(out$plot == g)
    expect_equal(length(unique(out$year_btw[idx])), 1L)
    expect_equal(length(unique(out$z_btw[idx])),    1L)
    # And equal the per-group mean.
    expect_equal(out$year_btw[idx[1]], mean(d$year[idx]))
    expect_equal(out$z_btw[idx[1]],    mean(d$z[idx]))
  }

  # Within-group deviations sum to zero per group.
  agg <- stats::aggregate(out$year_wtn, by = list(out$plot), FUN = sum)
  expect_true(all(abs(agg$x) < 1e-10))
})

test_that("within_between() recovers separate within and between slopes", {
  set.seed(2)
  n_plot  <- 30L
  n_year  <- 8L
  beta_btw <- 0.05   # cross-plot heterogeneity in baseline year
  beta_wtn <- 0.20   # within-plot temporal trend
  plot_id <- rep(seq_len(n_plot), each = n_year)
  # Each plot has a different mean year (cross-plot heterogeneity).
  plot_year_offset <- rep(stats::runif(n_plot, 0, 30), each = n_year)
  year_within <- rep(seq_len(n_year) - mean(seq_len(n_year)), times = n_plot)
  year <- 2000 + plot_year_offset + year_within

  # True data-generating model uses separate within and between coefficients.
  y <- beta_btw * (plot_year_offset + 2000) + beta_wtn * year_within +
       stats::rnorm(n_plot * n_year, sd = 0.05)

  d <- data.frame(plot = factor(plot_id), year = year, y = y)
  out <- within_between(d, group = "plot", vars = "year")

  fit <- stats::lm(y ~ year_btw + year_wtn, data = out)
  cf  <- stats::coef(fit)

  expect_lt(abs(cf["year_btw"] - beta_btw), 0.02)
  expect_lt(abs(cf["year_wtn"] - beta_wtn), 0.02)
})

test_that("within_between() validates inputs and refuses to overwrite", {
  d <- data.frame(g = c(1, 1, 2, 2), x = 1:4)

  expect_error(within_between(list(), "g", "x"),       "data frame")
  expect_error(within_between(d, c(),  "x"),           "character vector")
  expect_error(within_between(d, "g",  c()),           "character vector")
  expect_error(within_between(d, "missing", "x"),      "Group column")
  expect_error(within_between(d, "g",  "missing"),     "Variable column")

  d2 <- d; d2$x <- as.character(d2$x)
  expect_error(within_between(d2, "g", "x"),           "must be numeric")

  expect_error(within_between(d, "g", "x",
                              suffix = c("_a", "_a")), "distinct")

  d3 <- within_between(d, "g", "x")
  expect_error(within_between(d3, "g", "x"),           "already exist")
})

test_that("within_between() handles multi-column groups via interaction", {
  d <- data.frame(
    site   = rep(c("A", "B"), each = 4),
    season = rep(c("spring", "fall"), times = 4),
    y      = c(1, 2, 3, 4, 10, 20, 30, 40)
  )
  out <- within_between(d, group = c("site", "season"), vars = "y")

  # site=A, season=spring -> rows 1,3 -> mean 2
  # site=A, season=fall   -> rows 2,4 -> mean 3
  # site=B, season=spring -> rows 5,7 -> mean 20
  # site=B, season=fall   -> rows 6,8 -> mean 30
  expect_equal(out$y_btw, c(2, 3, 2, 3, 20, 30, 20, 30))
  expect_equal(out$y_btw + out$y_wtn, out$y)
})
