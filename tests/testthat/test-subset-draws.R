test_that(".tobs_subset_draws spaces rows over every stacked chain (#419)", {
  # Two chains of 1000 kept draws, pooled by rbind as the samplers do.
  chain <- rep(1:2, each = 1000L)
  draws <- cbind(chain = chain, b = seq_along(chain))
  sub <- tulpaObs:::.tobs_subset_draws(draws, 1000L)
  expect_lte(nrow(sub), 1000L)
  expect_gt(nrow(sub), 990L)
  expect_equal(as.numeric(table(sub[, "chain"])) / nrow(sub), c(0.5, 0.5),
               tolerance = 0.01)
  expect_equal(sub[[1, "b"]], 1L)
  expect_equal(sub[[nrow(sub), "b"]], 2000L)
  expect_identical(tulpaObs:::.tobs_subset_draws(draws, NULL), draws)
  expect_identical(tulpaObs:::.tobs_subset_draws(draws, 5000L), draws)
  one <- tulpaObs:::.tobs_subset_draws(draws[, 1, drop = FALSE], 10L)
  expect_identical(dim(one), c(10L, 1L))
})
