# helper mirroring the package internal so tests do not depend on ::: for it
.scan_is_true_test <- function(x) !is.na(x) & x

test_that("changepoint segmentation finds a mean shift", {
  rate <- c(rep(0.1, 6), rep(0.8, 6))
  seg <- tulpaObs:::.scan_changepoints(rate, min_seg = 1L, penalty = 0.05)
  expect_length(seg, length(rate))
  expect_equal(length(unique(seg)), 2L)
  # boundary between position 6 and 7
  expect_true(all(seg[1:6] == 1L))
  expect_true(all(seg[7:12] == 2L))
})

test_that("changepoint returns a single segment for a flat series", {
  rate <- rep(0.3, 8) + c(0, 0.01, -0.01, 0, 0.01, -0.01, 0, 0.01)
  seg <- tulpaObs:::.scan_changepoints(
    rate, min_seg = 1L,
    penalty = tulpaObs:::.scan_cp_penalty(rate, 1))
  expect_equal(length(unique(seg)), 1L)
})

test_that("constant-model negloglik is finite and minimised near the truth", {
  # one unit type: K = 3 occasions, simulate detected/undetected under psi,p
  set.seed(1)
  psi <- 0.5; p <- 0.4; K <- 3L; n <- 400L
  z <- rbinom(n, 1, psi)
  d <- vapply(seq_len(n), function(i) sum(rbinom(K, 1, z[i] * p)), integer(1))
  Kv <- rep(K, n)
  f0 <- tulpaObs:::.scan_occu_negll(c(qlogis(psi), qlogis(p)), Kv, d)
  expect_true(is.finite(f0))
  fit <- tulpaObs:::.scan_fit_const(Kv, d, max_condition = 1e6,
                                    optim_control = list(maxit = 200))
  expect_true(fit$identifiable)
  expect_lt(abs(fit$psi_hat - psi), 0.12)
  expect_lt(abs(fit$p_hat - p), 0.12)
})

# ---------------------------------------------------------------------------
# End-to-end: single-visit plot data, closure holds, pooling years rescues
# identifiability.
# ---------------------------------------------------------------------------

make_plot_data <- function(seed = 7) {
  set.seed(seed)
  grid <- expand.grid(x = 1:14, y = 1:14)
  n_plot <- nrow(grid)
  z <- rbinom(n_plot, 1, 0.5)          # occupancy fixed across years (closed)
  years <- 1:9
  p <- 0.4
  rows <- do.call(rbind, lapply(years, function(yr) {
    data.frame(plot = seq_len(n_plot), x = grid$x, y = grid$y, year = yr,
               y_det = rbinom(n_plot, 1, z * p))
  }))
  rows
}

test_that("scan flags single-visit (K=1) as non-identifiable, pooled years as identifiable", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    plot = "plot", cell_sizes = 0.9, block_lengths = c(1L, 3L),
    score = "info")

  expect_s3_class(res, "tobs_aggregation_scan")
  expect_equal(nrow(res$candidates), 2L)

  row_l1 <- res$candidates[res$candidates$block == "1yr", ]
  row_l3 <- res$candidates[res$candidates$block == "3yr", ]
  expect_false(.scan_is_true_test(row_l1$identifiable))
  expect_true(.scan_is_true_test(row_l3$identifiable))

  # recommended = the identifiable, least-pooling candidate (3yr here)
  expect_false(is.null(res$recommended))
  expect_equal(res$recommended$block, "3yr")
  expect_true(is.finite(res$recommended$se_p))
})

test_that("auto changepoint segmentation runs and yields one block under closure", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    cell_sizes = 0.9, block_lengths = NULL, score = "count")
  expect_false(is.null(res$segmentation))
  expect_equal(nrow(res$segmentation), length(unique(d$year)))
  expect_equal(length(unique(res$segmentation$segment)), 1L)
})

test_that("count mode reports structural conditions without fitting", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    cell_sizes = 0.9, block_lengths = c(1L, 3L), score = "count")
  expect_true(all(is.na(res$candidates$se_p)))
  row_l1 <- res$candidates[res$candidates$block == "1yr", ]
  expect_false(.scan_is_true_test(row_l1$identifiable))   # K=1, no replication
})

test_that("auto cell-size ladder is proposed when cell_sizes is NULL", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    cell_sizes = NULL, block_lengths = 3L, score = "count")
  expect_gt(length(res$cell_sizes), 1L)
  expect_equal(nrow(res$candidates), length(res$cell_sizes))
})

test_that("print.tobs_aggregation_scan() reports the scan summary (#276)", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    plot = "plot", cell_sizes = 0.9, block_lengths = c(1L, 3L),
    score = "info")
  expect_output(print(res), "Occupancy aggregation scan")
})

test_that("plot.tobs_aggregation_scan() draws the identifiability heatmap (#276)", {
  d <- make_plot_data()
  res <- occu_aggregation_scan(
    d, response = "y_det", coords = c("x", "y"), year = "year",
    plot = "plot", cell_sizes = 0.9, block_lengths = c(1L, 3L),
    score = "info")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(res))
})

test_that("input validation rejects bad response and missing columns", {
  d <- make_plot_data()
  d$bad <- d$y_det + 1L
  expect_error(
    occu_aggregation_scan(d, "bad", c("x", "y"), "year", cell_sizes = 1),
    "0/1")
  expect_error(
    occu_aggregation_scan(d, "nope", c("x", "y"), "year", cell_sizes = 1),
    "not found")
})
