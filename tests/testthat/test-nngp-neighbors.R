# The NNGP neighbour graph is built by `compute_nngp_neighbors()` and reached
# from every continuous-field term through `.tobs_build_nngp_block()`
# (`R/formula_terms.R`). The term doors admit exactly two coordinate columns and
# clamp the neighbour count to `n - 1`, so what those doors happen to reject is
# not evidence about the function. The tests here call it directly: its
# behaviour at N == 1, and the coordinate dimension it reads, are its own
# contract.
#
# The neighbour SELECTION and the neighbour COVARIANCE the compiled kernels
# build from it have to be computed under the same metric, and the kernels read
# every coordinate column.

test_that("compute_nngp_neighbors() builds the empty graph at N == 1", {
  one <- matrix(c(0.3, 0.7), nrow = 1L)
  g <- compute_nngp_neighbors(one, 4L)

  expect_identical(dim(g$nn_idx), c(1L, 4L))
  expect_true(all(g$nn_idx == 0L))
  expect_true(all(is.infinite(g$nn_dist)))
  expect_true(all(g$nn_neighbor_dist == 0))
  expect_identical(g$nn_order, 1L)
  expect_identical(g$nn_order_inv, 1L)
  expect_identical(g$k, 4L)

  # The term door clamps the neighbour count to what n points allow, so it
  # reaches the same function at k = 0. Both arities give a well-formed graph.
  b <- .tobs_build_nngp_block(one, 4L)
  expect_identical(b$nn, 0L)
  expect_length(b$nn_idx, 0L)
  expect_length(b$nn_dist, 0L)
  expect_length(b$nn_neighbor_dist, 0L)
  expect_identical(b$nn_order, 1L)
})


test_that("the block a term packs reads every coordinate column", {
  set.seed(3L)
  n <- 40L
  co <- cbind(stats::runif(n), stats::runif(n))
  b2 <- .tobs_build_nngp_block(co, 5L)

  # A constant third coordinate is the same domain, so the packed graph is
  # unchanged.
  bc <- .tobs_build_nngp_block(cbind(co, 7), 5L)
  expect_identical(bc$nn_idx, b2$nn_idx)
  expect_equal(bc$nn_dist, b2$nn_dist, tolerance = 0)
  expect_equal(bc$nn_neighbor_dist, b2$nn_neighbor_dist, tolerance = 0)

  # A third coordinate that varies moves the graph, so the block is built over
  # every column rather than the first two.
  b3 <- .tobs_build_nngp_block(cbind(co, stats::runif(n)), 5L)
  expect_false(identical(b3$nn_idx, b2$nn_idx))

  # One coordinate column builds a graph at all.
  b1 <- .tobs_build_nngp_block(matrix(sort(stats::runif(n)), ncol = 1L), 5L)
  expect_length(b1$nn_idx, n * 5L)
  expect_true(all(is.finite(b1$nn_neighbor_dist)))
})


test_that("the two-column graph is the lexicographic Euclidean one, exactly", {
  # The packing tests read neighbour ORDER, so pin it against the explicit
  # two-column formula. Exact, not toleranced.
  set.seed(11L)
  n <- 30L
  k <- 5L
  co <- cbind(stats::runif(n), stats::runif(n))
  g <- compute_nngp_neighbors(co, k)

  ord <- order(co[, 1], co[, 2])
  expect_identical(g$nn_order, ord)
  expect_identical(g$nn_order_inv, order(ord))

  cs <- co[ord, , drop = FALSE]
  for (i in c(2L, 6L, 20L, 30L)) {
    m <- min(i - 1L, k)
    d <- sqrt((cs[seq_len(i - 1L), 1] - cs[i, 1])^2 +
              (cs[seq_len(i - 1L), 2] - cs[i, 2])^2)
    o <- order(d)[seq_len(m)]
    expect_identical(g$nn_idx[i, seq_len(m)], as.integer(o))
    expect_equal(unname(g$nn_dist[i, seq_len(m)]), d[o], tolerance = 0)

    if (m > 1L) {
      nb <- cs[o, , drop = FALSE]
      pair <- outer(seq_len(m), seq_len(m), function(a, b) {
        sqrt((nb[a, 1] - nb[b, 1])^2 + (nb[a, 2] - nb[b, 2])^2)
      })
      expect_equal(unname(g$nn_neighbor_dist[i, seq_len(m), seq_len(m)]),
                   pair, tolerance = 0)
    }
  }

  # Rows beyond the neighbour count stay at the initialised sentinels.
  expect_true(all(g$nn_idx[1L, ] == 0L))
  expect_true(all(is.infinite(g$nn_dist[1L, ])))
})
