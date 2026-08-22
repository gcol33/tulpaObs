# The NNGP neighbour-pair distance array crosses into the compiled kernels as a
# flat vector, and the two sides disagreed about its ordering.
# `compute_nngp_neighbors()` returns [i, j1, j2]; the hmc_gp kernels read it
# row-major at i * k * k + j1 * k + j2, while the term packed it with
# `as.vector()`, which is column-major. Both orderings have N * k * k elements,
# so the kernel's bounds guard passed and nothing surfaced: the neighbour
# covariance went near-singular and the GP prior stopped constraining the field
# at all.
#
# `nn_order` / `nn_idx` are NOT part of this: they stay 1-based in the term and
# `populate_helpers.h` rebases them for the engine (which is why the twin below
# is handed `ord - 1L`, matching the SVC twin test). Only the pair-distance
# array is passed straight through, so only it has to be permuted in R.
#
# This pins the convention against the engine itself rather than against a
# recovery run, so a revert fails here in seconds instead of surfacing as a
# quietly wrong field.

.gpp_term <- function(N = 50L, k = 8L, seed = 7L) {
  set.seed(seed)
  co <- cbind(stats::runif(N), stats::runif(N))
  list(co = co, tm = .tobs_term_gp(coords = co, nn = k,
                                   prior_range = c(0.1, 0.05)))
}

.gpp_dense_ref <- function(co, ord, w, sigma2, phi) {
  # The kernels are range-parameterised, exp(-d / phi).
  D  <- as.matrix(stats::dist(co[ord, , drop = FALSE]))
  Sg <- sigma2 * exp(-D / phi)
  ch <- chol(Sg + diag(1e-8, nrow(Sg)))
  -0.5 * sum(backsolve(ch, w[ord], transpose = TRUE)^2) -
    sum(log(diag(ch))) - 0.5 * nrow(Sg) * log(2 * pi)
}


test_that("gp() packs the neighbour-pair distances in the kernel's order", {
  skip_if_not(exists("cpp_test_gp_nngp_twins", envir = asNamespace("tulpa"),
                     inherits = FALSE),
              "installed tulpa has no GP NNGP twin probe")
  N <- 50L; k <- 8L
  g  <- .gpp_term(N, k)
  tm <- g$tm
  nn_idx  <- matrix(as.integer(tm$nn_idx),  nrow = N, byrow = TRUE)
  nn_dist <- matrix(as.numeric(tm$nn_dist), nrow = N, byrow = TRUE)
  ord     <- as.integer(tm$nn_order)

  twin <- get("cpp_test_gp_nngp_twins", envir = asNamespace("tulpa"))
  set.seed(3)
  w <- as.numeric(scale(stats::rnorm(N)))
  sigma2 <- 1.4; phi <- 0.35

  # `populate_helpers.h` rebases nn_order for the engine; the pair-distance
  # array is handed over exactly as the term stores it.
  got <- twin(w, sigma2, phi, g$co, nn_idx, nn_dist,
              tm$nn_neighbor_dist, ord - 1L, order(ord) - 1L, 0L)

  # Both twins agree, and the density is finite -- the kernel signals a rejected
  # input by returning -1e10, which is what the column-major flattening produced.
  expect_equal(unname(got[["dbl"]]), unname(got[["ad"]]), tolerance = 1e-8)
  expect_true(is.finite(got[["dbl"]]))
  expect_gt(got[["dbl"]], -1e9)

  # An NNGP with k = 8 neighbours is a close approximation of the dense GP it
  # approximates: assert on the relative gap so this measures the ordering, not
  # the approximation.
  ref <- .gpp_dense_ref(g$co, ord, w, sigma2, phi)
  expect_lt(abs(got[["dbl"]] - ref) / abs(ref), 0.05)

  # The direction matters, so pin it: the old column-major flattening is
  # rejected outright by the same kernel at the same inputs.
  bad <- twin(w, sigma2, phi, g$co, nn_idx, nn_dist,
              as.vector(aperm(array(tm$nn_neighbor_dist, c(k, k, N)),
                              c(3L, 2L, 1L))),
              ord - 1L, order(ord) - 1L, 0L)
  expect_lt(bad[["ad"]], -1e9)
})


test_that("multiscale_gp() packs both scales in the kernel's order", {
  set.seed(9)
  N <- 40L; k <- 6L
  co <- cbind(stats::runif(N), stats::runif(N))
  tm <- .tobs_term_multiscale_gp(coords = co, nn_local = k, nn_regional = k)
  ref <- .tobs_nngp_pair_dist(compute_nngp_neighbors(co, k)$nn_neighbor_dist)
  expect_equal(tm$nn_neighbor_dist_local, ref)
  expect_equal(tm$nn_neighbor_dist_regional, ref)
  expect_length(tm$nn_neighbor_dist_local, N * k * k)
})


test_that("the pair-distance flattening is row-major in [i, j1, j2]", {
  # Independent of any kernel: element (i, j1, j2) must land at
  # i * k * k + j1 * k + j2, the stride hmc_gp_autodiff.h reads with.
  a <- array(0, dim = c(3L, 2L, 2L))
  for (i in 1:3) for (j1 in 1:2) for (j2 in 1:2) a[i, j1, j2] <- i * 100 + j1 * 10 + j2
  flat <- .tobs_nngp_pair_dist(a)
  expect_length(flat, 12L)
  for (i in 1:3) for (j1 in 1:2) for (j2 in 1:2) {
    expect_identical(flat[(i - 1L) * 4L + (j1 - 1L) * 2L + j2], a[i, j1, j2])
  }
  # the symmetric array a[i, j1, j2] == a[i, j2, j1] would hide a j1/j2 swap,
  # so the fixture above is deliberately asymmetric in the last two indices
  expect_false(identical(a[1, 1, 2], a[1, 2, 1]))
})
