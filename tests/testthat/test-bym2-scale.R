# The BYM2 Riebler scale factor (gcol33/tulpaObs#228). One implementation,
# reached through a dense door, a CSR door, and the fitter-argument resolver.
# The three used to be three separate implementations with two different
# eigenvalue filters.

lattice_adj <- function(n_row, n_col) {
  n <- n_row * n_col
  adj <- matrix(0, n, n)
  for (s in seq_len(n)) {
    r <- ((s - 1L) %/% n_col) + 1L
    cc <- ((s - 1L) %% n_col) + 1L
    if (r > 1L)     adj[s, s - n_col] <- 1
    if (r < n_row)  adj[s, s + n_col] <- 1
    if (cc > 1L)    adj[s, s - 1L] <- 1
    if (cc < n_col) adj[s, s + 1L] <- 1
  }
  adj
}

test_that("scale factor is the geometric mean of the non-zero ICAR eigenvalues", {
  for (dim in list(c(4, 4), c(5, 5), c(10, 10), c(4, 25))) {
    adj <- lattice_adj(dim[1], dim[2])
    Q <- diag(rowSums(adj)) - adj
    ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
    oracle <- exp(mean(log(ev[ev > 1e-10])))
    expect_equal(compute_bym2_scale(adj), oracle, tolerance = 1e-12)
  }
})

test_that("the dense and CSR doors return the same value to the bit", {
  for (dim in list(c(4, 4), c(6, 6), c(10, 10))) {
    adj <- lattice_adj(dim[1], dim[2])
    csr <- adjacency_to_csr(adj)
    expect_identical(
      compute_bym2_scale_csr(csr$row_ptr, csr$col_idx, nrow(adj)),
      compute_bym2_scale(adj))
  }
})

test_that("a disconnected graph drops one null eigenvalue per component", {
  b <- lattice_adj(5, 5)
  n <- nrow(b)
  adj <- matrix(0, 2 * n, 2 * n)
  adj[seq_len(n), seq_len(n)] <- b
  adj[n + seq_len(n), n + seq_len(n)] <- b
  # Two identical components: each contributes the same non-zero spectrum, so
  # the geometric mean over both is the single-component value.
  expect_equal(compute_bym2_scale(adj), compute_bym2_scale(b), tolerance = 1e-10)
})

test_that("an eigenvalue negative past the tolerance is an error, not a NaN", {
  # A precision that is genuinely indefinite used to come back NaN at two of
  # the three call sites and a silently-dropped term at the third, and the NaN
  # then multiplied into the BYM2 mixing weight.
  adj <- matrix(0, 4, 4)
  adj[1, 2] <- adj[2, 1] <- -3      # spectrum of D - adj is (2, 0, 0, -6)
  adj[3, 4] <- adj[4, 3] <- 1
  expect_error(compute_bym2_scale(adj), "negative eigenvalue")
})

test_that("a graph with no edges is an error, not a NaN", {
  expect_error(compute_bym2_scale(matrix(0, 5, 5)), "no non-zero eigenvalues")
})

test_that("a roundoff-scale negative eigenvalue is dropped, not rejected", {
  # The filter is a numerical-zero test, so a null eigenvalue returned with a
  # negative sign has to be dropped rather than reach log() or the error above.
  adj <- matrix(0, 4, 4)
  adj[1, 2] <- adj[2, 1] <- -1e-12  # spectrum of D - adj is (2, 0, 0, -2e-12)
  adj[3, 4] <- adj[4, 3] <- 1
  val <- compute_bym2_scale(adj)
  expect_true(is.finite(val) && val > 0)
  expect_equal(val, 2, tolerance = 1e-9)
})

test_that(".bym2_resolve_scale computes when absent and validates when given", {
  adj <- lattice_adj(5, 5)
  csr <- adjacency_to_csr(adj)
  expect_identical(
    .bym2_resolve_scale(NULL, csr$row_ptr, csr$col_idx, nrow(adj)),
    compute_bym2_scale(adj))
  expect_identical(.bym2_resolve_scale(2.5, csr$row_ptr, csr$col_idx, nrow(adj)), 2.5)
  for (bad in list(0, -1, NA_real_, NaN, c(1, 2), "1")) {
    expect_error(.bym2_resolve_scale(bad, csr$row_ptr, csr$col_idx, nrow(adj)),
                 "positive scalar")
  }
})

test_that("both areal count fitters treat a missing scale_factor the same way", {
  # removal_laplace_bym2() used to default it to 1 and lean on the caller to
  # supply one, while nmix_laplace_bym2() computed it; a missing argument meant
  # different things at the two doors.
  for (f in list(nmix_laplace_bym2, removal_laplace_bym2)) {
    fm <- formals(f)
    expect_true("scale_factor" %in% names(fm))
    expect_null(fm$scale_factor)
  }
})
