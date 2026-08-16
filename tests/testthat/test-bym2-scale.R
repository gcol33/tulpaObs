# The BYM2 Riebler scale factor (gcol33/tulpaObs#228, #232). One
# implementation, reached through a dense door, a Q door, a CSR door, and the
# fitter-argument resolver.
#
# #228 unified three implementations of one constant. #232 established that the
# constant they had unified on was the WRONG one -- the geometric mean of the
# eigenvalues of Q, where Riebler et al. (2016) is the geometric mean of the
# MARGINAL VARIANCES diag(Q^+), the quantity `INLA::inla.scale.model()` applies
# and the quantity `tulpa::compute_bym2_scale()` already reported. These assert
# the settled definition and both live spellings of the loading.

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

# diag(Q^+) the long way round, from an explicit pseudo-inverse.
icar_marginal_var <- function(adj) {
  Q <- diag(rowSums(adj)) - adj
  e <- eigen(Q, symmetric = TRUE)
  nz <- abs(e$values) > 1e-10
  V <- e$vectors[, nz, drop = FALSE]
  diag(V %*% (t(V) / e$values[nz]))
}

test_that("the scale factor is the geometric mean of the MARGINAL VARIANCES", {
  for (dim in list(c(4, 4), c(5, 5), c(10, 10), c(4, 25))) {
    adj <- lattice_adj(dim[1], dim[2])
    expect_equal(.bym2_scale(adj), exp(mean(log(icar_marginal_var(adj)))),
                 tolerance = 1e-10)
  }
})

test_that("it is NOT the geometric mean of the eigenvalues", {
  # The quantity this used to return. Kept as an explicit negative so the two
  # cannot be confused again -- and asserted on how they move with graph size,
  # not on their distance at one graph. On a lattice the two happen to pass
  # within 0.05 of each other around 10x10 and diverge either side of it, so a
  # fixed absolute gap is not the property that separates them.
  eig_mean <- function(adj) {
    ev <- eigen(diag(rowSums(adj)) - adj, symmetric = TRUE,
                only.values = TRUE)$values
    exp(mean(log(ev[ev > 1e-10])))
  }
  dims <- list(c(5, 5), c(10, 10), c(20, 20))
  adjs <- lapply(dims, function(d) lattice_adj(d[1], d[2]))
  s  <- vapply(adjs, .bym2_scale, numeric(1))
  em <- vapply(adjs, eig_mean, numeric(1))

  # Nowhere near equal, on any of the three, at any of the rescalings that
  # could plausibly have been meant.
  for (f in list(identity, function(x) 1 / x, sqrt, function(x) 1 / sqrt(x))) {
    expect_false(isTRUE(all.equal(s, f(em), tolerance = 1e-3)))
  }
  # And no fixed transform can relate them at all: the marginal-variance
  # constant RISES with graph size while every rescaling of the eigenvalue
  # constant falls, so the ratio is not even approximately constant.
  expect_true(all(diff(s) > 0))
  expect_true(all(diff(em) > 0))
  ratio <- s / (1 / sqrt(em))
  expect_gt(max(ratio) / min(ratio), 1.4)
})

test_that("it matches tulpa's implementation through the engine spelling", {
  # tulpa's compute_bym2_scale() returns the ENGINE loading 1 / sqrt(s) -- its
  # engine multiplies the structured block by it, where tulpaObs's own kernels
  # divide by sqrt(s). Two spellings of one loading, so the conversion has to
  # land exactly on the upstream value. This is the external cross-check on the
  # constant (gcol33/tulpaObs#232). Internal upstream, hence the triple colon.
  ref <- getFromNamespace("compute_bym2_scale", "tulpa")
  for (dim in list(c(4, 4), c(5, 5), c(10, 10), c(3, 12))) {
    adj <- lattice_adj(dim[1], dim[2])
    expect_equal(.bym2_engine_scale(.bym2_scale(adj)), ref(adj),
                 tolerance = 1e-10)
  }
})

test_that("the engine spelling is the reciprocal square root", {
  expect_equal(.bym2_engine_scale(0.25), 2)
  expect_equal(.bym2_engine_scale(4), 0.5)
  # A field loaded either way carries the same variance: sigma^2 rho / s.
  s <- .bym2_scale(lattice_adj(5, 5))
  expect_equal((1 / sqrt(s))^2, 1 / s)
})

test_that("the dense, Q and CSR doors return the same value to the bit", {
  for (dim in list(c(4, 4), c(6, 6), c(10, 10))) {
    adj <- lattice_adj(dim[1], dim[2])
    csr <- adjacency_to_csr(adj)
    expect_identical(.bym2_scale_csr(csr$row_ptr, csr$col_idx, nrow(adj)),
                     .bym2_scale(adj))
    expect_identical(.bym2_scale_from_Q(diag(rowSums(adj)) - adj),
                     .bym2_scale(adj))
  }
})

test_that("a disconnected graph drops one null eigenvalue per component", {
  b <- lattice_adj(5, 5)
  n <- nrow(b)
  adj <- matrix(0, 2 * n, 2 * n)
  adj[seq_len(n), seq_len(n)] <- b
  adj[n + seq_len(n), n + seq_len(n)] <- b
  # Q^+ of a block-diagonal Q is block-diagonal in the same blocks, so every
  # node keeps its single-component marginal variance and the geometric mean
  # over both components is the single-component value.
  expect_equal(.bym2_scale(adj), .bym2_scale(b), tolerance = 1e-8)
})

test_that("an eigenvalue negative past the tolerance is an error, not a NaN", {
  # A precision that is genuinely indefinite used to come back NaN at two of
  # the three call sites and a silently-dropped term at the third, and the NaN
  # then multiplied into the BYM2 mixing weight.
  adj <- matrix(0, 4, 4)
  adj[1, 2] <- adj[2, 1] <- -3      # spectrum of D - adj is (2, 0, 0, -6)
  adj[3, 4] <- adj[4, 3] <- 1
  expect_error(.bym2_scale(adj), "negative eigenvalue")
})

test_that("a graph with no edges is an error, not a NaN", {
  expect_error(.bym2_scale(matrix(0, 5, 5)), "no non-zero eigenvalues")
})

test_that("a roundoff-scale negative eigenvalue is dropped, not rejected", {
  # The filter is a numerical-zero test, so a null eigenvalue returned with a
  # negative sign has to be dropped rather than reach the error above. Built on
  # Q directly: perturbing an adjacency to produce one also disconnects the
  # node it touches, which is the separate error below.
  Q <- diag(c(1, 1, 0, 0)); Q[1, 2] <- Q[2, 1] <- -1
  Q[3, 4] <- Q[4, 3] <- -1e-13; Q[3, 3] <- Q[4, 4] <- 1e-13
  expect_error(.bym2_scale_from_Q(Q), "zero marginal variance")
  # The 2-node connected block alone: eigenvalues (2, 0), diag(Q^+) = 1/4 each.
  expect_equal(.bym2_scale_from_Q(matrix(c(1, -1, -1, 1), 2)), 0.25,
               tolerance = 1e-12)
})

test_that("an isolated node is an error, not a floored variance", {
  # It carries no ICAR prior at all, so its marginal variance is zero and the
  # geometric mean over the graph is zero. Flooring it (what the community
  # copy of this function did) turns that into a huge finite loading instead.
  adj <- lattice_adj(3, 3)
  adj <- rbind(cbind(adj, 0), 0)     # one extra node, no neighbours
  expect_error(.bym2_scale(adj), "zero marginal variance")
  expect_error(.bym2_scale(adj), "isolated")
})

test_that(".bym2_resolve_scale computes when absent and validates when given", {
  adj <- lattice_adj(5, 5)
  csr <- adjacency_to_csr(adj)
  expect_identical(
    .bym2_resolve_scale(NULL, csr$row_ptr, csr$col_idx, nrow(adj)),
    .bym2_scale(adj))
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

test_that("the term carries s and the engine boundary converts it", {
  # `bym2()` puts the Riebler constant on the term, because tulpaObs's own
  # kernels read it directly. Only the sites that hand it to the tulpa engine
  # convert (gcol33/tulpaObs#232).
  adj <- lattice_adj(4, 4)
  tm  <- .tobs_term_bym2(adj)
  expect_equal(tm$scale_factor, .bym2_scale(adj))
  expect_equal(.bym2_engine_scale(tm$scale_factor),
               getFromNamespace("compute_bym2_scale", "tulpa")(adj),
               tolerance = 1e-10)
  # A user-supplied value is taken as s, unconverted.
  expect_identical(.tobs_term_bym2(adj, scale_factor = 0.7)$scale_factor, 0.7)
})
