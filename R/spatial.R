# =============================================================================
# spatial.R — Internal precompute helpers for spatial formula terms
#
# Spatial structure is specified inside a tobs() formula via the term
# constructors in formula_terms.R (`icar()`, `bym2()`, `gp()`,
# `multiscale_gp()`, `spde()`). Those constructors call the helpers below to
# turn an adjacency matrix into CSR arrays, compute the BYM2 scale factor, and
# build the NNGP neighbor structure. The helpers are not user-facing.
# =============================================================================

# Compressed-sparse-row representation of a 0/1 adjacency matrix (0-based
# column indices, as the C++ ICAR/BYM2 likelihood expects).
adjacency_to_csr <- function(adj) {
  n <- nrow(adj)
  row_ptr <- integer(n + 1)
  col_idx <- integer(0)
  n_neighbors <- integer(n)

  for (i in seq_len(n)) {
    neighbors <- which(adj[i, ] > 0)
    n_neighbors[i] <- length(neighbors)
    col_idx <- c(col_idx, neighbors - 1L)  # 0-based
    row_ptr[i + 1] <- row_ptr[i] + n_neighbors[i]
  }

  list(
    row_ptr = as.integer(row_ptr),
    col_idx = as.integer(col_idx),
    n_neighbors = as.integer(n_neighbors)
  )
}

# Dense 0/1 adjacency from the CSR arrays, for the callers that carry only the
# compressed form. Inverse of adjacency_to_csr() for a binary graph.
csr_to_adjacency <- function(csr, n) {
  adj <- matrix(0, n, n)
  for (i in seq_len(n)) {
    a <- csr$row_ptr[i] + 1L
    b <- csr$row_ptr[i + 1L]
    if (b >= a) adj[i, csr$col_idx[a:b] + 1L] <- 1
  }
  adj
}

# log|Q(rho)| = log|D - rho W| per rho grid point (tau- and z-independent),
# precomputed from the dense graph via a Cholesky. Enters the outer-grid
# marginal as a weight on each hyperparameter node, so every family fitting a
# proper CAR field reads it from here.
.tobs_car_logdet_Q <- function(graph, rho_grid) {
  D <- diag(rowSums(graph))
  vapply(rho_grid, function(rho) {
    Q <- D - rho * graph
    ch <- tryCatch(chol(Q), error = function(e) NULL)
    if (is.null(ch)) return(-Inf)            # Q(rho) not PD -> skip grid point
    2 * sum(log(diag(ch)))
  }, numeric(1))
}

# BYM2 scaling factor s: the geometric mean of the MARGINAL VARIANCES of the
# intrinsic ICAR field, i.e. of diag(Q^+) for the unit-precision ICAR precision
# Q = D - W (Riebler et al. 2016, section 3.2). Dividing the structured
# component by sqrt(s) puts it on unit geometric-mean marginal variance, which
# is what makes the mixing weight rho mean "the structured share of the
# variance" independently of the graph, and what makes it comparable across
# implementations.
#
# This is the constant `INLA::inla.scale.model()` applies: it returns s * Q, the
# precision whose generalised inverse has geometric-mean diagonal 1. Verified
# against it directly on square lattices --
#
#   lattice     inla.scale.model     .bym2_scale()
#   5x5              0.516386             0.516386
#   10x10            0.644879             0.644879
#   20x20            0.765026             0.765027
#
# -- and NOT the geometric mean of the eigenvalues of Q, which this function
# used to return (2.646529 / 2.831882 / 2.984944 on the same three graphs) and
# which no reciprocal or square root of maps onto the reference.
#
# The engine and the R paths spell the loading differently; `.bym2_engine_scale()`
# below is the boundary between the two conventions.
#
# Q is positive semi-definite with one null direction per connected component,
# so the eigenvalue filter is a NUMERICAL-ZERO test, not a sign test: the null
# eigenvalues come back at roundoff scale with either sign, and everything
# surviving the test carries a real marginal variance. An eigenvalue that
# survives it and is still negative says the matrix handed in is not an ICAR
# precision, and it errors here rather than reaching log() and sending NaN on
# into the BYM2 mixing weight. 1e-10 is the historical floor; the
# n * eps * max|lambda| term is the scale the roundoff-zero eigenvalues actually
# sit at, and overtakes the floor only on graphs far larger than the floor was
# picked for.
.bym2_scale <- function(adj) {
  .bym2_scale_from_Q(diag(rowSums(adj)) - adj)
}

# The same constant from the ICAR precision directly, for the callers that
# already hold Q rather than the adjacency.
.bym2_scale_from_Q <- function(Q) {
  Q <- as.matrix(Q)
  n <- nrow(Q)
  e <- eigen(Q, symmetric = TRUE)
  tol <- max(1e-10, n * .Machine$double.eps * max(abs(e$values)))
  keep <- abs(e$values) > tol
  if (!any(keep)) {
    stop("ICAR precision has no non-zero eigenvalues (the graph carries no ",
         "edges); the BYM2 scale factor is undefined for it.", call. = FALSE)
  }
  if (any(e$values[keep] < 0)) {
    stop(sprintf(paste0(
      "ICAR precision has a negative eigenvalue (%.3e) past the numerical-zero ",
      "tolerance (%.3e); the BYM2 scale factor is undefined for it. The ",
      "adjacency matrix must be symmetric 0/1 with a zero diagonal."),
      min(e$values[keep]), tol), call. = FALSE)
  }
  # diag(Q^+), assembled from the surviving spectrum. A node that sits in the
  # null space contributes a zero marginal variance -- it carries no ICAR prior
  # at all -- and the geometric mean over the graph is then zero, so this is an
  # error rather than a floored value that would surface as a huge loading.
  V <- e$vectors[, keep, drop = FALSE]
  vdiag <- as.numeric((V^2) %*% (1 / e$values[keep]))
  if (any(vdiag <= tol)) {
    stop(sprintf(paste0(
      "%d graph node(s) have zero marginal variance under the ICAR prior (they ",
      "are isolated -- no neighbours); the BYM2 scale factor is undefined for ",
      "such a graph. Drop them or connect them."), sum(vdiag <= tol)),
      call. = FALSE)
  }
  exp(mean(log(vdiag)))
}

# Same constant for the fitters that carry only the compressed graph.
.bym2_scale_csr <- function(row_ptr, col_idx, n) {
  .bym2_scale(
    csr_to_adjacency(list(row_ptr = row_ptr, col_idx = col_idx), n))
}

# The BYM2 loading, in the two spellings that are live.
#
# The structured block is always built from the UNSCALED ICAR precision, so its
# realization v has covariance Q^+ and geometric-mean marginal variance s. To
# put the field on unit geometric-mean marginal variance the structured
# component has to be loaded at sigma * sqrt(rho / s).
#
# That is what tulpaObs's own R paths and its own C++ kernels write, so they
# take `s` and divide. The tulpa engine instead writes
# `sigma * sqrt(rho) * scale_factor` (src/latent_block.h states it as the block
# contract), so its `scale_factor` argument means 1 / sqrt(s). Both spellings
# describe the same loading; only the argument differs, and this is the one
# place the conversion happens, so a call site's convention is visible at the
# call rather than inferred from which kernel it reaches.
.bym2_engine_scale <- function(s) 1 / sqrt(s)

# Unit-scale realisations of the improper ICAR field on a graph: f ~ N(0, Q^+)
# on the sum-to-zero space, divided by sqrt(scale_q) so the geometric-mean
# marginal SD is 1 -- the Sorbye-Rue convention the fitters and every
# simulate_* entry share. Returns an [n_nodes x n_draw] matrix. `Q` and
# `scale_q` are passed in because every call site already holds them.
#
# The draw goes through the Cholesky of Q + 11'/N, NOT an eigenbasis of Q. Both
# give the same distribution, but a lattice Q has repeated eigenvalues in bulk
# (31 of 64 on an 8x8 grid, 39 of 81 on a 9x9), and inside a repeated block the
# eigenvector basis is fixed only up to an orthogonal rotation that no LAPACK
# build is obliged to resolve the same way. An eigen draw therefore makes the
# realisation, at a FIXED seed, a property of the linear-algebra library rather
# than of the seed: two exact eigendecompositions of the same Q give fields
# correlating at -0.01. That is what made three areal fixtures report different
# answers on the Linux runner than here, on nominally the same data (#279).
#
# chol() of a positive-definite matrix is unique, so this construction has no
# such freedom. Q + 11'/N is positive definite on a connected graph, and its
# inverse agrees with Q^+ on the sum-to-zero space that the centring projects
# onto. It is the construction the car_proper branch of
# simulate_ms_occu_cover_spatial() already used.
.tobs_draw_icar_unit <- function(Q, scale_q, n_draw = 1L) {
  N  <- nrow(Q)
  Rc <- chol(Q + matrix(1 / N, N, N))
  out <- vapply(seq_len(n_draw), function(i) {
    f <- backsolve(Rc, stats::rnorm(N))
    (f - mean(f)) / sqrt(scale_q)
  }, numeric(N))
  matrix(out, N, as.integer(n_draw))
}

# Resolve the scale factor at an areal-BYM2 fitter door: compute it from the
# graph when the caller left it out, validate whatever came in. Both count
# fitters (nmix, removal) go through here, so a missing argument means the same
# thing at each door. Returns `s` -- these doors reach tulpaObs's own kernels.
.bym2_resolve_scale <- function(scale_factor, row_ptr, col_idx, n) {
  if (is.null(scale_factor)) {
    scale_factor <- .bym2_scale_csr(row_ptr, col_idx, n)
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be a positive scalar.", call. = FALSE)
  }
  as.numeric(scale_factor)
}

# Flatten the [N, k, k] neighbour-pair distance array for the NNGP kernels.
#
# `compute_nngp_neighbors()` returns an R array indexed [i, j1, j2]; the hmc_gp
# kernels read it row-major, at i * k * k + j1 * k + j2 (tulpa
# src/hmc_gp_autodiff.h, and the convention is stated at
# inst/include/tulpa/sampler_model_data.h). `as.vector()` is column-major, so it
# hands every off-diagonal entry of the neighbour covariance a distance from an
# unrelated cell. Both orderings have N * k * k elements, so the kernel's bounds
# guard passes and the corruption is silent: the neighbour covariance goes
# near-singular, its conditional variance pins at the kernel's floor, and the
# field is left effectively unconstrained.
#
# `populate_helpers.h` rebases nn_order / nn_idx for the engine, but passes this
# array straight through, so the permutation has to happen here.
.tobs_nngp_pair_dist <- function(nn_neighbor_dist) {
  as.numeric(aperm(nn_neighbor_dist, c(3L, 2L, 1L)))
}

# Nearest-neighbor structure for the NNGP approximation: for each location
# (in a coordinate ordering) its `k` nearest predecessors, their distances,
# and the pairwise distances among those neighbors.
#
# Every step reads whatever coordinate dimension it is handed. The neighbour
# SELECTION here and the neighbour COVARIANCE the kernels build from it have to
# be computed under the same metric, and the kernels read every column, so a
# selection pinned to two columns would order by a projection of the domain the
# covariance does not use.
compute_nngp_neighbors <- function(coords, k) {
  N <- nrow(coords)
  d <- ncol(coords)

  # Lexicographic ordering over every coordinate column: a valid NNGP ordering
  # that conditions better than the raw input order.
  order_idx <- do.call(order, lapply(seq_len(d), function(j) coords[, j]))
  coords_ordered <- coords[order_idx, , drop = FALSE]

  nn_idx <- matrix(0L, nrow = N, ncol = k)
  nn_dist <- matrix(Inf, nrow = N, ncol = k)
  nn_neighbor_dist <- array(0, dim = c(N, k, k))

  # `seq_len(N)[-1]` is empty at N <= 1, where `2:N` counts DOWN to c(2, 1) and
  # indexes a row that does not exist.
  for (i in seq_len(N)[-1]) {
    n_candidates <- min(i - 1, k)
    if (n_candidates > 0) {
      # Euclidean distance from every earlier location to this one. The point is
      # laid out column-major by `each = `, the layout `prev` already has, so
      # this is one vectorised subtraction rather than a sweep.
      prev <- coords_ordered[1:(i-1), , drop = FALSE]
      dists <- sqrt(rowSums(
        (prev - rep(as.numeric(coords_ordered[i, ]), each = nrow(prev)))^2
      ))

      if (length(dists) <= k) {
        nn_order <- order(dists)
        nn_idx[i, seq_len(length(dists))] <- nn_order
        nn_dist[i, seq_len(length(dists))] <- dists[nn_order]
      } else {
        nn_order <- order(dists)[1:k]
        nn_idx[i, ] <- nn_order
        nn_dist[i, ] <- dists[nn_order]
      }

      n_neighbors <- sum(nn_idx[i, ] > 0)
      if (n_neighbors > 1) {
        neighbor_indices <- nn_idx[i, 1:n_neighbors]
        neighbor_coords <- coords_ordered[neighbor_indices, , drop = FALSE]
        # The whole pairwise block in one assignment, over every coordinate
        # column. The diagonal is exactly 0 by construction.
        nn_neighbor_dist[i, seq_len(n_neighbors), seq_len(n_neighbors)] <-
          as.matrix(stats::dist(neighbor_coords))
      }
    }
  }

  list(
    nn_idx = nn_idx,
    nn_dist = nn_dist,
    nn_neighbor_dist = nn_neighbor_dist,
    nn_order = order_idx,
    nn_order_inv = order(order_idx),
    k = k
  )
}
