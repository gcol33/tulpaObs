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

# BYM2 scaling factor: geometric mean of the non-zero generalized eigenvalues
# of the ICAR precision, so the mixing parameter has a graph-independent
# interpretation (Riebler et al. 2016).
compute_bym2_scale <- function(adj) {
  n <- nrow(adj)
  D <- diag(rowSums(adj))
  Q <- D - adj
  evals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  evals <- evals[evals > 1e-10]
  exp(mean(log(evals)))
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
compute_nngp_neighbors <- function(coords, k) {
  N <- nrow(coords)
  order_idx <- order(coords[, 1], coords[, 2])
  coords_ordered <- coords[order_idx, , drop = FALSE]

  nn_idx <- matrix(0L, nrow = N, ncol = k)
  nn_dist <- matrix(Inf, nrow = N, ncol = k)
  nn_neighbor_dist <- array(0, dim = c(N, k, k))

  for (i in 2:N) {
    n_candidates <- min(i - 1, k)
    if (n_candidates > 0) {
      dists <- sqrt(
        (coords_ordered[1:(i-1), 1] - coords_ordered[i, 1])^2 +
        (coords_ordered[1:(i-1), 2] - coords_ordered[i, 2])^2
      )

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
        for (j1 in 1:n_neighbors) {
          for (j2 in 1:n_neighbors) {
            if (j1 != j2) {
              nn_neighbor_dist[i, j1, j2] <- sqrt(
                (neighbor_coords[j1, 1] - neighbor_coords[j2, 1])^2 +
                (neighbor_coords[j1, 2] - neighbor_coords[j2, 2])^2
              )
            }
          }
        }
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
