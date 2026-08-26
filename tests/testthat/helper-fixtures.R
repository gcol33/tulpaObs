# Shared adjacency-graph fixtures for spatial recovery/unit tests.
#
# tulpaObs#276: a chain (path-graph) adjacency and a rook (4-neighbour grid)
# adjacency were each hand-rolled locally in dozens of test files, one of
# which (test-nested-laplace-families.R) built its chain as a double matrix
# instead of integer. ONE definition here, sourced automatically by
# testthat's helper- convention.

# Path graph on `n` nodes: node i adjacent to i+1 only.
chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  adj
}

# Rook (4-neighbour, N/S/E/W) adjacency on an n_row x n_col grid. Cells are
# numbered row-major: cell (r, c) -> (r - 1) * n_col + c. `n_col` defaults to
# `n_row` for the common square-grid call, `rook_adj(g)`.
rook_adj <- function(n_row, n_col = n_row) {
  n <- n_row * n_col
  A <- matrix(0L, n, n)
  idx <- function(r, c) (r - 1L) * n_col + c
  for (r in seq_len(n_row)) for (c in seq_len(n_col)) {
    i <- idx(r, c)
    if (r < n_row) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < n_col) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}
