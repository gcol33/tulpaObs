# Shared proper-CAR log-determinant (#182). The value is a weight on
# hyperparameter grid points, so a wrong one moves a reported rho posterior
# without changing any shape or raising an error. Cheap closed-form checks, no
# fits, so they run in every tier.

grid_adj <- function(nr, nc) {
  n <- nr * nc
  A <- matrix(0, n, n)
  id <- function(i, j) (j - 1L) * nr + i
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    if (i < nr) { A[id(i, j), id(i + 1L, j)] <- 1; A[id(i + 1L, j), id(i, j)] <- 1 }
    if (j < nc) { A[id(i, j), id(i, j + 1L)] <- 1; A[id(i, j + 1L), id(i, j)] <- 1 }
  }
  A
}

test_that(".tobs_car_logdet_Q equals log|D - rho W| by an independent route", {
  A <- grid_adj(5L, 5L)
  D <- diag(rowSums(A))
  rho_grid <- c(0.05, 0.3, 0.6, 0.9, 0.99)

  got <- .tobs_car_logdet_Q(A, rho_grid)
  ref <- vapply(rho_grid, function(rho) {
    as.numeric(determinant(D - rho * A, logarithm = TRUE)$modulus)
  }, numeric(1))

  expect_equal(got, ref, tolerance = 1e-10)
  expect_true(all(is.finite(got)))
  # log|Q(rho)| decreases in rho: the field prior loosens as rho -> 1.
  expect_true(all(diff(got) < 0))
})

test_that(".tobs_car_logdet_Q returns -Inf where Q(rho) is not positive definite", {
  # Two disconnected isolated nodes: Q is singular for every rho.
  A <- matrix(0, 3L, 3L)
  A[1L, 2L] <- 1; A[2L, 1L] <- 1
  expect_identical(.tobs_car_logdet_Q(A, 0.5), -Inf)
})

test_that("csr_to_adjacency inverts adjacency_to_csr", {
  set.seed(4)
  for (n in c(6L, 15L)) {
    A <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      if (stats::runif(1) < 0.3) { A[i, j] <- 1; A[j, i] <- 1 }
    }
    expect_identical(csr_to_adjacency(adjacency_to_csr(A), n), A)
  }
})

test_that("the CSR route and the dense route give the same log|Q(rho)|", {
  # ms_abun() carries only the CSR arrays; ms_occu()/ms_abun()'s EM path carry
  # the dense graph. Both must land on the same grid weights.
  A <- grid_adj(4L, 6L)
  rho_grid <- c(0.2, 0.5, 0.8)
  expect_identical(
    .tobs_car_logdet_Q(csr_to_adjacency(adjacency_to_csr(A), nrow(A)), rho_grid),
    .tobs_car_logdet_Q(A, rho_grid)
  )
})
