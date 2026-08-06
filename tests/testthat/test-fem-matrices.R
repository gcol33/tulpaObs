# =============================================================================
# test-fem-matrices.R - fem_matrices(), the SPDE mesh-assembly entry point
# re-exported from tulpaMesh (gcol33/tulpaObs#179).
#
# The Matern SPDE precision is built from these two matrices as
# Q = tau^2 (kappa^2 C + G)^2 / C, so a silently wrong C or G does not error --
# it produces a field whose prior is not the one the model claims. The
# reproducer in dev_notes/repro_tulpamesh_zero_triangles.R is exactly that
# failure: an all-zero C collapsed Q to a theta-independent ridge with no
# warning anywhere.
#
# So every assertion here is a mathematical identity of the FEM discretisation,
# checked against a quantity this file computes itself from the mesh geometry --
# never against a stored number, and never merely against a shape:
#
#   C_ij = int phi_i phi_j  =>  sum(C) = 1' C 1 = int 1 = the mesh area, which
#     the test recomputes from the vertices by the shoelace formula.
#   G_ij = int grad phi_i . grad phi_j  =>  G is symmetric PSD, G 1 = 0 (a
#     constant has no gradient), and for ANY linear field f = a x + b y + c
#     sampled at the nodes, f' G f = (a^2 + b^2) * area exactly, because the
#     P1 basis represents a linear function without error. That last identity
#     is what separates "some symmetric PSD matrix" from the actual stiffness
#     matrix.
#   A is barycentric, so its rows sum to 1, A %*% vertices returns the
#     observation coordinates, and A interpolates a linear field exactly.
# =============================================================================

.fem_mesh <- function() {
  co <- as.matrix(expand.grid(x = seq(0, 1, by = 0.25),
                              y = seq(0, 1, by = 0.25)))
  # max_edge is left unset: with it the triangulation collapses to zero
  # triangles at some settings (dev_notes/repro_tulpamesh_zero_triangles.R).
  tulpaMesh::tulpa_mesh(coords = co)
}

# Triangle areas by the shoelace formula, independent of the FEM assembly.
.fem_tri_areas <- function(mesh) {
  v <- mesh$vertices
  tv <- mesh$triangles
  a <- v[tv[, 1], , drop = FALSE]
  b <- v[tv[, 2], , drop = FALSE]
  cc <- v[tv[, 3], , drop = FALSE]
  0.5 * abs((b[, 1] - a[, 1]) * (cc[, 2] - a[, 2]) -
            (cc[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
}


test_that("the mass matrix integrates the mesh", {
  skip_if_no_tulpamesh()

  mesh <- .fem_mesh()
  area <- sum(.fem_tri_areas(mesh))
  fem <- fem_matrices(mesh)

  C <- as.matrix(fem$C)
  expect_equal(dim(C), c(mesh$n_vertices, mesh$n_vertices))
  expect_equal(fem$n_mesh, mesh$n_vertices)

  # sum(C) = 1' C 1 = int over the mesh of 1 = the total area.
  expect_equal(sum(C), area, tolerance = 1e-10)
  ones <- rep(1, nrow(C))
  expect_equal(as.numeric(t(ones) %*% C %*% ones), area, tolerance = 1e-10)

  # Every vertex owns a strictly positive share of the area, and C is a proper
  # inner-product matrix: symmetric and positive definite.
  expect_true(all(rowSums(C) > 0))
  expect_true(all(diag(C) > 0))
  expect_equal(C, t(C), tolerance = 1e-12)
  expect_gt(min(eigen(C, symmetric = TRUE)$values), 0)
})


test_that("the stiffness matrix is a Dirichlet form", {
  skip_if_no_tulpamesh()

  mesh <- .fem_mesh()
  area <- sum(.fem_tri_areas(mesh))
  G <- as.matrix(fem_matrices(mesh)$G)
  ev <- eigen(G, symmetric = TRUE)$values

  expect_equal(G, t(G), tolerance = 1e-12)

  # A constant field has zero gradient, so constants are in the null space.
  expect_lt(max(abs(G %*% rep(1, nrow(G)))), 1e-10 * max(abs(G)))

  # Positive semi-definite, with the constants as the ONLY null direction on a
  # connected mesh. A larger null space means disconnected components or
  # collapsed triangles; a negative eigenvalue means a broken assembly.
  expect_gt(min(ev), -1e-10 * max(ev))
  expect_equal(sum(ev < 1e-8 * max(ev)), 1L)

  # The identity that makes this the stiffness matrix and not just some
  # symmetric PSD matrix: P1 elements represent a linear field exactly, so for
  # f = a x + b y + c the discrete Dirichlet energy equals the analytic
  # integral of |grad f|^2 = a^2 + b^2 over the mesh.
  v <- mesh$vertices
  for (ab in list(c(2, -3), c(1, 0), c(-0.5, 0.25))) {
    f <- ab[1] * v[, 1] + ab[2] * v[, 2] + 7
    expect_equal(as.numeric(t(f) %*% G %*% f),
                 sum(ab^2) * area, tolerance = 1e-8)
  }
})


test_that("the projector is barycentric and reproduces coordinates", {
  skip_if_no_tulpamesh()

  mesh <- .fem_mesh()
  v <- mesh$vertices
  obs <- rbind(c(0.5, 0.25), c(0.13, 0.61), v[1, ], c(0.25, 0.25))
  fem <- fem_matrices(mesh, obs_coords = obs)
  A <- as.matrix(fem$A)

  expect_equal(dim(A), c(nrow(obs), mesh$n_vertices))

  # Barycentric weights: non-negative, summing to one, supported on at most the
  # three vertices of the containing triangle.
  expect_true(all(A >= -1e-12))
  expect_equal(rowSums(A), rep(1, nrow(obs)), tolerance = 1e-12)
  expect_true(all(rowSums(A != 0) <= 3L))

  # A point sitting on a vertex projects onto that vertex alone.
  expect_equal(sum(A[3, ] != 0), 1L)
  expect_equal(A[3, 1], 1, tolerance = 1e-12)

  # The strong statement: the weights reconstruct the observation coordinates,
  # and therefore interpolate any linear field exactly. A projector built
  # against the wrong triangle still has rows summing to one; it fails here.
  expect_equal(as.matrix(A %*% v), obs, tolerance = 1e-10,
               ignore_attr = TRUE)
  f <- 2 * v[, 1] - 3 * v[, 2] + 1
  expect_equal(as.numeric(A %*% f), 2 * obs[, 1] - 3 * obs[, 2] + 1,
               tolerance = 1e-10)

  # Supplying observations must not disturb the assembly itself.
  bare <- fem_matrices(mesh)
  expect_equal(as.matrix(fem$C), as.matrix(bare$C), tolerance = 1e-12)
  expect_equal(as.matrix(fem$G), as.matrix(bare$G), tolerance = 1e-12)

  # With no observations the projector is the identity on the mesh nodes.
  expect_equal(as.matrix(bare$A), diag(mesh$n_vertices), tolerance = 1e-12)
})


test_that("lumping conserves mass and reports the element areas", {
  skip_if_no_tulpamesh()

  mesh <- .fem_mesh()
  areas <- .fem_tri_areas(mesh)
  fem <- fem_matrices(mesh, lumped = TRUE)

  # Lumping adds the diagonal mass matrix and its parts; C and G are unchanged.
  expect_true(all(c("C0", "va", "ta") %in% names(fem)))
  expect_equal(fem$va, as.numeric(Matrix::rowSums(fem$C)), tolerance = 1e-12)
  expect_equal(sum(fem$va), sum(areas), tolerance = 1e-10)

  C0 <- as.matrix(fem$C0)
  expect_equal(max(abs(C0 - diag(diag(C0)))), 0)
  expect_equal(diag(C0), fem$va, tolerance = 1e-12, ignore_attr = TRUE)

  # The per-triangle areas the assembly used, against the shoelace formula.
  expect_equal(fem$ta, areas, tolerance = 1e-12)
  # Each vertex takes one third of each incident triangle, so the lumped total
  # is the mesh area -- lumping moves mass onto the diagonal without losing any.
  expect_equal(sum(diag(C0)), sum(as.matrix(fem$C)), tolerance = 1e-10)
})


test_that("parallel assembly does not change the matrices", {
  skip_if_no_tulpamesh()

  mesh <- .fem_mesh()
  serial <- fem_matrices(mesh)
  par <- fem_matrices(mesh, parallel = TRUE)

  # Per-triangle contributions are summed, so a different traversal order could
  # only move the last bits. Asserted with a tolerance and not with identity for
  # that reason; a real divergence is orders of magnitude larger.
  expect_equal(as.matrix(par$C), as.matrix(serial$C), tolerance = 1e-12)
  expect_equal(as.matrix(par$G), as.matrix(serial$G), tolerance = 1e-12)
})


test_that("a non-mesh input is refused", {
  skip_if_no_tulpamesh()

  expect_error(fem_matrices(list(vertices = cbind(0, 0))), "tulpa_mesh")

  # The zero-triangle mesh is the failure that used to return an all-zero C
  # (dev_notes/repro_tulpamesh_zero_triangles.R): it must error, not assemble.
  empty <- .fem_mesh()
  empty$triangles <- empty$triangles[0, , drop = FALSE]
  expect_error(fem_matrices(empty), "0 triangles")
})
