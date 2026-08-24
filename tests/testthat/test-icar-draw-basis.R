# The simulated ICAR field must be a property of the seed, not of the
# linear-algebra library underneath (#279).
#
# A lattice ICAR precision has repeated eigenvalues in bulk, and inside a
# repeated block the eigenvector basis is fixed only up to an orthogonal
# rotation that no LAPACK build is obliged to resolve the same way. Drawing the
# field as V diag(1/sqrt(lambda)) z therefore makes the realisation depend on
# which basis the library returned: three areal fixtures reported different
# answers on the Linux runner than on Windows, on nominally the same seed and
# the same data, because the data was in fact a different draw.

grid_adj_icar <- function(nr, nc) {
  g <- expand.grid(r = seq_len(nr), c = seq_len(nc))
  n <- nrow(g)
  a <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1L) a[i, j] <- 1
    }
  }
  a
}

# Re-express an eigendecomposition in a different, equally valid basis by
# rotating inside each repeated-eigenvalue block. This is exactly the freedom a
# LAPACK build has, and no more: the result still satisfies Q V = V L.
rotate_within_ties <- function(eig, tol = 1e-9) {
  V <- eig$vectors
  lam <- eig$values
  grp <- cumsum(c(TRUE, diff(lam) < -tol * max(abs(lam))))
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) < 2L) next
    A <- matrix(stats::rnorm(length(idx)^2), length(idx))
    V[, idx] <- V[, idx] %*% qr.Q(qr(A))
  }
  list(values = lam, vectors = V)
}

test_that("the lattice fixtures' ICAR precision is degenerate enough to matter", {
  # Not a property of the draw -- the premise the draw has to survive. If a
  # future graph fixture were non-degenerate this test would be vacuous.
  for (d in list(c(5L, 5L), c(8L, 8L), c(9L, 9L))) {
    adj <- grid_adj_icar(d[1L], d[2L])
    ev <- sort(eigen(tulpaObs:::.occu_cover_icar_Q(adj), symmetric = TRUE,
                     only.values = TRUE)$values)
    ties <- sum(diff(ev) < 1e-9 * max(ev))
    expect_gt(ties, length(ev) / 4,
              label = sprintf("%dx%d grid: repeated eigenvalues (%d of %d)",
                              d[1L], d[2L], ties, length(ev)))
  }
})

test_that(".tobs_draw_icar_unit draws through the Cholesky, not an eigenbasis", {
  adj <- grid_adj_icar(5L, 5L)
  Q <- tulpaObs:::.occu_cover_icar_Q(adj)
  scale_q <- tulpaObs:::.occu_cover_icar_scale(adj)
  N <- nrow(Q)

  # The construction, written out. chol() of a positive-definite matrix is
  # unique, so this reference carries no basis freedom either.
  set.seed(404)
  Rc <- chol(Q + matrix(1 / N, N, N))
  ref <- backsolve(Rc, stats::rnorm(N))
  ref <- (ref - mean(ref)) / sqrt(scale_q)

  set.seed(404)
  got <- as.numeric(tulpaObs:::.tobs_draw_icar_unit(Q, scale_q))
  expect_equal(got, ref, tolerance = 1e-12)

  # And the test is not vacuous: the eigen construction this replaced consumes a
  # different stream and lands somewhere else entirely.
  eig <- eigen(Q, symmetric = TRUE)
  keep <- eig$values > 1e-8
  set.seed(404)
  z <- stats::rnorm(sum(keep))
  old <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                      (z / sqrt(eig$values[keep])))
  old <- (old - mean(old)) / sqrt(scale_q)
  expect_gt(max(abs(got - old)), 1e-3)
})

test_that("the ICAR draw is invariant to the eigenbasis a LAPACK build returns", {
  adj <- grid_adj_icar(6L, 6L)
  Q <- tulpaObs:::.occu_cover_icar_Q(adj)
  scale_q <- tulpaObs:::.occu_cover_icar_scale(adj)
  N <- nrow(Q)

  eig <- eigen(Q, symmetric = TRUE)
  set.seed(11)
  eig_b <- rotate_within_ties(eig)
  # The rotated basis is an exact eigendecomposition of the same Q.
  expect_lt(max(abs(Q %*% eig_b$vectors -
                      eig_b$vectors %*% diag(eig_b$values))), 1e-10)

  # Rebuilding Q from either basis and drawing gives the same field, because the
  # draw never consults a basis.
  Qa <- eig$vectors   %*% (eig$values   * t(eig$vectors))
  Qb <- eig_b$vectors %*% (eig_b$values * t(eig_b$vectors))
  set.seed(7); fa <- as.numeric(tulpaObs:::.tobs_draw_icar_unit(Qa, scale_q))
  set.seed(7); fb <- as.numeric(tulpaObs:::.tobs_draw_icar_unit(Qb, scale_q))
  expect_gt(stats::cor(fa, fb), 1 - 1e-8)

  # The eigen construction, on the same two bases and the same seed, does not
  # agree -- which is the whole reason the draw moved off it.
  keep <- eig$values > 1e-8
  draw_eig <- function(E) {
    set.seed(7)
    zz <- stats::rnorm(sum(keep))
    v <- as.numeric(E$vectors[, keep, drop = FALSE] %*%
                      (zz / sqrt(E$values[keep])))
    v - mean(v)
  }
  expect_lt(abs(stats::cor(draw_eig(eig), draw_eig(eig_b))), 0.9)
})

test_that("the ICAR draw carries unit geometric-mean marginal SD", {
  skip_on_cran()
  adj <- grid_adj_icar(5L, 5L)
  Q <- tulpaObs:::.occu_cover_icar_Q(adj)
  scale_q <- tulpaObs:::.occu_cover_icar_scale(adj)
  set.seed(3)
  F <- tulpaObs:::.tobs_draw_icar_unit(Q, scale_q, 20000L)
  # Sorbye-Rue: the geometric mean of the marginal variances is 1.
  gm <- exp(mean(log(apply(F, 1L, stats::var))))
  expect_lt(abs(gm - 1), 0.05)
  expect_lt(max(abs(colMeans(F))), 1e-8)   # every draw sums to zero
})
