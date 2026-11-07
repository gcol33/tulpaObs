# Disconnected areal graphs: one sum-to-zero per connected component on every
# engine that fits an intrinsic field in-tree (#439), and BYM2 scaled per
# component (#438).

.two_island_adj <- function(n1 = 12L, n2 = 18L) {
  chain <- function(n) {
    A <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) A[i, i + 1L] <- A[i + 1L, i] <- 1
    A
  }
  n <- n1 + n2
  A <- matrix(0, n, n)
  A[seq_len(n1), seq_len(n1)] <- chain(n1)
  A[n1 + seq_len(n2), n1 + seq_len(n2)] <- chain(n2)
  A
}
.island_sums <- function(f, n1 = 12L) {
  c(sum(f[seq_len(n1)]), sum(f[-seq_len(n1)]))
}

test_that(".bym2_scaling scales each component to the reference variance (#438)", {
  A <- .two_island_adj()
  sc <- tulpaObs:::.bym2_scaling(A)
  expect_length(sc$node_prec, 30L)
  expect_length(unique(round(sc$node_prec, 12)), 2L)
  # Under diag(node_prec) Q every component's geometric-mean marginal variance
  # (generalised inverse) equals the reference scale.
  Q  <- diag(rowSums(A)) - A
  Qs <- sc$node_prec * Q
  for (cc in list(1:12, 13:30)) {
    expect_equal(tulpaObs:::.bym2_scale_from_Q(Qs[cc, cc]), sc$scale,
                 tolerance = 1e-8)
  }
  # The reference is the larger component, scaled as on its own.
  expect_equal(sc$scale, tulpaObs:::.bym2_scale(A[13:30, 13:30]))
  # A connected graph is unchanged.
  B <- .two_island_adj(); B[12, 13] <- B[13, 12] <- 1
  sb <- tulpaObs:::.bym2_scaling(B)
  expect_null(sb$node_prec)
  expect_identical(sb$scale, tulpaObs:::.bym2_scale(B))
  # The term carries it; a user scale.factor is one scale for the whole graph.
  t1 <- tulpaObs:::.tobs_term_bym2(graph = A)
  expect_equal(t1$node_prec, sc$node_prec)
  expect_null(tulpaObs:::.tobs_term_bym2(graph = A, scale.factor = 0.5)$node_prec)
  # Component labels for the C++ drivers: 0-based, empty when connected.
  expect_identical(tulpaObs:::.tobs_node_comp_arg(A), c(rep(0L, 12), rep(1L, 18)))
  expect_identical(tulpaObs:::.tobs_node_comp_arg(B), integer(0))
})

test_that("areal-BFGS blocks constrain per component (#439)", {
  A <- .two_island_adj()
  blk <- tulpaObs:::.areal_field_car(A, "icar", seq_len(30), 30L)
  expect_length(blk$constrain_sets, 2L)
  f <- blk$center(seq_len(30) + 0)
  expect_equal(.island_sums(f), c(0, 0))
  b2 <- tulpaObs:::.areal_field_bym2(A, 1, seq_len(30), 30L)
  g  <- b2$center(c(seq_len(30) + 0, rep(1, 30)))
  expect_equal(.island_sums(g[1:30]), c(0, 0))
  expect_equal(g[31:60], rep(1, 30))
  expect_length(tulpaObs:::.areal_field_car(A, "car_proper", seq_len(30), 30L)$constrain_sets, 0L)
})

test_that("C++ count, community and occupancy drivers centre each island (#439)", {
  skip_on_cran()
  A <- .two_island_adj()
  ctl <- list(verbose = FALSE, progress = FALSE)

  sa <- simulate_abun(N = 30, J = 3, n.abund.covs = 1, n.det.covs = 1,
                      beta.lambda = c(log(5), 0.3), beta.p = c(0.2, 0), seed = 5)
  fa <- tobs(~ abund_cov1 + icar(graph = A), data = sa$data, family = abun(),
             detection = ~ 1, y = sa$y, method = "nested_laplace", control = ctl)
  expect_equal(.island_sums(fa$spatial_field), c(0, 0), tolerance = 1e-6)
  expect_true(all(is.finite(fa$sds)))

  fb <- tobs(~ abund_cov1 + bym2(graph = A), data = sa$data, family = abun(),
             detection = ~ 1, y = sa$y, method = "nested_laplace", control = ctl)
  expect_true(all(is.finite(fb$means)))
  expect_true(all(is.finite(fb$sds)))

  sm <- simulate_ms_abun(n.species = 6, N = 30, J = 3, n.abund.covs = 1,
                         n.det.covs = 1, seed = 5)
  fm <- tobs(~ abund_cov1 + icar(graph = A), detection = ~ det_cov1,
             family = ms_abun(), data = sm$data, y = sm$y, species = sm$species,
             method = "nested_laplace", control = ctl)
  expect_equal(.island_sums(fm$spatial_field), c(0, 0), tolerance = 1e-6)

  so <- simulate_ms_occu(N = 30, J = 4, n.species = 6, seed = 5)
  fo <- tobs(~ x + icar(graph = A), data = so$data, family = ms_occu(),
             detection = ~ 1, y = so$y, species = paste0("sp", seq_len(6)),
             method = "nested_laplace", control = ctl)
  expect_equal(.island_sums(fo$spatial_field), c(0, 0), tolerance = 1e-6)
})
