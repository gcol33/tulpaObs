# test-ms-occu-cover-spatial.R - Stage-1 reduced-rank spatial-factor community
# occu_cover (gcol33/tulpa#67). Currently exercises the ground-truth simulator;
# the fitter recovery tests are added with the fitter increments.

.mscs_grid_adj <- function(nr, nc) {
  N <- nr * nc
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (c - 1L) * nr + r
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    k <- idx(r, c)
    if (r > 1L)  adj[k, idx(r - 1L, c)] <- 1L
    if (r < nr)  adj[k, idx(r + 1L, c)] <- 1L
    if (c > 1L)  adj[k, idx(r, c - 1L)] <- 1L
    if (c < nc)  adj[k, idx(r, c + 1L)] <- 1L
  }
  adj
}

test_that("simulate_ms_occu_cover_spatial returns well-formed K=1 community data", {
  adj <- .mscs_grid_adj(6L, 6L)        # N = 36 cells
  N <- nrow(adj); J <- 4L; S <- 8L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, J = J, seed = 11L)

  expect_identical(dim(sim$y),     c(N, J, S))
  expect_identical(dim(sim$y_pos), c(N, J, S))
  expect_length(sim$species, S)

  # Detections are 0/1; cover is present exactly at detected visits and positive
  # (lognormal support).
  expect_true(all(sim$y %in% c(0L, 1L)))
  expect_true(all(is.na(sim$y_pos) == (sim$y != 1L)))
  expect_true(all(sim$y_pos[!is.na(sim$y_pos)] > 0))

  # A detection implies presence: y == 1 only where the latent z == 1.
  for (s in seq_len(S)) {
    det_any <- rowSums(sim$y[, , s] == 1L) > 0
    expect_true(all(sim$truth$z[det_any, s] == 1L))
  }
})

test_that("the shared factor is unit-scaled and sign-anchored", {
  adj <- .mscs_grid_adj(7L, 7L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 6L, seed = 7L)
  w <- sim$truth$w

  expect_length(w, nrow(adj))
  expect_lt(abs(mean(w)), 1e-8)              # centred on the constrained space
  expect_gt(stats::sd(w), 0.3)               # genuine spatial amplitude
  # Sorbye-Rue unit-marginal scale: the field SD is O(1), not collapsed/blown up.
  expect_lt(stats::sd(w), 3)
  # Canonical sign anchor: reference species loading is positive.
  expect_gt(sim$truth$L[1L], 0)
})

test_that("loadings give per-species range heterogeneity (not the naive shared map)", {
  adj <- .mscs_grid_adj(8L, 8L)
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = 12L, sd_load = 1.2,
                                        seed = 99L)
  w <- sim$truth$w; L <- sim$truth$L

  # The spatial contribution to each species' occupancy predictor is L_s * w.
  # With K = 1 these are collinear in shape but differ in sign and amplitude:
  # species with opposite-sign loadings have ANTI-correlated spatial maps, which
  # a single shared field + intercept RE (naive structure) cannot produce.
  expect_gt(diff(range(L)), 0.5)             # loadings span a real range
  expect_true(any(L > 0) && any(L < 0))      # both signs present
  pos <- which(L > 0)[1L]; neg <- which(L < 0)[1L]
  expect_lt(stats::cor(L[pos] * w, L[neg] * w), -0.99)  # opposite maps

  # Per-cell psi genuinely varies across species (maps are not identical).
  cell_sd <- apply(sim$truth$psi, 1L, stats::sd)
  expect_gt(mean(cell_sd), 0.02)
})
