# test-ms-occu-cover-spatial-calibration.R - multi-seed identifiability +
# rare-species calibration for the Stage-1 reduced-rank spatial-factor community
# occu_cover (gcol33/tulpa#67). A single-seed recovery proves the fitter can
# work; these aggregate over seeds (recovery is not a lucky draw, the (L, w) sign
# anchor resolves consistently) and check that the per-species occupancy psi
# posterior intervals are calibrated -- the property that justifies the spatial
# class. The held-out variant blanks a random subset of cells and asks whether
# the shared ICAR field predicts the unobserved cells with nominal coverage,
# including for genuinely rare species.

.mcal_grid_adj <- function(nr, nc) {
  N <- nr * nc
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (c - 1L) * nr + r
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    k <- idx(r, c)
    if (r > 1L) adj[k, idx(r - 1L, c)] <- 1L
    if (r < nr) adj[k, idx(r + 1L, c)] <- 1L
    if (c > 1L) adj[k, idx(r, c - 1L)] <- 1L
    if (c < nc) adj[k, idx(r, c + 1L)] <- 1L
  }
  adj
}

.mcal_fit <- function(adj, sim, y = sim$y, y_pos = sim$y_pos) {
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = y, y_pos = y_pos, positive = "lognormal",
    species = sim$species, adj = adj)
  tulpaObs:::.tobs_fit_ms_occu_cover_spatial(model, sd_L = 1.2, max.em = 25L,
                                             tol = 1e-3)
}

test_that("recovery holds across seeds (field, loadings, sign anchor)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mcal_grid_adj(8L, 8L)
  seeds <- c(2024L, 4040L, 5151L, 77L, 303L)
  fcor <- lcor <- numeric(length(seeds))
  for (k in seq_along(seeds)) {
    sim <- simulate_ms_occu_cover_spatial(adj, n_species = 16L, J = 6L,
                                          sd_occ = 0.5, sd_load = 1.2,
                                          sigma_pos = 0.4, seed = seeds[k])
    fit <- .mcal_fit(adj, sim)
    # The sign anchor (sp1 loading made positive) must resolve the (L, w) ->
    # (-L, -w) symmetry the same way the simulator does, so the *signed*
    # correlations (not their absolute values) are the recovery measure.
    fcor[k] <- stats::cor(fit$w, sim$truth$w)
    lcor[k] <- stats::cor(fit$L, sim$truth$L)
  }
  # No seed collapses the latent factor or flips the anchor, and the typical
  # recovery clears the conditional-recovery milestone.
  expect_gt(min(fcor), 0.5)
  expect_gt(min(lcor), 0.8)
  expect_gt(stats::median(fcor), 0.75)
})

test_that("held-out cell psi intervals are calibrated, including rare species", {
  skip_on_cran()
  skip_if_fast()
  adj <- .mcal_grid_adj(9L, 9L); N <- nrow(adj)
  seeds <- c(2024L, 4040L, 77L, 909L)
  cov_held <- numeric(length(seeds))
  rare_hits <- rare_tot <- 0L

  for (k in seq_along(seeds)) {
    seed <- seeds[k]
    # Genuinely rare taxa: low occupancy intercept + wide community SD pushes a
    # handful of species to prevalence < 0.15.
    sim <- simulate_ms_occu_cover_spatial(adj, n_species = 18L, J = 6L,
                                          mu_occ = c(stats::qlogis(0.18), 0.7),
                                          sd_occ = 0.9, sd_load = 1.2,
                                          sigma_pos = 0.4, seed = seed)
    prev_s <- colMeans(sim$truth$z)

    # Blank all detection / cover observations at a random fifth of the cells;
    # their occupancy field is then pinned only by the ICAR neighbours.
    set.seed(seed * 13L + 5L)
    held <- sample.int(N, size = round(0.2 * N))
    y <- sim$y; y_pos <- sim$y_pos
    y[held, , ] <- NA; y_pos[held, , ] <- NA

    fit <- .mcal_fit(adj, sim, y = y, y_pos = y_pos)
    set.seed(seed * 7L + 1L)
    model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
      occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
      data = sim$data, y = y, y_pos = y_pos, positive = "lognormal",
      species = sim$species, adj = adj)
    psi <- tulpaObs:::.ms_ocs_psi_posterior(model, fit, n_draws = 400L)
    tr  <- sim$truth$psi
    S   <- dim(psi)[3L]

    cov_s <- numeric(S)
    for (s in seq_len(S)) {
      lo <- apply(psi[, , s], 2L, stats::quantile, 0.025)
      hi <- apply(psi[, , s], 2L, stats::quantile, 0.975)
      ins <- tr[, s] >= lo & tr[, s] <= hi
      cov_s[s] <- mean(ins[held])
      if (prev_s[s] <= 0.15) {
        rare_hits <- rare_hits + sum(ins[held])
        rare_tot  <- rare_tot  + length(held)
      }
    }
    cov_held[k] <- mean(cov_s)
    # Per-seed: held-out 95% intervals cover the latent psi at >= 0.80 (Laplace
    # is mildly under-dispersed for binary data, so a few points below 0.95 is
    # expected, not a defect).
    expect_gt(cov_held[k], 0.80)
  }

  # Aggregate held-out coverage is near nominal, and rare species (the regime
  # the shared field exists to rescue) are covered well above chance.
  expect_gt(mean(cov_held), 0.85)
  expect_gt(rare_hits / rare_tot, 0.75)
})
