# =============================================================================
# test-occu-cover-batch-fused.R - the bit-identity gate for the FUSED batched
# C++ driver. Distinct from test-occu-cover-batch.R, which gates the Stage-1
# looped backend behind the public tobs_batch API.
#
# This drives tulpa's batched joint nested-Laplace entry
# (cpp_nested_laplace_joint_multi_batch, via the tulpa_nl_joint_batch wrapper)
# on a multi-response occu_cover problem and asserts every species' per-cell
# modes + log-marginals are BIT-IDENTICAL to an independent single-species fit
# at the same outer grid (tulpa_nl_joint_single). The fused path only
# reorganises the work -- one design pass, B block-diagonal species solves --
# so the statistics must match exactly.
#
# Field-size coverage: the ICAR sum-to-zero penalty has a dense rank-1 (1 1')
# Hessian. The single-species sparse path stores it exactly for a field up to
# S2Z_DENSIFY_MAX = 256 nodes (densified block) and folds it in at solve time
# (Woodbury) above that. The per-species batched solve must be bit-identical to
# the independent fit in BOTH regimes, so the gate runs at a small densified
# field (N = 16) AND a large field that exercises the s2z rank-1 path on the
# single-species oracle (N = 300, n_x > 256).
# =============================================================================

# Build responses + an ICAR (multi-block) prior for one species of a tiny
# occu_cover problem, replicating the joint fitter's setup slice. Fixed
# sigma_pos (no phi grid); the batched path supplies it via phi_batch.
.batch_build_one <- function(cell_dat, adj, N, J, det_covs, y_det, y_pos) {
  vd_det <- tulpaObs:::.normalize_visits(det_covs, ~ det_cov1,
                                         n_sites = N, max_visits = J)
  vd_pos <- tulpaObs:::.normalize_visits(det_covs, ~ pos_cov1,
                                         n_sites = N, max_visits = J)
  model <- tulpaObs:::.tobs_build_occu_cover(
    occ_formula = ~ occ_cov1, det_formula = vd_det$det_formula,
    pos_formula = vd_pos$det_formula, data = cell_dat, y = y_det, y_pos = y_pos,
    positive = "lognormal",
    det_visit_formula = vd_det$det_visit_formula, det_visit_data = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula, pos_visit_data = vd_pos$visits)
  model$site_cell <- seq_len(N); model$n_cells <- N
  pos_vals  <- model$y_pos[model$valid & model$y == 1L]
  sigma_pos <- max(stats::sd(log(pos_vals)), 0.05) + 0.05
  arms_out  <- tulpaObs:::.occu_cover_build_joint_arms(
    model = model, sigma_pos_init = sigma_pos,
    alpha_axis = tulpaObs:::.tobs_alpha_axis(c(0, 1.0)),
    positive = "lognormal", multi = FALSE, n_cells = N,
    site_cell = seq_len(N), cover_aggregate = "none")
  responses <- arms_out$responses
  csr <- tulpaObs:::.occu_cover_adj_to_csr(adj)
  prior <- list(type = "icar", n_spatial_units = csr$n_spatial_units,
                adj_row_ptr = csr$adj_row_ptr, adj_col_idx = csr$adj_col_idx,
                n_neighbors = csr$n_neighbors, sigma_grid = c(0.5, 1.0),
                spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx)))
  list(responses = responses, prior = prior, sigma_pos = sigma_pos,
       n_arms = length(responses))
}

# One full bit-identity check at field size N: two distinct species, each fit
# independently (the oracle) and together through the batched driver; assert
# per-species modes + log-marginals + the stored inner-precision Q match
# cell-for-cell, and that the two species are genuinely distinct fits.
#
# `offset_arm`: when set to a coupled-arm index, thread a nonzero,
# heterogeneous per-observation exposure offset through that arm. The batched
# compute_eta_species must add it exactly as the single-species compute_eta
# does; with the offset dropped the batched fit disagrees with the oracle
# cell-for-cell.
.batch_fused_identity_at <- function(N, J = 4L, offset_arm = NULL) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) { if (s > 1L) adj[s, s-1L] <- 1L; if (s < N) adj[s, s+1L] <- 1L }

  mk <- function(seed) {
    sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", adj = adj,
                               sigma = 0.8, alpha = 1.0, seed = seed)
    long <- data.frame(site_id = rep(seq_len(N), each = J),
                       visit = rep(seq_len(J), times = N),
                       y = as.vector(t(sim$y)),
                       det_cov1 = sim$visit_data$det_cov1,
                       pos_cov1 = sim$visit_data$pos_cov1)
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    yp <- sim$y_pos; yp[is.na(yp)] <- 0
    list(od = od, sim = sim, y = od$y, y_pos = yp)
  }
  d1 <- mk(101L); d2 <- mk(202L)
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), d1$sim$data)

  b0 <- .batch_build_one(cell_dat, adj, N, J, d1$od$det.covs, d1$y,        d1$y_pos)
  b1 <- .batch_build_one(cell_dat, adj, N, J, d1$od$det.covs, d2$sim$y,
                         { yp <- d2$sim$y_pos; yp[is.na(yp)] <- 0; yp })

  # Exposure offset on a coupled arm. The batch design is species-invariant
  # (arms_list = b0$responses rides every species), so the same offset must sit
  # on b1's matching arm for the independent oracle of species 2 to see the same
  # design. Both species' arms share N (the design is response-invariant), so one
  # vector serves both.
  if (!is.null(offset_arm)) {
    nk  <- length(b0$responses[[offset_arm]]$y)
    off <- 0.6 * cos(seq_len(nk)) - 0.2
    b0$responses[[offset_arm]]$offset <- off
    b1$responses[[offset_arm]]$offset <- off
  }

  # Independent single-species fits at the dense grid (the oracle). store_Q so
  # the per-grid inner-covariance precision is returned and can be compared.
  s0 <- tulpa:::tulpa_nl_joint_single(b0$responses, b0$prior, copy = NULL,
          max_iter = 200L, tol = 1e-8, cell_coupling = "occu_cover_lognormal",
          store_Q = TRUE)
  s1 <- tulpa:::tulpa_nl_joint_single(b1$responses, b1$prior, copy = NULL,
          max_iter = 200L, tol = 1e-8, cell_coupling = "occu_cover_lognormal",
          store_Q = TRUE)

  # Batched fit: arms_list = species-0 responses (design); y_batch carries both.
  na <- b0$n_arms
  y_batch <- vector("list", na)
  y_batch[[2]] <- cbind(as.numeric(b0$responses[[2]]$y),
                        as.numeric(b1$responses[[2]]$y))
  y_batch[[3]] <- cbind(as.numeric(b0$responses[[3]]$y),
                        as.numeric(b1$responses[[3]]$y))
  phi_batch <- matrix(0, na, 2L)
  phi_batch[3, ] <- c(b0$sigma_pos, b1$sigma_pos)

  bat <- tulpa:::tulpa_nl_joint_batch(b0$responses, b0$prior, copy = NULL,
          n_batch = 2L, y_batch = y_batch, phi_batch = phi_batch,
          max_iter = 200L, tol = 1e-8, cell_coupling = "occu_cover_lognormal",
          store_Q = TRUE)

  expect_length(bat$per_species, 2L)

  # Per-species bit-identity (the fused path only reorganises the same math):
  # modes, log-marginals, AND the per-grid inner-covariance precision Q
  # (store_Q) must match the independent fits cell-for-cell.
  for (pair in list(list(s0, bat$per_species[[1]]),
                    list(s1, bat$per_species[[2]]))) {
    single <- pair[[1]]; bs <- pair[[2]]
    expect_equal(as.numeric(bs$modes), as.numeric(single$modes), tolerance = 1e-9)
    expect_equal(as.numeric(bs$log_marginal), as.numeric(single$log_marginal),
                 tolerance = 1e-9)
    expect_equal(bs$Q_csc_n, single$Q_csc_n)
    expect_equal(length(bs$Q_csc_x_per_grid), length(single$Q_csc_x_per_grid))
    for (k in seq_along(single$Q_csc_x_per_grid)) {
      expect_equal(as.numeric(bs$Q_csc_x_per_grid[[k]]),
                   as.numeric(single$Q_csc_x_per_grid[[k]]), tolerance = 1e-9)
      expect_equal(as.integer(bs$Q_csc_i_per_grid[[k]]),
                   as.integer(single$Q_csc_i_per_grid[[k]]))
    }
  }

  # Sanity: the two species are genuinely different fits.
  expect_gt(max(abs(as.numeric(bat$per_species[[1]]$modes) -
                    as.numeric(bat$per_species[[2]]$modes))), 0.1)
}


test_that("fused batched driver is per-species bit-identical to independent fits", {
  skip_on_cran()
  skip_if_fast()
  .batch_fused_identity_at(N = 16L)
})


test_that("fused batched driver threads a nonzero per-observation offset", {
  # Regression: the batched compute_eta_species built the per-arm linear
  # predictor without adding pa.offset, silently dropping any exposure/effort
  # offset on a coupled arm -- while the single-species compute_eta added it -- so
  # a joint fit with an offset disagreed cell-for-cell with the independent oracle
  # and violated the batch file's own bit-identity invariant. Threading a nonzero
  # heterogeneous offset on the positive (cover) arm and requiring the batched
  # per-species modes / log-marginals / inner precision to stay bit-identical to
  # the single-species fits (which carry the offset) pins the fix. Before the fix
  # this fails; after it, byte-identity holds.
  skip_on_cran()
  skip_if_fast()
  .batch_fused_identity_at(N = 16L, offset_arm = 3L)
})


test_that("fused batched driver is bit-identical at a large field (s2z rank-1 regime)", {
  # N = 300 > S2Z_DENSIFY_MAX (256), so the single-species sparse oracle folds
  # the ICAR sum-to-zero rank-1 in at solve time (Woodbury) rather than storing
  # the densified 1 1'. The batched per-species solve must still match the
  # independent fit cell-for-cell -- this is the regime flags as the one a
  # sparse-native batched driver must get right (the 2nd-species s2z correctness
  # contract).
  skip_on_cran()
  skip_if_fast()
  .batch_fused_identity_at(N = 300L)
})
