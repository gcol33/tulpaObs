# =============================================================================
# test-occu-cover-spatial.R - method / spatial-term pairing gates for
# occu_cover(): a spatial term needs method = "nested_laplace", and that
# method needs a spatial term. Recovery on the spatial path is gated in
# test-occu-cover-joint-coupled.R and test-occu-cover-coupling.R.
# =============================================================================


test_that("occu_cover() rejects laplace with spatial term", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 1, alpha = 1, seed = 1L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  expect_error(
    suppressWarnings(tobs(
      formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
      family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = ~ pos_cov1,
      y = od$y, y_pos = y_pos, visits = od$det.covs,
      method = "laplace", control = list(verbose = FALSE)
    )),
    "non-spatial"
  )

  expect_error(
    tobs(formula = ~ occ_cov1, data = cell_dat,
         family = occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = od$y, y_pos = y_pos, visits = od$det.covs,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "spatial term"
  )
})
