# =============================================================================
# test-joint-substrate-hyper-draws.R - gcol33/tulpaObs#359.
#
# `.tobs_joint_field_sd()` / `.tobs_joint_amp()` / `.tobs_joint_phi_at()` used
# to read a draw's sigma / alpha / phi_pos straight off `theta_grid[cells, ]`
# -- the raw outer-grid NODE coordinate, an atom rather than a marginal on a
# 3- to 15-node axis. `.tobs_joint_draws()` (and everything it feeds --
# predict(), the pointwise log-likelihood behind WAIC/LOO/CPO, ppc()) now
# reads `attr(tulpa::tulpa_posterior_draws(jf, ...), "theta")` instead: the
# SAME within-cell construction the fit's own `theta_ci_lo` / `theta_median` /
# `theta_ci_hi` are read from (`tulpa::tulpa_hyper_draws()`), so a draw's
# hyperparameter is continuized within its own cell rather than pinned to the
# cell's coordinate.
# =============================================================================

test_that("occu_cover joint draws continuize sigma/alpha/phi_pos, matching theta_ci", {
  skip_on_cran()
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "beta", phi = 25,
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 6789L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("beta"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L, engine = "joint")
  ))
  jf <- fit$joint_fit
  # This fit has no `share()` copy, so its outer grid is just the field
  # amplitude (`sigma`) and the beta arm's dispersion (`phi_pos`) -- no
  # `alpha` axis (the field is not coupled onto the cover arm).
  expect_true(all(c("sigma", "phi_pos") %in% colnames(jf$theta_grid)))
  n_grid <- nrow(jf$theta_grid)

  set.seed(1)
  D     <- tulpa::tulpa_posterior_draws(jf, idx = 1, n = 5000L)
  hyper <- attr(D, "theta")
  expect_true(is.matrix(hyper))
  expect_identical(colnames(hyper), colnames(jf$theta_grid))
  expect_equal(nrow(hyper), 5000L)

  for (ax in colnames(hyper)) {
    # A continuized draw takes far more distinct values than the grid has
    # nodes; the pre-fix atom read took at most n_grid distinct values.
    expect_gt(length(unique(hyper[, ax])), n_grid)
    # Reproduces the fit's own reported interval to Monte Carlo error (the
    # SAME within-cell construction, sampled rather than integrated).
    expect_equal(unname(median(hyper[, ax])), unname(jf$theta_median[[ax]]),
                 tolerance = 0.05)
    expect_equal(unname(stats::quantile(hyper[, ax], 0.025)),
                 unname(jf$theta_ci_lo[[ax]]), tolerance = 0.1)
    expect_equal(unname(stats::quantile(hyper[, ax], 0.975)),
                 unname(jf$theta_ci_hi[[ax]]), tolerance = 0.1)
  }

  # `.tobs_joint_draws()` (predict / WAIC / ppc's shared entry point) carries
  # the same continuized field amplitude and dispersion, not the grid atom.
  bundle <- tulpaObs:::.tobs_joint_draws(fit, n = 5000L)
  expect_gt(length(unique(bundle$blocks[[1L]]$amp_occ)), n_grid)
  expect_gt(length(unique(bundle$disp)), n_grid)
})
