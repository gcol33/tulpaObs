# Regression test for gcol33/tulpa#69.
#
# The intrinsic-ICAR sum-to-zero log-det in the coupled occu_cover joint Laplace
# was reconstructed via the matrix-determinant lemma from the Hessian A with the
# 1 1' penalty removed. A pins the field's constant direction only through the
# 1e-10 uniform ridge, so 1' A^{-1} 1 ~ 1 / ridge and log|A| ~ log(ridge); the
# lemma sum log|A| + log(1 + coef * 1' A^{-1} 1) is then a (-large)+(+large)
# cancellation that loses most of its digits once the spatial field is large.
# That biased the outer hyperparameter-integration weights (log_marginal) of any
# occu_cover fit whose field crosses the densify threshold (256 units); modes are
# unaffected (the Woodbury step is exact). tulpa's s2z_log_det_direct factors the
# well-conditioned H + coef 1 1' directly instead.
#
# This guards the rank-1 (> 256 node) log-det against the dense densify path at
# N = 300, where the two diverged by ~2.7 in log-marginal before the fix. The
# comparison is within a single simulated data set (only TULPA_S2Z_DENSIFY_MAX
# differs), so it is robust to engine/version drift. Generic single-arm fields do
# NOT reach this regime -- the cancellation needs the coupled occu_cover
# structure, where occupancy/detection zeros leave field units near-null-curvature.

test_that("occu_cover rank-1 s2z log-det matches densify above the 256-node field (gcol33/tulpa#69)", {
  skip_on_cran()
  old <- Sys.getenv("TULPA_S2Z_DENSIFY_MAX", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("TULPA_S2Z_DENSIFY_MAX")
          else Sys.setenv(TULPA_S2Z_DENSIFY_MAX = old), add = TRUE)

  N <- 300L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) { if (s > 1L) adj[s, s - 1L] <- 1L; if (s < N) adj[s, s + 1L] <- 1L }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", adj = adj,
                             sigma = 0.8, alpha = 1.0, seed = 101L)
  long <- data.frame(site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)), det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  yp <- sim$y_pos; yp[is.na(yp)] <- 0
  vd_det <- tulpaObs:::.normalize_visits(od$det.covs, ~ det_cov1, n_sites = N, max_visits = J)
  vd_pos <- tulpaObs:::.normalize_visits(od$det.covs, ~ pos_cov1, n_sites = N, max_visits = J)
  model <- tulpaObs:::.tobs_build_occu_cover(
    occ_formula = ~occ_cov1, det_formula = vd_det$det_formula, pos_formula = vd_pos$det_formula,
    data = cell_dat, y = od$y, y_pos = yp, positive = "lognormal",
    det_visit_formula = vd_det$det_visit_formula, det_visit_data = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula, pos_visit_data = vd_pos$visits)
  model$site_cell <- seq_len(N); model$n_cells <- N
  pv <- model$y_pos[model$valid & model$y == 1L]; spp <- max(stats::sd(log(pv)), 0.05) + 0.05
  arms_out <- tulpaObs:::.occu_cover_build_joint_arms(
    model = model, sigma_pos_init = spp, alpha_grid = c(0, 1.0), positive = "lognormal",
    multi = FALSE, n_cells = N, site_cell = seq_len(N), cover_aggregate = "none")
  responses <- arms_out$responses
  csr <- tulpaObs:::.occu_cover_adj_to_csr(adj)
  prior <- list(type = "icar", n_spatial_units = csr$n_spatial_units,
                adj_row_ptr = csr$adj_row_ptr, adj_col_idx = csr$adj_col_idx,
                n_neighbors = csr$n_neighbors, sigma_grid = c(0.5, 1.0),
                spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx)))

  fit_at <- function(densify_max) {
    Sys.setenv(TULPA_S2Z_DENSIFY_MAX = densify_max)
    tulpa:::tulpa_nl_joint_single(responses, prior, copy = NULL, max_iter = 200L,
                                  tol = 1e-8, cell_coupling = "occu_cover_lognormal",
                                  store_Q = TRUE)
  }
  dns <- fit_at("100000")  # store the full 1 1' (dense densify path)
  r1  <- fit_at("0")        # fold 1 1' in at solve time (the > 256-node path)

  # The field is in the rank-1 regime (n_x > S2Z_DENSIFY_MAX = 256).
  expect_gt(ncol(dns$modes), 256L)
  # Modes are unaffected by the log-det treatment (the Woodbury step is exact).
  expect_equal(max(abs(as.numeric(dns$modes) - as.numeric(r1$modes))), 0,
               tolerance = 1e-6)
  # With s2z_log_det_direct the rank-1 log-marginal matches densify per grid cell
  # (this gap was ~2.7 before the fix).
  expect_equal(r1$log_marginal, dns$log_marginal, tolerance = 1e-4)
})
