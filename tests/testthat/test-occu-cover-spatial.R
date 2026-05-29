# =============================================================================
# test-occu-cover-spatial.R - v2 nested-Laplace spatial path for occu_cover().
#
# What this test gates (v3 nested-Laplace, with z profiled out):
#   - Pipeline runs end-to-end across 10 seeds and produces finite estimates
#   - Slope recovery (psi / p / cover-arm) within 0.35 on average
#   - Dispersion recovery (log_sigma_pos) within 0.15
#   - Field SHAPE recovery: average cor(z_hat, f_true) > 0.85
#   - alpha recovery: within 0.50 (mean), CI coverage informal
#   - sigma recovery: within 0.40 (mean)
#
# v2's joint Laplace had a (z, alpha, sigma) ridge so couldn't gate alpha /
# sigma individually; v3 profiles z out per outer (alpha, sigma) candidate
# and breaks the ridge cleanly.
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


test_that("occu_cover() v3 nested-Laplace recovers slopes, hypers, field (10 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 10L
  N <- 100L; J <- 6L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }

  beta_occ_truth <- c(stats::qlogis(0.4),  0.7)
  beta_p_truth   <- c(0.0,                  0.8)
  beta_pos_truth <- c(log(0.20),           -0.4)
  sigma_pos_truth <- 0.35
  sigma_truth   <- 1.0
  alpha_truth   <- 1.0

  est_psi_x <- est_p_x <- est_pos_x <- est_sigma_pos <- field_cor <-
    est_alpha <- est_sigma <- rep(NA_real_, n_seeds)
  se_psi_x  <- se_p_x  <- se_pos_x  <- rep(NA_real_, n_seeds)
  conv <- logical(n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, sigma_pos = sigma_pos_truth,
      positive = "lognormal", adj = adj,
      sigma = sigma_truth, alpha = alpha_truth, seed = 4000L + s
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

    fit <- tryCatch(
      suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("lognormal"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(engine = "v3_nested",
                       verbose = FALSE, max.iter = 1500L)
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    conv[s] <- isTRUE(fit$convergence$converged)

    est_psi_x[s]     <- fit$means["psi_occ_cov1"]
    est_p_x[s]       <- fit$means["p_det_cov1"]
    est_pos_x[s]     <- fit$means["pos_pos_cov1"]
    est_sigma_pos[s] <- exp(fit$means["log_sigma_pos"])
    est_alpha[s]     <- fit$means["alpha"]
    est_sigma[s]     <- exp(fit$means["log_sigma"])
    se_psi_x[s]      <- fit$sds["psi_occ_cov1"]
    se_p_x[s]        <- fit$sds["p_det_cov1"]
    se_pos_x[s]      <- fit$sds["pos_pos_cov1"]
    field_cor[s]     <- abs(cor(fit$spatial_field, sim$truth$f))
  }

  # Inclusion: any fit that produced finite estimates. v2's joint Laplace on
  # 100+ params often hits BFGS max-iter without strict convergence even when
  # the recovery is fine; gate on usability, not on opt$convergence.
  ok <- is.finite(est_psi_x) & is.finite(est_pos_x)
  expect_gte(sum(ok), 7L)

  # Slope point recovery (v3 nested-Laplace gates).
  expect_lt(abs(mean(est_psi_x[ok])     - beta_occ_truth[2L]), 0.35)
  expect_lt(abs(mean(est_p_x[ok])       - beta_p_truth[2L]),   0.35)
  expect_lt(abs(mean(est_pos_x[ok])     - beta_pos_truth[2L]), 0.20)
  expect_lt(abs(mean(est_sigma_pos[ok]) - sigma_pos_truth),    0.15)

  # Hyperparameter recovery (v3 breaks the v2 ridge).
  expect_lt(abs(mean(est_alpha[ok]) - alpha_truth), 0.50)
  expect_lt(abs(mean(est_sigma[ok]) - sigma_truth), 0.40)

  # Field SHAPE recovery.
  expect_gt(mean(field_cor[ok]), 0.85)
})
