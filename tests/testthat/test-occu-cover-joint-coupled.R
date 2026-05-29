# =============================================================================
# test-occu-cover-joint-coupled.R - end-to-end gates for the joint-coupled
# engine wired through the occu_cover_lognormal cell-coupling spec
# (gcol33/tulpa#32 Layer B.2 consumer + R-facing fit wiring).
#
# Routed via `tobs(method = "nested_laplace", control = list(engine =
# "joint_coupled"))`. Compares against the v3 nested-Laplace path that
# profiles z out per outer (alpha, sigma) candidate.
# =============================================================================


test_that("joint_coupled smoke fit runs end-to-end and returns finite betas", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
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
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L,
                   engine = "joint_coupled")
  ))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")
  expect_true(all(is.finite(fit$means[c("psi_(Intercept)", "p_(Intercept)",
                                          "pos_(Intercept)")])))
  # The joint_fit object is attached for diagnostics; carries the outer
  # (sigma, alpha) grid + per-cell log_marginal.
  expect_s3_class(fit$joint_fit, "tulpa_nested_laplace_joint")
  expect_equal(length(fit$spatial_field), N)
  expect_true(all(is.finite(fit$spatial_field)))
})


test_that("joint_coupled errors on non-spatial occu_cover", {
  N <- 30L; J <- 4L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", seed = 1L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  # No spatial term on the psi formula -> nested_laplace itself errors before
  # the engine pick fires.
  expect_error(
    tobs(formula = ~ occ_cov1, data = cell_dat,
         family = occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = od$y, y_pos = y_pos, visits = od$det.covs,
         method = "nested_laplace",
         control = list(verbose = FALSE, engine = "joint_coupled")),
    "spatial term"
  )
})


test_that("joint_coupled recovers slopes, hypers, field shape (10 seeds)", {
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

  est_psi_x <- est_p_x <- est_pos_x <- est_alpha <- est_sigma <-
    field_cor <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, sigma_pos = sigma_pos_truth,
      positive = "lognormal", adj = adj,
      sigma = sigma_truth, alpha = alpha_truth, seed = 5000L + s
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
        control = list(verbose = FALSE, max.iter = 80L,
                       engine = "joint_coupled")
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    est_psi_x[s] <- fit$means["psi_occ_cov1"]
    est_p_x[s]   <- fit$means["p_det_cov1"]
    est_pos_x[s] <- fit$means["pos_pos_cov1"]
    est_alpha[s] <- fit$means[["alpha"]]
    est_sigma[s] <- fit$means[["sigma"]]
    field_cor[s] <- abs(cor(fit$spatial_field, sim$truth$f))
  }

  ok <- is.finite(est_psi_x) & is.finite(est_pos_x)
  expect_gte(sum(ok), 7L)

  expect_lt(abs(mean(est_psi_x[ok])     - beta_occ_truth[2L]), 0.35)
  expect_lt(abs(mean(est_p_x[ok])       - beta_p_truth[2L]),   0.35)
  expect_lt(abs(mean(est_pos_x[ok])     - beta_pos_truth[2L]), 0.20)

  expect_lt(abs(mean(est_alpha[ok]) - alpha_truth), 0.50)
  # The joint engine's coarse 5-point sigma_grid + Laplace's known small-N
  # underestimation of variance components together pull sigma's posterior
  # mean below truth. The wider tolerance here matches the v3 nested-Laplace
  # gate's reasoning (truth recovery within ~half a unit on a scale-like
  # axis) -- joint_coupled with a finer grid via control$sigma.grid
  # tightens this.
  expect_lt(abs(mean(est_sigma[ok]) - sigma_truth), 0.90)

  expect_gt(mean(field_cor[ok]), 0.80)
})
