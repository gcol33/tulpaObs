# =============================================================================
# test-occu-cover-trend.R - the joint_coupled engine with a SECOND coupled
# field: a spatially-varying temporal trend weighted by a per-cell covariate
# (gcol33/tulpaObs#15). The intercept field PLUS the trend field both couple
# onto the cover arm, each with its own scale (alpha, alpha_trend), via the
# multi-block nested-Laplace copy path.
#
# Requested with control = list(trend = list(weight = "time")), naming a
# per-cell numeric covariate in the cell data.
# =============================================================================


.trend_chain_adj <- function(N) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  adj
}

.trend_fit <- function(sim, N, J, max.iter = 100L) {
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = sim$adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = max.iter,
                   engine = "joint_coupled",
                   trend = list(weight = "time"))
  ))
}


test_that("occu_cover trend smoke fit runs end-to-end and exposes both fields", {
  N <- 30L; J <- 4L
  adj <- .trend_chain_adj(N)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal", adj = adj,
    sigma = 0.8, alpha = 1.0, trend = TRUE,
    sigma_trend = 0.7, alpha_trend = 0.9, seed = 31337L
  )
  expect_true("time" %in% names(sim$data))
  expect_length(sim$truth$f2, N)

  fit <- .trend_fit(sim, N, J, max.iter = 300L)
  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")

  # Both coupling scales surface as named hyperparameters.
  expect_true(all(c("sigma", "alpha", "sigma_trend", "alpha_trend") %in%
                    names(fit$means)))
  expect_true(all(is.finite(fit$means[c("psi_(Intercept)", "p_(Intercept)",
                                          "pos_(Intercept)", "alpha",
                                          "alpha_trend")])))

  # Both fields are returned, demeaned, length N.
  expect_length(fit$spatial_field, N)
  expect_length(fit$trend_field, N)
  expect_true(all(is.finite(fit$spatial_field)))
  expect_true(all(is.finite(fit$trend_field)))
  expect_lt(abs(mean(fit$spatial_field)), 1e-6)
  expect_lt(abs(mean(fit$trend_field)),   1e-6)

  # Issue-16 compat: the joint betas+field vcov covers BOTH fields (one
  # sum-to-zero constraint row per ICAR block), so it is (p_beta + 2N) square,
  # symmetric, and finite.
  p_beta <- length(fit$process_info[[1L]]$coef_names) +
            length(fit$process_info[[2L]]$coef_names) +
            length(fit$process_info[[3L]]$coef_names)
  Vj <- fit$joint_vcov
  expect_equal(dim(Vj), c(p_beta + 2L * N, p_beta + 2L * N))
  expect_true(isSymmetric(Vj))
  expect_true(all(is.finite(Vj)))
})


test_that("occu_cover trend recovers slopes, both couplings, both fields (10 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 10L
  N <- 100L; J <- 6L
  adj <- .trend_chain_adj(N)

  beta_occ_truth <- c(stats::qlogis(0.4), 0.7)
  beta_p_truth   <- c(0.0,                0.8)
  beta_pos_truth <- c(log(0.20),         -0.4)
  sigma_truth        <- 1.0
  alpha_truth        <- 1.0
  sigma_trend_truth  <- 0.9
  alpha_trend_truth  <- 0.9

  est_psi_x <- est_p_x <- est_pos_x <- est_alpha <- est_alpha_tr <-
    cor1 <- cor2 <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, sigma_pos = 0.35, positive = "lognormal",
      adj = adj, sigma = sigma_truth, alpha = alpha_truth, trend = TRUE,
      sigma_trend = sigma_trend_truth, alpha_trend = alpha_trend_truth,
      seed = 6000L + s
    )
    fit <- tryCatch(.trend_fit(sim, N, J, max.iter = 100L),
                    error = function(e) NULL)
    if (is.null(fit)) next

    est_psi_x[s]    <- fit$means["psi_occ_cov1"]
    est_p_x[s]      <- fit$means["p_det_cov1"]
    est_pos_x[s]    <- fit$means["pos_pos_cov1"]
    est_alpha[s]    <- fit$means[["alpha"]]
    est_alpha_tr[s] <- fit$means[["alpha_trend"]]
    cor1[s]         <- abs(cor(fit$spatial_field, sim$truth$f))
    cor2[s]         <- abs(cor(fit$trend_field,   sim$truth$f2))
  }

  ok <- is.finite(est_psi_x) & is.finite(est_pos_x)
  expect_gte(sum(ok), 7L)

  # Fixed-effect slopes recover within the issue's tolerance.
  expect_lt(abs(mean(est_psi_x[ok]) - beta_occ_truth[2L]), 0.35)
  expect_lt(abs(mean(est_p_x[ok])   - beta_p_truth[2L]),   0.35)
  expect_lt(abs(mean(est_pos_x[ok]) - beta_pos_truth[2L]), 0.35)

  # Both coupling scales recover within ~0.5 of truth.
  expect_lt(abs(mean(est_alpha[ok])    - alpha_truth),       0.50)
  expect_lt(abs(mean(est_alpha_tr[ok]) - alpha_trend_truth), 0.50)

  # Both field shapes recover (mean |cor| > 0.7). Gating the mean rather than
  # every seed tolerates the occasional draw where the intercept field is
  # weakly identified against the trend field on the same graph.
  expect_gt(mean(cor1[ok]), 0.70)
  expect_gt(mean(cor2[ok]), 0.70)
})
