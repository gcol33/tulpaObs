# =============================================================================
# test-occu-cover.R - joint occupancy-detection + cover hurdle (occu_cover())
#
# v1 scope: non-spatial Laplace via exact two-state marginal + optim. These
# tests cover:
#   - family constructor + tobs() smoke
#   - input shape validation
#   - point recovery and 95% Wald CI coverage. The family is now
#     `status = "working"` (gcol33/tulpaObs#96), so the gate is the working-family
#     standard per the recovery rubric in tulpaObs CLAUDE.md "Statistical Code
#     Needs Recovery Tests": pooled coverage (every coefficient x seed cell) at
#     the 0.85 finite-sample floor for a nominal-95% interval, plus a 0.80
#     per-coordinate floor for Monte-Carlo slack. Measured pooled coverage across
#     the non-spatial laplace / nuts and shared-field nested-Laplace paths, beta
#     and lognormal arms, is 0.92-0.96.
# =============================================================================


test_that("occu_cover() constructor returns a tobs_family", {
  f <- occu_cover()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "occu_cover")
  expect_equal(f$params$positive, "beta")
  expect_equal(occu_cover("lognormal")$params$positive, "lognormal")
})


test_that("occu_cover() rejects structured terms in v1", {
  sim <- simulate_occu_cover(N = 40L, J = 3L, positive = "lognormal", seed = 1L)
  od <- tobs_data(
    df = data.frame(
      site_id = rep(seq_len(nrow(sim$y)), each = ncol(sim$y)),
      visit   = rep(seq_len(ncol(sim$y)), times = nrow(sim$y)),
      y       = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1
    ),
    y = "y", site = "site_id", visit = "visit", det.covs = "det_cov1"
  )
  cell_dat <- data.frame(site_id = seq_len(nrow(sim$y)))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  # bym2() on the psi formula under method = "laplace" is rejected by the
  # dispatcher with a pointer to the v2 spatial path. (Isolated-node graph
  # would error if it ever reached the spatial fitter, but the dispatcher
  # rejects first.)
  adj_iso <- matrix(0L, 40, 40); for (i in seq_len(39)) adj_iso[i, i + 1L] <- adj_iso[i + 1L, i] <- 1L
  expect_error(
    suppressWarnings(tobs(
      formula   = ~ 1 + bym2(graph = adj_iso),
      data      = cell_dat,
      family    = occu_cover("lognormal"),
      detection = ~ det_cov1,
      positive  = ~ det_cov1,
      y         = od$y, y_pos = y_pos, visits = od$det.covs,
      method    = "laplace", control = list(verbose = FALSE)
    )),
    "non-spatial"
  )

})


test_that("occu_cover() nobs() counts the valid visit rows", {
  skip_on_cran()
  skip_if_fast()

  N <- 60L; J <- 4L
  sim <- simulate_occu_cover(N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L,
                             n_pos_covs = 1L, positive = "lognormal",
                             seed = 4242L)
  long <- data.frame(
    site_id  = rep(seq_len(N), each = J),
    visit    = rep(seq_len(J), times = N),
    y        = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1,
    pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- tobs(formula = ~ occ_cov1, data = cell_dat,
              family = occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = od$y, y_pos = y_pos,
              visits = od$det.covs, method = "laplace",
              control = list(verbose = FALSE, max.iter = 200L))

  # One observation per surveyed visit, read through the shared visit view so
  # the dense and compact layouts of the same data report the same count.
  expect_identical(nobs(fit), sum(!is.na(sim$y)))
  expect_gt(nobs(fit), 0L)
})


test_that("occu_cover() recovers parameters (lognormal positive, 20 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 20L
  N <- 300L; J <- 5L
  beta_occ_truth <- c(stats::qlogis(0.4), 0.9)
  beta_p_truth   <- c(0.0, 0.6)
  beta_pos_truth <- c(log(0.10), -0.4)
  sigma_pos_truth <- 0.4

  est <- list(
    psi_int = numeric(n_seeds), psi_x = numeric(n_seeds),
    p_int   = numeric(n_seeds), p_x   = numeric(n_seeds),
    pos_int = numeric(n_seeds), pos_x = numeric(n_seeds),
    sigma   = numeric(n_seeds)
  )
  se <- est
  conv <- logical(n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N         = N, J = J,
      n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
      beta_occ  = beta_occ_truth,
      beta_p    = beta_p_truth,
      beta_pos  = beta_pos_truth,
      sigma_pos = sigma_pos_truth,
      positive  = "lognormal",
      seed      = 9000L + s
    )

    long <- data.frame(
      site_id = rep(seq_len(N), each = J),
      visit   = rep(seq_len(J), times = N),
      y       = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1,
      pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

    fit <- tryCatch(
      tobs(formula   = ~ occ_cov1, data = cell_dat,
           family    = occu_cover("lognormal"),
           detection = ~ det_cov1,
           positive  = ~ pos_cov1,
           y         = od$y, y_pos = y_pos, visits = od$det.covs,
           method    = "laplace",
           control   = list(verbose = FALSE, max.iter = 500L)),
      error = function(e) NULL
    )
    if (is.null(fit)) { conv[s] <- FALSE; next }
    conv[s] <- isTRUE(fit$convergence$converged)

    est$psi_int[s] <- fit$means["psi_(Intercept)"]
    est$psi_x[s]   <- fit$means["psi_occ_cov1"]
    est$p_int[s]   <- fit$means["p_(Intercept)"]
    est$p_x[s]     <- fit$means["p_det_cov1"]
    est$pos_int[s] <- fit$means["pos_(Intercept)"]
    est$pos_x[s]   <- fit$means["pos_pos_cov1"]
    est$sigma[s]   <- exp(fit$means["log_sigma_pos"])

    se$psi_int[s] <- fit$sds["psi_(Intercept)"]
    se$psi_x[s]   <- fit$sds["psi_occ_cov1"]
    se$p_int[s]   <- fit$sds["p_(Intercept)"]
    se$p_x[s]     <- fit$sds["p_det_cov1"]
    se$pos_int[s] <- fit$sds["pos_(Intercept)"]
    se$pos_x[s]   <- fit$sds["pos_pos_cov1"]
    # delta-method SE on sigma = exp(log_sigma)
    se$sigma[s]   <- exp(fit$means["log_sigma_pos"]) * fit$sds["log_sigma_pos"]
  }

  expect_gte(mean(conv), 0.80)

  # Point recovery: bias |mean(est) - truth| < tol
  expect_lt(abs(mean(est$psi_int[conv]) - beta_occ_truth[1L]), 0.20)
  expect_lt(abs(mean(est$psi_x[conv])   - beta_occ_truth[2L]), 0.20)
  expect_lt(abs(mean(est$p_int[conv])   - beta_p_truth[1L]),   0.20)
  expect_lt(abs(mean(est$p_x[conv])     - beta_p_truth[2L]),   0.20)
  expect_lt(abs(mean(est$pos_int[conv]) - beta_pos_truth[1L]), 0.20)
  expect_lt(abs(mean(est$pos_x[conv])   - beta_pos_truth[2L]), 0.20)
  expect_lt(abs(mean(est$sigma[conv])   - sigma_pos_truth),    0.10)

  # 95% Wald CI coverage (working-family gate, gcol33/tulpaObs#96): pooled over
  # every coefficient x seed cell at the 0.85 floor, with a 0.80 per-coordinate
  # floor for Monte-Carlo slack. Measured pooled coverage ~0.95 at 30 seeds.
  cov_cells <- list(
    psi_int = abs(est$psi_int[conv] - beta_occ_truth[1L]) < 1.96 * se$psi_int[conv],
    psi_x   = abs(est$psi_x[conv]   - beta_occ_truth[2L]) < 1.96 * se$psi_x[conv],
    p_int   = abs(est$p_int[conv]   - beta_p_truth[1L])   < 1.96 * se$p_int[conv],
    p_x     = abs(est$p_x[conv]     - beta_p_truth[2L])   < 1.96 * se$p_x[conv],
    pos_int = abs(est$pos_int[conv] - beta_pos_truth[1L]) < 1.96 * se$pos_int[conv],
    pos_x   = abs(est$pos_x[conv]   - beta_pos_truth[2L]) < 1.96 * se$pos_x[conv]
  )
  per_coord <- vapply(cov_cells, mean, numeric(1))
  expect_gte(mean(unlist(cov_cells)), 0.85)   # pooled
  expect_gte(min(per_coord),          0.80)   # no coordinate collapses
})


test_that("occu_cover() recovers parameters (beta positive, 20 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 20L
  N <- 250L; J <- 5L
  beta_occ_truth <- c(stats::qlogis(0.5), 0.7)
  beta_p_truth   <- c(0.0, 0.5)
  beta_pos_truth <- c(stats::qlogis(0.3), -0.3)
  phi_truth      <- 30

  nm <- c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)", "p_det_cov1",
          "pos_(Intercept)", "pos_pos_cov1")
  truth <- c(beta_occ_truth, beta_p_truth, beta_pos_truth)
  est <- se <- matrix(NA_real_, n_seeds, 6L, dimnames = list(NULL, nm))
  phi_est <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, positive = "beta", phi = phi_truth,
      seed = 7000L + s
    )
    long <- data.frame(
      site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
      y = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1,
      pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
    # clip into (0, 1) for the beta arm.
    y_pos <- pmin(pmax(y_pos, 1e-6), 1 - 1e-6)

    fit <- tryCatch(
      tobs(formula = ~ occ_cov1, data = cell_dat,
           family = occu_cover("beta"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y = od$y, y_pos = y_pos, visits = od$det.covs,
           method = "laplace",
           control = list(verbose = FALSE, max.iter = 500L)),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    est[s, ] <- fit$means[nm]; se[s, ] <- fit$sds[nm]
    phi_est[s] <- exp(fit$means["log_phi"])
  }

  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 16L)

  # Point recovery of every coefficient (mean over seeds).
  bias <- colMeans(est[ok, , drop = FALSE]) - truth
  expect_true(all(abs(bias) < 0.30))
  expect_gt(mean(phi_est, na.rm = TRUE), 10)   # well above the 1 / 0 boundary

  # 95% Wald CI coverage (working-family gate, gcol33/tulpaObs#96): pooled at the
  # 0.85 floor, 0.80 per-coordinate floor. Measured pooled ~0.95 at 30 seeds.
  cov_cells <- abs(est[ok, , drop = FALSE] -
                   matrix(truth, sum(ok), 6L, byrow = TRUE)) <
               1.96 * se[ok, , drop = FALSE]
  expect_gte(mean(cov_cells),            0.85)   # pooled
  expect_gte(min(colMeans(cov_cells)),   0.80)   # no coordinate collapses
})
