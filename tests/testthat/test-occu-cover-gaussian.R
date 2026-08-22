# Tests for occu_cover(response = "gaussian") -- the identity-Gaussian cover arm
# on the joint occupancy-detection + cover hurdle. The non-spatial Laplace path
# evaluates the positive-arm density entirely in R (.occu_cover_pos_logdens), so
# this is the cheapest recovery route. mu = eta on the response scale (no log
# transform); the cover magnitude may be any real at detected visits.

test_that("occu_cover(response = 'gaussian') constructor is wired through", {
  f <- occu_cover("gaussian")
  expect_s3_class(f, "tobs_family")
  expect_equal(f$params$positive, "gaussian")
  expect_equal(f$observation, "detection_plus_gaussian")
})

test_that("occu_cover(gaussian): WAIC works, PPC gated, NUTS available", {
  N <- 120L; J <- 5L
  sim <- simulate_occu_cover(N = N, J = J, positive = "gaussian",
    beta_pos = c(2.0, -0.4), sigma_pos = 0.5, seed = 11L)
  long <- data.frame(site_id = rep(1:N, each = J), visit = rep(1:J, times = N),
    y = as.vector(t(sim$y)), det_cov1 = sim$visit_data$det_cov1,
    pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
    det.covs = c("det_cov1", "pos_cov1"))
  cell <- cbind(data.frame(site_id = 1:N), sim$data)
  yp <- sim$y_pos; yp[is.na(yp)] <- 0
  fit <- tobs(formula = ~ occ_cov1, data = cell, family = occu_cover("gaussian"),
    detection = ~ det_cov1, positive = ~ pos_cov1, y = od$y, y_pos = yp,
    visits = od$det.covs, method = "laplace", control = list(verbose = FALSE))
  w <- waic(fit)
  expect_true(is.finite(w$waic) && is.finite(w$p_waic))
  expect_error(ppc(fit), "not defined for occu_cover.*gaussian")
  # NUTS is wired for the gaussian arm; a short sample returns a fit rather
  # than routing gaussian through the lognormal dispatch.
  nut <- tobs(formula = ~ occ_cov1, data = cell, family = occu_cover("gaussian"),
    detection = ~ det_cov1, positive = ~ pos_cov1, y = od$y, y_pos = yp,
    visits = od$det.covs, method = "nuts",
    control = list(n.iter = 400L, n.warmup = 300L, verbose = FALSE))
  expect_equal(nut$method, "nuts")
})


test_that("occu_cover() recovers parameters (gaussian positive, 20 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 20L
  N <- 300L; J <- 5L
  beta_occ_truth <- c(stats::qlogis(0.4), 0.9)
  beta_p_truth   <- c(0.0, 0.6)
  beta_pos_truth <- c(2.0, -0.4)
  sigma_pos_truth <- 0.5

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
      positive  = "gaussian",
      seed      = 9200L + s
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
           family    = occu_cover("gaussian"),
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
    se$sigma[s]   <- exp(fit$means["log_sigma_pos"]) * fit$sds["log_sigma_pos"]
  }

  expect_gte(mean(conv), 0.80)

  expect_lt(abs(mean(est$psi_int[conv]) - beta_occ_truth[1L]), 0.20)
  expect_lt(abs(mean(est$psi_x[conv])   - beta_occ_truth[2L]), 0.20)
  expect_lt(abs(mean(est$p_int[conv])   - beta_p_truth[1L]),   0.20)
  expect_lt(abs(mean(est$p_x[conv])     - beta_p_truth[2L]),   0.20)
  expect_lt(abs(mean(est$pos_int[conv]) - beta_pos_truth[1L]), 0.15)
  expect_lt(abs(mean(est$pos_x[conv])   - beta_pos_truth[2L]), 0.15)
  expect_lt(abs(mean(est$sigma[conv])   - sigma_pos_truth),    0.10)

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
