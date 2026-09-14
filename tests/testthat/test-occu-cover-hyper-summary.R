# =============================================================================
# test-occu-cover-hyper-summary.R
# - `fit$hyper_axis` / `fit$hyper_summary` on an occu_cover joint fit.
#
# `means` / `sds` report a hyperparameter's grid-weighted moments, which read
# an SD of 0 wherever the outer weight sits on one cell. `hyper_summary` sets
# the engine's own per-axis read beside the mean: `theta_sd` (with the
# estimator `theta_sd_source` names) and the weighted quantiles
# `theta_ci_lo` / `theta_ci_hi`, taken off the axis `hyper_axis` records.
# =============================================================================

.hs_fixture <- function() {
  N <- 24L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 202L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim1$y)),
    det_cov1 = sim1$visit_data$det_cov1, pos_cov1 = sim1$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  yp1 <- sim1$y_pos; yp1[is.na(yp1)] <- 0
  yp2 <- sim2$y_pos; yp2[is.na(yp2)] <- 0
  list(adj = adj, od = od,
       cell_dat = cbind(data.frame(site_id = seq_len(N)), sim1$data),
       y = list(a = od$y, b = sim2$y), y_pos = list(yp1, yp2))
}

.hs_fit <- function(fx, y, y_pos, control) {
  adj <- fx$adj
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = fx$cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = y, y_pos = y_pos, visits = fx$od$det.covs,
    method = "nested_laplace", control = control))
}

.hs_ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint",
                 diagnose.k = FALSE, progress = FALSE)

# Each row against the engine read of its own axis, through the transform the
# row states.
.hs_expect_engine_read <- function(fit) {
  hs <- fit$hyper_summary
  jf <- fit$joint_fit
  n_hyper <- length(fit$means) - fit$n_fixed
  hyper_names <- utils::tail(names(fit$means), n_hyper)
  expect_s3_class(hs, "data.frame")
  expect_identical(names(hs), c("parameter", "mean", "sd", "lwr", "upr",
                                "sd_source", "axis"))
  expect_identical(hs$parameter, hyper_names)
  expect_identical(names(fit$hyper_axis), hyper_names)
  expect_identical(unname(fit$hyper_axis), hs$axis)
  expect_equal(hs$mean, unname(fit$means[hyper_names]), tolerance = 0)
  expect_true(all(is.na(hs$axis) | hs$axis %in% colnames(jf$theta_grid)))
  for (i in seq_len(nrow(hs))) {
    ax <- hs$axis[[i]]
    if (is.na(ax)) {
      expect_true(all(is.na(unlist(hs[i, c("sd", "lwr", "upr", "sd_source")]))))
      next
    }
    expect_identical(hs$sd_source[[i]], unname(jf$theta_sd_source[[ax]]))
    if (identical(hs$parameter[[i]], "phi_pos")) {
      # The lognormal cover arm's engine phi is the residual variance.
      expect_equal(hs$lwr[[i]], sqrt(jf$theta_ci_lo[[ax]]), tolerance = 1e-12)
      expect_equal(hs$upr[[i]], sqrt(jf$theta_ci_hi[[ax]]), tolerance = 1e-12)
      expect_equal(hs$sd[[i]],
                   jf$theta_sd[[ax]] / (2 * sqrt(jf$theta_mean[[ax]])),
                   tolerance = 1e-12)
    } else {
      expect_equal(hs$sd[[i]],  jf$theta_sd[[ax]],    tolerance = 0)
      expect_equal(hs$lwr[[i]], jf$theta_ci_lo[[ax]], tolerance = 0)
      expect_equal(hs$upr[[i]], jf$theta_ci_hi[[ax]], tolerance = 0)
    }
  }
}

test_that("hyper_summary reads each hyperparameter's own engine axis", {
  skip_on_cran()
  skip_if_fast()
  fx <- .hs_fixture()
  fit <- .hs_fit(fx, fx$y$a, fx$y_pos[[1L]], .hs_ctrl)
  expect_true(all(c("sigma", "phi_pos") %in% fit$hyper_summary$parameter))
  expect_identical(fit$hyper_axis[c("sigma", "phi_pos")],
                   c(sigma = "sigma", phi_pos = "phi_pos"))
  .hs_expect_engine_read(fit)
  expect_true(all(is.finite(fit$hyper_summary$sd)))
  expect_true(all(fit$hyper_summary$lwr <= fit$hyper_summary$upr))
})

test_that("a fused batch species carries the hyper_summary of its own fit", {
  skip_on_cran()
  skip_if_fast()
  fx <- .hs_fixture()
  set.seed(41L)
  batch <- .hs_fit(fx, fx$y, fx$y_pos, c(.hs_ctrl, list(batch.backend = "fused")))
  expect_identical(batch$backend, "fused")
  op <- options(tulpaObs.compress_nodet = FALSE)
  on.exit(options(op), add = TRUE)
  set.seed(41L)
  ind <- list(a = .hs_fit(fx, fx$y$a, fx$y_pos[[1L]], .hs_ctrl),
              b = .hs_fit(fx, fx$y$b, fx$y_pos[[2L]], .hs_ctrl))
  for (sp in c("a", "b")) {
    .hs_expect_engine_read(batch$fits[[sp]])
    expect_identical(batch$fits[[sp]]$hyper_axis, ind[[sp]]$hyper_axis)
    expect_equal(batch$fits[[sp]]$hyper_summary, ind[[sp]]$hyper_summary,
                 tolerance = 1e-9, info = sp)
  }
})

test_that("hyper_summary divides a rescaled row and leaves a derived row empty", {
  fit <- list(
    theta_grid      = matrix(0, 1L, 3L,
                             dimnames = list(NULL, c("b1.sigma", "b2.sigma", "b3.tau"))),
    theta_mean      = c(b1.sigma = 1.1, b2.sigma = 0.6, b3.tau = 4),
    theta_sd        = c(b1.sigma = 0.2, b2.sigma = 0.1, b3.tau = 1),
    theta_ci_lo     = c(b1.sigma = 0.8, b2.sigma = 0.4, b3.tau = 2),
    theta_ci_hi     = c(b1.sigma = 1.5, b2.sigma = 0.9, b3.tau = 7),
    theta_sd_source = c(b1.sigma = "stencil", b2.sigma = "moment", b3.tau = "stencil"))
  hs <- tulpaObs:::.occu_cover_hyper_summary(
    fit, c("sigma", "sigma_re_x", "sigma_pos_field"),
    c(sigma = 1.1, sigma_re_x = 0.3, sigma_pos_field = 0.5),
    c(sigma = "b1.sigma", sigma_re_x = "b2.sigma", sigma_pos_field = NA),
    scale = list(sigma_re_x = 2), family = list())
  expect_identical(hs$axis, c("b1.sigma", "b2.sigma", NA))
  expect_equal(hs$sd,  c(0.2, 0.05, NA))
  expect_equal(hs$lwr, c(0.8, 0.2,  NA))
  expect_equal(hs$upr, c(1.5, 0.45, NA))
  expect_identical(hs$sd_source, c("stencil", "moment", NA))
  expect_equal(hs$mean, c(1.1, 0.3, 0.5))
})
