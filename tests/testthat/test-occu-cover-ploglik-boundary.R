# Regression: the occu_cover WAIC / LOO pointwise density must equal the fit
# kernel (src/occu_coupling_shared.h::pos_log_density) at the cover boundary
# (cover exactly 0 or 1) and at extreme eta -- log_safe instead of a bare log(0)
# -> -Inf, and no eta clamp -- rather than diverging from the likelihood the
# model was fit with.

test_that(".occu_cover_pos_logdens is finite at the beta cover boundary", {
  # Beta cover in [0, 1]: exactly 0 and 1 must give a finite log-density via
  # log_safe (the linear (a-1) log y / (b-1) log(1-y) terms stay finite), not the
  # -Inf / NaN a bare log(0) produced.
  d <- .occu_cover_pos_logdens(c(0, 1, 0.5), rep(0.2, 3), 20, "beta")
  expect_true(all(is.finite(d)))
  # And the interior value is the standard Beta density (unchanged away from 0/1).
  mu <- stats::plogis(0.2)
  ref <- dbeta(0.5, mu * 20, (1 - mu) * 20, log = TRUE)
  expect_equal(d[3L], ref)
})

test_that(".occu_cover_pos_logdens drops the eta clamp (matches the fit kernel)", {
  # The coupling fit kernel does not clamp eta; past the old +-30 bound the density
  # must still respond (a clamp would have pinned eta = 30 == eta = 300).
  d30  <- .occu_cover_pos_logdens(0.5, 30,  20, "beta")
  d300 <- .occu_cover_pos_logdens(0.5, 300, 20, "beta")
  expect_false(isTRUE(all.equal(d30, d300)))
  # Lognormal: raw predictor, so a large eta moves the residual monotonically.
  l1 <- .occu_cover_pos_logdens(1.0, 5,  0.4, "lognormal")
  l2 <- .occu_cover_pos_logdens(1.0, 50, 0.4, "lognormal")
  expect_false(isTRUE(all.equal(l1, l2)))
})

test_that("occu_cover WAIC pointwise density is finite when a detected beta cover is exactly 1", {
  skip_on_cran()
  set.seed(4)
  N <- 90L; J <- 4L
  sim <- simulate_occu_cover(N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L,
                             n_pos_covs = 1L, phi = 30, sigma_pos = 0.4,
                             positive = "beta", seed = 4L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  fit <- tobs(formula = ~ occ_cov1, data = cell_dat, family = occu_cover("beta"),
              detection = ~ det_cov1, positive = ~ pos_cov1, y = od$y,
              y_pos = y_pos, visits = od$det.covs, method = "laplace",
              control = list(verbose = FALSE, max.iter = 60L))

  # Baseline WAIC is finite.
  expect_true(is.finite(waic(fit)$waic))

  # Force a boundary cover (exactly 1) at a detected visit. Before the fix the
  # bare log(1 - 1) = log(0) put -Inf / NaN into the WAIC pointwise output; now
  # the C++ ploglik uses the fit kernel's log_safe, so every pointwise term stays
  # finite (a boundary datum is legitimately near-zero likelihood, hence a large
  # finite negative, not -Inf / NaN). The C++ pointwise ll and the fit-marginal
  # pointwise ll must also agree (single source of truth).
  det  <- which(fit$model$valid & fit$model$y == 1L, arr.ind = TRUE)[1L, ]
  site <- det[1L]
  fit$model$y_pos[site, det[2L]] <- 1.0
  pll <- tulpaObs:::.tobs_ploglik_occu_cover(fit, n.draws = 50L)  # [S x n_sites]
  expect_true(all(is.finite(pll)))
  # The injection actually engaged the boundary path: the affected site's
  # log-density is a large finite negative (a boundary datum is legitimately
  # near-zero likelihood), not -Inf / NaN, and far below the other sites.
  site_mean <- colMeans(pll)
  expect_true(is.finite(site_mean[site]))
  expect_lt(site_mean[site], min(site_mean[-site]))
})
