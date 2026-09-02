# fitted() and the AIC / BIC penalty for occu_cover().
#
# Two defects this file pins:
#
#   * `fitted()` errored on EVERY occu_cover() fit (#291). The generic fallback
#     in `fitted.tobs_fit()` rebuilds eta from `model$X_processes`, which this
#     family never builds, so it died on `X_occ %*% beta_occ` with a NULL
#     design. `occu_cover` was the only family in the joint-cover group without
#     a handler; its siblings all had one.
# The compact (ragged) carrier is asserted in test-occu-cover-compact.R, which
# owns that fixture.
#
#   * `logLik()` reported `df = 0`, so AIC and BIC came back identical and
#     neither could rank models.

.fl_sim <- function(resp = "lognormal", seed = 9L, n = 6L, pos_field = FALSE) {
  adj <- rook_adj(n); N <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 5L, positive = resp, phi = 30,
    beta_occ = c(qlogis(0.7), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(if (identical(resp, "beta")) qlogis(0.3) else log(0.25), 0.2),
    sigma_pos = 0.3, adj = adj, sigma = 0.5, alpha = 1.0,
    pos_field = pos_field, seed = seed)
  if (is.null(sim$data$cell)) sim$data$cell <- seq_len(N)
  list(sim = sim, adj = adj, N = N)
}

.fl_check_shape <- function(fit, f, label) {
  v <- fitted(fit)
  expect_true(all(c("psi", "p", "cover", "z", "site_of_visit") %in% names(v)),
              info = label)
  expect_length(v$psi, f$N)
  expect_length(v$z, f$N)
  expect_true(all(is.finite(v$psi)) && all(v$psi > 0 & v$psi < 1), info = label)
  expect_true(all(is.finite(v$p))   && all(v$p   > 0 & v$p   < 1), info = label)
  expect_true(all(is.finite(v$cover)), info = label)
  # z is a probability, and is exactly 1 wherever the unit was ever detected.
  expect_true(all(v$z >= 0 & v$z <= 1), info = label)
  ever <- .occu_cover_visit_view(fit$model)$any_det
  expect_true(all(v$z[ever == 1L] == 1), info = label)
  expect_length(v$p, length(v$site_of_visit))
  expect_length(v$cover, length(v$site_of_visit))
  v
}

test_that("fitted() works on occu_cover, both engines and both cover families", {
  skip_if_fast(); skip_on_cran()
  for (resp in c("lognormal", "beta")) {
    f <- .fl_sim(resp)
    ns <- suppressWarnings(tobs(
      occurrence = ~ occ_cov1, detection = ~ 1, positive = ~ 1,
      family = occu_cover(response = resp),
      data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
      method = "laplace", control = list(progress = FALSE)))
    .fl_check_shape(ns, f, paste(resp, "laplace"))

    sp <- suppressWarnings(tobs(
      occurrence = ~ occ_cov1 + icar(graph = f$adj, group_var = "cell"),
      detection = ~ 1, positive = ~ 1 + share(spatial()),
      family = occu_cover(response = resp),
      data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
      method = "nested_laplace", control = list(progress = FALSE)))
    .fl_check_shape(sp, f, paste(resp, "nested_laplace"))
  }
})

test_that("fitted() carries the spatial field rather than scoring it at 0", {
  skip_if_fast(); skip_on_cran()
  # The cover arm has no covariate here, so its only source of between-site
  # variation is the copied field. A fit that scored the field at 0 would
  # return one distinct cover value.
  f <- .fl_sim("lognormal", seed = 4L)
  sp <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = f$adj, group_var = "cell"),
    detection = ~ 1, positive = ~ 1 + share(spatial()),
    family = occu_cover(response = "lognormal"),
    data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
    method = "nested_laplace", control = list(progress = FALSE)))
  expect_gt(length(unique(round(fitted(sp)$cover, 8))), 1L)
})

test_that("logLik() reports a real df, so AIC and BIC differ", {
  skip_if_fast(); skip_on_cran()
  f <- .fl_sim("lognormal", seed = 3L)
  ns <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1, detection = ~ 1, positive = ~ 1,
    family = occu_cover(response = "lognormal"),
    data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
    method = "laplace", control = list(progress = FALSE)))
  # psi (2) + p (1) + cover (1) + log dispersion (1).
  expect_equal(ns$n_fixed, 5L)
  expect_equal(attr(stats::logLik(ns), "df"), 5L)
  expect_false(isTRUE(all.equal(stats::AIC(ns), stats::BIC(ns))))

  sp <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = f$adj, group_var = "cell"),
    detection = ~ 1, positive = ~ 1 + share(spatial()),
    family = occu_cover(response = "lognormal"),
    data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
    method = "nested_laplace", control = list(progress = FALSE)))
  # The outer hyperparameters are integrated, not estimated, so they carry no
  # penalty: df counts the coefficient block alone, not `length(means)`.
  expect_lt(sp$n_fixed, length(sp$means))
  expect_equal(attr(stats::logLik(sp), "df"), sp$n_fixed)
  expect_false(isTRUE(all.equal(stats::AIC(sp), stats::BIC(sp))))
})
