# the joint nested_laplace occu_cover() path reported `sigma` as the raw
# amplitude against the unscaled intrinsic ICAR precision Q = D - W (tau = 1),
# while the #204 NUTS path reports `field_sd`, the geo-mean marginal SD
# (Sorbye-Rue) -- the only field-scale summary comparable across icar / bym2 /
# car_proper. The two differ by sqrt(scale_q), which is ~2.13 for the N = 30
# chain graph fixture below, so a fit compared against a simulator (or against
# the NUTS path) by reading `sigma` off the joint path was silently off by
# that factor. `field_sd` on the joint path should track `sigma_joint *
# sqrt(scale_q)`, and both backends' `field_sd` should land in the same
# convention on the same simulated field.

test_that("occu_cover() joint path exposes field_sd, not just the raw amplitude sigma", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L)
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
    formula = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L, engine = "joint")
  ))

  # Kept OUT of means/sds/vcov/draws: field_sd is a fixed deterministic
  # multiple of sigma (sqrt(scale_q), not fit-dependent), so folding it into
  # the sampled parameter vector makes the joint vcov singular (exact
  # collinearity) -- caught by chol(V) in test-occu-cover-joint.R on the first
  # cut of this fix. Reported instead on `fit$spatial`, alongside `sigma_mean`.
  expect_false("field_sd" %in% names(fit$means))
  expect_true(is.finite(fit$spatial$field_sd_mean))
  expect_true(is.finite(fit$spatial$field_sd_sd))

  scale_q <- tulpaObs:::.occu_cover_icar_scale(adj)
  expect_equal(fit$spatial$field_sd_mean,
               unname(fit$means[["sigma"]]) * sqrt(scale_q), tolerance = 1e-8)
  expect_equal(fit$spatial$field_sd_sd,
               unname(fit$sds[["sigma"]]) * sqrt(scale_q), tolerance = 1e-8)

  # field_sd, not raw sigma, is the number comparable to the simulator's own
  # sigma = 0.8 (geo-mean marginal SD convention, see ?simulate_occu_cover).
  expect_lt(abs(fit$spatial$field_sd_mean - 0.8), abs(fit$means[["sigma"]] - 0.8))
})

test_that("joint and NUTS occu_cover() backends report field_sd in the same convention", {
  skip_if_fast()
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  common <- list(
    formula = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs
  )
  fit_joint <- suppressWarnings(do.call(tobs, c(common, list(
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L, engine = "joint")
  ))))
  fit_nuts <- suppressWarnings(do.call(tobs, c(common, list(
    method = "nuts",
    control = list(verbose = FALSE, n.iter = 500L, n.warmup = 300L, seed = 1)
  ))))

  scale_q <- tulpaObs:::.occu_cover_icar_scale(adj)

  # The convention is a DETERMINISTIC relation on each backend -- field_sd is
  # sigma rescaled by a graph constant, not an estimate -- so it is asserted as
  # one, on every draw, the same way the joint path is asserted against
  # sigma_mean above. Comparing two point estimates instead cannot separate a
  # unit mismatch from ordinary estimator disagreement, and it is not what the
  # claim says.
  expect_equal(unname(fit_nuts$hyper_draws[, "field_sd"]),
               unname(fit_nuts$hyper_draws[, "sigma"]) * sqrt(scale_q),
               tolerance = 1e-12)

  fsd_joint <- fit_joint$spatial$field_sd_mean
  # Median, not mean. This posterior is strongly right-skewed -- on the
  # reference fixture the 95th percentile of field_sd sits an order of
  # magnitude above the median -- so its mean is carried by the tail and moves
  # with whatever draws a platform's arithmetic happens to produce. The mean
  # is what made this block fail on aarch64 and pass on x86_64 (issue #314);
  # the median is the summary the package already uses for variance-scale
  # posteriors elsewhere, for the same reason.
  fsd_nuts <- stats::median(fit_nuts$hyper_draws[, "field_sd"])

  # Two different estimators on 30 cells, so this is a gross-divergence guard
  # only. It is deliberately looser than the log(2) it replaces: the unit
  # question it used to stand in for is now carried by the exact assertion
  # above, which is the stricter instrument, so the slack here costs no signal.
  expect_lt(abs(log(fsd_joint / fsd_nuts)), log(4))
})
