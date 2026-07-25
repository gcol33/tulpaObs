# The continuous NNGP gp() surface: exposure + recovery (gcol33/tulpaObs#152).
#
# gp() was smoke-tested only (test-spatial-occ.R just checks the fit runs and
# `n_params > 2`). Once the NNGP neighbour-pair distance packing bug and the
# `gp_parameterization` draw-storage corruption were fixed (commit 318b682),
# `fit$spatial_field` and `fit$gp_layout` became readable for the first time,
# which is what makes an actual recovery test writable here.
#
# What that test then found: the surface DOES track a known truth's SHAPE, but
# its AMPLITUDE is severely attenuated -- sd(field_hat) sits at a few percent
# of sd(field_truth), reproducibly, even under near-perfect detection data and
# 5x longer warmup (ruling out both weak occupancy identifiability and
# insufficient adaptation as the explanation). Filed upstream as
# gcol33/tulpa#243: `compute_gp_spatial_prior` samples the field CENTERED (the
# NNGP density evaluated directly on the raw field), a textbook Neal's funnel
# for a hierarchical scale + high-dimensional field jointly sampled by NUTS.
# The engine already has a non-centered transform (`nngp_nc_forward` /
# `nngp_nc_backward` in hmc_gp_nc.h) but it is wired only into the post-hoc
# draw-storage path, never into the actual sampling gradient -- so every fit
# runs centered regardless of `gp_parameterization`. This mirrors what
# gcol33/tulpa#144 left open for svc() ("sigma_svc reads high... may be the
# funnel... untested either way") -- svc() shares the identical centered
# joint-hyperparameter-and-field architecture.
#
# Until gcol33/tulpa#243 lands, this test scores what the fit can actually
# deliver: the surface's SHAPE (correlation with truth) and that phi (range)
# recovers -- not the amplitude, which is the tracked, upstream-blocked gap.
# This is the same split test-occu-svc-nngp-recovery.R makes for svc(): a
# recovery test on cor(), a separate calibration test that does NOT assert on
# the surface's scale.

.gpr_sim <- function(N, J, seed, sigma_f = 1.2, phi_f = 0.25, p_det = 0.6,
                     b0 = 0.2) {
  set.seed(seed)
  lon <- stats::runif(N); lat <- stats::runif(N)
  D <- as.matrix(stats::dist(cbind(lon, lat)))
  K <- sigma_f^2 * exp(-D / phi_f)
  w <- as.numeric(t(chol(K + 1e-8 * diag(N))) %*% stats::rnorm(N))
  z <- stats::rbinom(N, 1, stats::plogis(b0 + w))
  y <- matrix(stats::rbinom(N * J, 1, p_det * rep(z, J)), N, J)
  list(data = data.frame(lon = lon, lat = lat), y = y, w = w,
       sigma_f = sigma_f, phi_f = phi_f)
}

.gpr_fit <- function(sim, seed, n_iter = 600L, n_warmup = 300L, nn = 10) {
  suppressWarnings(tobs(
    ~ gp(lon, lat, nn = nn, prior_range = c(0.1, 0.05)),
    data = sim$data,
    family = occu(), detection = ~ 1, y = sim$y,
    method = "nuts",
    control = list(n.iter = n_iter, n.warmup = n_warmup, seed = seed,
                   verbose = FALSE)))
}

test_that("the gp() field is exposed, named, and correctly shaped", {
  skip_on_cran()
  skip_if_fast()

  sim <- .gpr_sim(N = 40L, J = 4L, seed = 1L)
  fit <- .gpr_fit(sim, seed = 1L, n_iter = 120L, n_warmup = 60L)

  expect_false(is.null(fit$gp_layout))
  expect_identical(as.integer(fit$gp_layout$n_units), 40L)
  expect_false(isTRUE(fit$gp_layout$collapsed))

  expect_length(fit$spatial_field, 40L)
  expect_true(all(is.finite(fit$spatial_field)))

  cn <- colnames(fit$draws)
  expect_equal(sum(grepl("^gp_w\\[", cn)), 40L)
  expect_true("log_sigma2_gp" %in% cn)
  expect_true("log_phi_gp" %in% cn)
  expect_equal(sum(grepl("^param\\[", cn)), 0L)

  # fitted() reads the same offset the sampler put in the psi logit, not a
  # dropped/flat field (gcol33/tulpaObs#152's third reported bug).
  psi_hat <- fitted(fit)$psi
  expect_gt(length(unique(round(psi_hat, 8))), 1L)
})

test_that("occu() + gp() recovers the shape of a known GP surface", {
  skip_on_cran()
  skip_if_fast()

  seeds <- 1:4
  cr <- rep(NA_real_, length(seeds))
  for (i in seq_along(seeds)) {
    sim <- .gpr_sim(N = 40L, J = 6L, seed = seeds[i])
    fit <- try(.gpr_fit(sim, seed = seeds[i]), silent = TRUE)
    if (inherits(fit, "try-error")) next
    cr[i] <- stats::cor(fit$spatial_field, sim$w)
  }
  ok <- !is.na(cr)
  expect_gt(sum(ok), 2L)

  # Measured 0.046 / 0.720 / 0.736 / 0.699 (mean 0.55) at N = 40, J = 6, p =
  # 0.6 -- the same scale as gcol33/tulpaObs#152's own repro. Seed 1 is a weak
  # outlier; the floor stays well clear of it rather than being tuned to
  # exclude it. This asserts the surface carries real spatial signal; it does
  # NOT assert the fit is calibrated (the amplitude is known attenuated --
  # gcol33/tulpa#243).
  expect_gt(mean(cr[ok]), 0.35)
})

test_that("occu() + gp() identifies the NNGP range but not the amplitude", {
  skip_on_cran()
  skip_if_fast()

  sim <- .gpr_sim(N = 40L, J = 6L, seed = 1L, sigma_f = 1.2, phi_f = 0.25)
  fit <- .gpr_fit(sim, seed = 1L)

  # phi (range) recovers reasonably despite the amplitude funnel -- measured
  # 0.21 on this seed, 0.28-0.67 over several others (gcol33/tulpa#243) against
  # a truth of 0.25. This is the asymmetry that issue points at: the general
  # NNGP machinery (shared with icar/svc/spde) is sound, only the field/sigma2
  # coupling is not.
  phi_hat <- exp(as.numeric(fit$means[["log_phi_gp"]]))
  expect_gt(phi_hat, 0.05); expect_lt(phi_hat, 2.0)

  # No assertion on sigma2_gp or sd(fit$spatial_field) here. Measured
  # sd(field_hat)/sd(field_truth) at 1.7%-11.8% across seeds, unchanged by 5x
  # longer warmup and unrelated to divergence count (gcol33/tulpa#243) -- a
  # threshold here would either be so loose it asserts nothing or would record
  # the bug as expected behaviour. That claim belongs to the upstream issue,
  # not this test.
})
