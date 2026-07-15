# The continuous NNGP svc() surface: exposure + recovery (gcol33/tulpaObs#118).
#
# svc(lon, lat, indices=) was smoke-tested only: the fit ran and had the right
# shape, but the estimated surface was unnamed in the draws, so nothing could
# score it. `fit$svc_field` now exposes it (tulpa's ParamLayout carries the
# block offsets; the fitter returns them as `svc_layout`), which made an actual
# recovery test writable for the first time.
#
# What that test then found: the surface DOES track a known truth (cor 0.49-0.78
# over seeds), so the eta-assembly and the NNGP gradient are broadly right -- but
# the sampler is not healthy. ~75% of post-warmup draws diverge and the NNGP
# range phi lands at ~4 against a truth of 0.25, because tulpa priors phi as
# Uniform(lower, upper) (mean ~5 at the defaults) behind a hard -INFINITY
# rejection. See gcol33/tulpa#144. The calibration assertions below are
# skip()ped against that issue rather than loosened until they pass: a threshold
# tuned to a diverging sampler would record the bug as the expected behaviour.
#
# This is the continuous-NNGP flavour, NOT the areal weighted-bar SVC -- that
# arm arrives as a `spatial` term and is recovery-tested in
# test-occu-spatial-svc-recovery.R / test-occu-svc-joint-recovery.R. See the SVC
# note in CLAUDE.md.

# A known spatially-varying INTERCEPT surface drawn from the model the term
# assumes: an exponential-kernel GP over the site coordinates.
.svc_sim <- function(N, J, seed, sigma_f = 1.3, phi_f = 0.25, p_det = 0.6,
                     b0 = 0) {
  set.seed(seed)
  lon <- stats::runif(N); lat <- stats::runif(N)
  D <- as.matrix(stats::dist(cbind(lon, lat)))
  K <- sigma_f^2 * exp(-D / phi_f)
  f <- as.numeric(t(chol(K + 1e-8 * diag(N))) %*% stats::rnorm(N))
  # A varying intercept is confounded with the global intercept up to a
  # constant, so centre the truth and score the SHAPE -- the same reason the
  # areal SVC tests score fit$spatial_field with cor() rather than absolute
  # values.
  f <- f - mean(f)
  z <- stats::rbinom(N, 1, stats::plogis(b0 + f))
  y <- matrix(stats::rbinom(N * J, 1, p_det * rep(z, J)), N, J)
  list(data = data.frame(lon = lon, lat = lat), y = y, f = f)
}

# indices = 1L targets the occupancy intercept column. Deliberate: the intercept
# column is 1.0 in both the natural and the autoscaled design, and
# .unscale_fit_per_process() rewrites only the leading per-process beta slices,
# never the SVC block -- so for any other `indices` the returned weights are on
# the scaled-column scale and truth would need rescaling before scoring.
.svc_fit <- function(sim, seed, n_iter = 600L, n_warmup = 300L) {
  suppressWarnings(tobs(
    ~ svc(lon, lat, indices = 1L, nn = 12), data = sim$data,
    family = occu(), detection = ~ 1, y = sim$y,
    method = "nuts",
    control = list(n.iter = n_iter, n.warmup = n_warmup, seed = seed,
                   verbose = FALSE)))
}

test_that("the svc() surface is exposed, named, and correctly shaped", {
  skip_on_cran()
  skip_if_fast()

  sim <- .svc_sim(N = 60L, J = 4L, seed = 1L)
  fit <- .svc_fit(sim, seed = 1L, n_iter = 120L, n_warmup = 60L)

  # The block layout tulpa reports, and the surface sliced from it by position
  # (not by parsing names).
  expect_false(is.null(fit$svc_layout))
  expect_identical(as.integer(fit$svc_layout$n_obs), 60L)
  expect_identical(as.integer(fit$svc_layout$n_svc), 1L)

  # One weight per site, exposed as a bare vector for a single varying
  # coefficient (mirrors fit$spatial_field on the areal path).
  expect_length(as.numeric(fit$svc_field), 60L)
  expect_true(all(is.finite(as.numeric(fit$svc_field))))

  # The block is named rather than falling through to "param[k]" -- that
  # fallthrough is exactly why the surface used to be unreadable.
  cn <- colnames(fit$draws)
  expect_equal(sum(grepl("^svc_w", cn)), 60L)
  expect_equal(sum(grepl("^log_sigma2_svc", cn)), 1L)
  expect_equal(sum(grepl("^log_phi_svc", cn)), 1L)
  expect_equal(sum(grepl("^param", cn)), 0L)

  # The surface is the posterior mean of its own draw columns.
  d <- attr(fit$svc_field, "draws")
  expect_equal(dim(d), c(nrow(fit$draws), 60L))
  expect_equal(as.numeric(fit$svc_field), unname(colMeans(d)), tolerance = 1e-10)

  # fit$svc stays the input term; fit$svc_field is what was estimated.
  expect_s3_class(fit$svc, "tobs_svc")
})

test_that("occu() + svc() recovers the shape of a known varying-intercept surface", {
  skip_on_cran()
  skip_if_fast()

  seeds <- 1:3
  cr <- rep(NA_real_, length(seeds))
  for (i in seq_along(seeds)) {
    sim <- .svc_sim(N = 150L, J = 6L, seed = seeds[i])
    fit <- try(.svc_fit(sim, seed = seeds[i]), silent = TRUE)
    if (inherits(fit, "try-error")) next
    cr[i] <- stats::cor(as.numeric(fit$svc_field), sim$f)
  }
  ok <- !is.na(cr)
  expect_gt(sum(ok), 1L)

  # A deliberately loose floor: this asserts the surface carries real spatial
  # signal (the gradient is not wrong), NOT that the fit is calibrated. Measured
  # 0.49 / 0.72 / 0.78 at these settings, against a sampler that diverges ~75%
  # of the time (gcol33/tulpa#144). Tighten this once #144 lands -- do not
  # tighten it by tuning the simulation.
  expect_gt(mean(cr[ok]), 0.4)
})

test_that("occu() + svc() is calibrated and identifies its NNGP hyperparameters", {
  skip_on_cran()
  skip_if_fast()
  skip(paste0("blocked on gcol33/tulpa#144: svc() priors phi as Uniform(lower, ",
              "upper) behind a hard -INFINITY rejection, so NUTS diverges on ",
              "~75% of post-warmup draws and phi collapses onto the prior mean ",
              "(~4 vs a truth of 0.25). Enable once the range is ",
              "reparameterized."))

  sim <- .svc_sim(N = 150L, J = 6L, seed = 11L, sigma_f = 1.3, phi_f = 0.25)
  fit <- .svc_fit(sim, seed = 11L)

  # Divergences should be rare, not the norm.
  expect_lt(mean(fit$divergent), 0.05)

  # The marginal SD and the range should both find their truth.
  sigma_hat <- exp(0.5 * as.numeric(fit$means[["log_sigma2_svc[1]"]]))
  expect_gt(sigma_hat, 0.7); expect_lt(sigma_hat, 2.2)

  phi_hat <- exp(as.numeric(fit$means[["log_phi_svc[1]"]]))
  expect_gt(phi_hat, 0.1); expect_lt(phi_hat, 0.8)

  # And the surface should track the truth tightly once the geometry is sane.
  expect_gt(stats::cor(as.numeric(fit$svc_field), sim$f), 0.75)
})
