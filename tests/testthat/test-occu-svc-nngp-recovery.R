# The continuous NNGP svc() surface: exposure + recovery (gcol33/tulpaObs#118).
#
# svc(lon, lat, indices=) was smoke-tested only: the fit ran and had the right
# shape, but the estimated surface was unnamed in the draws, so nothing could
# score it. `fit$svc_field` now exposes it (tulpa's ParamLayout carries the
# block offsets; the fitter returns them as `svc_layout`), which made an actual
# recovery test writable for the first time.
#
# What that test then found: the surface DOES track a known truth (cor 0.49-0.78
# over seeds), so the eta-assembly and the NNGP gradient were broadly right --
# but the sampler was not healthy. ~75% of post-warmup draws diverged and the
# NNGP range phi landed at ~4 against a truth of 0.25. The calibration
# assertions were skip()ped against gcol33/tulpa#144 rather than loosened until
# they passed, because a threshold tuned to a diverging sampler records the bug
# as the expected behaviour.
#
# Both causes are now fixed upstream and measured here (gcol33/tulpaObs#119):
# the range prior (tulpa#144, a Uniform behind a hard -INFINITY rejection, which
# gives NUTS no gradient to recover from) and the marginal-SD prior (the SVC
# half-Cauchy was improper on the coordinate it is sampled on, so nothing
# bounded sigma from above). Over seeds 1/2/3/11 at these settings:
# divergences 72-83% -> 0%, phi ~4 -> 0.14-0.23 against a truth of 0.25.
# Surface correlation did NOT move (0.73 mean, vs 0.66 before) -- it is bounded
# by the information in the data, not by sampler health, which is why the
# calibration test no longer asserts on it.
#
# svc() now requires prior_range = c(r0, alpha): tulpa ships the PC range
# anchors unset and refuses rather than inventing a default.
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
    ~ svc(lon, lat, indices = 1L, nn = 12, prior_range = c(0.1, 0.05)),
    data = sim$data,
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

  # This asserts the surface carries real spatial signal, NOT that the fit is
  # calibrated -- the hyperparameters are scored in the test below. Measured
  # 0.76 / 0.60 / 0.81 (mean 0.73) with the range and marginal-SD priors fixed,
  # against 0.49 / 0.72 / 0.78 beforehand. The floor is raised off the old 0.4
  # on those numbers rather than by tuning the simulation, and stays clear of
  # the 0.60 seed: what moved was the sampler, and the surface correlation is
  # bounded by the information at N = 150, J = 6, p = 0.6 either way.
  expect_gt(mean(cr[ok]), 0.55)
})

test_that("occu() + svc() is calibrated and identifies its NNGP hyperparameters", {
  skip_on_cran()
  skip_if_fast()

  sim <- .svc_sim(N = 150L, J = 6L, seed = 11L, sigma_f = 1.3, phi_f = 0.25)
  fit <- .svc_fit(sim, seed = 11L)

  # Divergences are now rare rather than the norm: 0 of 300 post-warmup draws
  # on every seed measured (1, 2, 3, 11), against 72-83% before tulpa#144.
  expect_lt(mean(fit$divergent), 0.05)

  # The range finds its truth instead of the old Uniform prior's mean.
  # Measured 0.14-0.23 over those seeds against a truth of 0.25; it sat at ~4.
  phi_hat <- exp(as.numeric(fit$means[["log_phi_svc[1]"]]))
  expect_gt(phi_hat, 0.1); expect_lt(phi_hat, 0.8)

  # The marginal SD is the weakly identified end of the GP ridge -- it trades
  # off against the range, and spreads 1.06-2.31 over seeds against a truth of
  # 1.3. The band is wide because that spread belongs to the model at N = 150,
  # not to the sampler.
  sigma_hat <- exp(0.5 * as.numeric(fit$means[["log_sigma2_svc[1]"]]))
  expect_gt(sigma_hat, 0.7); expect_lt(sigma_hat, 2.2)

  # No assertion on cor(surface, truth) here. This test used to close with
  # `expect_gt(cor, 0.75)` on the reasoning that the surface would track the
  # truth tightly "once the geometry is sane". Measuring it once the geometry
  # WAS sane refuted that: correlation barely moved across either fix (seeds
  # 1/2/3 went 0.49/0.72/0.78 -> 0.76/0.60/0.81, and seed 11 sits at 0.73), so
  # sampler health and surface accuracy are separate axes -- the latter is set
  # by how much the data says at N = 150, J = 6, p = 0.6. The surface is scored
  # by the recovery test above, which is where that claim belongs.
})
