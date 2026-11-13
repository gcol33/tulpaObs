# The continuous NNGP svc() surface: exposure + recovery.
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
# assertions were skip()ped rather than loosened until they passed,
# because a threshold tuned to a diverging sampler records the bug as the
# expected behaviour.
#
# Both causes are now fixed upstream and measured here: the range prior (a
# Uniform behind a hard -INFINITY rejection, which gives NUTS no gradient to
# recover from) and the marginal-SD prior (the SVC half-Cauchy was improper on
# the coordinate it is sampled on, so nothing bounded sigma from above). Over
# seeds 1/2/3/11 at these settings: divergences 72-83% -> 0%, phi ~4 ->
# 0.14-0.23 against a truth of 0.25. Surface correlation did NOT move (0.73
# mean, vs 0.66 before) -- it is bounded by the information in the data, not by
# sampler health, which is why the calibration test no longer asserts on it.
#
# svc() now requires prior.range = c(r0, alpha): tulpa ships the PC range
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
#
# The marginal-SD prior is stated, at a unit median (P(sigma > 1) = 0.5). The
# default P(sigma > 1) = 0.01 has median ~0.15 and puts 0.25% of its mass above
# the simulated sigma_f = 1.3, so under it the fit reports the prior: seed 11
# landed at sigma 0.18 on Linux and 0.71 on Windows, with phi pushed up the
# ridge. A recovery test needs a prior that does not exclude its own truth.
.svc_fit <- function(sim, seed, n_iter = 600L, n_warmup = 300L) {
  suppressWarnings(tobs(
    ~ svc(lon, lat, indices = 1L, nn = 12, prior.range = c(0.1, 0.05),
          sigma2.prior.U = 1, sigma2.prior.alpha = 0.5),
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

# Seeds 1/2/3/11 at N = 150, J = 6, p = 0.6, fitted once and shared by the two
# blocks below.
.svc_cal_seeds <- c(1L, 2L, 3L, 11L)
.svc_cal <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- lapply(.svc_cal_seeds, function(s) {
        sim <- .svc_sim(N = 150L, J = 6L, seed = s)
        fit <- try(.svc_fit(sim, seed = s), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        sg <- exp(0.5 * fit$draws[, "log_sigma2_svc[1]"])
        ph <- exp(fit$draws[, "log_phi_svc[1]"])
        list(cor = stats::cor(as.numeric(fit$svc_field), sim$f),
             div = mean(fit$divergent),
             sigma_q = stats::quantile(sg, c(0.05, 0.95), names = FALSE),
             phi_q = stats::quantile(ph, c(0.05, 0.95), names = FALSE))
      })
      names(cache) <<- .svc_cal_seeds
    }
    cache
  }
})

test_that("occu() + svc() recovers the shape of a known varying-intercept surface", {
  skip_on_cran()
  skip_if_fast()

  res <- .svc_cal()[c("1", "2", "3")]
  cr <- vapply(res, function(r) if (is.null(r)) NA_real_ else r$cor, 0)
  ok <- !is.na(cr)
  expect_gt(sum(ok), 1L)

  # The surface carries real spatial signal; the hyperparameters are scored in
  # the block below. Measured 0.77 / 0.60 / 0.80 on seeds 1/2/3, and 0.33-0.81
  # (median 0.71) over seeds 1-20. The correlation is bounded by the
  # information at N = 150, J = 6, p = 0.6, not by the sampler.
  expect_gt(mean(cr[ok]), 0.55)
})

test_that("occu() + svc() is calibrated and identifies its NNGP hyperparameters", {
  skip_on_cran()
  skip_if_fast()

  res <- Filter(Negate(is.null), .svc_cal())
  expect_gte(length(res), 3L)

  # Divergences on every seed. Over seeds 1-20 the rate was 0 on 18 and at
  # most 1.3% on the other two.
  expect_true(all(vapply(res, function(r) r$div, 0) < 0.05))

  # sigma and phi sit on the GP ridge: a single seed's point estimate spreads
  # 0.30-3.37 (sigma, truth 1.3) and 0.20-1.18 (phi, truth 0.25) over seeds
  # 1-20, and which value a seed lands on moves with platform and BLAS. What
  # holds across seeds is the posterior interval: the 90% interval covered
  # sigma on 17 of 20 seeds and phi on 20 of 20. Requiring coverage on at least
  # 2 of the 4 seeds fails a calibrated fit with probability ~1% at 85%
  # per-seed coverage, and fails a fit whose SD collapses onto its prior on
  # every seed.
  covers <- function(q, truth) q[1L] <= truth && truth <= q[2L]
  expect_gte(sum(vapply(res, function(r) covers(r$sigma_q, 1.3), TRUE)), 2L)
  expect_gte(sum(vapply(res, function(r) covers(r$phi_q, 0.25), TRUE)), 2L)
})
