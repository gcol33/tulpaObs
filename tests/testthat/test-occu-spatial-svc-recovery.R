# =============================================================================
# test-occu-spatial-svc-recovery.R - the standalone occu() nested-Laplace path
# carrying a varying-coefficient spatial bar: a cell-indexed spatial intercept
# field PLUS a spatial trend field weighted by a per-site covariate, with NO
# cover arm. This is occu_cover() with the cover arm removed (the
# apples-to-apples occupancy-only match for spOccupancy svcTPGOcc).
#
# Mirrors the occu_cover trend / MCAR recovery tests (cover arm dropped):
# simulate occupancy with a KNOWN spatial intercept field f1, a KNOWN spatial
# trend field f2 (weighted by a per-site `time` covariate carrying the temporal
# trend), a known occupancy covariate slope and detection p, fit standalone
# occu() with the bar, and recover both field SDs (sigma, sigma_trend), the
# covariate slope, the field shapes, and the psi / p separation.
#
# The bar carries the time dependence through its trend field (the mean-zero
# spatial trend), exactly as the occu_cover trend test does -- it does not put a
# competing global slope on the same covariate as the field weight (that is
# collinear with the SVC field and not what svcTPGOcc identifies).
# =============================================================================


.svc_chain_adj <- function(N) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  adj
}

# A smooth mean-zero field over the chain (a low-frequency sine), scaled to a
# target SD. Deterministic given the phase so the truth is known.
.svc_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Simulate single-season occupancy with a cell-indexed intercept field f1, a
# trend field f2 weighted by a per-site `time` covariate, and a separate
# occupancy covariate xocc with a known slope. Sites are cell x rep (several
# sites share one cell field node), exercising the group_var path. The occupancy
# logit is b0 + b_xocc * xocc + f1[cell] + time * f2[cell]; detection is a
# constant p over J visits.
.svc_simulate <- function(n_cells, reps, J, b0, b_xocc, p,
                          sigma_truth, sigma_trend_truth, seed) {
  set.seed(seed)
  adj <- .svc_chain_adj(n_cells)
  f1  <- .svc_smooth_field(n_cells, sigma_truth,       phase = 0.7)
  f2  <- .svc_smooth_field(n_cells, sigma_trend_truth, phase = 2.3)

  n_sites <- n_cells * reps
  cell    <- rep(seq_len(n_cells), each = reps)
  time    <- as.numeric(scale(stats::rnorm(n_sites)))   # trend-field weight
  xocc    <- as.numeric(scale(stats::rnorm(n_sites)))   # plain occupancy covariate

  eta_psi <- b0 + b_xocc * xocc + f1[cell] + time * f2[cell]
  z       <- stats::rbinom(n_sites, 1L, plogis(eta_psi))

  y <- matrix(0L, n_sites, J)
  for (i in which(z == 1L)) y[i, ] <- stats::rbinom(J, 1L, p)

  list(adj = adj, y = y, f1 = f1, f2 = f2,
       data = data.frame(cell = cell, time = time, xocc = xocc))
}

.svc_fit <- function(sim, max.iter = 25L) {
  suppressWarnings(tobs(
    formula  = ~ xocc + spatial(~ 1 + time || cell, graph = sim$adj),
    data     = sim$data, family = occu(),
    detection = ~ 1, y = sim$y,
    method   = "nested_laplace",
    control  = list(verbose = FALSE, max.iter = max.iter, progress = FALSE)
  ))
}


test_that("occu() spatial bar fits end-to-end and exposes both field SDs", {
  skip_if_fast()
  n_cells <- 24L; reps <- 3L; J <- 5L
  sim <- .svc_simulate(
    n_cells = n_cells, reps = reps, J = J,
    b0 = stats::qlogis(0.45), b_xocc = 0.6, p = 0.45,
    sigma_truth = 0.9, sigma_trend_truth = 0.7, seed = 4242L
  )

  fit <- .svc_fit(sim, max.iter = 25L)
  expect_s3_class(fit, "tobs_fit")

  # Both field SDs surface as named hyperparameters; the betas and SDs are finite.
  expect_true(all(c("sigma", "sigma_trend") %in% names(fit$means)))
  expect_true(is.finite(fit$means[["sigma"]]))
  expect_true(is.finite(fit$means[["sigma_trend"]]))
  expect_gt(fit$means[["sigma"]], 0)
  expect_gt(fit$means[["sigma_trend"]], 0)

  # Both cell-indexed fields are returned, demeaned, length n_cells.
  expect_length(fit$spatial_field, n_cells)
  expect_length(fit$trend_field, n_cells)
  expect_true(all(is.finite(fit$spatial_field)))
  expect_true(all(is.finite(fit$trend_field)))
  expect_lt(abs(mean(fit$spatial_field)), 1e-6)
  expect_lt(abs(mean(fit$trend_field)),   1e-6)

  # The occupancy-only fit carries no cover (pos) arm.
  expect_false(any(grepl("^pos_", names(fit$means))))
})


test_that("occu() spatial bar recovers both fields, covariate slope, psi/p (6 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 6L
  n_cells <- 40L; reps <- 6L; J <- 8L
  b0_truth          <- stats::qlogis(0.5)
  b_xocc_truth      <- 0.6
  p_truth           <- 0.45
  sigma_truth       <- 0.9
  sigma_trend_truth <- 0.7

  est_sigma <- est_sigma_tr <- est_bxocc <- est_p0 <-
    cor1 <- cor2 <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- .svc_simulate(
      n_cells = n_cells, reps = reps, J = J,
      b0 = b0_truth, b_xocc = b_xocc_truth, p = p_truth,
      sigma_truth = sigma_truth, sigma_trend_truth = sigma_trend_truth,
      seed = 9000L + s
    )
    fit <- tryCatch(.svc_fit(sim, max.iter = 30L), error = function(e) NULL)
    if (is.null(fit)) next

    est_sigma[s]    <- fit$means[["sigma"]]
    est_sigma_tr[s] <- fit$means[["sigma_trend"]]
    est_bxocc[s]    <- fit$means[["psi_xocc"]]
    est_p0[s]       <- fit$means[["p_(Intercept)"]]
    cor1[s]         <- abs(stats::cor(fit$spatial_field, sim$f1))
    cor2[s]         <- abs(stats::cor(fit$trend_field,   sim$f2))
  }

  ok <- is.finite(est_sigma) & is.finite(est_bxocc)
  expect_gte(sum(ok), 5L)

  # Both field SHAPES recover with high fidelity (mean |cor| > 0.75). This is the
  # core acceptance for the varying-coefficient bar: the cell-indexed intercept
  # field AND the per-site-weighted trend field are each estimated from the
  # occupancy data alone, no cover arm. The trend field carries the per-site
  # `time` weight, so its recovery is the direct evidence the SVC weight is fit.
  expect_gt(mean(cor1[ok]), 0.75)
  expect_gt(mean(cor2[ok]), 0.75)

  # Both field SDs are finite, positive, and distinguished (the trend field's SD
  # is reported separately from the intercept field's). The absolute field SD is
  # attenuated toward zero relative to truth here: a latent areal field in a
  # binary occupancy GLM is weakly identified in amplitude at this n (the same
  # attenuation the existing single-field occu() icar path shows), so the
  # right-skewed sigma = 1/sqrt(tau) posterior sits below the truth. The shape
  # (above) is the recoverable part; the SDs are gated for sanity, not bias.
  expect_true(all(is.finite(est_sigma[ok])) && all(est_sigma[ok] > 0))
  expect_true(all(is.finite(est_sigma_tr[ok])) && all(est_sigma_tr[ok] > 0))
  expect_lt(mean(est_sigma[ok]),    sigma_truth * 1.3)
  expect_lt(mean(est_sigma_tr[ok]), sigma_trend_truth * 1.3)

  # The occupancy covariate slope recovers within ~0.5 (it carries the same mild
  # upward bias the existing occu() areal path shows; the SVC fields do not
  # worsen it).
  expect_lt(abs(mean(est_bxocc[ok]) - b_xocc_truth), 0.5)

  # psi / p separation holds: detection p is recovered (the marginalized state
  # likelihood separates psi from p, not the visit-stacked Bernoulli that
  # conflates them). The detection intercept is the cleanly-identified anchor of
  # the separation.
  expect_lt(abs(mean(est_p0[ok]) - stats::qlogis(p_truth)), 0.30)
})
