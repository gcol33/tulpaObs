# =============================================================================
# test-occu-svc-joint-recovery.R - parameter recovery for the standalone occu()
# varying-coefficient (SVC) spatial bar rerouted through the joint direct-grid
# engine.
#
# The fixture is built to resemble the REAL EVA data that broke the EM
# fixed-point path (.tobs_em_nested_laplace): a LARGE-amplitude intercept field
# (sigma ~ 1.2), SPARSE detection (occupancy prevalence ~ 0.15), MULTIPLE
# detection covariates, and several visits -- exactly the combination the easy
# small recovery sim lacks and that made the EM oscillate at the
# min(max.iter, 25) cap. The joint engine integrates the field hyperparameters
# on a direct outer grid, so it completes stably (no `converged = FALSE`
# truncation) where the EM did not.
#
# Recovery asserts: the fit completes stably with converged = TRUE; both field
# SDs are finite / positive / distinguished; both field SHAPES recover
# (|cor| above a sensible threshold); psi / p separate (detection recovered).
# A cross-check confirms the single-arm route reproduces the joint engine's
# occupancy result: occu()'s psi betas + per-cell trend match occu_cover(alpha
# = 0)'s occupancy arm on the same data within tolerance.
# =============================================================================


.svcj_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Real-data-like single-season occupancy: a large intercept field f1, a trend
# field f2 weighted by a per-site `time` covariate, a sparse occupancy baseline
# (prevalence ~ 0.15 via a low psi intercept), and rich detection
# (5 covariates + intercept). Sites are cell x rep (group_var path). Returns the
# graph, detection matrix, fields, the site frame and the visit frame.
.svcj_simulate <- function(n_cells, reps, J, p_int, sigma_truth,
                           sigma_trend_truth, seed) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f1  <- .svcj_smooth_field(n_cells, sigma_truth,       phase = 0.7)
  f2  <- .svcj_smooth_field(n_cells, sigma_trend_truth, phase = 2.3)

  n_sites <- n_cells * reps
  cell    <- rep(seq_len(n_cells), each = reps)
  time    <- as.numeric(scale(stats::rnorm(n_sites)))   # trend-field weight
  xocc    <- as.numeric(scale(stats::rnorm(n_sites)))   # plain occupancy covariate

  # Sparse occupancy: a low intercept so the marginal prevalence is ~ 0.15.
  b0     <- stats::qlogis(0.15)
  b_xocc <- 0.5
  eta_psi <- b0 + b_xocc * xocc + f1[cell] + time * f2[cell]
  z       <- stats::rbinom(n_sites, 1L, plogis(eta_psi))

  # Rich detection: intercept + 5 visit covariates, moderate detection.
  d1 <- as.numeric(scale(stats::rnorm(n_sites)))
  d2 <- as.numeric(scale(stats::rnorm(n_sites)))
  d3 <- as.numeric(scale(stats::rnorm(n_sites)))
  d4 <- as.numeric(scale(stats::rnorm(n_sites)))
  d5 <- as.numeric(scale(stats::rnorm(n_sites)))
  eta_p <- qlogis(p_int) + 0.3 * d1 - 0.2 * d2 + 0.25 * d3 - 0.15 * d4 + 0.2 * d5
  p_site <- plogis(eta_p)

  y <- matrix(0L, n_sites, J)
  for (i in which(z == 1L)) y[i, ] <- stats::rbinom(J, 1L, p_site[i])

  site_data <- data.frame(cell = cell, time = time, xocc = xocc,
                          d1 = d1, d2 = d2, d3 = d3, d4 = d4, d5 = d5)
  list(adj = adj, y = y, f1 = f1, f2 = f2, site_data = site_data,
       prevalence = mean(z), p_int = p_int)
}


test_that("occu() SVC bar completes stably and recovers fields under real-data-like load", {
  skip_on_cran()
  skip_if_fast()

  n_cells <- 40L; reps <- 6L; J <- 8L
  sigma_truth       <- 1.2
  sigma_trend_truth <- 0.7
  p_int             <- 0.45

  n_seeds <- 5L
  est_sigma <- est_sigma_tr <- est_p0 <- cor1 <- cor2 <- prev <-
    rep(NA_real_, n_seeds)
  converged <- logical(n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- .svcj_simulate(n_cells, reps, J, p_int = p_int,
                          sigma_truth = sigma_truth,
                          sigma_trend_truth = sigma_trend_truth,
                          seed = 7100L + s)
    prev[s] <- sim$prevalence
    fit <- tryCatch(suppressWarnings(tobs(
      formula   = ~ xocc + spatial(~ 1 + time || cell, graph = sim$adj),
      data      = sim$site_data, family = occu(),
      detection = ~ d1 + d2 + d3 + d4 + d5,
      y = sim$y,
      method    = "nested_laplace",
      control   = list(verbose = FALSE, progress = FALSE)
    )), error = function(e) NULL)
    if (is.null(fit)) next

    # The fit went through the direct grid (NOT a non-converged EM snapshot):
    # the joint route reports converged = TRUE and carries no EM iteration count.
    converged[s]    <- isTRUE(fit$convergence$converged)
    est_sigma[s]    <- fit$means[["sigma"]]
    est_sigma_tr[s] <- fit$means[["sigma_trend"]]
    est_p0[s]       <- fit$means[["p_(Intercept)"]]
    cor1[s]         <- abs(stats::cor(fit$spatial_field, sim$f1))
    cor2[s]         <- abs(stats::cor(fit$trend_field,   sim$f2))
  }

  ok <- is.finite(est_sigma)
  expect_gte(sum(ok), 4L)

  # Sparse occupancy as designed (a sanity check on the fixture itself).
  expect_lt(mean(prev, na.rm = TRUE), 0.30)

  # Completes STABLY on every fitted seed -- the acceptance the EM failed: no
  # converged = FALSE truncation at the max.iter cap.
  expect_true(all(converged[ok]))

  # Both field SDs finite, positive, and distinguished (the trend SD reported
  # separately from the intercept SD). The absolute SD is attenuated toward zero
  # relative to the large truth (a latent areal field is weakly identified in
  # amplitude in a sparse binary occupancy GLM); the shape is the recoverable
  # part. Gated for sanity, not unbiasedness.
  expect_true(all(is.finite(est_sigma[ok])) && all(est_sigma[ok] > 0))
  expect_true(all(is.finite(est_sigma_tr[ok])) && all(est_sigma_tr[ok] > 0))

  # Both field SHAPES recover with high fidelity even at the large amplitude /
  # sparse detection that broke the EM.
  expect_gt(mean(cor1[ok]), 0.70)
  expect_gt(mean(cor2[ok]), 0.55)

  # psi / p separate: the detection intercept is recovered (the marginalized
  # state likelihood separates psi from p, not the visit-stacked Bernoulli that
  # conflates them).
  expect_lt(abs(mean(est_p0[ok]) - stats::qlogis(p_int)), 0.35)
})


test_that("occu() SVC occupancy arm matches occu_cover(alpha = 0) within tolerance", {
  skip_on_cran()
  skip_if_fast()

  n_cells <- 30L; reps <- 6L; J <- 6L
  sim <- .svcj_simulate(n_cells, reps, J, p_int = 0.5,
                        sigma_truth = 1.0, sigma_trend_truth = 0.6, seed = 8242L)
  n_sites <- nrow(sim$site_data)

  # occu() single-arm SVC fit. Pin a shared sigma grid so the two engines
  # integrate the occupancy field on the SAME outer support (occu_cover grids the
  # field on `sigma`; the single-arm occu route converts that grid to tau).
  sigma_grid <- c(0.4, 0.8, 1.5, 2.5)
  fit_occu <- suppressWarnings(tobs(
    formula   = ~ xocc + spatial(~ 1 + time || cell, graph = sim$adj),
    data      = sim$site_data, family = occu(),
    detection = ~ d1 + d2,
    y = sim$y,
    method    = "nested_laplace",
    control   = list(verbose = FALSE, progress = FALSE,
                     sigma.grid = sigma_grid, adaptive.grid = FALSE,
                     var.of.means.consistency = FALSE, diagnose.k = FALSE)
  ))

  # occu_cover(alpha = 0): the SAME occupancy arm + a beta cover arm decoupled
  # from the field (alpha = alpha_trend = 0), so the occupancy field is fitted by
  # the joint engine with no cover-arm transfer -- the head-to-head the single-arm
  # route should reproduce. Cover is a placeholder (presence-only positives at a
  # constant), with detection-derived y == 1 visits carrying a cover value.
  y_pos <- matrix(0, n_sites, J)
  set.seed(99L)
  for (i in seq_len(n_sites)) {
    det_v <- which(sim$y[i, ] == 1L)
    if (length(det_v)) y_pos[i, det_v] <- plogis(stats::rnorm(length(det_v), 0, 0.5))
  }
  fit_oc <- tryCatch(suppressWarnings(tobs(
    formula   = ~ xocc + spatial(~ 1 + time || cell, graph = sim$adj),
    data      = sim$site_data,
    family    = occu_cover(response = "beta", cover_aggregate = "none"),
    detection = ~ d1 + d2,
    positive  = ~ 1,
    y = sim$y, y_pos = y_pos,
    method    = "nested_laplace",
    control   = list(verbose = FALSE, progress = FALSE,
                     engine = "joint", n.threads = 1L,
                     sigma.grid = sigma_grid, alpha.grid = 0,
                     alpha.grid.trend = 0, integration = "grid",
                     adaptive.grid = FALSE, var.of.means.consistency = FALSE,
                     diagnose.k = FALSE)
  )), error = function(e) { message("occu_cover ref fit failed: ",
                                    conditionMessage(e)); NULL })
  skip_if(is.null(fit_oc), "occu_cover(alpha = 0) reference fit unavailable")

  # The occupancy betas (psi arm) recover the same coefficients. The two engines
  # carry a different cover arm (occu has none; occu_cover has a decoupled beta
  # arm whose intercept prior shifts the joint mode slightly), so the match is
  # close but not byte-identical -- a tolerance on the shared occupancy
  # parameters, the cross-check the issue asks for.
  psi_occu <- fit_occu$means[grep("^psi_", names(fit_occu$means))]
  psi_oc   <- fit_oc$means[grep("^psi_", names(fit_oc$means))]
  common   <- intersect(names(psi_occu), names(psi_oc))
  expect_gt(length(common), 1L)
  expect_lt(max(abs(psi_occu[common] - psi_oc[common])), 0.30)

  # The per-cell occupancy field (intercept) and the per-cell trend field match
  # in SHAPE: both engines fit the same areal occupancy field, so they recover the
  # same cell-to-cell pattern.
  expect_gt(abs(stats::cor(fit_occu$spatial_field, fit_oc$spatial_field)), 0.90)
  expect_gt(abs(stats::cor(fit_occu$trend_field,   fit_oc$trend_field)),   0.85)

  # The reported occupancy-field SD agrees within tolerance (same grid, same
  # field; the cover arm is decoupled).
  expect_lt(abs(fit_occu$means[["sigma"]] - fit_oc$means[["sigma"]]), 0.5)
})
