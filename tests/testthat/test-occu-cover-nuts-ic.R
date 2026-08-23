# =============================================================================
# test-occu-cover-nuts-ic.R - the occu_cover() information criteria score the
# structured terms a sampled fit carries.
#
# A `(1 | g)` on the detection / positive-cover formula samples under method =
# "nuts" and the coupled areal field samples with its hypers. Both are OFFSETS
# rather than coefficients -- the random effect per (site, visit), the field per
# site -- so .tobs_occu_cover_components() returns them alongside the
# coefficient draws and every diagnostic built on it folds them in: the
# pointwise log-likelihood (WAIC / LOO / CPO), the posterior predictive check,
# and the PIT / LOO-PIT. Tests:
#   - the per-visit offsets equal the fitter's own per-(site, visit) predictor
#     construction, on both observation arms
#   - the criteria MOVE with the random effect, and a zeroed offset reproduces
#     the no-offset score to the bit (the "RE at zero variance" reduction)
#   - byte-identity: a fit carrying neither term scores exactly as it did with
#     no offset argument at all, on the log-likelihood, the PPC and the PIT
#   - the sampled coupled field reaches the same criteria, per site
# =============================================================================


.ocic_sim <- function(seed, arm = "p", N = 100L, J = 5L, n_g = 8L,
                      sigma_re = 0.9) {
  simulate_occu_cover(
    N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", sigma_pos = 0.4,
    re_det_groups = if (identical(arm, "p")) n_g else NULL, sigma_re_p = sigma_re,
    re_pos_groups = if (identical(arm, "pos")) n_g else NULL,
    sigma_re_pos = sigma_re, seed = seed)
}

.ocic_fit <- function(sim, detection = ~ det_cov1, positive = ~ pos_cov1,
                      n.iter = 250L, seed = 3L) {
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  tobs(formula = ~ occ_cov1, data = sim$data, family = occu_cover("lognormal"),
       detection = detection, positive = positive, y = sim$y, y_pos = y_pos,
       visits = sim$visit_data, method = "nuts",
       control = list(verbose = FALSE, progress = FALSE, n.iter = n.iter,
                      n.warmup = n.iter, n.chains = 1L, seed = seed))
}

# lppd, the WAIC penalty and elpd from an [S x n_sites] pointwise log-likelihood.
.ocic_score <- function(ll) {
  lse  <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
  lppd <- sum(apply(ll, 2L, lse))
  pw   <- sum(apply(ll, 2L, stats::var))
  c(lppd = lppd, p_waic = pw, elpd = lppd - pw)
}

# The pointwise log-likelihood of a fit's own draws, with the per-visit offsets
# passed, suppressed, or replaced by explicit zeros.
.ocic_ll <- function(fit, c0, off = c("as_fitted", "absent", "zero")) {
  off <- match.arg(off)
  V   <- tulpaObs:::.occu_cover_visit_view(fit$model)$V
  S   <- nrow(c0$b_occ)
  args <- switch(off,
    as_fitted = list(off_det = c0$off_det, off_pos = c0$off_pos),
    absent    = list(),
    zero      = list(off_det = matrix(0, V, S), off_pos = matrix(0, V, S)))
  do.call(tulpaObs:::.occu_cover_ploglik_core,
          c(list(fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp,
                 c0$field_occ, c0$field_pos), args))
}


test_that("the per-visit RE offsets equal the fitter's own predictor", {
  skip_on_cran()
  skip_if_fast()
  for (arm in c("p", "pos")) {
    sim <- .ocic_sim(if (identical(arm, "p")) 11L else 12L, arm = arm)
    det <- if (identical(arm, "p")) ~ det_cov1 + (1 | habitat) else ~ det_cov1
    pos <- if (identical(arm, "pos")) ~ pos_cov1 + (1 | habitat) else ~ pos_cov1
    fit <- .ocic_fit(sim, det, pos)

    S  <- 60L
    c0 <- tulpaObs:::.tobs_occu_cover_components(fit, S)
    expect_identical(is.null(c0$off_det), !identical(arm, "p"))
    expect_identical(is.null(c0$off_pos), !identical(arm, "pos"))
    got <- if (identical(arm, "p")) c0$off_det else c0$off_pos
    expect_identical(dim(got), c(sum(fit$model$valid), S))

    # The fitter builds the same offset on the padded [n_sites x max_visits]
    # grid (.occu_cover_nuts_re_offsets, the predictor its own log-likelihood is
    # evaluated at). Flattening that grid site-major, visit-ascending and taking
    # the valid cells is exactly the view's per-visit rows, so the two
    # constructions have to agree cell for cell at every draw.
    blocks <- tulpaObs:::.occu_cover_nuts_re_blocks(fit$model)
    vw     <- tulpaObs:::.occu_cover_visit_view(fit$model)
    N <- fit$model$n_sites; J <- fit$model$max_visits
    for (d in c(1L, 17L, S)) {
      k <- 0L
      b_list <- lapply(blocks, function(b) {
        idx <- k + seq_len(b$n_groups); k <<- k + b$n_groups + 1L
        exp(fit$re_draws[d, k]) * fit$re_draws[d, idx]
      })
      ref <- tulpaObs:::.occu_cover_nuts_re_offsets(blocks, b_list, N, J)
      ref <- as.numeric(t(ref[[arm]]))[vw$flat_idx]
      expect_equal(got[, d], ref, tolerance = 1e-12)
    }
  }
})


test_that("the criteria move with the random effect and reduce to it at zero", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocic_sim(11L, arm = "p")
  fit <- .ocic_fit(sim, detection = ~ det_cov1 + (1 | habitat))
  expect_gt(fit$re[[1L]]$sigma_median, 0.4)

  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, 250L)
  on  <- .ocic_score(.ocic_ll(fit, c0, "as_fitted"))
  off <- .ocic_score(.ocic_ll(fit, c0, "absent"))

  # A detection random intercept the sampler put at sigma ~ 1 explains a large
  # share of the detection heterogeneity, so dropping it from the score costs
  # many nats of pointwise log-likelihood -- the criteria describe a different
  # model with and without it.
  expect_gt(on[["lppd"]] - off[["lppd"]], 5)
  expect_gt(on[["elpd"]] - off[["elpd"]], 5)

  # A random effect held at zero is the no-random-effect model, to the bit.
  expect_identical(.ocic_ll(fit, c0, "zero"), .ocic_ll(fit, c0, "absent"))

  # waic() is the scored version, not the zeroed one.
  w <- waic(fit, n.draws = 250L)
  expect_equal(w$elpd_waic, unname(on[["elpd"]]), tolerance = 1e-8)
  expect_gt(w$elpd_waic, off[["elpd"]] + 5)
})


test_that("a cover-arm random effect is scored on the cover arm", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocic_sim(12L, arm = "pos")
  fit <- .ocic_fit(sim, positive = ~ pos_cov1 + (1 | habitat))
  c0  <- tulpaObs:::.tobs_occu_cover_components(fit, 250L)
  expect_null(c0$off_det)
  expect_false(is.null(c0$off_pos))
  on  <- .ocic_score(.ocic_ll(fit, c0, "as_fitted"))
  off <- .ocic_score(.ocic_ll(fit, c0, "absent"))
  expect_gt(on[["elpd"]] - off[["elpd"]], 5)
  expect_identical(.ocic_ll(fit, c0, "zero"), .ocic_ll(fit, c0, "absent"))
})


test_that("a fit carrying neither term scores bit for bit as with no offset", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocic_sim(11L, arm = "p")
  fit <- .ocic_fit(sim)                       # habitat left out of the formula
  c0  <- tulpaObs:::.tobs_occu_cover_components(fit, 200L)
  expect_null(c0$off_det)
  expect_null(c0$off_pos)
  expect_true(all(c0$field_occ == 0) && all(c0$field_pos == 0))

  expect_identical(.ocic_ll(fit, c0, "as_fitted"), .ocic_ll(fit, c0, "absent"))
  expect_identical(.ocic_ll(fit, c0, "zero"),      .ocic_ll(fit, c0, "absent"))

  vw <- tulpaObs:::.occu_cover_visit_view(fit$model)
  S  <- nrow(c0$b_occ)
  zero <- matrix(0, vw$V, S)
  none <- matrix(0, vw$V, 0L)
  vd <- function(X) tulpaObs:::.occu_cover_visit_design(X, vw$V)
  ppc_kernel <- function(od, op) {
    set.seed(101)
    tulpaObs:::cpp_occu_cover_ppc(
      X_occ = fit$model$X_occ, X_det_site = fit$model$X_det_site,
      X_pos_site = fit$model$X_pos_site,
      X_det_visit = vd(vw$X_det_visit), X_pos_visit = vd(vw$X_pos_visit),
      site_of_visit = vw$site_of_visit, y_det_visit = vw$y_det_visit,
      y_pos_visit = vw$y_pos_visit, b_occ = c0$b_occ, b_det = c0$b_det,
      b_pos = c0$b_pos, disp = c0$disp, field_occ = c0$field_occ,
      field_pos = c0$field_pos, off_det = od, off_pos = op,
      any_det = vw$any_det, n_valid = as.integer(vw$n_valid),
      positive = 0L, eta_bound = tulpaObs:::.TOBS_ETA_BOUND, freeman = TRUE)
  }
  expect_identical(ppc_kernel(zero, zero), ppc_kernel(none, none))

  lim <- function(od) tulpaObs:::cpp_occu_cover_cdf_limits(
    X_occ = fit$model$X_occ, X_det_site = fit$model$X_det_site,
    X_det_visit = vd(vw$X_det_visit), site_of_visit = vw$site_of_visit,
    b_occ = c0$b_occ, b_det = c0$b_det, field_occ = c0$field_occ, off_det = od,
    any_det = vw$any_det, eta_bound = tulpaObs:::.TOBS_ETA_BOUND, n_threads = 1L)
  expect_identical(lim(zero), lim(none))

  # The whole criteria stack still runs on a fit with no structured term.
  expect_true(is.finite(waic(fit, n.draws = 200L)$waic))
  expect_true(is.finite(cpo(fit, n.draws = 200L)$elpd_loo))
  set.seed(7); expect_true(is.finite(ppc(fit, n.samples = 100L)$bayesian.p))
})


test_that("the PPC and the PIT see the random-effect offsets too", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocic_sim(11L, arm = "p")
  fit <- .ocic_fit(sim, detection = ~ det_cov1 + (1 | habitat))
  c0  <- tulpaObs:::.tobs_occu_cover_components(fit, 150L)
  vw  <- tulpaObs:::.occu_cover_visit_view(fit$model)
  vd  <- function(X) tulpaObs:::.occu_cover_visit_design(X, vw$V)
  none <- matrix(0, vw$V, 0L)

  ppc_kernel <- function(od) {
    set.seed(55)
    tulpaObs:::cpp_occu_cover_ppc(
      X_occ = fit$model$X_occ, X_det_site = fit$model$X_det_site,
      X_pos_site = fit$model$X_pos_site,
      X_det_visit = vd(vw$X_det_visit), X_pos_visit = vd(vw$X_pos_visit),
      site_of_visit = vw$site_of_visit, y_det_visit = vw$y_det_visit,
      y_pos_visit = vw$y_pos_visit, b_occ = c0$b_occ, b_det = c0$b_det,
      b_pos = c0$b_pos, disp = c0$disp, field_occ = c0$field_occ,
      field_pos = c0$field_pos, off_det = od, off_pos = none,
      any_det = vw$any_det, n_valid = as.integer(vw$n_valid),
      positive = 0L, eta_bound = tulpaObs:::.TOBS_ETA_BOUND, freeman = TRUE)
  }
  expect_false(isTRUE(all.equal(ppc_kernel(c0$off_det)$fit.y, ppc_kernel(none)$fit.y)))

  lim <- function(od) tulpaObs:::cpp_occu_cover_cdf_limits(
    X_occ = fit$model$X_occ, X_det_site = fit$model$X_det_site,
    X_det_visit = vd(vw$X_det_visit), site_of_visit = vw$site_of_visit,
    b_occ = c0$b_occ, b_det = c0$b_det, field_occ = c0$field_occ, off_det = od,
    any_det = vw$any_det, eta_bound = tulpaObs:::.TOBS_ETA_BOUND, n_threads = 1L)
  # A per-site detection summary is more probable once the group's detection
  # offset is folded in, so the CDF limits move.
  expect_false(isTRUE(all.equal(lim(c0$off_det)$cdf_upper, lim(none)$cdf_upper)))

  set.seed(8); r <- pit_residuals(fit, n.samples = 150L)
  expect_length(r, fit$model$n_sites)
  expect_true(all(is.finite(r)))
})


test_that("a sampled coupled field is scored per site", {
  skip_on_cran()
  skip_if_fast()
  side <- 6L; J <- 4L; N <- side * side
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (cc in seq_len(side)) {
    i <- idx(r, cc)
    if (r > 1L)   adj[i, idx(r - 1L, cc)] <- 1L
    if (r < side) adj[i, idx(r + 1L, cc)] <- 1L
    if (cc > 1L)  adj[i, idx(r, cc - 1L)] <- 1L
    if (cc < side) adj[i, idx(r, cc + 1L)] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal", beta_occ = c(stats::qlogis(0.5), 0.8),
    beta_p = c(0.3, 0.5), beta_pos = c(log(0.12), -0.4), sigma_pos = 0.4,
    adj = adj, sigma = 0.7, alpha = 1.0, seed = 11L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)), det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  # copy(spatial()) is what puts the field on the cover arm as well, which is
  # what the simulated alpha = 1 does and what keeps the field_pos assertion
  # below non-vacuous: without it the amplitude is 0 and both sides of that
  # comparison are zero matrices.
  fit <- tobs(formula = ~ occ_cov1 + car_proper(graph = adj),
              data = cbind(data.frame(site_id = seq_len(N)), sim$data),
              family = occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1 + copy(spatial()), y = od$y, y_pos = y_pos,
              visits = od$det.covs, method = "nuts",
              control = list(verbose = FALSE, progress = FALSE, n.iter = 300L,
                             n.warmup = 300L, n.chains = 1L, seed = 4L))

  S  <- 300L
  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, S)
  expect_identical(dim(c0$field_occ), c(N, S))
  expect_false(all(c0$field_occ == 0))
  # Per site, the field draws are the sampled per-cell field at that site's cell,
  # and the cover arm carries the same field scaled by that draw's alpha.
  sc <- fit$model$site_cell %||% seq_len(N)
  expect_equal(rowMeans(c0$field_occ), unname(fit$spatial_field[sc]),
               tolerance = 1e-12)
  expect_equal(c0$field_pos,
               sweep(c0$field_occ, 2L, fit$hyper_draws[seq_len(S), "alpha"], "*"),
               tolerance = 1e-12)
  expect_false(all(c0$field_pos == 0))

  zero <- matrix(0, N, S)
  core <- tulpaObs:::.occu_cover_ploglik_core
  on  <- .ocic_score(core(fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp,
                          c0$field_occ, c0$field_pos))
  off <- .ocic_score(core(fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp,
                          zero, zero))
  expect_gt(on[["elpd"]] - off[["elpd"]], 2)
  expect_equal(waic(fit, n.draws = S)$elpd_waic, unname(on[["elpd"]]),
               tolerance = 1e-8)
  expect_true(is.finite(cpo(fit, n.draws = S)$elpd_loo))
  set.seed(9)
  expect_true(is.finite(ppc(fit, n.samples = 150L)$bayesian.p))
})
