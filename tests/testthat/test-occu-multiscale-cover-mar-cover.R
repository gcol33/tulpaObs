# =============================================================================
# test-occu-multiscale-cover-mar-cover.R - missing-at-random cover on
# occu_multiscale_cover() (tulpaObs #262).
#
# A detected visit (y == 1) may carry a missing cover (y_pos = NA). The cell,
# plot and detection arms keep the visit; only that visit's cover factor drops
# out -- the rule the two-level twin occu_cover() applies, and the one the
# family's own pointwise-loglik kernel already applied while its three fit
# kernels did not. Coverage:
#   (1) the cell-coupling density drops the missing-cover visit and zeroes its
#       cover-arm derivatives, leaving every other block byte-identical
#       (beta + lognormal);
#   (2) the NUTS target does the same and stays byte-exact against the R oracle
#       -- which is the Laplace path's own marginal, so one check covers both;
#   (3) the builder admits NA at a detected visit, carries it as the sentinel
#       and keeps the detection; a build with no missing cover gains no NA;
#   (4) a fit runs with NA cover present. The cover factor is ADDITIVE in this
#       marginal (a detected plot factorises), so dropping cover leaves the
#       occupancy / availability / detection estimates where the full-data fit
#       puts them, and moves the cover arm alone.
# =============================================================================

# Self-contained cover log-densities, so the mechanism test needs no engine
# internals.
.mmc_beta_ld <- function(y, eta, phi) {
  mu <- plogis(eta); a <- mu * phi; b <- (1 - mu) * phi
  lgamma(phi) - lgamma(a) - lgamma(b) + (a - 1) * log(y) + (b - 1) * log(1 - y)
}
.mmc_lnorm_ld <- function(y, eta, sigma) {
  -log(y) - log(sigma) - 0.5 * log(2 * pi) - 0.5 * ((log(y) - eta) / sigma)^2
}

# One cell, two plots of two visits each: plot 1 detected at both visits (so it
# carries cover), plot 2 undetected (the within-plot availability mixture). The
# cell has a detection, so this is the any-detection branch.
.mmc_cell <- list(
  eta_psi   = 0.2,
  eta_theta = c(0.3, -0.1),
  eta_p     = c(-0.3, 0.4, 0.1, -0.2),
  eta_pos   = c(0.1, -0.2, 0.0, 0.0),
  y_det     = c(1L, 1L, 0L, 0L),
  sizes     = c(2L, 2L)
)


# --- (1) the cell-coupling density ----------------------------------------

test_that("multiscale coupling drops a detected visit with missing cover (beta)", {
  cc  <- .mmc_cell
  phi <- 20
  ref <- cpp_eval_occu_multiscale_cover_beta_cell(
    cc$eta_psi, cc$eta_theta, cc$eta_p, cc$eta_pos, cc$y_det,
    c(0.30, 0.55, 0, 0), cc$sizes, phi, "observed")
  mar <- cpp_eval_occu_multiscale_cover_beta_cell(
    cc$eta_psi, cc$eta_theta, cc$eta_p, cc$eta_pos, cc$y_det,
    c(0.30, NA,   0, 0), cc$sizes, phi, "observed")

  # The only change is visit 2's cover factor dropping out.
  f2 <- .mmc_beta_ld(0.55, cc$eta_pos[2], phi)
  expect_equal(mar$cell_ll, ref$cell_ll - f2, tolerance = 1e-10)

  # Visit 2's cover grad / Hessian are zeroed; every other block is untouched.
  expect_equal(mar$grad_pos[2], 0)
  expect_equal(mar$neg_hess_pos[2], 0)
  expect_equal(mar$grad_pos[1],  ref$grad_pos[1],  tolerance = 1e-12)
  expect_equal(mar$grad_p,       ref$grad_p,       tolerance = 1e-12)
  expect_equal(mar$grad_theta,   ref$grad_theta,   tolerance = 1e-12)
  expect_equal(mar$grad_psi,     ref$grad_psi,     tolerance = 1e-12)
  expect_equal(mar$cross_theta_p, ref$cross_theta_p, tolerance = 1e-12)
  expect_equal(mar$cross_p_p,     ref$cross_p_p,     tolerance = 1e-12)
})


test_that("multiscale coupling drops a detected visit with missing cover (lognormal)", {
  cc    <- .mmc_cell
  sigma <- 0.35
  eta_pos <- c(1.1, 0.9, 0, 0)
  ref <- cpp_eval_occu_multiscale_cover_lognormal_cell(
    cc$eta_psi, cc$eta_theta, cc$eta_p, eta_pos, cc$y_det,
    c(2.1, 3.0, 0, 0), cc$sizes, sigma, "observed")
  mar <- cpp_eval_occu_multiscale_cover_lognormal_cell(
    cc$eta_psi, cc$eta_theta, cc$eta_p, eta_pos, cc$y_det,
    c(2.1, NA,  0, 0), cc$sizes, sigma, "observed")

  f2 <- .mmc_lnorm_ld(3.0, eta_pos[2], sigma)
  expect_equal(mar$cell_ll, ref$cell_ll - f2, tolerance = 1e-10)
  expect_equal(mar$grad_pos[2], 0)
  expect_equal(mar$neg_hess_pos[2], 0)
  expect_equal(mar$grad_pos[1], ref$grad_pos[1], tolerance = 1e-12)
  expect_equal(mar$grad_p,     ref$grad_p,     tolerance = 1e-12)
  expect_equal(mar$grad_theta, ref$grad_theta, tolerance = 1e-12)
  expect_equal(mar$grad_psi,   ref$grad_psi,   tolerance = 1e-12)

  # Every cover missing on a detected plot: the plot keeps its detection terms
  # and contributes no cover factor at all.
  none <- cpp_eval_occu_multiscale_cover_lognormal_cell(
    cc$eta_psi, cc$eta_theta, cc$eta_p, eta_pos, cc$y_det,
    c(NA, NA, 0, 0), cc$sizes, sigma, "observed")
  f1 <- .mmc_lnorm_ld(2.1, eta_pos[1], sigma)
  expect_equal(none$cell_ll, ref$cell_ll - f1 - f2, tolerance = 1e-10)
  expect_equal(none$grad_pos, rep(0, length(cc$eta_p)))
  expect_equal(none$grad_p,   ref$grad_p, tolerance = 1e-12)
})


# --- (2) the NUTS target, against the shared R marginal --------------------

test_that("multiscale NUTS target drops missing cover and matches the R oracle", {
  for (positive in c("lognormal", "beta")) {
    set.seed(if (positive == "beta") 22L else 11L)
    sim <- simulate_occu_multiscale_cover(
      n_cells = 25L, plots_per_cell = 4L, visits_per_plot = 3L,
      beta_pos = if (positive == "beta") c(stats::qlogis(0.3), -0.3)
                 else c(log(0.10), -0.4),
      positive = positive, phi = if (positive == "beta") 12 else 0.35,
      sigma = 0, alpha = 0, seed = if (positive == "beta") 22L else 11L)

    # Drop a third of the recorded covers at random.
    det <- which(!is.na(sim$y) & sim$y == 1L)
    expect_gt(length(det), 20L)
    set.seed(99L)
    drop <- sample(det, floor(length(det) / 3))
    y_pos_mar <- sim$y_pos
    y_pos_mar[drop] <- NA_real_

    model <- .omc_bind_model(sim, positive, y_pos = y_pos_mar)
    lay   <- tulpaObs:::.tobs_occu_mscale_cover_nuts_layout(model)
    spec  <- tulpaObs:::.tobs_occu_mscale_cover_nuts_spec(model)

    set.seed(3)
    theta <- stats::rnorm(lay$total) * 0.3
    theta[lay$disp] <- log(if (positive == "beta") 12 else 0.35)
    sb <- 5

    # The oracle runs .occu_mscale_cover_nonspatial_ll, the Laplace path's own
    # marginal, so agreement here pins all three kernels to one gate.
    lp_R <- tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(theta, model, lay, sb)
    cpp  <- cpp_occu_mscale_cover_nuts_joint_logpost(spec, theta, sb, sb)
    expect_true(is.finite(lp_R))
    expect_equal(cpp$lp, lp_R, tolerance = 1e-9)

    num <- vapply(seq_along(theta), function(j) {
      tp <- theta; tp[j] <- tp[j] + 1e-6
      tm <- theta; tm[j] <- tm[j] - 1e-6
      (tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(tp, model, lay, sb) -
       tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(tm, model, lay, sb)) / 2e-6
    }, 0)
    expect_lt(max(abs(cpp$grad - num)), 1e-5)

    # The dropped covers, and only those, leave the target: the gap to the
    # full-data target is exactly the sum of their densities.
    full <- .omc_bind_model(sim, positive)
    lp_full <- tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(theta, full, lay, sb)
    eta_pos <- tulpaObs:::.occu_ms_eta_visit(
      model$X_pos_site, theta[lay$pos_site],
      model$X_pos_visit, theta[lay$pos_visit], nrow(sim$y), ncol(sim$y))
    disp <- exp(theta[lay$disp])
    gap  <- if (positive == "beta")
      sum(.mmc_beta_ld(sim$y_pos[drop], eta_pos[drop], disp))
    else
      sum(.mmc_lnorm_ld(sim$y_pos[drop], eta_pos[drop], disp))
    expect_equal(lp_full - lp_R, gap, tolerance = 1e-8)
  }
})


# --- (3) the builder -------------------------------------------------------

test_that("multiscale builder carries missing cover as the sentinel", {
  sim <- simulate_occu_multiscale_cover(
    n_cells = 20L, plots_per_cell = 3L, visits_per_plot = 2L,
    positive = "lognormal", sigma = 0, alpha = 0, seed = 7L)

  full <- .omc_bind_model(sim, "lognormal")
  expect_false(anyNA(full$y_pos))

  det <- which(!is.na(sim$y) & sim$y == 1L)
  y_pos_mar <- sim$y_pos
  y_pos_mar[det[1L]] <- NA_real_

  model <- expect_silent(.omc_bind_model(sim, "lognormal", y_pos = y_pos_mar))
  # The visit keeps its detection; only the cover is the sentinel.
  expect_true(is.na(model$y_pos[det[1L]]))
  expect_identical(model$y[det[1L]], 1L)
  expect_true(model$valid[det[1L]])
  expect_equal(sum(is.na(model$y_pos)), 1L)
  # Undetected visits stay zero-filled, as before.
  expect_true(all(model$y_pos[!(model$valid & model$y == 1L)] == 0))

  # A cover outside the arm's support is still refused where it IS recorded.
  y_bad <- sim$y_pos
  y_bad[det[1L]] <- -1
  expect_error(.omc_bind_model(sim, "lognormal", y_pos = y_bad),
               "y_pos > 0")
})


# --- (4) a fit with missing cover present ----------------------------------

test_that("multiscale fit with missing cover moves the cover arm alone", {
  skip_on_cran()

  sim <- simulate_occu_multiscale_cover(
    n_cells = 60L, plots_per_cell = 4L, visits_per_plot = 3L,
    beta_pos = c(log(0.10), -0.4), positive = "lognormal", phi = 0.35,
    sigma = 0, alpha = 0, seed = 31L)

  det <- which(!is.na(sim$y) & sim$y == 1L)
  set.seed(31L)
  y_pos_mar <- sim$y_pos
  y_pos_mar[sample(det, floor(length(det) / 3))] <- NA_real_

  fit_one <- function(y_pos) suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = y_pos, method = "laplace",
    control = list(verbose = FALSE)))

  full <- fit_one(sim$y_pos)
  mar  <- fit_one(y_pos_mar)

  expect_s3_class(mar, "tobs_fit")
  expect_true(all(is.finite(mar$means)))
  expect_true(all(is.finite(mar$sds)))

  # A detected plot's density factorises, so the cover factor is additive in
  # this marginal and the other three arms optimise the same objective either
  # way: their estimates and SEs are the full-data ones up to the optimiser.
  keep <- grep("^(psi|theta|p)_", names(full$means), value = TRUE)
  expect_gte(length(keep), 6L)
  expect_equal(mar$means[keep], full$means[keep], tolerance = 1e-3)
  expect_equal(mar$sds[keep],   full$sds[keep],   tolerance = 1e-3)

  # The cover arm keeps recovering its truth off the covers that remain, with a
  # wider interval for the ones that left.
  expect_lt(abs(unname(mar$means["pos_x_cov"]) - (-0.4)), 0.25)
  expect_gt(unname(mar$sds["pos_x_cov"]), unname(full$sds["pos_x_cov"]))

  # The scoring doors read the same gate the fit did: the family's pointwise
  # loglik already dropped a missing cover, so nothing here goes non-finite.
  ll <- tulpaObs:::.tobs_ploglik_occu_multiscale_cover(mar, n.draws = 50L)
  expect_equal(sum(!is.finite(ll)), 0L)
  expect_true(is.finite(suppressWarnings(waic(mar))$elpd_waic))
  expect_true(is.finite(
    suppressWarnings(loo(mar))$estimates["elpd_loo", "Estimate"]))
  expect_false(anyNA(unlist(fitted(mar))))
})


# --- (5) the shared-field engine, which is the family's deliverable --------

test_that("multiscale shared-field fit takes missing cover", {
  skip_on_cran()
  skip_if_fast()

  sim <- simulate_occu_multiscale_cover(
    n_cells = 40L, plots_per_cell = 4L, visits_per_plot = 3L,
    beta_pos = c(log(0.10), -0.4), positive = "lognormal", phi = 0.35,
    sigma = 0.5, alpha = 1, seed = 31L)

  det <- which(!is.na(sim$y) & sim$y == 1L)
  set.seed(31L)
  y_pos_mar <- sim$y_pos
  y_pos_mar[sample(det, floor(length(det) / 3))] <- NA_real_

  fit_one <- function(y_pos) suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = y_pos, method = "nested_laplace",
    control = list(sigma.grid = c(0.3, 0.6, 1.0), alpha.grid = c(0, 0.5, 1, 2),
                   diagnose.k = FALSE, max.iter = 500L, verbose = FALSE)))

  full <- fit_one(sim$y_pos)
  mar  <- fit_one(y_pos_mar)

  expect_s3_class(mar, "tobs_fit")
  expect_true(all(is.finite(mar$means)))
  expect_true(all(is.finite(mar$sds)))

  # The occupancy field is shared across the arms here, so the other arms are
  # no longer exactly the full-data ones; they stay close, and the field the
  # deliverable reports is the same surface.
  keep <- grep("^(psi|theta|p)_", names(full$means), value = TRUE)
  expect_lt(max(abs(mar$means[keep] - full$means[keep])), 0.1)
  expect_gt(stats::cor(full$spatial_field, mar$spatial_field), 0.95)

  ll <- tulpaObs:::.tobs_ploglik_occu_multiscale_cover(mar, n.draws = 50L)
  expect_equal(sum(!is.finite(ll)), 0L)
})
