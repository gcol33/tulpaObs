# =============================================================================
# test-cover-phi-scale.R
# - the joint nested-Laplace path reports the cover dispersion on cover()'s own
#   surface, against a SIMULATED truth (#309, #310).
#
# The arm's `phi` slot is the ENGINE convention, which for the gaussian arm that
# `lognormal` and `gaussian` both compile to is the residual VARIANCE, while
# every consumer -- `fit$means`, the per-draw `$disp` the diagnostics read, and
# `.occu_mscale_cover_sigma_pos()` feeding the lognormal conditional mean --
# takes it for the SD. #309 was that conversion missing on both ends of both
# joint families, at every arm coefficient's estimate and every recovery
# assertion in the suite unchanged, because the mean model does not see the
# dispersion scale.
#
# WHAT IS ASSERTED, and why not the estimate against its truth. A single small
# fixture puts the estimate anywhere from 0.33 to 0.49 for a truth of 0.40, so a
# band tight enough to exclude the squared value (0.16) is tight enough to flake.
# Two reads that do NOT depend on the per-seed level are used instead:
#
#   - the RATIO across two truths at one seed. Everything but the dispersion is
#     shared, so the level cancels: measured 1.757 / 1.775 / 1.806 against
#     0.70 / 0.40 = 1.75, where the variance scale would give 3.06.
#   - the estimate is nearer its truth than its truth SQUARED, which needs no
#     tolerance at all and is exactly the claim.
#
# The beta arm is the control: its phi is a precision on BOTH levels and must
# pass through unconverted, so applying the gaussian conversion to it would be
# the mirror defect. Its ratio is measured 3.888 / 4.331 / 4.894 against a true
# 20 / 5 = 4, where a square root would give 2.
#
# Cheap on purpose: every fit here is a 20-cell / 3-visit joint fit taking well
# under a second, so this runs at the smoke tier rather than only in the tier
# that has completed once.
# =============================================================================

.phi_scale_sig_grid <- exp(seq(log(0.05), log(1.5),  length.out = 9L))
.phi_scale_ln_grid  <- exp(seq(log(0.08), log(1.4),  length.out = 13L))
.phi_scale_be_grid  <- exp(seq(log(1.5),  log(60),   length.out = 13L))

# One coupled joint fit whose truth carries essentially no field (`sigma`
# 0.02), so the cover-arm spread IS the dispersion and nothing can absorb it.
.phi_scale_fit <- function(positive, truth, seed, N = 20L, J = 3L) {
  adj <- chain_adj(N)
  args <- list(N = N, J = J, positive = positive, adj = adj,
               beta_occ = c(0.2, 0.6), beta_p = c(0.4, -0.5),
               sigma = 0.02, alpha = 1.0, seed = seed)
  args$beta_pos <- if (positive == "beta") c(-0.6, 0.3) else c(log(0.25), 0.3)
  if (positive == "beta") args$phi <- truth else args$sigma_pos <- truth
  sim <- do.call(simulate_occu_cover, args)

  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = adj),
    data = cbind(data.frame(site_id = seq_len(N)), sim$data),
    family = occu_cover(positive), detection = ~ det_cov1,
    positive = ~ pos_cov1 + share(spatial()),
    y = od$y, y_pos = y_pos, visits = od$det.covs, method = "nested_laplace",
    control = list(verbose = FALSE, engine = "joint", progress = FALSE,
                   sigma.grid = .phi_scale_sig_grid,
                   phi.grid.pos = if (positive == "beta") .phi_scale_be_grid
                                  else .phi_scale_ln_grid)))
}

.phi_of <- function(fit) unname(fit$means[["phi_pos"]])

# Nearer its truth than its truth on the wrong scale. `wrong` is the value the
# unconverted axis would report: the square for a variance-convention arm, the
# square root for a precision one.
expect_on_scale <- function(est, truth, wrong, label) {
  expect_lt(abs(est - truth), abs(est - wrong),
            label = sprintf("%s: %.4f is nearer %.4f (wrong scale) than %.4f",
                            label, est, wrong, truth))
}


test_that("occu_cover() joint reports the lognormal cover SD, not its square", {
  lo <- .phi_of(.phi_scale_fit("lognormal", 0.40, 101L))
  hi <- .phi_of(.phi_scale_fit("lognormal", 0.70, 101L))

  expect_on_scale(lo, 0.40, 0.40^2, "lognormal 0.40")
  expect_on_scale(hi, 0.70, 0.70^2, "lognormal 0.70")

  # 1.75 on the SD scale, 3.06 on the variance scale.
  expect_gt(hi / lo, 1.45)
  expect_lt(hi / lo, 2.30)
})


test_that("occu_cover() joint reports the gaussian cover SD, not its square", {
  lo <- .phi_of(.phi_scale_fit("gaussian", 0.40, 102L))
  hi <- .phi_of(.phi_scale_fit("gaussian", 0.70, 102L))

  expect_on_scale(lo, 0.40, 0.40^2, "gaussian 0.40")
  expect_on_scale(hi, 0.70, 0.70^2, "gaussian 0.70")
  expect_gt(hi / lo, 1.45)
  expect_lt(hi / lo, 2.30)
})


test_that("occu_cover() joint leaves the beta precision unconverted", {
  lo <- .phi_of(.phi_scale_fit("beta", 5, 101L))
  hi <- .phi_of(.phi_scale_fit("beta", 20, 101L))

  # The mirror defect: a precision is an SD on both levels, so the gaussian
  # conversion applied here would report its square root.
  expect_on_scale(lo, 5, sqrt(5), "beta phi = 5")
  expect_on_scale(hi, 20, sqrt(20), "beta phi = 20")

  # 4 unconverted, 2 if square-rooted.
  expect_gt(hi / lo, 2.80)
  expect_lt(hi / lo, 6.50)
})


test_that("occu_cover() joint agrees between its reported and per-draw dispersion", {
  # The tree passed through a state where the per-draw axis was converted and
  # the reported one was not, and nothing said so: `$disp` is what the SBC
  # simulator, the WAIC/LOO/CPO pointwise log-likelihood, the PPC and
  # `predict()` read, `fit$means` is what a user reads, and a fit is only
  # coherent when they are the same quantity.
  fit <- .phi_scale_fit("lognormal", 0.40, 103L)
  drw <- mean(tulpaObs:::.tobs_joint_draws(fit, n = 2000L)$disp)
  expect_equal(drw, .phi_of(fit), tolerance = 0.05)
})


# -----------------------------------------------------------------------------
# occu_multiscale_cover(): the same axis, one consumer further on.
#
# `.occu_mscale_cover_sigma_pos()` reads the reported entry to build the
# lognormal conditional mean, so the dispersion scale reaches `fitted()` and
# `predict()` here rather than stopping at what is printed -- exp(sigma^2 / 2)
# evaluated at the square of the dispersion was about 9% on every predicted
# cover value.
# -----------------------------------------------------------------------------

.phi_scale_ms_fit <- function(truth, seed, n_cells = 20L) {
  sim <- simulate_occu_multiscale_cover(n_cells = n_cells, plots_per_cell = 3L,
                                        visits_per_plot = 2L, phi = truth,
                                        sigma = 0.02, seed = seed)
  suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot,
    positive = ~ x_cov + share(spatial(), alpha = grid(c(0, 0.5, 1, 2))),
    y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE,
                   sigma.grid = .phi_scale_sig_grid,
                   phi.grid.pos = .phi_scale_ln_grid)))
}


test_that("occu_multiscale_cover() joint reports the cover SD, not its square", {
  lo <- .phi_of(.phi_scale_ms_fit(0.40, 101L))
  hi <- .phi_of(.phi_scale_ms_fit(0.70, 101L))

  expect_on_scale(lo, 0.40, 0.40^2, "multiscale 0.40")
  expect_on_scale(hi, 0.70, 0.70^2, "multiscale 0.70")

  # Measured 1.713 / 1.754 / 1.756 across seeds, against 3.06 on the variance
  # scale.
  expect_gt(hi / lo, 1.45)
  expect_lt(hi / lo, 2.30)
})


test_that("occu_multiscale_cover() agrees across all three of its dispersion reads", {
  # The reported entry, the helper that turns it into a conditional mean, and
  # the per-draw axis the diagnostics read are one quantity or the fit is not
  # coherent. They disagreed in the tree once, between two commits of #309, and
  # nothing failed: the reported one was raw while the per-draw one was
  # converted, so `predict()` and `waic()` described different dispersions.
  fit <- .phi_scale_ms_fit(0.40, 103L)
  reported <- .phi_of(fit)
  helper   <- tulpaObs:::.occu_mscale_cover_sigma_pos(fit$means)
  drawn    <- mean(tulpaObs:::.tobs_joint_draws(fit, n = 2000L)$disp)

  expect_identical(helper, reported)
  expect_equal(drawn, reported, tolerance = 0.05)
})
