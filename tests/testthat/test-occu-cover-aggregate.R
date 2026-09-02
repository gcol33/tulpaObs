# =============================================================================
# test-occu-cover-aggregate.R
# - cell-aggregated cover.
#
# occu_cover(cover_aggregate = "mean"/"median") collapses the cover arm to one
# observation per occupancy unit (the mean / median cover over that unit's
# detected visits) so the cover arm contributes to the shared field at the cell
# scale rather than the per-visit scale. Covers: aggregation resolution +
# fall-back + error gates (fast), the arm shrinks to one row per detected cell,
# and multi-seed recovery of the field / occupancy / cover coefficients.
# =============================================================================

# group_var occu_cover with a CELL-level positive covariate + shared field.
.agg_sim <- function(seed, n_cells = 30L, n_per = 5L, J = 12L,
                     b_pos = c(stats::qlogis(0.3), 0.6), b_occ1 = 0.7,
                     sigma = 0.8, alpha = 1.2, phi = 30) {
  set.seed(seed)
  adj <- chain_adj(n_cells); n_sites <- n_cells * n_per
  Q  <- tulpaObs:::.occu_cover_icar_Q(adj)
  sq <- tulpaObs:::.occu_cover_icar_scale(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  zk <- stats::rnorm(sum(keep))
  f  <- as.numeric(eig$vectors[, keep, drop = FALSE] %*% (zk / sqrt(eig$values[keep])))
  f  <- (f - mean(f)) / sqrt(sq)
  cell_idx <- rep(seq_len(n_cells), times = n_per)
  xcell <- as.numeric(scale(stats::rnorm(n_cells)))
  xocc  <- as.numeric(scale(stats::rnorm(n_cells)))
  site <- data.frame(cell_idx = cell_idx, xpos = xcell[cell_idx], xocc = xocc[cell_idx])
  eta_psi <- 0.3 + b_occ1 * xocc[cell_idx] + sigma * f[cell_idx]
  z <- stats::rbinom(n_sites, 1L, stats::plogis(eta_psi))
  det_cov <- matrix(stats::rnorm(n_sites * J), n_sites, J)
  Y <- matrix(0L, n_sites, J); Ypos <- matrix(0, n_sites, J)
  for (i in seq_len(n_sites)) for (j in seq_len(J)) if (z[i] == 1L) {
    d <- stats::rbinom(1L, 1L, stats::plogis(stats::qlogis(0.7) + 0.4 * det_cov[i, j]))
    Y[i, j] <- d
    if (d == 1L) {
      mu <- stats::plogis(b_pos[1] + b_pos[2] * xcell[cell_idx[i]] +
                          alpha * sigma * f[cell_idx[i]])
      Ypos[i, j] <- stats::rbeta(1L, mu * phi, (1 - mu) * phi)
    }
  }
  Ypos[Ypos <= 0] <- 0; Ypos[Ypos >= 1] <- 1 - 1e-6
  list(site = site, Y = Y, Ypos = Ypos,
       vd = data.frame(det_cov = as.vector(t(det_cov))),
       adj = adj, f = f, b_pos = b_pos, b_occ1 = b_occ1, n_cells = n_cells)
}

.agg_fit <- function(sim, cover_aggregate = "mean", max.iter = 300L) {
  suppressWarnings(tobs(
    formula = ~ xocc + icar(graph = sim$adj, group_var = "cell_idx"),
    data = sim$site,
    family = occu_cover("beta", cover_aggregate = cover_aggregate),
    detection = ~ det_cov,
    positive = ~ xpos + share(spatial(), alpha = grid(c(0, 0.8, 1.5))),
    y = sim$Y, y_pos = sim$Ypos, visits = sim$vd,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = max.iter, engine = "joint",
                   sigma.grid = exp(seq(log(0.4), log(1.6), length.out = 4)),
                   adaptive.grid = FALSE,
                   diagnose.k = FALSE)))
}


test_that("aggregation resolution, fall-back, and error gates", {
  # The family carries the cover_aggregate choice (NULL until set).
  expect_null(occu_cover("beta")$params$cover_aggregate)
  expect_identical(occu_cover("beta", cover_aggregate = "mean")$params$cover_aggregate,
                   "mean")
  expect_identical(occu_cover("beta", cover_aggregate = "median")$params$cover_aggregate,
                   "median")
  expect_error(occu_cover("beta", cover_aggregate = "nonsense"), "should be one of")

  sim <- .agg_sim(seed = 11L, n_cells = 12L, n_per = 3L, J = 6L)

  # Default on the spatial path with a CELL-level positive design -> "mean".
  fit_def <- .agg_fit(structure(sim), cover_aggregate = NULL)
  expect_identical(fit_def$model$cover_aggregate, "mean")

  # A VISIT-level positive covariate with the default quietly falls back to
  # per-visit cover (no opt-out needed).
  sim2 <- sim
  sim2$vd <- data.frame(det_cov = as.vector(t(matrix(rnorm(nrow(sim$site)*6L),
                                                     nrow(sim$site), 6L))),
                        pcov = as.vector(t(matrix(rnorm(nrow(sim$site)*6L),
                                                  nrow(sim$site), 6L))))
  fit_fallback <- suppressWarnings(tobs(
    formula = ~ icar(graph = sim$adj, group_var = "cell_idx"),
    data = sim$site, family = occu_cover("beta"),     # default
    detection = ~ det_cov,                            # visit-level pos covariate
    positive = ~ pcov + share(spatial(), alpha = grid(c(0, 1.0))),
    y = sim$Y, y_pos = sim$Ypos, visits = sim2$vd,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 200L, engine = "joint",
                   sigma.grid = exp(seq(log(0.5), log(1.5), length.out = 3)),
                   adaptive.grid = FALSE, diagnose.k = FALSE)))
  expect_identical(fit_fallback$model$cover_aggregate, "none")

  # Explicit aggregation + a visit-level positive covariate -> error.
  expect_error(
    tobs(formula = ~ icar(graph = sim$adj, group_var = "cell_idx"),
         data = sim$site, family = occu_cover("beta", cover_aggregate = "mean"),
         detection = ~ det_cov, positive = ~ pcov,
         y = sim$Y, y_pos = sim$Ypos, visits = sim2$vd,
         method = "nested_laplace",
         control = list(verbose = FALSE, engine = "joint")),
    "cell-level positive design")

  # Explicit aggregation on the non-spatial laplace path -> error.
  expect_error(
    tobs(formula = ~ xocc, data = sim$site,
         family = occu_cover("beta", cover_aggregate = "mean"),
         detection = ~ det_cov, positive = ~ xpos,
         y = sim$Y, y_pos = sim$Ypos, visits = sim$vd, method = "laplace"),
    "shared-field spatial path")

  # Explicit aggregation on a v2/v3 escape hatch -> error.
  expect_error(
    tobs(formula = ~ icar(graph = sim$adj, group_var = "cell_idx"),
         data = sim$site, family = occu_cover("beta", cover_aggregate = "mean"),
         detection = ~ det_cov, positive = ~ xpos,
         y = sim$Y, y_pos = sim$Ypos, visits = sim$vd, method = "nested_laplace",
         control = list(engine = "v3_nested")),
    "joint engine")
})


test_that("aggregated cover arm holds one row per detected occupancy unit", {
  sim <- .agg_sim(seed = 21L, n_cells = 20L, n_per = 4L, J = 8L)
  n_detected_sites <- sum(rowSums(sim$Y == 1L) > 0L)
  n_valid_visits   <- sum(!is.na(sim$Y))   # per-visit arm spans every valid visit
  expect_gt(n_valid_visits, n_detected_sites)   # genuinely many visits/site

  # Build the model the dispatcher would, with a cell-level positive design, then
  # build the joint-coupled arms in each aggregation mode and inspect the cover
  # (pos) arm directly -- no fit needed.
  model <- tulpaObs:::.tobs_build_occu_cover(
    occ_formula = ~ xocc, det_formula = ~ 1, pos_formula = ~ xpos,
    data = sim$site, y = sim$Y, y_pos = sim$Ypos, positive = "beta")
  site_cell <- as.integer(sim$site$cell_idx)

  arms_mean <- tulpaObs:::.occu_cover_build_joint_arms(
    model, sigma_pos_init = 10, positive = "beta",
    alpha_axis = tulpaObs:::.tobs_alpha_axis(c(0, 1)),
    multi = FALSE, n_cells = sim$n_cells, site_cell = site_cell,
    cover_aggregate = "mean")
  arms_none <- tulpaObs:::.occu_cover_build_joint_arms(
    model, sigma_pos_init = 10, positive = "beta",
    alpha_axis = tulpaObs:::.tobs_alpha_axis(c(0, 1)),
    multi = FALSE, n_cells = sim$n_cells, site_cell = site_cell,
    cover_aggregate = "none")

  # Cover arm: one row per detected site (aggregated) vs one per valid visit
  # (per-visit; the spec gates the cover term on detection internally).
  expect_equal(length(arms_mean$responses$pos$y), n_detected_sites)
  expect_equal(length(arms_none$responses$pos$y), n_valid_visits)
  expect_equal(arms_mean$n_pos_rows, n_detected_sites)
  # Each aggregated cover row maps to a distinct occupancy unit (site).
  expect_equal(anyDuplicated(arms_mean$responses$pos$cell_obs_map), 0L)
  # The aggregated value is the per-site mean cover over its detected visits.
  det_mat <- model$valid & (model$y == 1L)
  sw <- which(rowSums(det_mat) > 0L)
  expect_equal(arms_mean$responses$pos$y,
               vapply(sw, function(i) mean(model$y_pos[i, det_mat[i, ]]), numeric(1)))
  # Occupancy and detection arms are identical across modes.
  expect_identical(arms_mean$responses$psi$y, arms_none$responses$psi$y)
  expect_identical(arms_mean$responses$p$y,   arms_none$responses$p$y)
})


test_that("aggregated cover recovers field + occupancy + cover coefficients (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  pos1 <- occ1 <- fcor <- numeric(0)
  for (sd in seeds) {
    sim <- .agg_sim(seed = sd)
    fit <- .agg_fit(sim)
    pos1 <- c(pos1, unname(fit$means[["pos_xpos"]]))
    occ1 <- c(occ1, unname(fit$means[["psi_xocc"]]))
    fcor <- c(fcor, cor(fit$spatial_field, sim$f))
  }
  # Field shape reconstructs (the cover arm still informs it at the cell scale).
  expect_gt(median(fcor), 0.85)
  # Occupancy slope recovers (truth 0.7).
  expect_gt(median(occ1), 0.4)
  # Cover slope recovers sign and rough magnitude (truth 0.6). Aggregation trades
  # cover-arm precision for field balance, so it attenuates but stays clearly
  # positive across seeds.
  expect_gt(median(pos1), 0.25)
  expect_lt(median(pos1), 0.95)
})


test_that("median aggregation runs and tracks the mean on symmetric cover", {
  skip_on_cran()
  skip_if_fast()
  sim <- .agg_sim(seed = 717L)
  fit_mean   <- .agg_fit(sim, cover_aggregate = "mean")
  fit_median <- .agg_fit(sim, cover_aggregate = "median")
  expect_s3_class(fit_median, "tobs_fit")
  expect_identical(fit_median$model$cover_aggregate, "median")
  # Mean and median of a near-symmetric per-cell cover sit close; the cover
  # intercept should agree to within a modest tolerance.
  expect_equal(unname(fit_median$means[["pos_(Intercept)"]]),
               unname(fit_mean$means[["pos_(Intercept)"]]), tolerance = 0.4)
})


test_that("WAIC scores aggregated cover at the unit scale, not per visit (#34)", {
  skip_on_cran()
  skip_if_fast()
  # An aggregated fit collapses the cover arm to one observation per occupancy
  # unit, so its pointwise log-likelihood must score the cover term once per
  # unit. The pre-fix GOF summed the per-visit cover density at the fitted
  # (aggregated, tight) dispersion, so p_waic grew super-linearly in the
  # visits-per-site J. Post-fix p_waic stays on the unit scale and flat in J.
  pw <- vapply(c(12L, 60L), function(J) {
    sim <- .agg_sim(seed = 909L, J = J)
    fit <- .agg_fit(sim, cover_aggregate = "mean")
    waic(fit, n.draws = 200L)$p_waic
  }, numeric(1))
  n_sites <- 30L * 5L
  # Effective parameter count stays well below the site count (it scaled past it
  # before the fix) and does not balloon with J (a 5x J increase).
  expect_lt(pw[1], n_sites)
  expect_lt(pw[2], n_sites)
  expect_lt(pw[2], 3 * pw[1])

  # The pointwise log-likelihood reads the dispersion the spec held fixed (a beta
  # precision well away from 1), not the bare unit default of the pre-fix path.
  sim <- .agg_sim(seed = 909L)
  fit <- .agg_fit(sim, cover_aggregate = "mean")
  expect_gt(tulpaObs:::.tobs_joint_draws(fit, n = 20L)$disp[1], 2)
})


test_that("tobs_ppc scores aggregated cover at the unit scale, not per visit (#34)", {
  skip_on_cran()
  skip_if_fast()
  # The aggregated fit's dispersion is fit on the per-unit cover (tight, high phi).
  # The pre-fix PPC replicated that phi per detected visit, so the positive-part
  # discrepancy ran off scale and drove the Bayesian p-value to the boundary,
  # worse as the visits-per-site J grew. Post-fix the cover discrepancy is one
  # term per detected unit, so the p-value stays interior and does not collapse
  # with J.
  bp <- vapply(c(12L, 60L), function(J) {
    sim <- .agg_sim(seed = 707L, J = J)
    fit <- .agg_fit(sim, cover_aggregate = "mean")
    set.seed(1L)
    ppc(fit, n.samples = 200L)$bayesian.p
  }, numeric(1))
  expect_true(all(bp > 0.02 & bp < 0.98))
})
