# =============================================================================
# test-occu-cover-nl-ic.R - the grid-integrated occu_cover() information criteria
# score every random effect the fit carries (gcol33/tulpaObs#215).
#
# A `(1 | g)` on the detection / positive-cover formula rides the nested-Laplace
# joint engine as an iid prior block (gcol33/tulpaObs#102, #103), and one on the
# occurrence formula alongside the shared field does the same on the occupancy
# arm (gcol33/tulpaObs#56). .tobs_joint_draws() samples every such latent from
# the same outer-grid posterior as the coefficients, so they reach the criteria
# as offsets by the pathway the sampled route already uses
# (gcol33/tulpaObs#211): an observation-arm term enters the shared `Arms`
# predictor view per visit, an occupancy-arm term joins the per-site occupancy
# offset the shared field loads. The pointwise log-likelihood, the posterior
# predictive check and the PIT / LOO-PIT all read both.
# Tests:
#   - the components carry the per-visit offset on the arm that holds the term,
#     and its per-group posterior mean agrees with the fit's own BLUPs
#   - the criteria MOVE with the random effect on either observation arm, and a
#     zeroed offset reproduces the population-mean score to the bit
#   - a fit carrying no random effect leaves the offsets absent, so its score is
#     the arithmetic it had before
#   - an occupancy-arm random effect is scored: its per-site offset carries the
#     fit's own BLUPs, and removing it moves the criteria
# =============================================================================


.ocnl_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

.ocnl_sim <- function(seed, arm = "p", side = 6L, J = 6L, n_g = 6L,
                      sigma_re = 0.9) {
  adj <- .ocnl_grid_adj(side)
  simulate_occu_cover(
    N = nrow(adj), J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
    re_det_groups = if (identical(arm, "p")) n_g else NULL,
    sigma_re_p = sigma_re,
    re_pos_groups = if (identical(arm, "pos")) n_g else NULL,
    sigma_re_pos = sigma_re, seed = seed)
}

.ocnl_fit <- function(sim, adj, arm = "p") {
  det <- if (identical(arm, "p")) ~ det_cov1 + (1 | habitat) else ~ det_cov1
  pos <- if (identical(arm, "pos")) ~ pos_cov1 + (1 | habitat) else ~ pos_cov1
  tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
       family = occu_cover("lognormal"), detection = det, positive = pos,
       y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

# lppd, the WAIC penalty and elpd from an [S x n_sites] pointwise log-likelihood.
.ocnl_score <- function(ll) {
  lse  <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
  lppd <- sum(apply(ll, 2L, lse))
  pw   <- sum(apply(ll, 2L, stats::var))
  c(lppd = lppd, p_waic = pw, elpd = lppd - pw)
}

# The pointwise log-likelihood of a fit's own draws, with the per-visit offsets
# passed, suppressed (the population-mean score), or replaced by explicit zeros.
.ocnl_ll <- function(fit, c0, off = c("as_fitted", "absent", "zero")) {
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


test_that("the grid-integrated components carry the per-visit RE offset", {
  skip_on_cran()
  skip_if_fast()
  for (arm in c("p", "pos")) {
    sim <- .ocnl_sim(if (identical(arm, "p")) 1L else 3L, arm = arm)
    adj <- .ocnl_grid_adj(6L)
    fit <- .ocnl_fit(sim, adj, arm)

    S  <- 400L
    set.seed(11)
    c0 <- tulpaObs:::.tobs_occu_cover_components(fit, S)
    expect_identical(is.null(c0$off_det), !identical(arm, "p"))
    expect_identical(is.null(c0$off_pos), !identical(arm, "pos"))
    got <- if (identical(arm, "p")) c0$off_det else c0$off_pos
    vw  <- tulpaObs:::.occu_cover_visit_view(fit$model)
    expect_identical(dim(got), c(vw$V, S))

    # Every visit in a group carries that group's offset, so the per-visit
    # posterior means collapse to one value per level. The fit reports the same
    # quantity independently as its BLUPs, centred over groups; the draws are
    # not, so the two agree up to that common level.
    codes <- vw$re[[arm]][[1L]]$codes
    grp   <- tapply(rowMeans(got), codes, mean)
    expect_length(grp, 6L)
    blup  <- fit$re[[1L]]$blup
    expect_equal(as.numeric(grp) - mean(grp), blup, tolerance = 0.05)
  }
})


test_that("the criteria move with the detection RE and reduce to it at zero", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocnl_sim(1L, arm = "p")
  adj <- .ocnl_grid_adj(6L)
  fit <- .ocnl_fit(sim, adj, "p")
  expect_gt(fit$means[["sigma_re_p"]], 0.4)

  set.seed(11)
  c0  <- tulpaObs:::.tobs_occu_cover_components(fit, 400L)
  on  <- .ocnl_score(.ocnl_ll(fit, c0, "as_fitted"))
  off <- .ocnl_score(.ocnl_ll(fit, c0, "absent"))

  # A detection random intercept the grid put at sigma ~ 1.1 explains a large
  # share of the detection heterogeneity, so scoring it at the population mean
  # costs many nats: the criteria describe a different model with and without it.
  expect_gt(on[["lppd"]] - off[["lppd"]], 5)
  expect_gt(on[["elpd"]] - off[["elpd"]], 5)

  # Zeroing the offset is the "RE at zero variance" reduction: same arithmetic,
  # to the bit, as passing none at all.
  expect_identical(.ocnl_ll(fit, c0, "zero"), .ocnl_ll(fit, c0, "absent"))

  # The public criteria read the same components, so WAIC and LOO land on the
  # scored model rather than the population-mean one.
  w <- tobs_waic(fit, n.draws = 400L)
  expect_gt(w$elpd_waic, off[["elpd"]] + 3)
  cp <- tobs_cpo(fit, n.draws = 400L)
  expect_true(is.finite(cp$elpd_loo))
  expect_gt(cp$elpd_loo, off[["elpd"]] + 3)
})


test_that("the criteria move with the cover-arm RE", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocnl_sim(3L, arm = "pos", side = 7L, J = 8L)
  adj <- .ocnl_grid_adj(7L)
  fit <- .ocnl_fit(sim, adj, "pos")
  expect_gt(fit$means[["sigma_re_pos"]], 0.4)

  set.seed(11)
  c0  <- tulpaObs:::.tobs_occu_cover_components(fit, 400L)
  on  <- .ocnl_score(.ocnl_ll(fit, c0, "as_fitted"))
  off <- .ocnl_score(.ocnl_ll(fit, c0, "absent"))
  expect_gt(on[["lppd"]] - off[["lppd"]], 5)
  expect_gt(on[["elpd"]] - off[["elpd"]], 5)
  expect_identical(.ocnl_ll(fit, c0, "zero"), .ocnl_ll(fit, c0, "absent"))

  # The posterior predictive check and the PIT read the same offsets, so both
  # stay finite on a fit that carries one.
  pp <- tobs_ppc(fit, n.samples = 40L)
  expect_true(is.finite(pp$bayesian.p))
  pit <- tobs_pit_residuals(fit, n.samples = 40L)
  expect_true(all(is.finite(as.numeric(unlist(pit)))))
})


test_that("a fit with no random effect keeps the no-offset arithmetic", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocnl_sim(2L, arm = "none")
  adj <- .ocnl_grid_adj(6L)
  fit <- .ocnl_fit(sim, adj, "none")
  expect_null(fit$re)
  expect_null(fit$model$re_psi)

  set.seed(11)
  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, 200L)
  expect_null(c0$off_det)
  expect_null(c0$off_pos)
  # No occupancy-arm term either, so the occupancy predictor keeps the field
  # alone -- the offset builder returns the "none" signal the components read.
  expect_null(tulpaObs:::.occu_cover_re_site_offset(
    fit$model, tulpaObs:::.tobs_joint_draws(fit, n = 20L)$re,
    fit$model$n_sites, 20L))
  # No offset -> the kernels take the zero-column "arm carries none" signal, so
  # the score is the arithmetic of the no-offset path, bit for bit.
  expect_identical(.ocnl_ll(fit, c0, "zero"), .ocnl_ll(fit, c0, "absent"))
  set.seed(7); w1 <- tobs_waic(fit, n.draws = 200L)
  set.seed(7); w2 <- tobs_waic(fit, n.draws = 200L)
  expect_identical(w1$elpd_waic, w2$elpd_waic)
  expect_true(is.finite(w1$elpd_waic))
})


# An occupancy-arm random intercept alongside the shared field: a smoothed ICAR
# surface on psi (copied onto cover by alpha) plus a per-region deviation.
.ocnl_psi_sim <- function(seed, side = 7L, J = 5L, n_g = 7L, sigma_u = 1.6) {
  set.seed(seed)
  adj <- .ocnl_grid_adj(side); n <- nrow(adj)
  phi <- as.numeric(scale(stats::rnorm(n)))
  for (r in 1:4) { pn <- phi; for (i in seq_len(n)) {
    nb <- which(adj[i, ] == 1L); pn[i] <- 0.5*phi[i] + 0.5*mean(phi[nb]) }; phi <- pn }
  phi <- 0.6 * as.numeric(scale(phi)); phi <- phi - mean(phi)
  grp <- factor(paste0("reg", rep(seq_len(n_g), length.out = n)))
  u <- stats::rnorm(n_g, 0, sigma_u); u <- u - mean(u)
  x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, stats::plogis(0.7*x + phi + u[as.integer(grp)]))
  p <- stats::plogis(0.4)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) y[i, ] <- if (z[i] == 1L) stats::rbinom(J, 1, p) else 0L
  y_pos <- matrix(0, n, J)
  mu <- log(25) + 0.5 * phi
  for (i in seq_len(n)) for (j in seq_len(J))
    if (y[i, j] == 1L) y_pos[i, j] <- exp(stats::rnorm(1, mu[i], 0.4))
  list(adj = adj, data = data.frame(x = x, region = grp), y = y, y_pos = y_pos)
}


test_that("the criteria score an occupancy-arm random effect", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocnl_psi_sim(2L)
  fit <- tobs(occurrence = ~ x + icar(graph = sim$adj) + (1 | region),
              data = sim$data, family = occu_cover("lognormal"),
              detection = ~ 1, positive = ~ 1,
              y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$re$psi$arm, "psi")
  expect_gt(fit$means[["sigma_re"]], 0.4)
  # The grouping the fit ran on is carried on the model, one code per occupancy
  # unit, and the term reports its variable and levels like an observation arm.
  expect_identical(fit$re$psi$var, "region")
  expect_length(fit$re$psi$levels, 7L)
  expect_length(fit$model$re_psi$group_idx, fit$model$n_sites)

  S <- 400L
  set.seed(11)
  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, S)
  # The components sample the grid-integrated posterior first, so re-seeding
  # reproduces the SAME draw and isolates the offset that went into field_occ.
  set.seed(11)
  bd <- tulpaObs:::.tobs_joint_draws(fit, n = S)
  expect_identical(c0$b_occ, bd$b$occ)
  so <- tulpaObs:::.occu_cover_re_site_offset(fit$model, bd$re,
                                              fit$model$n_sites, S)
  expect_identical(dim(so), c(fit$model$n_sites, S))
  # A per-site term has no per-visit offset; it rides the occupancy predictor.
  expect_null(c0$off_det)
  expect_null(c0$off_pos)

  base <- list(fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp)
  on  <- .ocnl_score(do.call(tulpaObs:::.occu_cover_ploglik_core,
                             c(base, list(c0$field_occ, c0$field_pos))))
  off <- .ocnl_score(do.call(tulpaObs:::.occu_cover_ploglik_core,
                             c(base, list(c0$field_occ - so, c0$field_pos))))
  # Removing the per-site offset is the population-mean score: the fit describes
  # its own data materially worse, and the criteria land elsewhere.
  expect_gt(on[["lppd"]] - off[["lppd"]], 3)
  expect_gt(abs(on[["elpd"]] - off[["elpd"]]), 3)

  # Every visit of a site carries that site's group offset, so the per-site
  # posterior means collapse to one value per level; the fit reports the same
  # quantity independently as its BLUPs, centred over groups.
  grp <- tapply(rowMeans(so), fit$model$re_psi$group_idx, mean)
  expect_lt(max(abs(as.numeric(grp) - mean(grp) - fit$re$psi$blup)), 0.1)

  # The public criteria read the same components, so they score the term rather
  # than the population-mean model -- and say nothing about dropping it.
  expect_warning(w <- tobs_waic(fit, n.draws = S), NA)
  expect_lt(abs(w$elpd_waic - on[["elpd"]]),
            abs(w$elpd_waic - off[["elpd"]]))
  cp <- tobs_cpo(fit, n.draws = S)
  expect_true(is.finite(cp$elpd_loo))
  expect_lt(abs(cp$elpd_loo - on[["elpd"]]), abs(cp$elpd_loo - off[["elpd"]]))
  pp <- tobs_ppc(fit, n.samples = 40L)
  expect_true(is.finite(pp$bayesian.p))
})
