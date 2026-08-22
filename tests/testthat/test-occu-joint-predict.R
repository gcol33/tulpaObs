# =============================================================================
# test-occu-joint-predict.R - the per-arm predictor the joint predict routes
# assemble (R/occu_cover_predict.R). One state builder serves the occu_cover() /
# cover() arms and the occupancy-only occu() SVC route, so the arms a fit carries
# drive which quantities are computed, and every arm's linear predictor picks up
# that arm's RE BLUP offset: the occupancy arm reads the fit's "psi" terms, the
# detection arm its "p" terms.
#
# The bundle-level blocks drive the state builder directly from a hand-built draw
# bundle, so they are engine-independent -- the same footing as the RE design
# builder block in test-occu-cover-obs-re.R.
# =============================================================================

# An occu()-shaped joint fit: the occupancy / detection formulas and the fitted
# designs the arm-slot resolver reads widths off.
.ojp_fit <- function(det_formula = ~ 1, p_det = 1L) {
  model <- list(
    formulas    = list(occ = ~ x, det = det_formula),
    X_processes = list(matrix(0, 1L, 2L), matrix(0, 1L, p_det)))
  structure(list(model = model, occu_only_joint = TRUE),
            class = c("tobs_fit", "tulpa_fit"))
}

# Draw bundle in the shape .tobs_joint_draws_occu() returns: an occupancy and a
# detection coefficient block, and one field block at zero amplitude on both arms
# so the field drops out of the accumulator and the linear predictor is the
# betas plus the RE offset exactly. `re` carries whatever RE terms are pinned.
.ojp_bundle <- function(n_draws = 5L, n_cells = 6L, re = NULL, p_det = 1L) {
  list(n = n_draws, positive = NA_character_, cells = seq_len(n_draws),
       disp = rep(1, n_draws),
       b = list(occ = cbind(rep(0.5, n_draws), rep(-0.25, n_draws)),
                det = matrix(0.3, n_draws, p_det),
                pos = NULL),
       blocks = list(list(z = matrix(0, n_draws, n_cells),
                          amp_occ = rep(0, n_draws),
                          amp_pos = rep(0, n_draws),
                          weight  = NULL)),
       n_cells = n_cells, re = re)
}

# Random intercept on `arm`: one latent column per group, constant across draws
# so the offset each row picks up is exactly its group's value.
.ojp_re_intercept <- function(arm, n_draws, blup) {
  ng <- length(blup)
  list(arm = arm, var = "habitat", levels = paste0("h", seq_len(ng)),
       n_coefs = 1L, n_groups = ng,
       draws = matrix(rep(blup, each = n_draws), n_draws, ng))
}

# Random slope on `arm`: coefficient-major latent draws, the intercept block
# followed by the `x` slope block.
.ojp_re_slope <- function(arm, n_draws, b0, b1) {
  ng <- length(b0)
  list(arm = arm, var = "habitat", levels = paste0("h", seq_len(ng)),
       n_coefs = 2L, n_groups = ng, coef_names = c("(Intercept)", "x"),
       draws = cbind(matrix(rep(b0, each = n_draws), n_draws, ng),
                     matrix(rep(b1, each = n_draws), n_draws, ng)))
}

.ojp_newdata <- function() {
  data.frame(x       = c(0, 0.5, -0.5, 1, -1, 0.25),
             habitat = factor(paste0("h", c(1L, 2L, 3L, 1L, 2L, 3L))))
}

.ojp_chain_adj <- function(N) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  adj
}


test_that("occupancy-only joint predict adds the psi-arm RE offset", {
  fit  <- .ojp_fit()
  n    <- 5L
  blup <- c(-1.5, 0.8, 1.5)
  bun  <- .ojp_bundle(n_draws = n, re = list(.ojp_re_intercept("psi", n, blup)))
  nd   <- .ojp_newdata()
  cell <- seq_len(nrow(nd))
  g    <- c(1L, 2L, 3L, 1L, 2L, 3L)

  with_g <- tulpaObs:::.tobs_joint_arm_states(fit, bun, nd, cell, "occ")$p
  no_g   <- tulpaObs:::.tobs_joint_arm_states(
    fit, bun, nd[, "x", drop = FALSE], cell, "occ")$p

  # The two predictions differ on every row: the group's latent draws enter psi
  # when `newdata` supplies the grouping column, and only then.
  expect_false(isTRUE(all.equal(with_g, no_g)))
  expect_equal(stats::qlogis(with_g) - stats::qlogis(no_g),
               matrix(blup[g], nrow(nd), n), ignore_attr = TRUE)

  # Without the grouping column the row sits at the population mean -- the
  # occupancy betas alone, since the field block carries zero amplitude.
  expect_equal(no_g, matrix(stats::plogis(0.5 - 0.25 * nd$x), nrow(nd), n),
               ignore_attr = TRUE)

  # An unseen level shrinks the term the same way a missing column does.
  nd_new <- nd
  nd_new$habitat <- factor(rep("newHab", nrow(nd)))
  expect_equal(tulpaObs:::.tobs_joint_arm_states(fit, bun, nd_new, cell, "occ")$p,
               no_g)
})


test_that("occupancy-only joint predict adds the detection-arm RE offset", {
  fit  <- .ojp_fit()
  n    <- 5L
  blup <- c(-1.0, 0.6, 1.2)
  nd   <- .ojp_newdata()
  cell <- seq_len(nrow(nd))
  g    <- c(1L, 2L, 3L, 1L, 2L, 3L)

  bun_p <- .ojp_bundle(n_draws = n, re = list(.ojp_re_intercept("p", n, blup)))
  with_g <- tulpaObs:::.tobs_joint_arm_states(fit, bun_p, nd, cell, "det")$p_det
  no_g   <- tulpaObs:::.tobs_joint_arm_states(
    fit, bun_p, nd[, "x", drop = FALSE], cell, "det")$p_det

  expect_false(isTRUE(all.equal(with_g, no_g)))
  expect_equal(stats::qlogis(with_g) - stats::qlogis(no_g),
               matrix(blup[g], nrow(nd), n), ignore_attr = TRUE)
  expect_equal(no_g, matrix(stats::plogis(0.3), nrow(nd), n),
               ignore_attr = TRUE)

  # The arm keys are separate: an occupancy-arm term leaves detection at the
  # population mean, and vice versa.
  bun_psi <- .ojp_bundle(n_draws = n, re = list(.ojp_re_intercept("psi", n, blup)))
  expect_equal(
    tulpaObs:::.tobs_joint_arm_states(fit, bun_psi, nd, cell, "det")$p_det, no_g)
  expect_equal(
    tulpaObs:::.tobs_joint_arm_states(fit, bun_p, nd, cell, "occ")$p,
    matrix(stats::plogis(0.5 - 0.25 * nd$x), nrow(nd), n), ignore_attr = TRUE)
})


test_that("occupancy-only joint predict weights a random slope by its covariate", {
  fit  <- .ojp_fit()
  n    <- 4L
  b0   <- c(-0.6, 0.2, 0.9)
  b1   <- c(0.4, -0.7, 1.1)
  bun  <- .ojp_bundle(n_draws = n, re = list(.ojp_re_slope("psi", n, b0, b1)))
  nd   <- .ojp_newdata()
  cell <- seq_len(nrow(nd))
  g    <- c(1L, 2L, 3L, 1L, 2L, 3L)

  with_g <- tulpaObs:::.tobs_joint_arm_states(fit, bun, nd, cell, "occ")$p
  no_g   <- tulpaObs:::.tobs_joint_arm_states(
    fit, bun, nd[, "x", drop = FALSE], cell, "occ")$p

  # The row offset is the group intercept plus the slope coefficient weighted by
  # the covariate value in `newdata` (row 1 sits at x = 0, so it carries the
  # group intercept alone).
  expect_equal(stats::qlogis(with_g) - stats::qlogis(no_g),
               matrix(b0[g] + nd$x * b1[g], nrow(nd), n), ignore_attr = TRUE)
  expect_equal(stats::qlogis(with_g[1L, ]) - stats::qlogis(no_g[1L, ]),
               rep(b0[1L], n))
})


test_that("joint predict designs only the arms the caller reads", {
  # The detection formula names a covariate `newdata` does not carry, so a
  # detection prediction errors while an occupancy prediction is unaffected.
  fit  <- .ojp_fit(det_formula = ~ dcov, p_det = 2L)
  bun  <- .ojp_bundle(n_draws = 3L, p_det = 2L)
  nd   <- .ojp_newdata()
  cell <- seq_len(nrow(nd))

  psi <- tulpaObs:::.tobs_joint_arm_states(fit, bun, nd, cell, "occ")
  expect_equal(dim(psi$p), c(nrow(nd), 3L))
  expect_null(psi$p_det)
  expect_error(tulpaObs:::.tobs_joint_arm_states(fit, bun, nd, cell, "det"))

  # A fit with no cover arm reports no cover quantities even when the caller
  # asks for every arm the fit could carry.
  st <- tulpaObs:::.tobs_joint_arm_states(.ojp_fit(), .ojp_bundle(n_draws = 3L),
                                          nd, cell)
  expect_named(st, c("p", "p_det"))
})


test_that("occu() SVC joint fit predicts occupancy, detection and change", {
  skip_on_cran()
  skip_if_fast()

  set.seed(246L)
  N   <- 20L
  adj <- .ojp_chain_adj(N)
  cell_dat <- data.frame(site_id = seq_len(N), x = rnorm(N))
  long <- data.frame(site_id = rep(seq_len(N), each = 3L),
                     visit = rep(1:3, N), y = rbinom(3L * N, 1L, 0.4),
                     w = rnorm(3L * N))
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = "w")
  fit <- suppressWarnings(tobs(
    formula = ~ x + icar(graph = adj, group_var = "site_id") +
                icar(graph = adj, weight = x, group_var = "site_id"),
    data = cell_dat, family = occu(), detection = ~ w,
    y = od$y, visits = od$det.covs,
    method = "nested_laplace", control = list(verbose = FALSE, max.iter = 50L)))

  nd <- data.frame(cell = seq_len(N), x = 0)
  po <- predict(fit, newdata = nd, type = "occupancy", nsim = 200, draws = FALSE)
  expect_s3_class(po, "tobs_prediction")
  expect_equal(nrow(po), N)
  expect_true(all(po$mean > 0 & po$mean < 1))

  pb <- predict(fit, newdata = nd, type = "both", nsim = 200, draws = FALSE)
  expect_named(pb, c("occupancy", "detection"))
  expect_true(all(pb$detection$mean > 0 & pb$detection$mean < 1))

  ch <- as.data.frame(predict(fit, newdata = nd, type = "change",
                              times = c(-1, 1), time_col = "x",
                              nsim = 200, draws = FALSE))
  expect_true(all(c("psi_T1", "psi_T2", "delta_psi") %in% names(ch)))
  expect_equal(ch$delta_psi, ch$psi_T2 - ch$psi_T1, tolerance = 1e-8)
})
