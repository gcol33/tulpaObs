# =============================================================================
# test-int-occu-areal-recovery.R - int_occu() + an areal field under
# nested_laplace: RECOVERY, not plumbing.
#
# The companion to test-dyn-occu-areal-recovery.R. When the dynamic path's field
# defect was fixed, int_occu was recorded as "inherits the fix, untested". It did
# not inherit it -- it shared the disease and none of the cure. Three separate
# defects, none visible to the only test of this combination (a shape smoke test
# in test-nested-laplace-families.R asserting class, a type string, and
# length(fit$spatial_field) on a fixture containing NO field):
#
#   1. build_integrated_callbacks' m_step_encode had no latent_prior branch. It
#      read `if (is.null(spatial_occ)) M <- 1000L else M <- 4L`, and a nested
#      areal block arrives with spatial_occ = NULL (the prior is attached to
#      occ$prior upstream, not to the spatial slot) -- so the nested case took the
#      M = 1000 branch. A state row carries ONE binary occupancy observation;
#      M = 1000 overstates its information ~1000x and swamps the ICAR prior.
#   2. .tobs_laplace_nested() never passed latent_prior to
#      build_integrated_callbacks, so the field-aware E-step that the single-season
#      and dynamic arms use was dead code on this path.
#   3. use_louis was gated on model_type %in% c("single", "dynamic"), so the
#      integrated state block fell through to .se_from_laplace_fit(), which finds
#      no H_beta on a nested-Laplace fit and returns NA -- every integrated nested
#      fit reported NA state SEs while the per-source detection arms were finite.
#
# The Louis identity is not family-specific: the state arm's complete-data score
# is x_i (z_i - psi_i) here too, and integrated shares one psi across sources, so
# em_result$weights is already the per-site w = E[z_i | y] the identity needs (no
# column slice, unlike dynamic's [n_sites x n_seasons] matrix).
#
# Measured after the fixes (8 seeds, medians, dev_notes/probe_int_occu_areal.R;
# truth slope 0.5, two sources at J = 4 / 3, 40 cells x 6 sites):
#
#                     slope    SEs    field sd (truth)   coverage
#   null field        0.4888   8/8    0.1783  (0)        0.88
#   field             0.5830   8/8    0.9758  (1.0)      0.88
#
# Thresholds are set from those MEDIANS with margin, never from one seed.
# =============================================================================

.iar_chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

.iar_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Sites are cell x rep so the icar term's group_var maps sites > field nodes.
# Two sources with different visit budgets share one psi -- the point of the
# family, and the arm the field enters.
.iar_simulate <- function(seed, field, n_cells = 40L, reps = 6L, J = c(4L, 3L),
                          sigma_f = 1.0, b_x = 0.5) {
  set.seed(seed)
  adj <- .iar_chain_adj(n_cells)
  f0 <- if (field) .iar_smooth_field(n_cells, sigma_f, phase = 0.7) else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L,
                     plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  y <- lapply(J, function(j) {
    m <- matrix(0L, n_sites, j)
    for (i in seq_len(n_sites)) m[i, ] <- stats::rbinom(j, 1L, z[i] * plogis(0.4))
    m
  })
  names(y) <- c("s1", "s2")
  list(adj = adj, y = y, f0 = f0, data = data.frame(x = x, cell = cell))
}

.iar_fit <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"),
       detection = ~ 1, data = s$data, family = int_occu(), y = s$y,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

.iar_sweep <- function(field, n_seed = 8L) {
  est <- se <- fsd <- fcor <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- .iar_simulate(i, field)
    f <- .iar_fit(s)
    est[i] <- f$means[["psi_x"]]
    se[i]  <- if ("psi_x" %in% names(f$sds)) f$sds[["psi_x"]] else NA_real_
    ff <- as.numeric(f$spatial_field)
    fsd[i] <- stats::sd(ff)
    if (field && length(ff) == length(s$f0)) fcor[i] <- stats::cor(ff, s$f0)
  }
  list(est = est, se = se, fsd = fsd, fcor = fcor)
}


test_that("int_occu() + icar() shrinks a field that is not there", {
  skip_on_cran()
  skip_if_fast()

  # Truth: no field. The between-cell variation is pure binomial noise (6
  # sites/cell) and the ICAR must shrink it. Measured median 0.1783.
  #
  # The 0.6 guard is set from that median with margin. This fixture's pre-fix
  # value was never measured -- the M-step branch was missing rather than
  # mis-set, so the fix landed before a pre-fix sweep -- so the threshold is NOT
  # calibrated against a known broken value here. The analogous dynamic arm on
  # its own fixture returned 1.35 at M = 1000 and 0.67 at M = 4 against the same
  # truth of 0, which is the order of magnitude 0.6 is chosen to exclude.
  #
  # This asserts shrinkage, NOT identification: sigma = 0 is outside the outer
  # grid (see test-occu-areal-recovery.R), so the marginal can only pile onto the
  # smallest field the grid can express. The interior-field arm below is what
  # tests that the precision is actually identified.
  r <- .iar_sweep(field = FALSE)
  expect_true(all(is.finite(r$fsd)))
  expect_lt(stats::median(r$fsd), 0.6)
})


test_that("int_occu() + icar() recovers a field that is there", {
  skip_on_cran()
  skip_if_fast()

  # Truth: field sd 1.0, which the outer grid can represent in its interior.
  # Measured median sd 0.9758.
  r <- .iar_sweep(field = TRUE)
  expect_true(all(is.finite(r$fsd)))
  expect_gt(stats::median(r$fsd), 0.65)
  expect_lt(stats::median(r$fsd), 1.5)

  # The surface itself, not just its scale.
  expect_gt(stats::median(r$fcor), 0.8)
})


test_that("int_occu() + icar() recovers the shared slope and reports its SE", {
  skip_on_cran()
  skip_if_fast()

  r <- .iar_sweep(field = FALSE)

  # The SE guard is the regression test for the use_louis gate: before the fix
  # EVERY integrated nested fit reported an NA state SE, which no shape test
  # could see. This is a single assertion that would have caught it.
  expect_true(all(is.finite(r$se)), info = "state SEs must not be NA")

  # Slope recovery on the MEDIAN over seeds; a single seed's slope has sd ~0.15.
  expect_lt(abs(stats::median(r$est) - 0.5), 0.3)

  # Interval calibration. Measured 0.88 over 8 seeds against a nominal 0.95.
  # Asserted at 0.6 because 8 seeds put large binomial noise on the coverage
  # estimate itself; it still fails hard on the pre-fix state, where every
  # interval was NA and coverage undefined.
  hit <- abs(r$est - 0.5) <= 1.96 * r$se
  expect_gt(mean(hit), 0.6)
})
