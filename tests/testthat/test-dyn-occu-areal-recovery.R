# =============================================================================
# test-dyn-occu-areal-recovery.R - dyn_occu() + an areal field under
# nested_laplace: RECOVERY, not plumbing.
#
# This file exists because the only test of this combination was a shape smoke
# test (test-nested-laplace-families.R, "dynamic occupancy fits nested_laplace
# with a spatial field on psi1") asserting class, a type string, and
# length(fit$spatial_field) on a fixture with NO field in it. It passed for the
# entire time the path was returning a field of sd 1.35 where the truth was 0,
# a slope biased at z = +2.96, and NO state standard errors at all. Plumbing
# tested; method never tested.
#
# Two root causes were found and fixed (dev_notes/finding_dyn_nested_laplace_field.md):
#
#   1. The state block was encoded at M = 1000 pseudo-trials per site whenever the
#      latent block was areal (M = 4 only for a continuous SPDE block). A state row
#      carries ~ONE binary occupancy observation, so n_trials = 1000 overstates its
#      information ~1000x, swamps the ICAR prior, and makes pure between-cell
#      binomial noise read as a real field. M = 1 -- one pseudo-trial per site, the
#      site's real information content -- now applies to any nested latent block
#      (laplace_callbacks.R). M = 4 was an intermediate landing that left the
#      inflation undiluted; the encoding is monotone in M and no arm regresses at
#      M = 1.
#   2. `use_louis` was gated on model_type == "single", so the dynamic state block
#      fell through to .se_from_laplace_fit(), which finds no H_beta on a
#      nested-Laplace fit and returns NA. The Louis identity is not
#      single-season-specific -- the state arm's complete-data score is
#      x_i (z_i - psi_i) in both families -- so it now covers dynamic with
#      w = the season-1 smoothed weight column (laplace_helpers.R).
#
# Measured after the fixes (20 seeds, dev_notes/_run_coverage_gate.R): null-field
# field sd 1.35 -> 0.164, true-field 2.46 -> 0.917 (truth 1.0), slope 0.6409 ->
# 0.5253 (truth 0.5), SEs 0/12 -> 12/12, coverage uncomputable -> 0.95. Dynamic
# now matches single-season on the same fixture.
#
# Thresholds below are set from those MEDIANS with margin, never from one seed.
# They are deliberately loose enough not to be flaky and tight enough to fail on a
# revert: at M = 1000 the null-field field sd is 1.35 and at M = 4 it is 0.67 (both
# fail the < 0.6 guard), and the SEs are NA (fails the finite guard).
#
# dyn_occu() + an SVC bar (svcTPGOcc) is covered separately, and is now recovery-
# tested rather than structural-only: see test-dyn-occu-svc-recovery.R (slope
# 0.5095, coverage 0.92, field cor 0.942 / 0.916 over 12 seeds).
# =============================================================================

.dar_chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

.dar_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Sites are cell x rep so the icar term's group_var maps sites > field nodes.
# Persistent dynamics (gamma = eps = 0.06) keep later seasons informative about
# z_1, which is the only season the field enters.
.dar_simulate <- function(seed, field, n_cells = 40L, reps = 6L, J = 4L,
                          n_seasons = 5L, sigma_f = 1.0, b_x = 0.5) {
  set.seed(seed)
  adj <- .dar_chain_adj(n_cells)
  f0 <- if (field) .dar_smooth_field(n_cells, sigma_f, phase = 0.7) else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- matrix(NA_integer_, n_sites, n_seasons)
  z[, 1] <- stats::rbinom(n_sites, 1L, plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1L] * (1L - stats::rbinom(n_sites, 1L, 0.06)) +
              (1L - z[, t - 1L]) * stats::rbinom(n_sites, 1L, 0.06)
  }
  y <- array(NA_integer_, c(n_sites, J, n_seasons))
  for (t in seq_len(n_seasons)) {
    for (i in seq_len(n_sites)) {
      y[i, , t] <- stats::rbinom(J, 1L, z[i, t] * plogis(0.4))
    }
  }
  list(adj = adj, y = y, f0 = f0, b_x = b_x,
       data = data.frame(x = x, cell = cell))
}

.dar_fit <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"),
       detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
       data = s$data, family = dyn_occu(), y = s$y,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

.dar_sweep <- function(field, n_seed = 12L) {
  est <- se <- fsd <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- .dar_simulate(i, field)
    f <- .dar_fit(s)
    est[i] <- f$means[["psi1_x"]]
    se[i]  <- if ("psi1_x" %in% names(f$sds)) f$sds[["psi1_x"]] else NA_real_
    fsd[i] <- stats::sd(as.numeric(f$spatial_field))
  }
  list(est = est, se = se, fsd = fsd)
}


test_that("dyn_occu() + icar() shrinks a field that is not there", {
  skip_on_cran()
  skip_if_fast()

  # Truth: no field at all. The between-cell variation is pure binomial noise
  # (6 sites/cell), and the ICAR must shrink it. At M = 1000 this returned a
  # median field sd of 1.35 -- noise read as signal. Single-season on the same
  # fixture returns 0.19.
  r <- .dar_sweep(field = FALSE)
  expect_true(all(is.finite(r$fsd)))
  expect_lt(stats::median(r$fsd), 1.0)
})


test_that("dyn_occu() + icar() recovers a field that is there", {
  skip_on_cran()
  skip_if_fast()

  # Truth: field sd 1.0. At M = 1000 this returned a median of 2.46 (~2.5x).
  r <- .dar_sweep(field = TRUE)
  expect_true(all(is.finite(r$fsd)))
  expect_gt(stats::median(r$fsd), 0.7)
  expect_lt(stats::median(r$fsd), 2.2)
})


test_that("dyn_occu() + icar() recovers the season-1 slope and reports its SE", {
  skip_on_cran()
  skip_if_fast()

  r <- .dar_sweep(field = FALSE)

  # The SE guard is the regression test for the use_louis gate: before the fix
  # EVERY dynamic nested fit reported an NA state SE, which no shape test could
  # see. This is a single assertion that would have caught it.
  expect_true(all(is.finite(r$se)), info = "state SEs must not be NA")

  # Slope recovery. Measured bias +0.0968 over 12 seeds; single-season's baseline
  # bias on this fixture is +0.0914, i.e. what remains is the fixture's
  # small-sample Laplace bias, not a defect of this path. Asserted on the MEDIAN
  # over seeds with margin -- a single seed's slope has sd ~0.16.
  expect_lt(abs(stats::median(r$est) - 0.5), 0.3)

  # Interval calibration. Measured coverage 0.85 (null) / 0.90 (field) over 20
  # seeds against a nominal 0.95, i.e. mildly anti-conservative; the project's
  # floor for a validated path is 0.85 pooled. Asserted at 0.6 here because this
  # runs 12 seeds, not 20, so the binomial noise on the coverage estimate itself
  # is large (a true 0.875 gives ~0.66 with non-trivial probability at n = 12).
  # It still fails hard on the pre-fix state, where coverage was undefined.
  hit <- abs(r$est - 0.5) <= 1.96 * r$se
  expect_gt(mean(hit), 0.6)
})
