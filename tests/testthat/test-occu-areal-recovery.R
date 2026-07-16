# =============================================================================
# test-occu-areal-recovery.R - single-season occu() + an areal field under
# nested_laplace: RECOVERY, not plumbing.
#
# The generic nested-Laplace areal path had no recovery test at all. Its only
# tests (test-nested-laplace-occu.R) say so in their own header: "these tests
# only assert that the pipeline runs and returns a sensibly-shaped tobs_fit.
# Recovery is tracked under a follow-up." This is that follow-up.
#
# WHY THIS FIXTURE USES A FIELD, AND NOT A NULL ONE
#
# Every "single-season is healthy" claim in the dynamic-field investigation
# (dev_notes/finding_dyn_nested_laplace_field.md) rested on null-field fixtures,
# and a null-field fixture CANNOT test this path. The outer grid is one axis,
# b1.tau, 9 log-spaced cells (dev_notes/_probe_grid_dump.R):
#
#   tau  : 0.3000 0.5335 0.9487 1.6870 3.0000 5.3348 9.4868 16.8702 30.0000
#   sigma: 1.8257 1.3691 1.0267 0.7699 0.5774 0.4330 0.3247  0.2435  0.1826
#
# sigma = 0 is OUTSIDE that grid -- 0.1826 is the smallest field those cells can
# express. Against a null field the marginal therefore piles onto the LAST cell
# (measured weight 0.603 on tau = 30) and reports sigma ~= 0.22 because that is
# the grid's floor, not because it identified anything. "Correctly shrunk" and
# "pinned at the grid ceiling" produce the same number, so the test would pass
# either way -- which is exactly how the dynamic path shipped broken.
#
# So the fixture uses a field the grid can represent in its INTERIOR, and the
# test asserts the peak cell lands off BOTH boundaries. Measured on a truth of
# marginal field sd 1.0: peak cell 6 of 9, weight 0.309, with the weight spread
# over cells 4-9 -- an interior mode, i.e. the precision is identified rather
# than pinned.
#
# Note tau is the ICAR CONDITIONAL precision on increments, not 1/sd^2 of the
# realised surface: a marginal field sd of 1.0 on this chain graph fits at
# sigma ~= 0.40. So recovery is asserted on the field surface (its sd and its
# correlation with the truth), never on the reported `sigma` hyperparameter.
#
# Measured (8 seeds, medians, dev_notes/_probe_single_areal_interior.R; truth
# slope 0.5, 40 cells x 6 sites, J = 4):
#
#                    slope    SEs    field sd (truth)   field cor   coverage
#   null field       0.4775   8/8    0.1658  (0)        -           0.88
#   interior field   0.5563   8/8    0.9189  (1.0)      0.943       0.88
#
# Thresholds are set from those MEDIANS with margin, never from one seed.
# =============================================================================

.sar_chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

.sar_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Sites are cell x rep so the icar term's group_var maps sites > field nodes.
.sar_simulate <- function(seed, sigma_f, n_cells = 40L, reps = 6L, J = 4L,
                          b_x = 0.5) {
  set.seed(seed)
  adj <- .sar_chain_adj(n_cells)
  f0 <- if (sigma_f > 0) .sar_smooth_field(n_cells, sigma_f, phase = 0.7)
        else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L,
                     plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) y[i, ] <- stats::rbinom(J, 1L, z[i] * plogis(0.4))
  list(adj = adj, y = y, f0 = f0, data = data.frame(x = x, cell = cell))
}

.sar_fit <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"),
       detection = ~ 1, data = s$data, family = occu(), y = s$y,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

.sar_sweep <- function(sigma_f, n_seed = 8L) {
  est <- se <- fsd <- fcor <- peak <- ncell <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- .sar_simulate(i, sigma_f)
    f <- .sar_fit(s)
    est[i] <- f$means[["psi_x"]]
    se[i]  <- if ("psi_x" %in% names(f$sds)) f$sds[["psi_x"]] else NA_real_
    ff <- as.numeric(f$spatial_field)
    fsd[i] <- stats::sd(ff)
    if (sigma_f > 0 && length(ff) == length(s$f0)) fcor[i] <- stats::cor(ff, s$f0)
    w <- f$nested_laplace$occ_fit$weights
    if (!is.null(w)) { peak[i] <- which.max(w); ncell[i] <- length(w) }
  }
  list(est = est, se = se, fsd = fsd, fcor = fcor, peak = peak, ncell = ncell)
}


test_that("occu() + icar() recovers an interior field surface", {
  skip_on_cran()
  skip_if_fast()

  # Truth: marginal field sd 1.0. Measured median sd 0.9189, cor 0.943
  # (min 0.9147 over the 8 seeds).
  r <- .sar_sweep(sigma_f = 1.0)
  expect_true(all(is.finite(r$fsd)))
  expect_gt(stats::median(r$fsd), 0.65)
  expect_lt(stats::median(r$fsd), 1.4)

  # The surface itself, not just its scale. Asserted on the median AND the
  # minimum: a median-only guard passes while individual seeds collapse.
  expect_gt(stats::median(r$fcor), 0.85)
  expect_gt(min(r$fcor), 0.7)
})


test_that("occu() + icar() identifies the field precision rather than pinning it", {
  skip_on_cran()
  skip_if_fast()

  # THE assertion this file exists for. Against a null field the outer grid's
  # weight piles onto the last cell because sigma = 0 is off the grid, so a
  # null-field fixture cannot distinguish an identified precision from one
  # pinned at a boundary. On a field the grid can represent in its interior the
  # mode must sit strictly inside: measured peak cell 6 of 9 (weight 0.309),
  # never 1 or 9 across the seeds.
  r <- .sar_sweep(sigma_f = 1.0)
  expect_true(all(is.finite(r$peak)))
  expect_true(all(r$peak > 1), info = "grid mode pinned at the smallest tau")
  expect_true(all(r$peak < r$ncell), info = "grid mode pinned at the largest tau")
})


test_that("occu() + icar() shrinks a field that is not there", {
  skip_on_cran()
  skip_if_fast()

  # Truth: no field; the between-cell variation is pure binomial noise. Measured
  # median 0.1658.
  #
  # This asserts SHRINKAGE, not identification -- see the file header: the mode
  # sits on the grid's last cell here by construction, so this arm says only that
  # the path does not invent a field. It is the weak half of the pair, and it is
  # the arm the dynamic path passed while broken. Deliberately not asserting
  # peak == 9: that would lock in the grid's floor as expected behaviour and fail
  # if the tau axis is later widened, which would be an improvement.
  r <- .sar_sweep(sigma_f = 0)
  expect_true(all(is.finite(r$fsd)))
  expect_lt(stats::median(r$fsd), 0.6)
})


test_that("occu() + icar() recovers the slope and reports its SE", {
  skip_on_cran()
  skip_if_fast()

  r <- .sar_sweep(sigma_f = 1.0)

  expect_true(all(is.finite(r$se)), info = "state SEs must not be NA")

  # Slope recovery on the MEDIAN over seeds; a single seed's slope has sd ~0.15.
  expect_lt(abs(stats::median(r$est) - 0.5), 0.3)

  # Interval calibration. Measured 0.88 over 8 seeds against a nominal 0.95;
  # asserted at 0.6 because 8 seeds put large binomial noise on the coverage
  # estimate itself.
  hit <- abs(r$est - 0.5) <= 1.96 * r$se
  expect_gt(mean(hit), 0.6)
})
