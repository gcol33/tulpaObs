# =============================================================================
# test-dyn-occu-svc-recovery.R - multi-season occupancy with an areal
# varying-coefficient (SVC) bar: the spOccupancy `svcTPGOcc` analogue. RECOVERY.
#
#   dyn_occu() + spatial(~ 1 + w || cell, graph)
#
# An intercept field plus one covariate-weighted trend field, both on the same
# graph, both entering season-1 occupancy. Recovery is asserted for the slope,
# BOTH field surfaces, and the interval coverage.
#
# This file only became writable after the M encoding was fixed. Before it, the
# same fixture gave slope 0.9735 against a truth of 0.5 (z = 4.84), coverage 0.58,
# and field correlations of 0.518 / 0.461 -- a coin flip. The state block was
# encoded at M = 1000 pseudo-trials per site while a site carries ONE binary
# occupancy observation, so the field prior was swamped ~1000x, both fields
# inflated ~2x, and the slope inflated with them through the logistic
# conditional-vs-marginal factor sqrt(1 + 0.346 sigma^2). laplace_callbacks.R now
# encodes any nested latent block at M = 1.
#
# Measured after the fix, 12 seeds (dev_notes/_run_dyn_svc_seeds2.R):
#   slope    mean 0.5095, bias 0.0095, mcse 0.0496, z = 0.19   (truth 0.5)
#   coverage 0.92                                              (floor 0.85)
#   f0 cor   median 0.942, min 0.915                           (truth sd 1.0)
#   f1 cor   median 0.916, min 0.814                           (truth sd 0.8)
#
# Thresholds are set from the 12-seed MINIMA with margin, never from one seed. An
# earlier draft of the sibling file asserted cor > 0.6 tuned on seed 1 (0.866 /
# 0.791) when the medians were 0.425 / 0.303; it passed by luck and was not
# shipped. The guards below fail hard on the pre-fix state (field cor ~0.5, slope
# ~0.97, coverage 0.58).
# =============================================================================

.svcr_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Sites are cell x rep, so the bar's `cell` grouping exercises the sites >
# field-nodes map. Persistent dynamics (gamma = eps = 0.06) keep the later
# seasons informative about z_1 -- the only season the fields enter.
.svcr_simulate <- function(seed, n_cells = 40L, reps = 6L, J = 4L,
                           n_seasons = 5L, sigma_f0 = 1.0, sigma_f1 = 0.8,
                           b_x = 0.5) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f0 <- .svcr_smooth_field(n_cells, sigma_f0, phase = 0.7)
  f1 <- .svcr_smooth_field(n_cells, sigma_f1, phase = 2.3)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  w <- as.numeric(scale(stats::rnorm(n_sites)))
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  eta1 <- stats::qlogis(0.35) + b_x * x + f0[cell] + w * f1[cell]
  z <- matrix(NA_integer_, n_sites, n_seasons)
  z[, 1] <- stats::rbinom(n_sites, 1L, plogis(eta1))
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
  list(adj = adj, y = y, f0 = f0, f1 = f1, b_x = b_x,
       data = data.frame(x = x, w = w, cell = cell))
}

.svcr_sweep <- function(n_seed = 12L) {
  est <- se <- c0 <- c1 <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- .svcr_simulate(i)
    fit <- tobs(
      ~ x + spatial(~ 1 + w || cell, graph = s$adj),
      detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
      data = s$data, family = dyn_occu(), y = s$y,
      method = "nested_laplace",
      control = list(verbose = FALSE, progress = FALSE)
    )
    tf <- fit$trend_field %||% fit$trend_fields[[1L]]
    est[i] <- fit$means[["psi1_x"]]
    se[i]  <- if ("psi1_x" %in% names(fit$sds)) fit$sds[["psi1_x"]] else NA_real_
    c0[i]  <- stats::cor(as.numeric(fit$spatial_field), s$f0)
    c1[i]  <- stats::cor(as.numeric(tf), s$f1)
  }
  list(est = est, se = se, c0 = c0, c1 = c1)
}


test_that("dyn_occu() + an SVC bar recovers both field surfaces", {
  skip_on_cran()
  skip_if_fast()

  r <- .svcr_sweep()

  # Measured medians 0.942 / 0.916, minima 0.915 / 0.814 over 12 seeds. Asserted
  # on the MEDIAN with margin; the pre-fix state medians 0.518 / 0.461.
  expect_gt(stats::median(r$c0), 0.80)
  expect_gt(stats::median(r$c1), 0.70)

  # Every seed recovers, not just the typical one -- the guard against a
  # median carried by a few lucky draws.
  expect_gt(min(r$c0), 0.70)
  expect_gt(min(r$c1), 0.60)
})


test_that("dyn_occu() + an SVC bar recovers the season-1 slope with calibrated intervals", {
  skip_on_cran()
  skip_if_fast()

  r <- .svcr_sweep()
  expect_true(all(is.finite(r$se)), info = "state SEs must not be NA")

  # Measured: mean 0.5095, bias 0.0095, z = 0.19 over 12 seeds (truth 0.5).
  # Pre-fix: 0.9735, z = 4.84.
  expect_lt(abs(stats::median(r$est) - 0.5), 0.2)

  # Measured coverage 0.92 (floor 0.85, nominal 0.95); pre-fix 0.58. Asserted at
  # 0.7 for the binomial noise on a 12-seed coverage estimate, not as a lowered
  # bar -- it still fails hard on the pre-fix state.
  hit <- abs(r$est - 0.5) <= 1.96 * r$se
  expect_gt(mean(hit), 0.7)
})
