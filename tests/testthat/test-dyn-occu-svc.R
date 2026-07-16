# =============================================================================
# test-dyn-occu-svc.R - multi-season occupancy with an areal varying-coefficient
# (SVC) spatial bar: the spOccupancy `svcTPGOcc` analogue.
#
#   dyn_occu() + spatial(~ 1 + w || cell, graph)
#
# The state formula carries an intercept field plus one covariate-weighted trend
# field, both on the same graph, entering season-1 occupancy. The nested-Laplace
# driver already supported model_type "dynamic" and its multifield expansion was
# already generic; what the dynamic binder lacked was the site-level `data` slot
# the bar resolves its node index against.
#
# SCOPE: structural. The RECOVERY tests for the SVC bar live in
# test-dyn-occu-svc-recovery.R; this file guards the wiring (the binder slot the
# bar resolves against, and that two distinct fields come back on the documented
# slots).
#
# History worth keeping, because it cost a day. An earlier draft of this file
# asserted cor > 0.6 on both fields, tuned on seed 1 (0.866 / 0.791) when the
# 12-seed MEDIANS were 0.425 / 0.303 -- the thresholds passed by luck. That
# cherry-pick was caught by running seeds, and the underlying weakness turned out
# NOT to be the fixture (the theory at the time) but a real defect: the state
# block was encoded at M = 1000 pseudo-trials per site when a site carries ONE
# binary occupancy observation, so the field prior was swamped ~1000x and both
# fields inflated. Fixed in laplace_callbacks.R (M = 1 for any nested latent
# block); the same fixture now recovers both fields at cor ~0.94 / ~0.92 with a
# 12-seed MINIMUM of 0.915 / 0.814 -- better than the number that was once
# cherry-picked. See dev_notes/finding_dyn_nested_laplace_field.md.
# =============================================================================

.svct_chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

.svct_smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

# Multi-season occupancy on a chain graph. Sites are cell x rep, so the bar's
# `cell` grouping exercises the sites > field-nodes map. Both fields enter
# season-1 occupancy; colonization / extinction are field-free.
.svct_simulate <- function(n_cells = 30L, reps = 2L, J = 4L, n_seasons = 4L,
                           sigma_f0 = 1.0, sigma_f1 = 0.8, seed = 1L) {
  set.seed(seed)
  adj <- .svct_chain_adj(n_cells)
  f0  <- .svct_smooth_field(n_cells, sigma_f0, phase = 0.7)
  f1  <- .svct_smooth_field(n_cells, sigma_f1, phase = 2.3)

  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  w    <- as.numeric(scale(stats::rnorm(n_sites)))
  x    <- as.numeric(scale(stats::rnorm(n_sites)))

  b0 <- stats::qlogis(0.35); b_x <- 0.5
  eta1 <- b0 + b_x * x + f0[cell] + w * f1[cell]

  z <- matrix(NA_integer_, n_sites, n_seasons)
  z[, 1] <- stats::rbinom(n_sites, 1L, plogis(eta1))
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1L] * (1L - stats::rbinom(n_sites, 1L, 0.15)) +
              (1L - z[, t - 1L]) * stats::rbinom(n_sites, 1L, 0.25)
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


test_that("the dynamic binder carries the site frame structured terms resolve against", {
  # A structured term resolves its node index / group_var against model$data, and
  # those columns need not appear in any arm's fixed-effect formula -- `cell`
  # below is referenced only inside the spatial term. The dynamic and integrated
  # binders omitted this slot while the single-season one carried it, which is
  # what blocked the multi-season SVC path. Guards the slot directly so a future
  # binder change surfaces here rather than as a downstream "argument is of
  # length zero" from the bar expansion.
  s <- .svct_simulate(n_cells = 6L, reps = 2L, J = 2L, n_seasons = 2L, seed = 3L)
  bind <- .tobs_bind_formulas(
    list(psi1 = ~ x, p = ~ 1, gamma = ~ 1, epsilon = ~ 1), s$data)
  model <- .tobs_build_dynamic(
    occ_formula = bind$fe$psi1, det_formula = bind$fe$p, data = s$data,
    y = s$y, col_formula = bind$fe$gamma, ext_formula = bind$fe$epsilon)

  expect_false(is.null(model$data))
  expect_true("cell" %in% names(model$data))
  expect_identical(nrow(model$data), model$n_sites)
})


test_that("dyn_occu() + an SVC spatial bar fits and returns two distinct fields", {
  skip_on_cran()
  skip_if_fast()

  # STRUCTURAL ONLY -- see the file header. Asserts the path runs and reports
  # two separate surfaces on the documented slots. The checks below are shape
  # guards (finite, right length, not the same surface reported twice), NOT
  # recovery: correlation against truth is deliberately not asserted, because
  # on this fixture it is a coin flip across seeds.
  s <- .svct_simulate(seed = 1L)
  fit <- tobs(
    ~ x + spatial(~ 1 + w || cell, graph = s$adj),
    detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
    data = s$data, family = dyn_occu(), y = s$y,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE)
  )

  f0_hat <- as.numeric(fit$spatial_field)
  f1_hat <- as.numeric(fit$trend_field %||% fit$trend_fields[[1L]])

  expect_length(f0_hat, length(s$f0))
  expect_length(f1_hat, length(s$f1))
  expect_true(all(is.finite(f0_hat)), info = "intercept field is finite")
  expect_true(all(is.finite(f1_hat)), info = "trend field is finite")

  # Two fields were fitted, not one surface reported on both slots.
  expect_lt(abs(stats::cor(f0_hat, f1_hat)), 0.95)

  # The season-1 slope is reported under the name the rest of the suite uses.
  # coef() on the dynamic family reports the four per-process intercepts as
  # probabilities; the betas live in fit$means. Its VALUE is not asserted --
  # sd across seeds is 3.65 on this fixture.
  expect_true("psi1_x" %in% names(fit$means))
  expect_true(is.finite(fit$means[["psi1_x"]]))
})


test_that("dyn_occu() + a plain areal field still fits (no SVC regression)", {
  skip_on_cran()
  skip_if_fast()

  s <- .svct_simulate(seed = 2L)
  fit <- tobs(
    ~ x + icar(graph = s$adj, group_var = "cell"),
    detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
    data = s$data, family = dyn_occu(), y = s$y,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE)
  )
  f0_hat <- as.numeric(fit$spatial_field)
  expect_length(f0_hat, length(s$f0))
  expect_true(all(is.finite(f0_hat)))
})
