# The outer-grid cost controls the engine implements, reached through
# `occu_cover()`: the cheap-pass screen (`prune` / `prune.tol`) and the two
# placement knobs (`auto.recenter`, `recenter.pilot`).
#
# Each of these is a knob on `tulpa_nested_laplace_joint()`'s control list, and
# `occu_cover()` builds that list itself rather than forwarding `control`
# wholesale, so a key reaches the engine only if the family declares it AND the
# builder names it. Both halves are asserted here: the fits below would pass a
# value-only test even if the control were dropped on the floor, so each one
# also reads the provenance field the engine writes ONLY when the knob fired.
#
# The pruned fit is compared against the unpruned one because the engine's
# safety gate re-solves the full grid whenever the cheap ranking is unreliable
# -- a prune that reaches the engine either matches the full-grid answer or
# falls back to it, and either outcome is a pass on the numbers.

oc_grid_fixture <- function(seed = 4242L, N = 30L, J = 4L) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = seed
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(adj = adj, od = od, y_pos = y_pos,
       cell_dat = cbind(data.frame(site_id = seq_len(N)), sim$data))
}

oc_grid_fit <- function(fx, ctrl = list()) {
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = fx$adj), data = fx$cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = fx$od$y, y_pos = fx$y_pos, visits = fx$od$det.covs,
    method = "nested_laplace",
    control = c(list(verbose = FALSE, max.iter = 500L, engine = "joint"), ctrl)
  ))
}

oc_betas <- function(fit) {
  fit$means[c("psi_(Intercept)", "psi_occ_cov1",
              "p_(Intercept)", "pos_(Intercept)")]
}


test_that("occu_cover() declares the outer-grid cost controls", {
  keys <- occu_cover("lognormal")$control_keys
  expect_true(all(c("prune", "prune.tol",
                    "auto.recenter", "recenter.pilot") %in% keys))
})


test_that("prune reaches the engine and does not move the answer", {
  fx   <- oc_grid_fixture()
  full <- oc_grid_fit(fx)
  cut  <- oc_grid_fit(fx, list(prune = TRUE, prune.tol = 1e-3))

  # The screen ran: the kernel emits prune_mask only when it actually screened,
  # and the safety gate replaces the result with a full-grid one carrying
  # prune_fallback_triggered when the cheap ranking was unreliable.
  jf <- cut$joint_fit
  expect_true(!is.null(jf$prune_mask) || isTRUE(jf$prune_fallback_triggered))

  expect_equal(oc_betas(cut), oc_betas(full), tolerance = 1e-4)
})


test_that("prune.tol alone is not read as prune", {
  # `prune` is a unique prefix of `prune.tol`, so a `$` read of the control
  # list returns the TOLERANCE on a fit that sets only the tolerance, which
  # switches screening on where the caller asked only to size it.
  fx  <- oc_grid_fixture()
  fit <- oc_grid_fit(fx, list(prune.tol = 1e-3))
  expect_null(fit$joint_fit$prune_mask)
  expect_null(fit$joint_fit$prune_fallback_triggered)
})


test_that("recenter.pilot reaches the engine's placement pass", {
  fx  <- oc_grid_fixture()
  fit <- oc_grid_fit(fx, list(recenter.pilot = TRUE))
  jf  <- fit$joint_fit
  expect_true(!is.null(jf$outer_grid_pilot) ||
                !is.null(jf$outer_grid_pilot_declined))
})


test_that("auto.recenter = FALSE reaches the engine and is recorded", {
  fx  <- oc_grid_fixture()
  fit <- oc_grid_fit(fx, list(auto.recenter = FALSE))
  expect_identical(fit$joint_fit$outer_grid_recenter_declined,
                   "auto_recenter_disabled")
})


test_that("the joint path refuses a per-axis auto.recenter policy", {
  fx <- oc_grid_fixture()
  expect_error(oc_grid_fit(fx, list(auto.recenter = "rail")),
               "per-axis placement policy")
})
