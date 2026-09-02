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


test_that("a pruned cell is not reported as a non-converged one", {
  # A cell the cheap-pass screen dropped is never solved, so its log-marginal is
  # non-finite for the same reason a failed inner Newton's is. `prune_mask`
  # separates them; without it a pruned fit tells its caller that most of its
  # grid failed, which on a 25 km run is most of the grid.
  fit <- list(log_marginal = c(-1, -Inf, -Inf, -2),
              weights      = c(0.5, 0, 0, 0.5),
              prune_mask   = c(FALSE, TRUE, TRUE, FALSE))
  expect_no_warning(oc <- .tobs_joint_ok_cells(fit, "test route"))
  expect_identical(oc$ok_cells, c(1L, 4L))

  # A genuine failure alongside pruned cells is still reported, counted over
  # the cells that were actually solved.
  bad <- list(log_marginal = c(-1, -Inf, -Inf, NaN),
              weights      = c(1, 0, 0, 0),
              prune_mask   = c(FALSE, TRUE, TRUE, FALSE))
  expect_warning(.tobs_joint_ok_cells(bad, "test route"),
                 "dropping 1 / 2 outer-grid")

  # The screen runs on the base tensor and the refinement passes append their
  # cells after it, so a mask shorter than the grid indexes the base cells and
  # the appended ones were solved.
  short <- list(log_marginal = c(-1, -Inf, -2),
                weights      = c(0.5, 0, 0.5),
                prune_mask   = c(FALSE, TRUE))
  expect_no_warning(.tobs_joint_ok_cells(short, "test route"))

  # Without a mask every non-finite cell is a failure, which is what an
  # unpruned fit means by one.
  expect_warning(.tobs_joint_ok_cells(
    list(log_marginal = c(-1, -Inf), weights = c(1, 0)), "test route"),
    "dropping 1 / 2 outer-grid")
})


test_that("a pruned occu_cover fit does not warn about convergence", {
  # The package-default axes give this fixture a 5-cell grid, which the screen
  # keeps whole; the declared axes below are what make the screen drop cells,
  # so the assertion is about a fit that actually pruned. Measured: 5 of 20
  # base cells pruned, and the refinement pass appends 4 more afterwards, so
  # the mask is shorter than the grid it indexes.
  fx <- oc_grid_fixture()
  warns <- character(0)
  fit <- withCallingHandlers(
    tobs(formula = ~ occ_cov1 + bym2(graph = fx$adj), data = fx$cell_dat,
         family = occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = fx$od$y, y_pos = fx$y_pos, visits = fx$od$det.covs,
         method = "nested_laplace",
         control = list(verbose = FALSE, max.iter = 500L, engine = "joint",
                        prune = TRUE, prune.tol = 1e-3,
                        sigma.grid = c(0.3, 0.5, 0.8, 1.2, 1.8),
                        alpha.grid = c(0, 0.5, 1, 1.5))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  jf <- fit$joint_fit
  expect_gt(sum(jf$prune_mask), 0L)
  expect_gt(sum(!is.finite(jf$log_marginal)), 0L)
  expect_false(any(grepl("did not converge", warns)))
})


# --- the grid_adaptive lattice-builder knobs --------------------------------
#
# `integration = "grid_adaptive"` evaluates a strict subset of the same tensor
# lattice. Its four tuning knobs are a DIFFERENT mechanism from the three
# post-integration refinement knobs that share the `adaptive.grid` prefix, and
# all seven resolve in one place (`.tobs_adaptive_grid_control()`).

test_that("the builder knobs resolve unset so the engine owns their defaults", {
  r <- .tobs_adaptive_grid_control(list())
  expect_true(r$adaptive_grid)
  expect_equal(r$adaptive_grid_edge_thresh, 0.02)
  expect_equal(r$adaptive_grid_max_passes, 1L)
  for (k in c("adaptive_grid_cutoff", "adaptive_grid_stride",
              "adaptive_grid_max_frac", "adaptive_grid_min_cells")) {
    expect_true(k %in% names(r))
    expect_null(r[[k]])
  }
})

test_that("a sub-knob alone does not become the master flag's value", {
  # Every name in this family has `adaptive.grid` as a prefix, so a `$` read of
  # the flag returns a sub-knob's VALUE whenever exactly one sub-knob is set and
  # the flag is not -- `list(adaptive.grid.edge.thresh = 0.05)$adaptive.grid`
  # is 0.05, not NULL. The resolver reads exactly.
  for (k in c("adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
              "adaptive.grid.cutoff", "adaptive.grid.stride",
              "adaptive.grid.max.frac", "adaptive.grid.min.cells")) {
    ctl <- stats::setNames(list(0.05), k)
    expect_identical(.tobs_adaptive_grid_control(ctl)$adaptive_grid, TRUE,
                     info = k)
  }
  expect_false(.tobs_adaptive_grid_control(list(adaptive.grid = FALSE))$adaptive_grid)
})

test_that("every knob passes control validation on the families that host one", {
  keys <- c("adaptive.grid.cutoff", "adaptive.grid.stride",
            "adaptive.grid.max.frac", "adaptive.grid.min.cells")
  for (fam in list(occu_cover(), occu_multiscale_cover(), cover())) {
    expect_true(all(keys %in% fam$control_keys), info = fam$name)
  }
  # occu()'s spatial (SVC) reroute reaches the same engine through the shared
  # route group rather than its own family key list.
  expect_true(all(keys %in% .tobs_control_groups$nested_laplace_joint))
})
