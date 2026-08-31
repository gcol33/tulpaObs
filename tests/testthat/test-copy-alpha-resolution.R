# The copy coefficient's outer axis is set by STATING its nodes or by naming a
# RESOLUTION, and the two are different requests.
#
# The axis carries prior structure -- a point mass at alpha = 0 ("no coupling")
# and a log-spaced slab above it -- so `control$alpha.grid` (and
# `copy(alpha = grid(...))`) restate that structure along with the nodes, which
# is why a fit that only wants the axis integrated more finely cannot use them.
# `control$alpha.n` names the number of slab nodes and the ENGINE re-reads its
# own axis at that resolution, atom and slab bounds unchanged.
#
# It is the only way to raise this axis: it does not densify when the donor
# `sigma.grid` does, so on informative data the outer quadrature's effective
# sample size saturates on the copy amplitude while every other axis tracks the
# request.

engine_alpha_axis <- function(n = NULL) {
  as.numeric(tulpa:::.nl_grid_axis("copy_alpha", n = n))
}

# The copy-axis coordinates a fit actually integrated. `theta_grid` labels its
# axes "alpha" on the single-block path and "b<k>.alpha" on the multi-block one.
# These are the nodes the fit ENDED with: the engine's var-of-means consistency
# pass adds cells to a sharply peaked axis after integrating, so the declared
# placement is a subset of them, not the whole set.
alpha_nodes_of <- function(fit, name = "alpha") {
  tg <- fit$joint$theta_grid
  cn <- colnames(tg) %||% fit$joint$theta_names
  sort(unique(as.numeric(tg[, which(cn == name)])))
}

# Single shared ICAR field, presence-anchored and copied onto the cover arm.
simulate_copy_axis_cover <- function(seed = 4242L, N = 300L, n_s = 25L,
                                     sigma_true = 0.6, alpha_true = 1.2,
                                     phi_true = 30,
                                     beta_occ = c(-0.3, 0.7),
                                     beta_pos = c(0.4, -0.5)) {
  set.seed(seed)
  s_idx <- sample.int(n_s, N, replace = TRUE)
  w_s   <- sigma_true * as.numeric(scale(rnorm(n_s)))
  x     <- rnorm(N)
  occur <- rbinom(N, 1, plogis(beta_occ[1] + beta_occ[2] * x + w_s[s_idx]))
  mu    <- plogis(beta_pos[1] + beta_pos[2] * x + alpha_true * w_s[s_idx])
  y     <- numeric(N)
  pos   <- occur == 1L
  y[pos] <- rbeta(sum(pos), mu[pos] * phi_true, (1 - mu[pos]) * phi_true)
  list(data = data.frame(x = x, region = factor(s_idx)),
       y    = pmin(pmax(y, 0), 1 - 1e-6))
}

fit_copy_axis_cover <- function(sim, adj, ctrl) {
  suppressWarnings(tobs(
    formula = ~ x + icar(graph = adj, group_var = "region"),
    data = sim$data, family = cover("beta"), y = sim$y,
    method = "nested_laplace",
    control = c(list(sigma.grid = c(0.3, 0.6), phi.grid = c(12, 40),
                     adaptive.grid = FALSE, var.of.means.consistency = FALSE),
                ctrl)))
}


# ---- the resolver ---------------------------------------------------------

test_that("the axis is EITHER stated nodes OR a resolution", {
  # Unset: the engine's own axis, marked as ours, and no resolution asked for.
  ax <- tulpaObs:::.tobs_alpha_axis()
  expect_equal(as.numeric(ax$alpha_grid), engine_alpha_axis())
  expect_null(ax$alpha_n)
  expect_true(tulpa::is_auto_grid(ax$alpha_grid))

  # A resolution hands the engine no nodes at all -- it resolves its own axis,
  # so the structure is never restated here.
  ax_n <- tulpaObs:::.tobs_alpha_axis(n = 9)
  expect_null(ax_n$alpha_grid)
  expect_identical(ax_n$alpha_n, 9L)

  # Stated nodes win over a resolution: a field block with no copy() arrives
  # pinned (grid = 0) and has no axis left to resolve.
  ax_pin <- tulpaObs:::.tobs_alpha_axis(grid = 0, n = 9)
  expect_equal(as.numeric(ax_pin$alpha_grid), 0)
  expect_null(ax_pin$alpha_n)
})

test_that("a resolution reads the engine's axis, atom and bounds unchanged", {
  nodes <- tulpaObs:::.tobs_alpha_nodes(tulpaObs:::.tobs_alpha_axis(n = 9))
  expect_equal(as.numeric(nodes), engine_alpha_axis(9))
  # `n` counts the SLAB, so the axis comes back one node longer than asked.
  expect_length(nodes, 10L)
  # Same atom, same slab bounds -- only the spacing between them moves.
  dflt <- engine_alpha_axis()
  expect_equal(min(nodes), min(dflt))
  expect_equal(max(nodes), max(dflt))
  expect_true(any(nodes == 0))
  # The default resolution reproduces the default axis exactly, so the two
  # spellings agree where they overlap.
  expect_equal(as.numeric(tulpaObs:::.tobs_alpha_nodes(
    tulpaObs:::.tobs_alpha_axis(n = length(dflt) - 1L))), dflt)
})

test_that("a resolution must be a single positive integer", {
  expect_error(tulpaObs:::.tobs_alpha_axis(n = 0), "integer >= 1")
  expect_error(tulpaObs:::.tobs_alpha_axis(n = c(4, 8)), "integer >= 1")
  expect_error(tulpaObs:::.tobs_alpha_axis(n = "many"), "integer >= 1")
})

test_that("the trend block inherits the base axis unless it names its own", {
  base <- tulpaObs:::.tobs_alpha_axis(n = 9)
  expect_identical(tulpaObs:::.tobs_alpha_axis_trend(list(), base), base)
  expect_identical(
    tulpaObs:::.tobs_alpha_axis_trend(list(alpha.n.trend = 4), base)$alpha_n, 4L)
  expect_equal(as.numeric(tulpaObs:::.tobs_alpha_axis_trend(
    list(alpha.grid.trend = c(0, 1)), base)$alpha_grid), c(0, 1))
})

test_that("the control reader is exact, not prefix-matching", {
  # `$` partial-matches on a list, so a fit setting only the trend knob would
  # otherwise read it as the intercept block's.
  ax <- tulpaObs:::.tobs_alpha_axis_base(list(alpha.n.trend = 12))
  expect_null(ax$alpha_n)
  expect_equal(as.numeric(ax$alpha_grid), engine_alpha_axis())
})


# ---- the two spellings are refused together -------------------------------

test_that("stating nodes and naming a resolution on one block is an error", {
  expect_error(
    tulpaObs:::.tobs_check_alpha_control(
      list(alpha.grid = c(0, 1), alpha.n = 9), "cover()"),
    "not both")
  expect_error(
    tulpaObs:::.tobs_check_alpha_control(
      list(alpha.grid.trend = c(0, 1), alpha.n.trend = 9), "cover()"),
    "not both")
  # Different blocks: each takes its own spelling.
  expect_true(tulpaObs:::.tobs_check_alpha_control(
    list(alpha.grid = c(0, 1), alpha.n.trend = 9), "cover()"))
})

test_that("a copy() stating nodes and a resolution are refused together", {
  states <- tulpaObs:::.tobs_copy_states_nodes
  expect_true(states(list(alpha_grid = c(0, 1), alpha_integrate = TRUE)))
  expect_true(states(list(alpha_grid = 0.8, alpha_integrate = FALSE)))
  # copy(spatial()) with no amplitude asks for the DEFAULT axis, so it composes
  # with a resolution instead of colliding with it.
  expect_false(states(list(alpha_grid = NULL, alpha_integrate = NA)))

  expect_error(
    tulpaObs:::.tobs_check_alpha_copy(TRUE, list(alpha.n = 9), "occu_cover()"),
    "states the copy axis's nodes")
  expect_true(
    tulpaObs:::.tobs_check_alpha_copy(FALSE, list(alpha.n = 9), "occu_cover()"))
})

test_that("the resolution is a declared control key on the joint families", {
  for (fam in list(cover("beta"), occu_cover(), occu_multiscale_cover())) {
    expect_true(all(c("alpha.n", "alpha.n.trend") %in% fam$control_keys))
  }
  sim <- simulate_copy_axis_cover(N = 40L, n_s = 4L)
  expect_error(
    tobs(formula = ~ x + icar(graph = chain_adj(4L), group_var = "region"),
         data = sim$data, family = cover("beta"), y = sim$y,
         method = "nested_laplace",
         control = list(alpha.grid = c(0, 1), alpha.n = 9)),
    "not both")
})


# ---- it reaches the engine ------------------------------------------------

test_that("the resolution places the copy axis the fit integrates", {
  skip_on_cran()
  skip_if_fast()
  adj <- chain_adj(25L)
  sim <- simulate_copy_axis_cover()

  dflt <- fit_copy_axis_cover(sim, adj, list())
  fine <- fit_copy_axis_cover(sim, adj, list(alpha.n = 11))

  expect_true(all(engine_alpha_axis()   %in% alpha_nodes_of(dflt)))
  expect_true(all(engine_alpha_axis(11) %in% alpha_nodes_of(fine)))
  # The declared axis is the one that moved: the default's interior slab nodes
  # are not the fine axis's, so this is a re-placement, not an extension.
  expect_false(all(engine_alpha_axis() %in% engine_alpha_axis(11)))
  expect_gt(length(alpha_nodes_of(fine)), length(alpha_nodes_of(dflt)))

  # Naming the engine's OWN resolution reproduces the default fit, so the
  # pass-through resolves the same axis rather than a parallel one.
  same <- fit_copy_axis_cover(
    sim, adj, list(alpha.n = length(engine_alpha_axis()) - 1L))
  expect_equal(alpha_nodes_of(same), alpha_nodes_of(dflt))
  expect_equal(same$joint$log_marginal, dflt$joint$log_marginal)

  # And the finer axis is integrated, not merely declared: more nodes over the
  # same span is a different quadrature of the same posterior.
  expect_gt(length(fine$joint$log_marginal), length(dflt$joint$log_marginal))
  expect_equal(unname(coef(fine)$occ), unname(coef(dflt)$occ), tolerance = 0.05)
})

test_that("the trend block carries its own resolution", {
  skip_on_cran()
  skip_if_fast()
  adj <- chain_adj(25L)
  sim <- simulate_copy_axis_cover()
  sim$data$time.sc <- as.numeric(scale(seq_len(nrow(sim$data))))

  fit <- suppressWarnings(tobs(
    formula = ~ x + icar(graph = adj, group_var = "region") +
                icar(graph = adj, weight = time.sc, group_var = "region"),
    data = sim$data, family = cover("beta"), y = sim$y,
    method = "nested_laplace",
    control = list(sigma.grid = c(0.4, 0.8), phi.grid = c(12, 40),
                   alpha.n = 4, alpha.n.trend = 2,
                   adaptive.grid = FALSE, integration = "grid")))

  # Each block reads its own resolution off its own knob.
  expect_true(all(engine_alpha_axis(4) %in% alpha_nodes_of(fit, "b1.alpha")))
  expect_true(all(engine_alpha_axis(2) %in% alpha_nodes_of(fit, "b2.alpha")))
})

test_that("a bare copy() takes the resolution; a copy() with nodes refuses it", {
  skip_on_cran()
  skip_if_fast()
  N <- 30L; J <- 4L
  adj <- chain_adj(N)
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", adj = adj,
                             sigma = 0.8, alpha = 1.0, seed = 31337L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit_at <- function(pos_formula, ctrl) {
    suppressWarnings(tobs(
      formula = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
      family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = pos_formula,
      y = od$y, y_pos = y_pos, visits = od$det.covs,
      method = "nested_laplace",
      control = c(list(verbose = FALSE, engine = "joint", max.iter = 100L),
                  ctrl)))
  }

  # `copy(spatial())` names no amplitude: it asks for the engine's own axis, and
  # the resolution says how finely to read it. The outer grid this fit
  # integrated is `fit$joint_fit`, and auto placement may EXPAND the axis, so
  # what is asserted is the node count it carries rather than the node values.
  alpha_axis_size <- function(fit) {
    tg <- fit$joint_fit$theta_grid
    length(unique(as.numeric(tg[, which(colnames(tg) == "alpha")])))
  }
  dflt <- fit_at(~ pos_cov1 + copy(spatial()), list())
  fine <- fit_at(~ pos_cov1 + copy(spatial()), list(alpha.n = 11))
  expect_gt(alpha_axis_size(fine), alpha_axis_size(dflt))

  # At `alpha.n = 1` the slab collapses to its lower bound, so a truth of
  # alpha = 1.0 is off the axis and cannot be reached -- the knob reaches the
  # likelihood, not just the grid's shape.
  narrow <- fit_at(~ pos_cov1 + copy(spatial()), list(alpha.n = 1))
  expect_lte(narrow$means[["alpha"]], max(engine_alpha_axis(1)))
  expect_gt(dflt$means[["alpha"]], narrow$means[["alpha"]])

  # `copy(alpha = grid(...))` states the nodes, so the two are the same axis
  # written twice.
  expect_error(
    fit_at(~ pos_cov1 + copy(spatial(), alpha = grid(c(0, 0.5, 1))),
           list(alpha.n = 11)),
    "states the copy axis's nodes")
})
