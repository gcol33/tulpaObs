# The copy coefficient's outer axis is set by STATING its nodes or by naming a
# RESOLUTION, and the two are different requests.
#
# The axis carries prior structure -- a point mass at alpha = 0 ("no coupling")
# and a log-spaced slab above it -- so `share(alpha = grid(...))` (and
# `control$alpha.grid`) restate that structure along with the nodes, which is
# why a fit that only wants the axis integrated more finely cannot use them.
# `share(alpha = grid(n = ))` (and `control$alpha.n`) names the number of slab
# nodes and the ENGINE re-reads its own axis at that resolution, atom and slab
# bounds unchanged.
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

  # Stated nodes win over a resolution: a field block with no share() arrives
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

test_that("a share() stating an amplitude and a resolution knob are refused together", {
  states <- tulpaObs:::.tobs_copy_states_amplitude
  expect_true(states(list(alpha_grid = c(0, 1), alpha_integrate = TRUE)))
  expect_true(states(list(alpha_grid = 0.8, alpha_integrate = FALSE)))
  # grid(n = ) states the axis too -- as a resolution rather than as nodes --
  # so it collides with the control spelling of the same request.
  expect_true(states(list(alpha_grid = NULL, alpha_n = 9L,
                          alpha_integrate = TRUE)))
  # share(spatial()) with no amplitude asks for the DEFAULT axis, so it composes
  # with a resolution instead of colliding with it.
  expect_false(states(list(alpha_grid = NULL, alpha_integrate = NA)))

  expect_error(
    tulpaObs:::.tobs_check_alpha_copy(TRUE, list(alpha.n = 9), "occu_cover()"),
    "states the copy axis")
  expect_true(
    tulpaObs:::.tobs_check_alpha_copy(FALSE, list(alpha.n = 9), "occu_cover()"))
})


# ---- the resolution is written in the formula -----------------------------

test_that("grid() states nodes OR a resolution, never both and never neither", {
  g <- tulpaObs:::.tobs_terms$grid
  expect_equal(g(c(0.25, 0.5, 1))$values, c(0.25, 0.5, 1))
  expect_null(g(c(0.25, 0.5, 1))$n)
  expect_identical(g(n = 9)$n, 9L)
  expect_null(g(n = 9)$values)
  expect_error(g(c(1, 2), n = 3), "not both")
  expect_error(g(), "grid\\(n = 9\\)")
  expect_error(g(n = 0), "integer >= 1")
})

test_that("share(alpha = grid(n = )) reaches the engine's own axis at that n", {
  reg <- list2env(tulpaObs:::.tobs_terms, parent = environment())
  cp  <- function(e) eval(e, reg)
  fv  <- list(fields = list(list(component = "intercept")),
              group_var = "cell_idx")
  app <- tulpaObs:::.occu_cover_apply_copy_coupling

  ctrl <- app(list(cp(quote(share(spatial(), alpha = grid(n = 9))))), fv, list())
  expect_identical(ctrl[["alpha.n"]], 9L)
  expect_null(ctrl[["alpha.grid"]])
  # And it resolves to the engine's axis at that resolution, not to nodes of
  # our own: the axis a fit integrates is the one `alpha.n = 9` names.
  expect_equal(
    as.numeric(tulpaObs:::.tobs_alpha_nodes(
      tulpaObs:::.tobs_alpha_axis_base(ctrl))),
    engine_alpha_axis(9))

  # Nodes and a resolution are per BLOCK, so one field can take one of each.
  fv2 <- list(fields = list(list(component = "intercept"),
                            list(component = "time.sc")),
              group_var = "cell_idx")
  ctrl2 <- app(list(cp(quote(share(spatial(), terms = list(
    intercept = grid(n = 7), time.sc = grid(c(0, 1))))))), fv2, list())
  expect_identical(ctrl2[["alpha.n"]], 7L)
  expect_null(ctrl2[["alpha.grid"]])
  expect_equal(as.numeric(ctrl2[["alpha.grid.trend"]]), c(0, 1))
  expect_null(ctrl2[["alpha.n.trend"]])
})

test_that("a bare share() still composes with the control resolution", {
  # The formula spelling does not overwrite what it did not state: a share()
  # naming no amplitude asks for the engine's default axis, which is what the
  # resolution knob then reads more finely.
  reg <- list2env(tulpaObs:::.tobs_terms, parent = environment())
  fv  <- list(fields = list(list(component = "intercept")),
              group_var = "cell_idx")
  ctrl <- tulpaObs:::.occu_cover_apply_copy_coupling(
    list(eval(quote(share(spatial())), reg)), fv, list(alpha.n = 5))
  expect_equal(ctrl[["alpha.n"]], 5)
  expect_null(ctrl[["alpha.grid"]])
})


# ---- the prior on the copy coefficient ------------------------------------

test_that("share(prior = ) carries the hyperprior, and is checked where written", {
  reg <- list2env(tulpaObs:::.tobs_terms, parent = environment())
  cp  <- function(e) eval(e, reg)
  pc  <- list("pc.prec", c(4, 0.01))
  fv  <- list(fields = list(list(component = "intercept")),
              group_var = "cell_idx")
  fv2 <- list(fields = list(list(component = "intercept"),
                            list(component = "time.sc")),
              group_var = "cell_idx")
  app <- tulpaObs:::.occu_cover_apply_copy_coupling

  expect_identical(
    app(list(cp(quote(share(spatial(), prior = pc)))), fv, list())[["prior.alpha"]],
    pc)
  # Shape is checked at the call site that writes it; the families and their
  # parameter ranges belong to the engine.
  expect_error(cp(quote(share(spatial(), prior = c(4, 0.01)))), "list\\(<family>")
  # One prior reaches the engine per fit, so the two spellings are exclusive.
  expect_error(
    app(list(cp(quote(share(spatial(), prior = pc)))), fv,
        list(prior.alpha = list("half_normal", 2))),
    "not both")
  # And the engine bakes it on the FIRST copied block (gcol33/tulpa#655), so a
  # fit copying two is refused rather than given it on one of them.
  expect_error(app(list(cp(quote(share(spatial(), prior = pc)))), fv2, list()),
               "copies 2")
  # A decoupled fit has no coefficient to regularize.
  expect_null(app(list(), fv, list())[["prior.alpha"]])
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

test_that("a bare share() takes the resolution; a share() with nodes refuses it", {
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

  # `share(spatial())` names no amplitude: it asks for the engine's own axis, and
  # the resolution says how finely to read it. The outer grid this fit
  # integrated is `fit$joint_fit`, and auto placement may EXPAND the axis, so
  # what is asserted is the node count it carries rather than the node values.
  alpha_axis_size <- function(fit) {
    tg <- fit$joint_fit$theta_grid
    length(unique(as.numeric(tg[, which(colnames(tg) == "alpha")])))
  }
  dflt <- fit_at(~ pos_cov1 + share(spatial()), list())
  fine <- fit_at(~ pos_cov1 + share(spatial()), list(alpha.n = 11))
  expect_gt(alpha_axis_size(fine), alpha_axis_size(dflt))

  # At `alpha.n = 1` the slab collapses to its lower bound, so a truth of
  # alpha = 1.0 is off the axis and cannot be reached -- the knob reaches the
  # likelihood, not just the grid's shape.
  narrow <- fit_at(~ pos_cov1 + share(spatial()), list(alpha.n = 1))
  expect_lte(narrow$means[["alpha"]], max(engine_alpha_axis(1)))
  expect_gt(dflt$means[["alpha"]], narrow$means[["alpha"]])

  # `share(alpha = grid(...))` states the nodes, so the two are the same axis
  # written twice.
  expect_error(
    fit_at(~ pos_cov1 + share(spatial(), alpha = grid(c(0, 0.5, 1))),
           list(alpha.n = 11)),
    "share(alpha = ) states the copy axis", fixed = TRUE)
})


# --- `alpha` exists only when there is a copy ---------------------------------
#
# `.occu_cover_apply_copy_coupling()` spells "no share() named this block" as an
# amplitude axis stating the single node 0, which is the decoupled model. That
# is not a parameter, and reporting it as one left an all-zero row and column in
# vcov() (rank-deficient, so solve() and any joint Wald failed) and one too many
# in n_params.
#
# The two engine paths remove it differently, and the difference is deliberate:
# the single-block backend takes a numeric `field_coef = 0` on the cover arm and
# builds no amplitude axis at all, while the multi-block driver offers the SD
# parameterization ONLY to a copied block, so a decoupled field there keeps its
# copy spec (pinned at 0) and the axis is suppressed at reporting instead.

.aq_adj <- function(n = 8L) rook_adj(n)

.aq_sim <- function(seed = 5L, pos_field = FALSE, alpha = 1.0) {
  adj <- .aq_adj(); N <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 6L, positive = "lognormal",
    beta_occ = c(qlogis(0.7), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3,
    adj = adj, sigma = 0.5, alpha = alpha,
    pos_field = pos_field, sigma_pos_int = 0.6, sigma_pos_trend = 0.0,
    seed = seed)
  if (is.null(sim$data$cell)) sim$data$cell <- seq_len(N)
  list(sim = sim, adj = adj)
}

.aq_fit <- function(f, positive, ...) {
  suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = f$adj, group_var = "cell"),
    detection = ~ 1, positive = positive,
    family = occu_cover(response = "lognormal"),
    data = f$sim$data, y = f$sim$y, y_pos = f$sim$y_pos,
    method = "nested_laplace",
    control = c(list(progress = FALSE), list(...))))
}

test_that("a fit with no share() reports no alpha and has a full-rank vcov", {
  skip_if_fast(); skip_on_cran()
  f   <- .aq_sim()
  fit <- .aq_fit(f, ~ 1)

  expect_false("alpha" %in% names(fit$means))
  expect_false("alpha" %in% colnames(fit$vcov))
  expect_false("alpha" %in% names(fit$sds))
  expect_equal(fit$n_params, length(fit$means))
  # The defect this guards: an all-zero row/column left vcov() singular.
  expect_equal(qr(fit$vcov)$rank, ncol(fit$vcov))
  expect_no_error(solve(fit$vcov))
  # Single-block path: the amplitude axis is not built at all.
  expect_false("alpha" %in% fit$joint_fit$theta_names)
  # The cover arm genuinely carries none of the field.
  expect_true(all(fit$field_pos == 0))
})

test_that("a share() restores the amplitude and it is estimated", {
  skip_if_fast(); skip_on_cran()
  f   <- .aq_sim()
  fit <- .aq_fit(f, ~ 1 + share(spatial()))

  expect_true("alpha" %in% names(fit$means))
  expect_gt(fit$means[["alpha"]], 0)
  expect_equal(qr(fit$vcov)$rank, ncol(fit$vcov))
})

test_that("a decoupled multi-block fit keeps the field's SD axis", {
  skip_if_fast(); skip_on_cran()
  # REGRESSION GUARD. Dropping the copy spec here would move the shared field
  # from `b1.sigma` onto its precision axis `b1.tau`. Those are the same implied
  # SDs but a different prior measure on the field's scale, and it is not
  # cosmetic: it pulled the measured `sigma_pos_field` recovery from a median of
  # ~0.55 to ~0.42 against a truth of 0.6. The copy stays; only the reporting of
  # its pinned axis is dropped.
  f   <- .aq_sim(pos_field = TRUE, alpha = 0.0)
  fit <- .aq_fit(f, ~ 1 + spatial(~ 1 || cell, graph = f$adj),
                 integration = "ccd")

  tn <- fit$joint_fit$theta_names
  expect_true("b1.sigma" %in% tn)
  expect_false("b1.tau" %in% tn)
  # The pinned amplitude axis is still integrated, and still not reported.
  expect_true("b1.alpha" %in% tn)
  expect_false("alpha" %in% names(fit$means))
  expect_equal(qr(fit$vcov)$rank, ncol(fit$vcov))
})
