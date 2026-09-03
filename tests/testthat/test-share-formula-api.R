# =============================================================================
# test-share-formula-api.R - the formula-native cross-arm coupling API.
#
# The occu_cover() positive arm carries a scaled copy of the occurrence arm's
# spatial effect, selected structurally (no name needed) by a constructor:
#
#   occurrence = ~ ... + spatial(~ 1 + time || cell_idx, graph = adj)
#   positive   = ~ ... + share(spatial(), alpha = grid(c(...)))
#   control    = list(engine = "joint")
#
# replacing the control-based path (formula =, control$alpha.grid[.trend],
# engine = "joint"), which #295 retired as user surface: those keys are now the
# wire format the formula compiles into, and a fit that sets one is refused.
# =============================================================================


# Build a small occu_cover dataset (optionally with a per-cell trend covariate)
# in the shape the joint path consumes.
.cfa_data <- function(N = 30L, J = 4L, trend = FALSE, seed = 12345L) {
  adj <- chain_adj(N)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal", adj = adj,
    sigma = 0.8, alpha = 1.0,
    trend = trend, sigma_trend = 0.7, alpha_trend = 0.9, seed = seed)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N), cell_idx = seq_len(N)),
                    sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(adj = adj, od = od, cell_dat = cell_dat, y_pos = y_pos, N = N)
}

.cfa_fit <- function(args) suppressWarnings(suppressMessages(do.call(tobs, args)))

# Compare the load-bearing pieces of two fits for byte-identity.
.cfa_expect_identical_fit <- function(a, b) {
  testthat::expect_identical(a$means, b$means)
  testthat::expect_identical(a$sds, b$sds)
  testthat::expect_identical(a$spatial_field, b$spatial_field)
  testthat::expect_identical(a$joint_fit$log_marginal, b$joint_fit$log_marginal)
  if (!is.null(a$trend_field) || !is.null(b$trend_field)) {
    testthat::expect_identical(a$trend_field, b$trend_field)
  }
}


# ---------------------------------------------------------------------------
# Parser-level: the new tokens are recognized and translated correctly.
# ---------------------------------------------------------------------------

test_that("grid() marks a coupling axis for integration; scalar fixes it", {
  g <- tulpaObs:::.tobs_term_grid(c(0.25, 0.5, 1.0))
  expect_s3_class(g, "tobs_alpha_grid")
  expect_identical(g$values, c(0.25, 0.5, 1.0))

  res_g <- tulpaObs:::.tobs_resolve_copy_alpha(g)
  expect_true(res_g$integrate)
  expect_identical(res_g$grid, c(0.25, 0.5, 1.0))

  res_s <- tulpaObs:::.tobs_resolve_copy_alpha(0.5)
  expect_false(res_s$integrate)
  expect_identical(res_s$grid, 0.5)

  # A bare numeric vector must be wrapped in grid() to integrate.
  expect_error(tulpaObs:::.tobs_resolve_copy_alpha(c(0.25, 0.5)), "grid\\(")
})

test_that("share(spatial()) selects the spatial effect with no name", {
  cp <- tulpaObs:::.tobs_term_copy(
    spatial(), alpha = tulpaObs:::.tobs_term_grid(c(0.25, 0.5)))
  expect_s3_class(cp, "tobs_copy")
  expect_identical(cp$selector_type, "spatial")
  expect_null(cp$selector_group)
  expect_null(cp$ref)
  expect_true(cp$alpha_integrate)
  expect_identical(cp$alpha_grid, c(0.25, 0.5))
})

test_that("share(spatial(grp)) carries the grouping-variable discriminator", {
  cp <- tulpaObs:::.tobs_term_copy(spatial(cell_idx), alpha = 0.5)
  expect_identical(cp$selector_type, "spatial")
  expect_identical(cp$selector_group, "cell_idx")
  expect_false(cp$alpha_integrate)
})

test_that("share(spatial(), terms =) resolves a per-component amplitude map", {
  cp <- tulpaObs:::.tobs_term_copy(
    spatial(),
    terms = list(intercept = tulpaObs:::.tobs_term_grid(c(0.25, 0.5)),
                 time      = tulpaObs:::.tobs_term_grid(c(0.5, 1.0))))
  expect_named(cp$copy_terms, c("intercept", "time"))
  expect_true(cp$copy_terms$intercept$integrate)
  expect_identical(cp$copy_terms$time$grid, c(0.5, 1.0))
})

test_that("share(spatial(1)) integer position is rejected", {
  expect_error(tulpaObs:::.tobs_term_copy(spatial(1), alpha = 0.5),
               "integer position")
})

test_that("share(temporal()) non-spatial selector is rejected", {
  expect_error(tulpaObs:::.tobs_term_copy(temporal(), alpha = 0.5),
               "spatial")
})

test_that("share() with both alpha and terms errors", {
  expect_error(
    tulpaObs:::.tobs_term_copy(spatial(), alpha = 0.5,
                               terms = list(intercept = 0.5)),
    "not both")
})

test_that("share(\"name\") explicit string form still parses (back-compat)", {
  cp <- tulpaObs:::.tobs_term_copy(
    "occ_space", alpha = tulpaObs:::.tobs_term_grid(c(0.25, 0.5)))
  expect_identical(cp$ref, "occ_space")
  expect_null(cp$selector_type)
  expect_null(cp$component)
  expect_true(cp$alpha_integrate)

  cp2 <- tulpaObs:::.tobs_term_copy("occ_space.trend", alpha = 0.7)
  expect_identical(cp2$ref, "occ_space")
  expect_identical(cp2$component, "trend")
  expect_false(cp2$alpha_integrate)
})

test_that("share(spatial()) inside a positive formula parses without a name", {
  pf <- ~ pos_cov1 + share(spatial(), alpha = grid(c(0.25, 0.5, 1.0)))
  parsed <- tulpaObs:::.tobs_parse_formula(pf, data = NULL)
  expect_identical(deparse(parsed$fe_formula), "~pos_cov1")
  cp <- Filter(function(t) inherits(t, "tobs_copy"), parsed$terms)[[1L]]
  expect_identical(cp$selector_type, "spatial")
  expect_true(cp$alpha_integrate)
  expect_identical(cp$alpha_grid, c(0.25, 0.5, 1.0))
})

test_that("spatial(name =) is optional; a bar and a single term still take it", {
  adj <- chain_adj(4L)
  dat <- data.frame(cell_idx = 1:4, time.sc = c(-1, 0, 1, 2))
  of <- ~ spatial(~ 1 + time.sc || cell_idx, graph = adj)
  op <- tulpaObs:::.tobs_parse_formula(of, data = dat)
  sp <- Filter(function(t) inherits(t, "tobs_spatial"), op$terms)[[1L]]
  expect_true(isTRUE(sp$is_bar))
  expect_null(sp$field_name)

  on <- ~ spatial(~ 1 + time.sc || cell_idx, graph = adj, name = "occ_space")
  on_p <- tulpaObs:::.tobs_parse_formula(on, data = dat)
  spn <- Filter(function(t) inherits(t, "tobs_spatial"), on_p$terms)[[1L]]
  expect_identical(spn$field_name, "occ_space")
})


# ---------------------------------------------------------------------------
# What the formula compiles to.
# ---------------------------------------------------------------------------

# Each of these fitted the SAME model twice -- once through
# control$alpha.grid[.trend], once through share() -- and asserted the two doors
# agreed. That door is retired (#295), so what is asserted now is the half that
# still carries the claim: the formula compiles to the wire values the control
# key used to state. Read at the compile boundary rather than through a fit,
# which is both exact and free.

.cfa_compiled <- function(pos_formula, occ_formula, data) {
  copies <- tulpaObs:::.occu_cover_extract_pos_copies(pos_formula)$copies
  info   <- tulpaObs:::.occu_cover_spatial_fields(occ_formula, data)
  tulpaObs:::.occu_cover_apply_copy_coupling(copies, info, list())
}

test_that("share(spatial(), alpha = grid()) compiles to the intercept axis", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  g   <- c(0.25, 0.5, 1.0, 1.5)
  ctl <- .cfa_compiled(~ pos_cov1 + share(spatial(), alpha = grid(g)),
                       ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
                       d$cell_dat)
  expect_equal(as.numeric(ctl[["alpha.grid"]]), g)
})

test_that("omitting share() compiles to a pinned, decoupled axis", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  ctl <- .cfa_compiled(~ pos_cov1,
                       ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
                       d$cell_dat)
  expect_equal(as.numeric(ctl[["alpha.grid"]]), 0)
})

test_that("a scalar alpha compiles to a one-node axis", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  ctl <- .cfa_compiled(~ pos_cov1 + share(spatial(), alpha = 0.5),
                       ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
                       d$cell_dat)
  expect_equal(as.numeric(ctl[["alpha.grid"]]), 0.5)
})

test_that("per-component terms= compile to one axis per block", {
  d   <- .cfa_data(N = 20L, trend = TRUE, seed = 31337L)
  adj <- d$adj
  gi  <- c(0.25, 0.5, 1.0)
  gt  <- c(0.5, 1.0, 1.5)
  ctl <- .cfa_compiled(
    ~ pos_cov1 + share(spatial(),
                       terms = list(intercept = grid(gi), time = grid(gt))),
    ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
    d$cell_dat)
  expect_equal(as.numeric(ctl[["alpha.grid"]]), gi)
  expect_equal(as.numeric(ctl[["alpha.grid.trend"]]), gt)
})

test_that("whole-field share() compiles the same amplitude onto every block", {
  d   <- .cfa_data(N = 20L, trend = TRUE, seed = 31337L)
  adj <- d$adj
  gw  <- c(0.25, 0.5, 1.0)
  ctl <- .cfa_compiled(~ pos_cov1 + share(spatial(), alpha = grid(gw)),
                       ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
                       d$cell_dat)
  expect_equal(as.numeric(ctl[["alpha.grid"]]), gw)
  expect_equal(as.numeric(ctl[["alpha.grid.trend"]]), gw)
})

# One end-to-end fit still runs the compiled axis through the engine, so the
# compile assertions above are anchored to a model that actually fits.
test_that("the compiled axis is the axis the fit integrates", {
  skip_if_fast()
  d   <- .cfa_data()
  adj <- d$adj
  g   <- c(0.25, 0.5, 1.0, 1.5)
  fit <- .cfa_fit(list(
    occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
    positive   = ~ pos_cov1 + share(spatial(), alpha = grid(g)),
    data = d$cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
    visits = d$od$det.covs, method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L, engine = "joint",
                   adaptive.grid = FALSE,
                   var.of.means.consistency = FALSE)))
  # A genuinely fixed outer grid takes BOTH knobs: refinement reaches a stated
  # axis from the adaptive passes AND from the var-of-means consistency pass,
  # and the second runs whatever `adaptive.grid` is. With only the first off,
  # this fit integrates c(0.25, 0.5, 0.7303, 1, 1.5) -- inside the stated range,
  # so the axis is still bounded by what the formula wrote, but not equal to it.
  # theta_grid labels the copy axis "alpha" on a single-block fit and
  # "b<k>.alpha" on a multi-block one.
  tg    <- fit$joint_fit$theta_grid %||% fit$joint$theta_grid
  col   <- if ("b1.alpha" %in% colnames(tg)) "b1.alpha" else "alpha"
  nodes <- sort(unique(as.numeric(tg[, col])))
  expect_equal(nodes, sort(g))
})


# ---------------------------------------------------------------------------
# Guard rails.
# ---------------------------------------------------------------------------

test_that("control$alpha.grid is refused: it is not user surface (#295)", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + share(spatial(), alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint", alpha.grid = c(0.5, 1))))),
    "not user surface")
})

test_that("share(spatial(grp)) on a mismatched grouping variable errors", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + share(spatial(region_idx), alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "grouped on")
})

test_that("share(spatial(), terms =) must address every field block", {
  d   <- .cfa_data(N = 20L, trend = TRUE)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + share(spatial(), terms = list(intercept = grid(c(0.5, 1)))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "every field block")
})

test_that("share(string) on the positive arm errors (selector required)", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + share("occ_space", alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "not a string")
})

test_that("the positive-arm strip evaluates only share(), leaving the rest (#296)", {
  f  <- tulpaObs:::.occu_cover_extract_pos_copies
  dp <- function(x) paste(deparse(x$formula), collapse = "")

  # A structured term naming a data column used to die inside the strip
  # ("object year not found"), because the whole formula was parsed against
  # `data = NULL`. It has to survive the strip so the ARM rejector can name the
  # unsupported feature instead.
  r <- f(~ x_cov + temporal(year))
  expect_length(r$copies, 0L)
  expect_true(grepl("temporal(year)", dp(r), fixed = TRUE))

  r2 <- f(~ x_cov + re(cell))
  expect_length(r2$copies, 0L)
  expect_true(grepl("re(cell)", dp(r2), fixed = TRUE))

  # share() is still stripped, and composes with a term left in place.
  r3 <- f(~ x_cov + share(spatial()) + re(cell))
  expect_length(r3$copies, 1L)
  expect_s3_class(r3$copies[[1L]], "tobs_copy")
  expect_true(grepl("re(cell)", dp(r3), fixed = TRUE))
  expect_false(grepl("share(", dp(r3), fixed = TRUE))

  # The rebuilt fixed-effects formula keeps the intercept term and bar syntax.
  expect_true(grepl("- 1", dp(f(~ 0 + x_cov)), fixed = TRUE))
  expect_true(grepl("(1 | g)", dp(f(~ x_cov + (1 | g))), fixed = TRUE))
})

test_that("a structured term on the positive arm reaches the arm rejector (#296)", {
  skip_on_cran()
  skip_if_fast()
  d <- .cfa_data(N = 20L)
  # `visit` is not a column of the cell-level frame, so the pre-#296 strip
  # reported it as a missing object and pointed the user at their data.
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + icar(graph = d$adj), data = d$cell_dat,
      family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + temporal(visit),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "not supported on the positive cover arm")
})

test_that("a stated alpha grid states NODES, and the span reaches past them (#301)", {
  skip_on_cran()
  skip_if_fast()
  d <- .cfa_data(N = 24L)
  fit_at <- function(nodes, ...) suppressWarnings(suppressMessages(tobs(
    occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = d$adj),
    data = d$cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1,
    positive = eval(bquote(
      ~ pos_cov1 + share(spatial(), alpha = grid(.(nodes))))),
    y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
    method = "nested_laplace",
    control = c(list(verbose = FALSE, engine = "joint",
                     sigma.grid = c(0.5, 1)), list(...)))))
  span <- function(nodes, ...) {
    b <- tulpaObs:::.occu_cover_nuts_hyper_bounds(fit_at(nodes, ...), "icar")
    b$alpha
  }
  ratio <- function(s, lo, hi) (log(s[2L]) - log(s[1L])) / (log(hi) - log(lo))
  nodes9 <- exp(seq(log(0.2), log(0.5), length.out = 9))

  # With refinement OFF the stated nodes ARE the axis, and the only widening is
  # the half node step each end: each node represents the cell around it, so the
  # weights sum over a whole axis rather than stopping on the last node. For k
  # equally log-spaced nodes that is an EXACT identity in the axis's own log
  # coordinate -- (k - 1) steps between the nodes plus a half step at each end
  # -- so the span is k / (k - 1) times the stated range: 2x at two nodes,
  # 1.125x at nine. Asserted as the identity rather than to a tolerance, since
  # it is arithmetic and not an estimate.
  s2 <- span(c(0.2, 0.5), adaptive.grid = FALSE)
  expect_lt(s2[1L], 0.2); expect_gt(s2[2L], 0.5)
  expect_equal(ratio(s2, 0.2, 0.5), 2 / 1, tolerance = 1e-6)
  s9 <- span(nodes9, adaptive.grid = FALSE)
  expect_equal(ratio(s9, 0.2, 0.5), 9 / 8, tolerance = 1e-6)
  expect_lt(ratio(s9, 0.2, 0.5), ratio(s2, 0.2, 0.5))

  # With refinement ON -- the DEFAULT -- the stated nodes BOUND the axis: cells
  # are placed inside the range and never past its ends, so the node range is
  # what the formula said, to the bit. What refinement moves is the density.
  nodes_of <- function(...) {
    tg <- fit_at(...)$joint_fit$theta_grid
    sort(unique(as.numeric(tg[, match("alpha", colnames(tg))])))
  }
  n2a <- nodes_of(c(0.2, 0.5))
  expect_equal(range(n2a), c(0.2, 0.5))
  expect_gt(length(n2a), 2L)

  # Denser cells make the outer half step -- half of the outermost interval --
  # smaller, so the span tightens TOWARD the stated range from both sides. It
  # never reaches it: a node always represents the cell around it, so some
  # overhang survives however fine the grid gets.
  s2a <- span(c(0.2, 0.5))
  expect_lt(ratio(s2a, 0.2, 0.5), ratio(s2, 0.2, 0.5))
  expect_gt(ratio(s2a, 0.2, 0.5), 1)
  expect_gt(s2a[1L], s2[1L])
  expect_lt(s2a[2L], s2[2L])
})

test_that("the retired copy() spelling names its replacement", {
  d <- data.frame(y = c(0, 1, 1, 0), x = rnorm(4), cell_idx = 1:4)
  # Every scanner that resolves a call head against the term registry: the
  # shared parser, and the per-arm scanners the cover families run.
  expect_error(tulpaObs:::.tobs_parse_formula(~ x + copy(spatial()), data = d),
               "copy\\(\\) is now share\\(\\)")
  expect_error(tulpaObs:::.tobs_check_retired_term("copy"), "share\\(spatial")
  # A live term and an ordinary fixed effect are untouched by the check.
  expect_null(tulpaObs:::.tobs_check_retired_term("share"))
  expect_null(tulpaObs:::.tobs_check_retired_term(NA_character_))
  expect_null(tulpaObs:::.tobs_check_retired_term("log"))
})

test_that("giving both occurrence and formula errors", {
  expect_error(
    tobs(formula = ~ a, occurrence = ~ b, family = occu_cover("lognormal"),
         y = 1, data = data.frame(a = 1, b = 1)),
    "not both")
})
