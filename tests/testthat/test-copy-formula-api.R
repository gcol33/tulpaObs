# =============================================================================
# test-copy-formula-api.R - the formula-native cross-arm coupling API.
#
# The occu_cover() positive arm carries a scaled copy of the occurrence arm's
# spatial effect, selected structurally (no name needed) by a constructor:
#
#   occurrence = ~ ... + spatial(~ 1 + time || cell_idx, graph = adj)
#   positive   = ~ ... + copy(spatial(), alpha = grid(c(...)))
#   control    = list(engine = "joint")
#
# replacing the control-based path (formula =, control$alpha.grid[.trend],
# engine = "joint"). A pure API change moves no number: the equivalence
# tests below fit the SAME model both ways and assert byte-identical results.
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

test_that("copy(spatial()) selects the spatial effect with no name", {
  cp <- tulpaObs:::.tobs_term_copy(
    spatial(), alpha = tulpaObs:::.tobs_term_grid(c(0.25, 0.5)))
  expect_s3_class(cp, "tobs_copy")
  expect_identical(cp$selector_type, "spatial")
  expect_null(cp$selector_group)
  expect_null(cp$ref)
  expect_true(cp$alpha_integrate)
  expect_identical(cp$alpha_grid, c(0.25, 0.5))
})

test_that("copy(spatial(grp)) carries the grouping-variable discriminator", {
  cp <- tulpaObs:::.tobs_term_copy(spatial(cell_idx), alpha = 0.5)
  expect_identical(cp$selector_type, "spatial")
  expect_identical(cp$selector_group, "cell_idx")
  expect_false(cp$alpha_integrate)
})

test_that("copy(spatial(), terms =) resolves a per-component amplitude map", {
  cp <- tulpaObs:::.tobs_term_copy(
    spatial(),
    terms = list(intercept = tulpaObs:::.tobs_term_grid(c(0.25, 0.5)),
                 time      = tulpaObs:::.tobs_term_grid(c(0.5, 1.0))))
  expect_named(cp$copy_terms, c("intercept", "time"))
  expect_true(cp$copy_terms$intercept$integrate)
  expect_identical(cp$copy_terms$time$grid, c(0.5, 1.0))
})

test_that("copy(spatial(1)) integer position is rejected", {
  expect_error(tulpaObs:::.tobs_term_copy(spatial(1), alpha = 0.5),
               "integer position")
})

test_that("copy(temporal()) non-spatial selector is rejected", {
  expect_error(tulpaObs:::.tobs_term_copy(temporal(), alpha = 0.5),
               "spatial")
})

test_that("copy() with both alpha and terms errors", {
  expect_error(
    tulpaObs:::.tobs_term_copy(spatial(), alpha = 0.5,
                               terms = list(intercept = 0.5)),
    "not both")
})

test_that("copy(\"name\") explicit string form still parses (back-compat)", {
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

test_that("copy(spatial()) inside a positive formula parses without a name", {
  pf <- ~ pos_cov1 + copy(spatial(), alpha = grid(c(0.25, 0.5, 1.0)))
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
# Equivalence: OLD control path == NEW formula-native path (no number moves).
# ---------------------------------------------------------------------------

test_that("OLD (control alpha.grid) == NEW (copy(spatial()) grid) intercept field", {
  skip_if_fast()
  d   <- .cfa_data()
  adj <- d$adj
  g   <- c(0.25, 0.5, 1.0, 1.5)
  base <- list(data = d$cell_dat, family = occu_cover("lognormal"),
               detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
               visits = d$od$det.covs, method = "nested_laplace")

  old <- .cfa_fit(c(base, list(
    formula  = ~ occ_cov1 + icar(graph = adj),
    positive = ~ pos_cov1,
    control  = list(verbose = FALSE, max.iter = 500L,
                    engine = "joint", alpha.grid = g))))
  new <- .cfa_fit(c(base, list(
    occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
    positive   = ~ pos_cov1 + copy(spatial(), alpha = grid(g)),
    control    = list(verbose = FALSE, max.iter = 500L, engine = "joint"))))

  expect_identical(names(old$means), names(new$means))
  .cfa_expect_identical_fit(old, new)
})

test_that("decouple: OLD alpha.grid = 0 == NEW omitting copy()", {
  skip_if_fast()
  d   <- .cfa_data()
  adj <- d$adj
  base <- list(data = d$cell_dat, family = occu_cover("lognormal"),
               detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
               visits = d$od$det.covs, method = "nested_laplace")

  old <- .cfa_fit(c(base, list(
    formula  = ~ occ_cov1 + icar(graph = adj),
    positive = ~ pos_cov1,
    control  = list(verbose = FALSE, max.iter = 500L,
                    engine = "joint", alpha.grid = 0))))
  new <- .cfa_fit(c(base, list(
    occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
    positive   = ~ pos_cov1,   # no copy() => decoupled
    control    = list(verbose = FALSE, max.iter = 500L, engine = "joint"))))

  .cfa_expect_identical_fit(old, new)
})

test_that("scalar alpha fixes the amplitude (OLD alpha.grid = 0.5 == NEW copy(spatial(), alpha = 0.5))", {
  skip_if_fast()
  d   <- .cfa_data()
  adj <- d$adj
  base <- list(data = d$cell_dat, family = occu_cover("lognormal"),
               detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
               visits = d$od$det.covs, method = "nested_laplace")

  old <- .cfa_fit(c(base, list(
    formula  = ~ occ_cov1 + icar(graph = adj),
    positive = ~ pos_cov1,
    control  = list(verbose = FALSE, max.iter = 500L,
                    engine = "joint", alpha.grid = 0.5))))
  new <- .cfa_fit(c(base, list(
    occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
    positive   = ~ pos_cov1 + copy(spatial(), alpha = 0.5),
    control    = list(verbose = FALSE, max.iter = 500L, engine = "joint"))))

  .cfa_expect_identical_fit(old, new)
})

test_that("per-component: OLD alpha.grid vs alpha.grid.trend == NEW copy(spatial(), terms =)", {
  skip_if_fast()
  d   <- .cfa_data(trend = TRUE, seed = 31337L)
  adj <- d$adj
  gi  <- c(0.25, 0.5, 1.0)
  gt  <- c(0.5, 1.0, 1.5)
  base <- list(data = d$cell_dat, family = occu_cover("lognormal"),
               detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
               visits = d$od$det.covs, method = "nested_laplace")

  old <- .cfa_fit(c(base, list(
    formula  = ~ occ_cov1 + icar(graph = adj) + icar(graph = adj, weight = time),
    positive = ~ pos_cov1,
    control  = list(verbose = FALSE, max.iter = 300L, engine = "joint",
                    alpha.grid = gi, alpha.grid.trend = gt))))
  new <- .cfa_fit(c(base, list(
    occurrence = ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
    positive   = ~ pos_cov1 +
      copy(spatial(), terms = list(intercept = grid(gi), time = grid(gt))),
    control    = list(verbose = FALSE, max.iter = 300L, engine = "joint"))))

  .cfa_expect_identical_fit(old, new)
})

test_that("whole-field copy(spatial()) scales every block with one amplitude", {
  skip_if_fast()
  d   <- .cfa_data(trend = TRUE, seed = 31337L)
  adj <- d$adj
  gw  <- c(0.25, 0.5, 1.0)
  base <- list(data = d$cell_dat, family = occu_cover("lognormal"),
               detection = ~ det_cov1, y = d$od$y, y_pos = d$y_pos,
               visits = d$od$det.covs, method = "nested_laplace")

  old <- .cfa_fit(c(base, list(
    formula  = ~ occ_cov1 + icar(graph = adj) + icar(graph = adj, weight = time),
    positive = ~ pos_cov1,
    control  = list(verbose = FALSE, max.iter = 300L, engine = "joint",
                    alpha.grid = gw))))   # alpha.grid.trend defaults to alpha.grid
  new <- .cfa_fit(c(base, list(
    occurrence = ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
    positive   = ~ pos_cov1 + copy(spatial(), alpha = grid(gw)),
    control    = list(verbose = FALSE, max.iter = 300L, engine = "joint"))))

  .cfa_expect_identical_fit(old, new)
})


# ---------------------------------------------------------------------------
# Guard rails.
# ---------------------------------------------------------------------------

test_that("copy() and control$alpha.grid together is an error", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + copy(spatial(), alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint", alpha.grid = c(0.5, 1))))),
    "not both")
})

test_that("copy(spatial(grp)) on a mismatched grouping variable errors", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + copy(spatial(region_idx), alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "grouped on")
})

test_that("copy(spatial(), terms =) must address every field block", {
  d   <- .cfa_data(N = 20L, trend = TRUE)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 + time || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + copy(spatial(), terms = list(intercept = grid(c(0.5, 1)))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "every field block")
})

test_that("copy(string) on the positive arm errors (selector required)", {
  d   <- .cfa_data(N = 20L)
  adj <- d$adj
  expect_error(
    suppressWarnings(suppressMessages(tobs(
      occurrence = ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj),
      data = d$cell_dat, family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1 + copy("occ_space", alpha = grid(c(0.5, 1))),
      y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, engine = "joint")))),
    "not a string")
})

test_that("giving both occurrence and formula errors", {
  expect_error(
    tobs(formula = ~ a, occurrence = ~ b, family = occu_cover("lognormal"),
         y = 1, data = data.frame(a = 1, b = 1)),
    "not both")
})
