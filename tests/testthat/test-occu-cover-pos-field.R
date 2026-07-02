# =============================================================================
# test-occu-cover-pos-field.R - arm-specific spatial field on the cover
# (positive) arm of occu_cover() (gcol33/tulpaObs#110).
#
# `spatial(~ 1 + w || cell, graph = adj, to = "positive")` on the occurrence
# formula adds an INDEPENDENT, non-copied areal field on the cover arm alone --
# decoupled from the occupancy field's alpha copy. This is the opt-in for a
# spatially-structured cover trend that is not a scalar multiple of the occupancy
# field (the alpha copy collapses to a global slope when the two shapes differ,
# flattening delta_cover_cond). The occupancy field still drives psi and, via the
# alpha copy, delta_cover_exp.
# =============================================================================

.pf_grid_adj <- function(side) {
  N <- side * side
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    if (r > 1L)   adj[idx(r, c), idx(r - 1L, c)] <- 1L
    if (r < side) adj[idx(r, c), idx(r + 1L, c)] <- 1L
    if (c > 1L)   adj[idx(r, c), idx(r, c - 1L)] <- 1L
    if (c < side) adj[idx(r, c), idx(r, c + 1L)] <- 1L
  }
  adj
}

.pf_fit <- function(sim, occurrence, integration = "ccd") {
  suppressWarnings(tobs(
    occurrence = occurrence,
    detection  = ~ 1,
    family     = occu_cover(response = "lognormal"),
    data       = sim$data, y = sim$y, y_pos = sim$y_pos,
    method     = "nested_laplace",
    control    = list(progress = FALSE, integration = integration)))
}


# --- parse-time: the to = "positive" bar splits out as an arm-specific field ---

test_that("to = \"positive\" bar resolves to an arm-specific cover field", {
  adj  <- .pf_grid_adj(4L)
  n    <- nrow(adj)
  data <- data.frame(cell = seq_len(n), occ_cov1 = rnorm(n),
                     time = as.numeric(scale(rnorm(n))))
  f <- psi ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
       spatial(~ 1 + time || cell, graph = adj, to = "positive")
  si <- .occu_cover_spatial_fields(f, data)

  expect_length(si$fields, 1L)                       # the occupancy field only
  expect_false(is.null(si$pos_armspec))
  expect_identical(si$pos_armspec$type, "icar")
  expect_length(si$pos_armspec$fields, 2L)           # intercept + time trend
  cols <- vapply(si$pos_armspec$fields, function(x) x$column_name, "")
  expect_true("time" %in% cols)
  expect_true(any(vapply(si$pos_armspec$fields,
                         function(x) isTRUE(x$is_intercept), logical(1))))
})

test_that("a single-arm to = \"presence\" spatial bar is rejected", {
  adj  <- .pf_grid_adj(4L)
  n    <- nrow(adj)
  data <- data.frame(cell = seq_len(n), occ_cov1 = rnorm(n))
  f <- psi ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
       spatial(~ 1 || cell, graph = adj, to = "presence")
  expect_error(.occu_cover_spatial_fields(f, data),
               "no separate presence arm")
})

test_that("a detection-arm spatial bar resolves onto the detection (p) arm", {
  adj  <- .pf_grid_adj(4L)
  n    <- nrow(adj)
  data <- data.frame(cell = seq_len(n), occ_cov1 = rnorm(n))
  f <- psi ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
       spatial(~ 1 || cell, graph = adj, to = "detection")
  si <- .occu_cover_spatial_fields(f, data)
  expect_false(is.null(si$armspec[["p"]]))
  expect_true(is.null(si$armspec[["pos"]]))
})

test_that("detection-arm field recovers once the substrate scatters onto p", {
  # The parse -> block -> per-arm-sigma plumbing is arm-generic; the detection arm
  # carries the non-copied field block with field_coef = 1 (the shared field is
  # kept off detection by the spatial_idx = 0 sentinel), the same mechanism the
  # detection RE uses (gcol33/tulpa#140, gcol33/tulpaObs#102).
  skip_if_fast()
  skip_on_cran()
  adj <- .pf_grid_adj(8L); N <- nrow(adj); truth <- 0.7
  rec <- vapply(1:6, function(s) {
    sim <- simulate_occu_cover(
      N = N, J = 8L, positive = "lognormal",
      beta_occ = c(qlogis(0.7), 0.3), beta_p = c(qlogis(0.6), 0.1),
      beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3, adj = adj, sigma = 0.5,
      alpha = 0.0, det_field = TRUE, sigma_p_int = 0.0, sigma_p_trend = truth,
      seed = s)
    fit <- suppressWarnings(tobs(
      occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell"),
      detection  = ~ 1 + spatial(~ 0 + time || cell, graph = adj),
      positive   = ~ 1, family = occu_cover(response = "lognormal"),
      data = sim$data, y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
      control = list(progress = FALSE, integration = "ccd")))
    nm <- grep("^sigma_p_field", names(fit$means), value = TRUE)[1L]
    fit$means[[nm]]
  }, numeric(1))
  expect_lt(abs(stats::median(rec) - truth), 0.25)
})

test_that("an arm-specific cover field does not compose with the `|` MCAR field", {
  adj  <- .pf_grid_adj(4L)
  n    <- nrow(adj)
  data <- data.frame(cell = seq_len(n), occ_cov1 = rnorm(n),
                     time = as.numeric(scale(rnorm(n))))
  # A correlated `|` bar must be the only spatial term, so an extra to="positive"
  # bar alongside it errors at parse.
  f <- psi ~ occ_cov1 + spatial(~ 1 + time | cell, graph = adj) +
       spatial(~ 1 || cell, graph = adj, to = "positive")
  expect_error(.occu_cover_spatial_fields(f, data))
})


# --- smoke fit: structure + a NON-CONSTANT delta_cover_cond ------------------

test_that("occu_cover cover-arm field fits and yields a non-constant delta_cover_cond", {
  skip_if_fast()
  skip_on_cran()

  adj <- .pf_grid_adj(8L)
  N   <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 5L, positive = "lognormal",
    beta_occ = c(qlogis(0.6), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3,
    adj = adj, sigma = 0.5, alpha = 0.0,
    pos_field = TRUE, sigma_pos_int = 0.0, sigma_pos_trend = 0.7, seed = 42L)

  fit <- .pf_fit(sim, ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
                        spatial(~ 0 + time || cell, graph = adj, to = "positive"))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")
  # The occupancy field (sigma, alpha) plus the independent cover-arm field SD.
  expect_true(all(c("sigma", "alpha") %in% names(fit$means)))
  sig_nm <- grep("^sigma_pos_field", names(fit$means), value = TRUE)
  expect_length(sig_nm, 1L)
  expect_true(is.finite(fit$means[[sig_nm]]) && fit$means[[sig_nm]] > 0)

  # The cover-arm field posterior is surfaced separately from the occupancy field.
  expect_length(fit$pos_field, N)
  expect_true(all(is.finite(fit$pos_field)))
  expect_s3_class(fit$pos_field_table, "data.frame")

  # A change over the time covariate: delta_cover_cond must VARY across cells (the
  # core acceptance criterion -- an alpha-copy-only fit gives one distinct value).
  nd <- sim$data; nd$time <- 0
  pr <- predict(fit, newdata = nd, type = "change", times = c(-1, 1),
                time_col = "time", nsim = 400L)
  dcc <- pr$delta_cover_cond
  expect_gt(length(unique(round(dcc, 6))), N %/% 2L)
  expect_gt(stats::sd(dcc), 0.02)
})


# --- recovery: the cover-arm trend field SD is recovered ---------------------

test_that("occu_cover cover-arm trend field SD recovers across seeds", {
  skip_if_fast()
  skip_on_cran()

  adj <- .pf_grid_adj(8L)
  N   <- nrow(adj)
  truth <- 0.7
  seeds <- 1:6
  rec <- vapply(seeds, function(s) {
    sim <- simulate_occu_cover(
      N = N, J = 6L, positive = "lognormal",
      beta_occ = c(qlogis(0.6), 0.3), beta_p = c(qlogis(0.65), 0.1),
      beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3,
      adj = adj, sigma = 0.5, alpha = 0.0,
      pos_field = TRUE, sigma_pos_int = 0.0, sigma_pos_trend = truth, seed = s)
    fit <- .pf_fit(sim, ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
                          spatial(~ 0 + time || cell, graph = adj,
                                  to = "positive"))
    nm <- grep("^sigma_pos_field", names(fit$means), value = TRUE)[1L]
    fit$means[[nm]]
  }, numeric(1))

  expect_true(all(is.finite(rec)))
  # Grid-weighted marginal SD; median over seeds tracks truth (small-N spread).
  expect_lt(abs(stats::median(rec) - truth), 0.25)
  # No collapse to zero and no runaway.
  expect_true(all(rec > 0.2 & rec < 1.5))
})


# --- placement: a field in the positive formula == to = "positive" -----------

test_that("a spatial field in the positive formula equals the to = \"positive\" spelling", {
  skip_if_fast()
  skip_on_cran()

  adj <- .pf_grid_adj(8L)
  N   <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 5L, positive = "lognormal",
    beta_occ = c(qlogis(0.6), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3,
    adj = adj, sigma = 0.5, alpha = 0.0,
    pos_field = TRUE, sigma_pos_int = 0.0, sigma_pos_trend = 0.7, seed = 42L)
  ctrl <- list(progress = FALSE, integration = "ccd")

  fit_to <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
                   spatial(~ 0 + time || cell, graph = adj, to = "positive"),
    detection = ~ 1, positive = ~ time,
    family = occu_cover(response = "lognormal"),
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    method = "nested_laplace", control = ctrl))
  fit_place <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell"),
    detection = ~ 1, positive = ~ time + spatial(~ 0 + time || cell, graph = adj),
    family = occu_cover(response = "lognormal"),
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    method = "nested_laplace", control = ctrl))

  expect_identical(names(fit_place$means), names(fit_to$means))
  expect_equal(unname(unlist(fit_place$means)), unname(unlist(fit_to$means)),
               tolerance = 1e-10)
})


# --- control: the cover-field SD rides its own grid --------------------------

test_that("control$sigma.grid.pos.field is accepted and sets the cover-field grid", {
  skip_if_fast()
  skip_on_cran()

  adj <- .pf_grid_adj(8L)
  N   <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 5L, positive = "lognormal",
    beta_occ = c(qlogis(0.6), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3,
    adj = adj, sigma = 0.5, alpha = 0.0,
    pos_field = TRUE, sigma_pos_int = 0.0, sigma_pos_trend = 0.7, seed = 7L)

  # A one-point field grid pins the cover-field SD at that value (the marginal is
  # the weighted mean over a single grid node). Reaching the fit at all confirms
  # the control passes validation; the exact value confirms the grid was applied.
  fit <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell") +
                   spatial(~ 0 + time || cell, graph = adj, to = "positive"),
    detection  = ~ 1, positive = ~ time,
    family     = occu_cover(response = "lognormal"),
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    method  = "nested_laplace",
    control = list(progress = FALSE, integration = "ccd",
                   sigma.grid.pos.field = 0.5)))

  nm <- grep("^sigma_pos_field", names(fit$means), value = TRUE)[1L]
  expect_equal(unname(fit$means[[nm]]), 0.5, tolerance = 1e-6)
})
