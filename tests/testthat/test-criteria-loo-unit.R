# =============================================================================
# test-criteria-loo-unit.R -- loo.unit = c("obs", "cell") on waic() / loo() /
# cpo(). The convenience wrapper that auto-supplies the fit's per-observation cell
# map as tulpa_criteria(group =) for leave-one-group- out cross-validation
# (LOGO-CV), so cover() / occu_cover() report plot/site-level (default) AND
# cell-level LOO without the caller hand-building the cell map. waic() and cpo()
# reach the fold through tulpa_criteria(group =); loo() builds a psis_loo through
# loo::loo(), which scores whatever columns it is handed, so it applies the same
# fold to the pointwise matrix first (.tobs_loglik_fold_group).
#
# The first half is structural / dispatch unit tests on the family cell-map
# plumbing (.tobs_loo_cell_map), the front-door group resolution
# (.tobs_criteria_group) and the column fold, driven by lightweight mock fits and
# synthetic matrices -- no model fitting, so they run in every tier (CRAN
# included). The second half fits small cover() models to pin the loo::loo /
# loo::loo_compare dispatch end to end and is gated on the slow tiers. The
# numeric equivalence of cell-level LOO to an explicit group and to a hand
# aggregation on the waic() / cpo() side lives with the cover() / occu_cover()
# fixtures in test-cover-hurdle-aggregate-pos.R and test-occu-cover-ic-fullmodel.R.
# =============================================================================

test_that(".tobs_loo_cell_map reads the family-specific per-observation cell map", {
  # cover(): the map is spi_full, the per-row spatial node (cell) index.
  cover_fit <- structure(list(spi_full = c(3L, 3L, 1L, 2L, 1L)),
                          class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  expect_identical(tulpaObs:::.tobs_loo_cell_map(cover_fit),
                   c(3L, 3L, 1L, 2L, 1L))

  # occu_cover(): the map is site_cell, the per-site field cell.
  oc_fit <- structure(
    list(model = list(model_type = "occu_cover", n_sites = 4L,
                      site_cell = c(1L, 1L, 2L, 2L))),
    class = c("tobs_fit", "tulpa_fit"))
  expect_identical(tulpaObs:::.tobs_loo_cell_map(oc_fit), c(1L, 1L, 2L, 2L))

  # occu_cover() without a group_var: site_cell defaults to the identity (each
  # site is its own cell -> cell-level LOO equals site-level there).
  oc_id <- structure(
    list(model = list(model_type = "occu_cover", n_sites = 3L)),
    class = c("tobs_fit", "tulpa_fit"))
  expect_identical(tulpaObs:::.tobs_loo_cell_map(oc_id), 1:3)
})

test_that(".tobs_loo_cell_map returns NULL when the fit carries no cell structure", {
  # A non-spatial cover() fit has no areal field, so no spi_full.
  cover_nospat <- structure(list(spi_full = NULL),
                            class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  expect_null(tulpaObs:::.tobs_loo_cell_map(cover_nospat))

  # A family with no cell-level unit (e.g. single-season occupancy).
  other <- structure(list(model = list(model_type = "single")),
                     class = c("tobs_fit", "tulpa_fit"))
  expect_null(tulpaObs:::.tobs_loo_cell_map(other))
})

test_that(".tobs_criteria_group: loo.unit = 'obs' leaves the criteria call untouched", {
  cover_fit <- structure(list(spi_full = c(1L, 1L, 2L)),
                         class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  # "obs" must NOT inject a group -- byte-identical to the ungrouped call -- and
  # must leave any pass-through arg (e.g. an explicit group) exactly as supplied.
  expect_identical(tulpaObs:::.tobs_criteria_group(cover_fit, "obs", list()),
                   list())
  expect_identical(
    tulpaObs:::.tobs_criteria_group(cover_fit, "obs", list(chunk_size = 7L)),
    list(chunk_size = 7L))
  expect_identical(
    tulpaObs:::.tobs_criteria_group(cover_fit, "obs", list(group = c(1, 1, 2))),
    list(group = c(1, 1, 2)))
})

test_that(".tobs_criteria_group: loo.unit = 'cell' injects the auto cell map", {
  cover_fit <- structure(list(spi_full = c(3L, 3L, 1L, 2L)),
                         class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  got <- tulpaObs:::.tobs_criteria_group(cover_fit, "cell", list(chunk_size = 5L))
  expect_identical(got$group, c(3L, 3L, 1L, 2L))
  expect_identical(got$chunk_size, 5L)
})

test_that(".tobs_criteria_group errors on a fit with no cell map", {
  cover_nospat <- structure(list(spi_full = NULL),
                            class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  expect_error(
    tulpaObs:::.tobs_criteria_group(cover_nospat, "cell", list()),
    "needs a per-observation cell map")
})

test_that(".tobs_criteria_group errors when loo.unit = 'cell' collides with an explicit group", {
  cover_fit <- structure(list(spi_full = c(1L, 1L, 2L)),
                         class = c("cover_fit", "tobs_fit", "tulpa_fit"))
  expect_error(
    tulpaObs:::.tobs_criteria_group(cover_fit, "cell", list(group = c(1, 2, 3))),
    "either .* or an explicit")
})

test_that("waic() / loo() / cpo() reject an unknown loo.unit", {
  fit <- structure(list(model = list(model_type = "single")),
                   class = c("tobs_fit", "tulpa_fit"))
  expect_error(waic(fit, loo.unit = "plot"), "should be one of")
  expect_error(loo(fit,  loo.unit = "plot"), "should be one of")
  expect_error(cpo(fit,  loo.unit = "plot"), "should be one of")
})

test_that("dic() rejects loo.unit instead of forwarding it to tulpa_criteria()", {
  fit <- structure(list(model = list(model_type = "single")),
                   class = c("tobs_fit", "tulpa_fit"))
  # DIC is a plug-in deviance over every observation, so no fold applies to it
  # and tulpa_criteria() has no `loo.unit` formal to take.
  expect_error(dic(fit, loo.unit = "cell"), "no cross-validation unit")
  expect_error(dic(fit, loo.unit = "obs"),  "no cross-validation unit")
  expect_false("loo.unit" %in% names(formals(tulpa::tulpa_criteria)))
})

test_that("loo.tobs_fit declares loo.unit rather than absorbing it into ...", {
  fm <- formals(tulpaObs:::loo.tobs_fit)
  expect_true("loo.unit" %in% names(fm))
  expect_identical(eval(fm$loo.unit), c("obs", "cell"))
  # The three cross-validated doors offer the same unit.
  expect_identical(eval(formals(tulpaObs:::waic.tobs_fit)$loo.unit),
                   eval(fm$loo.unit))
  expect_identical(eval(formals(tulpaObs:::cpo.tobs_fit)$loo.unit),
                   eval(fm$loo.unit))
})

test_that(".tobs_loglik_fold_group folds columns the way tulpa_criteria(group =) does", {
  set.seed(7L)
  S <- 200L; N <- 12L
  ll <- matrix(stats::rnorm(S * N, -1, 0.5), S, N)
  g  <- c(3L, 3L, 1L, 2L, 1L, 2L, 3L, 1L, 2L, 2L, 1L, 3L)

  folded <- tulpaObs:::.tobs_loglik_fold_group(ll, g)
  cells  <- sort(unique(g))
  hand   <- vapply(cells, function(k) rowSums(ll[, g == k, drop = FALSE]),
                   numeric(S))
  expect_equal(folded, hand)
  expect_equal(dim(folded), c(S, length(cells)))
  expect_null(dimnames(folded))

  # Fold for fold, the same reduction tulpa_criteria() applies internally: PSIS
  # on the folded matrix must reproduce the grouped call per-fold elpd in the
  # same order, or loo() and waic() / cpo() would be scoring different folds.
  grouped <- tulpa::tulpa_criteria(ll, criteria = "loo", pointwise = TRUE,
                                   group = g)
  direct  <- tulpa::tulpa_criteria(folded, criteria = "loo", pointwise = TRUE)
  expect_equal(grouped$pointwise$elpd_loo, direct$pointwise$elpd_loo)
  expect_equal(grouped$pointwise$pareto_k, direct$pointwise$pareto_k)
  expect_equal(grouped$elpd_loo, direct$elpd_loo)
})

test_that(".tobs_loglik_fold_group leaves an ungrouped matrix untouched", {
  ll <- matrix(stats::rnorm(20L), 5L, 4L)
  expect_identical(tulpaObs:::.tobs_loglik_fold_group(ll, NULL), ll)
  expect_identical(tulpaObs:::.tobs_loglik_fold_group(ll), ll)
  expect_error(tulpaObs:::.tobs_loglik_fold_group(ll, c(1L, 1L, 2L)),
               "must have length n_obs = 4")
})


# ---------------------------------------------------------------------------
# End-to-end: the loo::loo / loo::loo_compare dispatch on a fitted tobs_fit.
# ---------------------------------------------------------------------------

.clu_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  adj
}

# Non-spatial cover() hurdle: x drives both arms, z is noise. Both fits score the
# same N plots, which is what loo_compare() requires of its members.
.clu_cover_data <- function(seed = 501L, N = 150L) {
  set.seed(seed)
  x <- stats::rnorm(N)
  z <- stats::rnorm(N)
  occur <- stats::rbinom(N, 1L, stats::plogis(-0.2 + 0.6 * x))
  y <- ifelse(occur == 1L,
              pmin(exp(stats::rnorm(N, 0.3 - 0.4 * x, 0.4)), 1 - 1e-6), 0)
  list(data = data.frame(x = x, z = z), y = y, N = N)
}

.clu_cover_fit <- function(d, f) {
  tobs(formula = f, data = d$data, family = cover("lognormal"), y = d$y,
       method = "laplace")
}

# Spatial cover() fit whose plots are grouped into cells via group_var, so
# spi_full is a genuine many-to-one map (n_per plots per cell) and the cell fold
# is not the identity.
.clu_cover_cell_fit <- function(seed = 502L, n_s = 8L, n_per = 15L) {
  set.seed(seed)
  adj  <- .clu_chain_adj(n_s)
  cell <- rep(seq_len(n_s), each = n_per)
  N    <- length(cell)
  w_s  <- 0.6 * (sqrt(0.7) * stats::rnorm(n_s) + sqrt(0.3) * stats::rnorm(n_s))
  x    <- stats::rnorm(N)
  occur <- stats::rbinom(N, 1L, stats::plogis(-0.3 + 0.7 * x + w_s[cell]))
  y <- ifelse(occur == 1L,
              pmin(exp(stats::rnorm(N, 0.4 - 0.5 * x + w_s[cell], 0.4)),
                   1 - 1e-6), 0)
  dat <- data.frame(x = x, region = factor(cell))
  tobs(formula = ~ x + bym2(graph = adj, group_var = "region"),
       data = dat, family = cover("lognormal"), y = y,
       method = "nested_laplace",
       control = list(verbose = FALSE, sigma.grid = c(0.4, 0.8),
                      rho.grid = c(0.5, 0.9)))
}

test_that("loo() on a tobs_fit returns a genuine psis_loo", {
  skip_on_cran()
  skip_if_fast()
  d   <- .clu_cover_data()
  fit <- .clu_cover_fit(d, ~ x)
  l   <- suppressWarnings(loo(fit, n.draws = 300L))

  expect_s3_class(l, "psis_loo")
  expect_s3_class(l, "loo")
  expect_true(loo::is.loo(l))
  expect_true(is.finite(l$estimates["elpd_loo", "Estimate"]))
  expect_true(is.finite(l$estimates["elpd_loo", "SE"]))
  # One PSIS fold per plot at the default unit, each with its own Pareto k.
  expect_equal(nrow(l$pointwise), d$N)
  expect_length(loo::pareto_k_values(l), d$N)
  expect_equal(unname(l$estimates["looic", "Estimate"]),
               -2 * unname(l$estimates["elpd_loo", "Estimate"]))
  expect_identical(attr(l, "loo_unit"), "obs")
})

test_that("loo::loo_compare() reads two tobs_fit objects", {
  skip_on_cran()
  skip_if_fast()
  d  <- .clu_cover_data()
  l1 <- suppressWarnings(loo(.clu_cover_fit(d, ~ x), n.draws = 300L))
  l2 <- suppressWarnings(loo(.clu_cover_fit(d, ~ z), n.draws = 300L))

  # loo_compare() only accepts members scoring the same observations, so this
  # also pins that two fits to one dataset produce aligned pointwise matrices.
  cmp <- loo::loo_compare(list(with_x = l1, noise = l2))
  expect_s3_class(cmp, "compare.loo")
  expect_equal(nrow(cmp), 2L)
  # loo >= 2.10 returns a data.frame carrying the member names in a `model`
  # column; before that the comparison was a matrix naming them in its rows.
  # DESCRIPTION admits both, so read the names from whichever the installed
  # version populates.
  member <- if ("model" %in% colnames(cmp)) cmp$model else rownames(cmp)
  expect_setequal(member, c("with_x", "noise"))
  expect_true(all(c("elpd_diff", "se_diff", "elpd_loo") %in% colnames(cmp)))
  expect_true(all(is.finite(cmp[, "elpd_diff"])))
  # loo_compare() sorts the best member first, at a zero difference to itself.
  expect_equal(unname(cmp[1L, "elpd_diff"]), 0)
  expect_lte(unname(cmp[2L, "elpd_diff"]), 0)
})

test_that("loo(): the cell unit folds the pointwise matrix to the cell folds", {
  skip_on_cran()
  skip_if_fast()
  fit <- .clu_cover_cell_fit()
  map <- tulpaObs:::.tobs_loo_cell_map(fit)
  expect_equal(length(unique(map)), 8L)
  expect_lt(length(unique(map)), length(map))

  set.seed(21L); l_cell <- suppressWarnings(loo(fit, n.draws = 200L,
                                                loo.unit = "cell"))
  set.seed(21L); l_obs  <- suppressWarnings(loo(fit, n.draws = 200L))

  expect_s3_class(l_cell, "psis_loo")
  expect_identical(attr(l_cell, "loo_unit"), "cell")
  expect_identical(attr(l_obs,  "loo_unit"), "obs")
  # The unit is one PSIS fold: a cell at "cell", a plot at "obs".
  expect_equal(nrow(l_cell$pointwise), 8L)
  expect_equal(nrow(l_obs$pointwise), length(map))
  expect_length(loo::pareto_k_values(l_cell), 8L)
  expect_true(is.finite(l_cell$estimates["elpd_loo", "Estimate"]))

  # Identical to hand-aggregating a cell columns into its joint conditional
  # log-likelihood per draw and running the same PSIS door on the result.
  set.seed(21L)
  ll     <- tulpaObs:::.tobs_pointwise_loglik(fit, n.draws = 200L)
  cells  <- sort(unique(map))
  ll_agg <- vapply(cells, function(k) rowSums(ll[, map == k, drop = FALSE]),
                   numeric(nrow(ll)))
  hand <- suppressWarnings(tulpaObs:::.tobs_loo_one(ll_agg, fit$chain_id))
  expect_equal(l_cell$estimates["elpd_loo", "Estimate"],
               hand$estimates["elpd_loo", "Estimate"])
  expect_equal(loo::pareto_k_values(l_cell), loo::pareto_k_values(hand))

  # Cell-level is a coarser estimand than the per-plot default; the numbers move.
  expect_false(isTRUE(all.equal(l_cell$estimates["elpd_loo", "Estimate"],
                                l_obs$estimates["elpd_loo", "Estimate"])))
})

test_that("loo(): the cell unit errors on a fit with no cell map", {
  skip_on_cran()
  skip_if_fast()
  d   <- .clu_cover_data()
  fit <- .clu_cover_fit(d, ~ x)
  expect_error(loo(fit, n.draws = 200L, loo.unit = "cell"),
               "needs a per-observation cell map")
})
