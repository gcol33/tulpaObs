# =============================================================================
# test-criteria-loo-unit.R -- loo.unit = c("obs", "cell") on cpo() /
# waic() (tulpaObs#105). The convenience wrapper that auto-supplies the
# fit's per-observation cell map as tulpa_criteria(group =) for leave-one-group-
# out cross-validation (LOGO-CV), so cover() / occu_cover() report plot/site-level
# (default) AND cell-level LOO without the caller hand-building the cell map.
#
# These are structural / dispatch unit tests on the family cell-map plumbing
# (.tobs_loo_cell_map) and the front-door group resolution (.tobs_criteria_group),
# driven by lightweight mock fits -- no model fitting, so they run in every tier
# (CRAN included). The end-to-end numeric equivalence on real fits (cell-level LOO
# == explicit group == hand-aggregated) lives with the cover() / occu_cover()
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

test_that("waic() / cpo() reject an unknown loo.unit", {
  fit <- structure(list(model = list(model_type = "single")),
                   class = c("tobs_fit", "tulpa_fit"))
  expect_error(waic(fit, loo.unit = "plot"), "should be one of")
  expect_error(cpo(fit, loo.unit = "plot"),  "should be one of")
})
