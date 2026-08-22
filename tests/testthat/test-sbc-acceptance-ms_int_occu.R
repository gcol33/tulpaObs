# =============================================================================
# test-sbc-acceptance-ms_int_occu.R
# -- posterior SBC acceptance for ms_int_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_int_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6l) /. SPECIES-COUNT SCOPED the same way
  # ms_occu was: a direct probe at three seeds on a small fixture found
  # `sp3_p2_(Intercept)` stuck at p_unif ~0.0029-0.0030 (reproducible,
  # ms_occu's exact failure signature). At S=14 (the
  # `.SBC_REG_FIXTURES$ms_int_occu()` fixture below, matching this family's
  # own recovery-test fixture size), the plain Laplace-EM -- no debiasing --
  # calibrates cleanly: 5 seeds during development (0-4), min p_unif range
  # 0.0013-0.052, 0 quantities below 1e-3 out of 43 possible across all 5
  # runs, no reproducible failing coefficient. Do NOT shrink the fixture's
  # N/species count for speed. Measured (seed = 0): posterior min p_unif
  # 0.052, narrow max <1e-3.
  fit <- .SBC_REG_FIXTURES$ms_int_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_int_occu$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
