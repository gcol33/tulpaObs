# =============================================================================
# test-sbc-acceptance-ms_count.R
# -- posterior SBC acceptance for ms_count()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_count posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6n) /. SPECIES-COUNT SCOPED the same way
  # ms_occu/ms_int_occu were: a multi-seed (0, 1, 2) probe on a small
  # fixture originally found `sp3_mu_(Intercept)` pinned at p_unif
  # ~9.6e-7-9.9e-7 every time (worse than ms_occu/ ms_int_occu's own
  # small-fixture failures, plausibly because this family has no detection
  # arm to dilute it). At S=20 (the `.SBC_REG_FIXTURES$ms_count()` fixture
  # below, matching ms_occu's own resolved scale), the plain Laplace-EM --
  # no debiasing -- calibrates cleanly: 5 seeds during development (0-4),
  # min p_unif range 0.0016-0.086, 0 quantities below 1e-3 out of 41
  # possible across all 5 runs, no reproducible failing coefficient. Do NOT
  # shrink the fixture's N/species count for speed. Measured (seed = 0):
  # posterior min p_unif 0.034, narrow max <1e-3.
  fit <- .SBC_REG_FIXTURES$ms_count()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_count$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
