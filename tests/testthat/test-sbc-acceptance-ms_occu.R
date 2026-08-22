# =============================================================================
# test-sbc-acceptance-ms_occu.R
# -- posterior SBC acceptance for ms_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6j) /. SPECIES-COUNT SCOPED: at S=5 this
  # family's Laplace-EM posterior is measurably non-Gaussian (validated
  # against method="nuts", the exact reference posterior: Rhat 1.011, ESS
  # 523, 0 divergences -- Vf/Cinv came back 3-15x too narrow and point
  # estimates measurably off), and posterior SBC fails hard (p_unif as low
  # as 0 on some coefficients). At S=20 (the `.SBC_REG_FIXTURES$ms_occu()`
  # fixture below), the SAME plain Laplace-EM -- no AGHQ debiasing --
  # calibrates cleanly: 5 seeds during development (0-4), min p_unif range
  # 0.0017-0.032, 0 quantities below 1e-3 out of 81 possible across all 5
  # runs, no reproducible failing coefficient. Do NOT shrink the fixture's
  # N/species count for speed -- that resurrects the S=5 failure this test
  # exists to guard against. Measured (seed = 0): posterior min p_unif
  # 0.010, narrow max <1e-3.
  fit <- .SBC_REG_FIXTURES$ms_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_occu$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
