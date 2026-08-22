# =============================================================================
# test-sbc-acceptance-cover.R
# -- posterior SBC acceptance for cover()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multiarm-S3 group). Same two-independent-block shape as occu_categorical
  # (presence, positive); positive = "lognormal" only for v1, dispersion held
  # fixed (no SE anywhere in the package for it). Checked at both N=200 and
  # N=600 during development -- consistent, no anomaly like the ms_occu
  # near-miss (#226). Measured (seed = 0, N=200): posterior min p_unif 0.124,
  # narrow max 1.8e-8.
  fit <- .SBC_REG_FIXTURES$cover()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$cover$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
