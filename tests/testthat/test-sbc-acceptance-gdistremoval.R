# =============================================================================
# test-sbc-acceptance-gdistremoval.R
# -- posterior SBC acceptance for gdistremoval()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("gdistremoval posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multi-response group). bad.factor = 2.0 -- the single time-step,
  # two-response-matrix shape needed a bigger tick than the other
  # registered families to clear 1e-3 on both arms cleanly. Measured
  # (seed = 0): posterior min p_unif 0.449, narrow max 2.0e-9.
  fit <- .SBC_REG_FIXTURES$gdistremoval()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 2.0, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
