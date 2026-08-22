# =============================================================================
# test-sbc-acceptance-occu.R
# -- posterior SBC acceptance for occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  fit <- .SBC_REG_FIXTURES$occu()
  # `bad.factor` is 1.5 rather than the default 1.25 because 100 simulations on
  # four coefficients do not resolve a 20% mis-scale: measured over two seeds
  # the 1.25 control landed at 1.8e-3 and 2.1e-2, so an assertion on it would
  # be reporting the seed. At 1.5 it lands at 7.7e-7 and 2.2e-6.
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.5, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  # The gate is the smallest uniformity p-value over the four coefficients, at
  # a threshold well below any per-quantity band level: each band holds at 0.95
  # SEPARATELY, so requiring four at once fails on a calibrated algorithm about
  # a fifth of the time. This is a multiplicity-aware regression guard.
  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)

  # A deliberately mis-scaled posterior has to fail the same read, or the band
  # is not measuring anything. Same simulations, same fits, narrower report.
  expect_lt(min(pu("narrow")), 1e-3)
})
