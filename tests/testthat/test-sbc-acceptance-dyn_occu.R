# =============================================================================
# test-sbc-acceptance-dyn_occu.R
# -- posterior SBC acceptance for dyn_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("dyn_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multi-season group). bad.factor = 1.75 -- the extra tick past
  # double_observer's 1.5 that epsilon's weak identification at 5 seasons
  # needs to clear 1e-3 cleanly. Measured (seed = 0): posterior min p_unif
  # 0.090 (epsilon), narrow max 1.5e-4 (epsilon).
  fit <- .SBC_REG_FIXTURES$dyn_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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
