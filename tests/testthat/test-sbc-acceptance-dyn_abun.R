# =============================================================================
# test-sbc-acceptance-dyn_abun.R
# -- posterior SBC acceptance for dyn_abun()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("dyn_abun posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multi-season group). dyn_abun shares dyn_occu's 3D response and
  # site-axis pooling but has its own working simulate(), so the replicate is
  # the shared simple-family route, not a bespoke forward simulator. Measured
  # (seed = 0): posterior min p_unif 0.087, narrow max 1.6e-7.
  fit <- .SBC_REG_FIXTURES$dyn_abun()
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

  ok_dyn_abun <- pu("posterior")
  expect_length(ok_dyn_abun, length(qs))
  expect_gt(min(ok_dyn_abun), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
