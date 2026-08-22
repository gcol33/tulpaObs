# =============================================================================
# test-sbc-acceptance-dyn_int_occu.R
# -- posterior SBC acceptance for dyn_int_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("dyn_int_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (the multi-season x multi-source product shape, section 6h). A named
  # list of S per-source 3D arrays, pooled with `.tobs_sbc_pool_named_3d`;
  # simulate() wraps the family's own `.tobs_simulate_dyn_int_occu()`
  # handler. Measured (seed = 0): posterior min p_unif 0.033, narrow max
  # 2.5e-5.
  fit <- .SBC_REG_FIXTURES$dyn_int_occu()
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
