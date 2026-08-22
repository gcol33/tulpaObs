# =============================================================================
# test-sbc-acceptance-occu_categorical.R
# -- posterior SBC acceptance for occu_categorical()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("occu_categorical posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multiarm-S3 group). occu_categorical has no fit$means/fit$draws -- two
  # independent Laplace-Gaussian blocks (presence, class) instead of a joint
  # MVN draw matrix -- so the scored quantity names come off the registry's
  # own draws() rather than fit$means. Measured (seed = 0): posterior min
  # p_unif 0.057, narrow max 2.8e-5.
  fit <- .SBC_REG_FIXTURES$occu_categorical()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$occu_categorical$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
