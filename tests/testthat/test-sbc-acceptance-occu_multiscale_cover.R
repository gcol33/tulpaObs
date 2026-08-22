# =============================================================================
# test-sbc-acceptance-occu_multiscale_cover.R
# -- posterior SBC acceptance for occu_multiscale_cover()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("occu_multiscale_cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multiarm-S3 group, section 6m). The standard single-block fit shape
  # (unlike cover()); the exchangeable unit is the CELL, not the plot, so
  # pooling/site track cell indices. Checked at three configurations
  # (n_cells=40 seeds 0/1, n_cells=120 seed 0) before registering -- all
  # consistent, no anomaly like the ms_occu/ms_int_occu near-misses (#226).
  # Measured (n_cells=40, seed=0): posterior min p_unif 0.093
  # (psi_(Intercept)), narrow max 5.5e-6.
  fit <- .SBC_REG_FIXTURES$occu_multiscale_cover()
  res <- suppressMessages(sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L))

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs2 <- names(fit$means)
  pu2 <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs2, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok2 <- pu2("posterior")
  expect_length(ok2, length(qs2))
  expect_gt(min(ok2), 1e-3)
  expect_lt(min(pu2("narrow")), 1e-3)
})
