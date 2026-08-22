# =============================================================================
# test-sbc-acceptance-ms_distance.R
# -- posterior SBC acceptance for ms_distance()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_distance posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6p). Multi-seed reproducibility probe (seeds
  # 0-1) on this same fixture found a DIFFERENT coefficient dipping
  # moderately low each time (seed 0: sp15_sigma_(Intercept) at 6.3e-4; seed
  # 1: sp19_lambda_abund_cov1 at 2.7e-4), consistent w/ the ~6% expected
  # false-positive rate across 60 tested quantities per run, NOT a
  # reproducible calibration bug (that signature is the SAME coefficient
  # pinned at ~1e-6-1e-7 on every seed). bad.factor=1.75's posterior arm
  # matched bad.factor=3's exactly (min_post=0.0006338732 at
  # sp15_sigma_(Intercept), 1/60 below 1e-3) -- expected, bad.factor only
  # rescales the narrow control arm, not the well-specified posterior.
  # Measured (seed = 0, bad.factor = 3): posterior min p_unif 6.3e-4 (1/60
  # below 1e-3), narrow max p_unif 3.0e-15.
  fit <- .SBC_REG_FIXTURES$ms_distance()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 3.0, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_distance$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-4)
  expect_lt(min(pu("narrow")), 1e-3)
})
