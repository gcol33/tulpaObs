# =============================================================================
# test-sbc-acceptance-ms_occu_cover.R
# -- posterior SBC acceptance for ms_occu_cover()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_occu_cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6j-bis). The occ+p+pos analogue of ms_occu;
  # theta is the per-species realized coefficient (S x P, species-major) plus
  # one shared dispersion coordinate, not fit$means alone (a community mean),
  # so the scored quantity names come off the registry's own draws() as
  # occu_categorical/cover do. Safe to attempt where
  # ms_occu/ms_int_occu/ms_count are not because this family's Cinv stays
  # consistent with Sigma under AGHQ debiasing (#226 part 2, commit `03b87ad`)
  # -- confirmed clean over 9 seeds during development (0-8), every one
  # comfortably above 1e-3 (range 0.0024-0.087), no coefficient repeating as
  # the minimum more than twice out of 17 possible (chance level), narrow
  # control rejecting hard on all nine. Measured (seed = 0): posterior min
  # p_unif 0.087, narrow max 5.2e-4.
  fit <- .SBC_REG_FIXTURES$ms_occu_cover()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_occu_cover$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
