# =============================================================================
# test-sbc-acceptance-jsdm.R
# -- posterior SBC acceptance for jsdm()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("jsdm posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6o). jsdm() shares ms_count()'s exact
  # community Laplace-EM, so the same S=20 fixture that resolved
  # ms_count's #226 failure mode applies directly -- 3 seeds (0, 1, 2)
  # during development, min p_unif range 0.0093-0.0108, 0 quantities below
  # 1e-3 out of 40 possible across all 3 runs, no reproducible failing
  # coefficient. Do NOT shrink the fixture's N/species count for speed.
  #
  # bad.factor = 3.0, not the 1.75 every sibling community family uses: a
  # single refit here is cheap (~0.3s, ~35s for the whole n.sim=100 run), so
  # unlike distsamp_open the standard n.sim did not need reducing -- but at
  # 1.75 the narrow arm only reached max p_unif ~0.0009-0.0012 (borderline,
  # not the comprehensive failure every other family's narrow arm shows).
  # bad.factor=3.0 (same n.sim, ~35s either way) pushed every quantity to
  # p_unif 0 across all 3 seeds while leaving the posterior arm's numbers
  # unchanged (bad.factor only touches the narrow arm).
  fit <- .SBC_REG_FIXTURES$jsdm()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 3.0, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$jsdm$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
