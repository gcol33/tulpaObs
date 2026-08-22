# =============================================================================
# test-sbc-acceptance-int_occu.R
# -- posterior SBC acceptance for int_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("int_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (via #220). Root-caused to two compounding bugs, not an SBC adapter
  # issue: (1) int_occu()'s per-source detection design matrix lost its
  # column names when padded to full-site width (R/occu.R), so the autoscale
  # unscale step couldn't find the intercept column and left the detection
  # intercept in standardized-covariate units while correctly unscaling the
  # slope; (2) model_type == "integrated" never got the exact-marginal
  # Newton debiasing single-season occu() and dyn_occu() already have
  # (R/int_occu_marginal.R). Verified: int_occu() with one source now
  # matches occu() on the same data to full optim precision (all four
  # coefficients and both SDs bit-identical). bad.factor
  # = 1.75, the same tuning dyn_occu() needed. Measured (seed = 0): posterior
  # min p_unif 0.171, narrow max 1.0e-5.
  fit <- .SBC_REG_FIXTURES$int_occu()
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

  ok_int_occu <- pu("posterior")
  expect_length(ok_int_occu, length(qs))
  expect_gt(min(ok_int_occu), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
