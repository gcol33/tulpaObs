# =============================================================================
# test-sbc-acceptance-ms_dyn_occu.R
# -- posterior SBC acceptance for ms_dyn_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_dyn_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6q). Shared gamma/eps globals condition each
  # species' psi1/p draw on the FULL (mu, global) vector via the (P+G) x P
  # `Bf[[s]]` cross-Hessian block (verified by inspecting dim(Bf[[1]])
  # directly, not assumed square). Measured (seed = 0): posterior min
  # p_unif 4.1e-3 (0/42 below 1e-3) at BOTH bad.factor = 1.75 and 3
  # (bad.factor rescales the narrow control only); narrow max p_unif 6.5e-6
  # (bad.factor=1.75) / 1.2e-15 (bad.factor=3).
  fit <- .SBC_REG_FIXTURES$ms_dyn_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_dyn_occu$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
