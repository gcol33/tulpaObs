# =============================================================================
# test-sbc-acceptance-distsamp_open.R
# -- posterior SBC acceptance for distsamp_open()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("distsamp_open posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (multi-season group). Constant-dynamics, Poisson only for v1 -- shares
  # dyn_abun's 3D response/pooling and fit$means/fit$draws shape.
  #
  # n.sim = 15, not the standard 100 every sibling family uses: a single
  # refit on this fixture measured ~100-130s, and the full posterior-SBC
  # experiment (pooled/augmented refit) measured 758.5 sec/sim end to end
  # (4-sim pilot). At n.sim=100 that is ~21h -- what the first acceptance
  # attempt actually ran into before it was killed with no output. n.sim=15
  # keeps a real run under ~4.3h.
  #
  # bad.factor = 4.0, not the 1.75 every sibling family uses: at n.sim=15,
  # 1.75 only pushed 3/6 quantities into clear failure (narrow p_unif for
  # lambda_(Intercept) and sigma_det_cov1 stayed at 0.256 / 0.642 -- no
  # signal at all), because the same perturbation needs to be stronger to be
  # detectable with fewer simulations. bad.factor=4.0 (same n.sim, same
  # cost) pushed every quantity below 1e-5.
  fit <- .SBC_REG_FIXTURES$distsamp_open()
  res <- sbc(fit, n.sim = 15L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 4.0, seed = 0L)

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
