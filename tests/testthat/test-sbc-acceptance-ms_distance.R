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

  # (community group, section 6g). Simulation i draws theta at seed + i, so
  # runs at neighbouring seeds share all but one replicate and are NOT
  # independent evidence about calibration. A single replicate whose refit
  # stalls one species far from its truth pins that coefficient's rank at the
  # extreme and fails the posterior arm on its own: the community EM's joint
  # Newton step accepted a direction that descended a species whose block was
  # indefinite at the start, and the species then never moved. The engine's
  # modified Newton step (absolute-eigenvalue modification of each species
  # block) removes that. bad.factor only rescales the narrow control arm, so
  # it leaves the posterior arm unchanged.
  fit <- .SBC_REG_FIXTURES$ms_distance()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 3.0, seed = 0L)

  expect_s3_class(res, "tulpa_sbc")
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
