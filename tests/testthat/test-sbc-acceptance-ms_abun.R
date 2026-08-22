# =============================================================================
# test-sbc-acceptance-ms_abun.R
# -- posterior SBC acceptance for ms_abun()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("ms_abun posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (community group, section 6r). ms_abun() needed TWO upstream tulpa
  # engine changes before it could even be attempted -- (blup_cross) exposed
  # the mode/theta cross-Hessian, but that alone was not enough:
  # tulpa_re_aghq() also only ever exposed the per-RE-TERM DIAGONAL of a
  # group's posterior covariance (blup_var). This family's lambda
  # (abundance) and p (detection) arms share a per-species grouping factor
  # with real cross-arm posterior covariance -- the same lambda/p
  # identifiability ridge tobs()'s penalized-EM exists to break for
  # occupancy psi/p -- so drawing a species' b_lambda_s and b_p_s
  # independently would repeat the #226 bug one level deeper (inside a
  # species instead of between the community mean and a species). Fixed by
  # adding blup_cov_g/blup_cross_g (the full per-group joint covariance/
  # cross-Hessian, unsliced by term), validated against a closed-form
  # joint-Hessian construction on a toy model with deliberately collinear RE
  # terms (tulpa test-re-aghq-cross-hessian.R).
  #
  # Only reachable via the AGHQ/joint_fd engine (control = list(optimizer =
  # "joint_fd", n.quad > 1)) -- the DEFAULT n.quad = 1 Laplace-EM
  # (cpp_nmix_community_em(), a different, faster engine entirely) does not
  # expose Cinv/Bf at all; .tobs_sbc_reject_ms_abun_scope() errors with a
  # pointer rather than silently drawing an independent (mu, b_s).
  #
  # Registered at REDUCED n.sim = 15 (distsamp_open's own precedent for the
  # identical problem, NEWS.md 0.0.227): a n.sim = 3 probe at reduced
  # n.draws/n.ref measured ~449s/sim, so a full n.sim = 100 at standard
  # n.draws/n.ref was projected at >12h; the actual n.sim = 15 run at
  # standard n.draws = 1000/n.ref = 200 completed in 39.9 min (faster than
  # projected). Measured (seed = 0, bad.factor = 3.0): posterior min p_unif
  # 3.9e-3 (0/40 below 1e-3), narrow min p_unif 2.9e-12 (the pass criterion
  # every sibling family uses). ONE narrow-arm quantity did not fail as
  # comprehensively as siblings' (max_narrow 1.3e-2, vs ~1e-6 to 1e-15
  # elsewhere) -- not re-investigated given cost (a re-run is another ~40
  # min), but does not weaken the calibration signal: the posterior arm
  # (the primary criterion) calibrates cleanly, and the narrow arm's role is
  # only to demonstrate SBC has power to detect a genuinely mis-scaled
  # model, which the min narrow p_unif already does at 2.9e-12.
  fit <- .SBC_REG_FIXTURES$ms_abun()
  res <- sbc(fit, n.sim = 15L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 3.0, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$ms_abun$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
