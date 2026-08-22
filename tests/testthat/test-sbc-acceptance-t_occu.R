# =============================================================================
# test-sbc-acceptance-t_occu.R
# -- posterior SBC acceptance for t_occu()
#
# One family per file. An acceptance block is a full SBC run -- n.sim fits of
# simulated data, plus the reference draws -- so a shard carries one family's
# measurement and an overrun costs that family's evidence alone. Fixtures and
# the means accessor are in helper-sbc-registry.R; the registry contract is in
# test-sbc-registry.R.
# =============================================================================

test_that("t_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # (section 6i). A pg_gibbs family whose fit$draws is already the real
  # pooled posterior sample; loglik_many is a Laplace approximation to the
  # AR1 year effect's marginal (FD-validated gradient and Hessian,
  # brute-force-cross-checked against a dense grid at T = 2). `qs` only
  # scores the reported psi/p/AR1-hyperparameter coefficients (matching
  # every other family's test); log_lik's own presence/finiteness is
  # covered by the generic cross-family CONTRACT test. Measured (seed = 0):
  # posterior min p_unif 0.081, narrow max 1.2e-4.
  fit <- .SBC_REG_FIXTURES$t_occu()
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

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
