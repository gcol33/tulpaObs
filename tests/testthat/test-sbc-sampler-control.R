# =============================================================================
# test-sbc-sampler-control.R
# -- a NUTS fit's own sampler sizing seeds its sbc() refit control (#358)
#
# Before this, every registered spec defaulted a NUTS refit's `control` from
# `.tobs_sbc_default_control()` (verbose/progress only) merged with the
# caller's `fit.control`, never from the observed fit's own resolved
# n.iter/n.warmup/n.chains/.... A fit built with a short chain (for speed, or
# because that is what converges at a given N) got refit hundreds of times at
# the SAMPLER'S defaults (1000/1000) instead -- sized for one real fit, not a
# cheap SBC replicate.
#
# These tests are STRUCTURE-tier: they check that a NUTS fit stamps
# `fit$sampler_control` and that the registered specs read it back, not the
# full calibration measurement (that is test-sbc-acceptance-*.R).
# =============================================================================

test_that("a NUTS fit records the sampler control it actually resolved", {
  skip_on_cran(); skip_if_fast()
  set.seed(1)
  sim <- simulate_abun(N = 15L, J = 3L, n_abund_covs = 1L, seed = 1L)
  fit <- suppressWarnings(tobs(
    ~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
    detection = ~ det_cov1, method = "nuts",
    control = list(n.iter = 40L, n.warmup = 40L, n.chains = 1L, seed = 1L,
                   verbose = FALSE, progress = FALSE)))
  expect_identical(fit$method, "nuts")
  expect_type(fit$sampler_control, "list")
  expect_identical(fit$sampler_control$n.iter, 40L)
  expect_identical(fit$sampler_control$n.warmup, 40L)
  expect_identical(fit$sampler_control$n.chains, 1L)
  expect_identical(fit$sampler_control$seed, 1L)
  # Not the (engine, family) table's defaults (1000/1000) -- what a refit
  # would have been sized at before this fit.
  expect_false(isTRUE(fit$sampler_control$n.iter == 1000L))
})

test_that("a laplace fit carries no sampler_control (nothing to default from)", {
  skip_on_cran(); skip_if_fast()
  sim <- simulate_abun(N = 15L, J = 3L, n_abund_covs = 1L, seed = 2L)
  fit <- suppressWarnings(tobs(
    ~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
    detection = ~ det_cov1, method = "laplace",
    control = list(verbose = FALSE)))
  expect_null(fit$sampler_control)
})

test_that(".tobs_sbc_control() defaults a NUTS spec's control from fit$sampler_control", {
  skip_on_cran(); skip_if_fast()
  sim <- simulate_abun(N = 15L, J = 3L, n_abund_covs = 1L, seed = 3L)
  fit <- suppressWarnings(tobs(
    ~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
    detection = ~ det_cov1, method = "nuts",
    control = list(n.iter = 30L, n.warmup = 25L, n.chains = 1L, seed = 4L,
                   verbose = FALSE, progress = FALSE)))

  spec <- tulpaObs:::.TOBS_SBC_REGISTRY$abun$spec(fit, fit.control = list())
  # Sized like the observed fit, not the sampler's (engine, family) defaults.
  expect_identical(spec$control$n.iter, 30L)
  expect_identical(spec$control$n.warmup, 25L)
  # The default control's own knobs (verbose/progress) still apply where the
  # fit's own sampler_control doesn't carry them.
  expect_identical(spec$control$verbose, FALSE)

  # An explicit fit.control still wins over the recorded sampler_control.
  spec2 <- tulpaObs:::.TOBS_SBC_REGISTRY$abun$spec(
    fit, fit.control = list(n.iter = 999L))
  expect_identical(spec2$control$n.iter, 999L)
})

test_that(".tobs_sbc_control() falls back to engine defaults for a non-NUTS fit", {
  skip_on_cran(); skip_if_fast()
  sim <- simulate_abun(N = 15L, J = 3L, n_abund_covs = 1L, seed = 5L)
  fit <- suppressWarnings(tobs(
    ~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
    detection = ~ det_cov1, method = "laplace",
    control = list(verbose = FALSE)))
  spec <- tulpaObs:::.TOBS_SBC_REGISTRY$abun$spec(fit, fit.control = list())
  expect_null(spec$control$n.iter)
  expect_identical(spec$control$verbose, FALSE)
  expect_identical(spec$control$progress, TRUE)
})

test_that("every NUTS fit-assembly writer forwards sampler_control to fit", {
  skip_on_cran(); skip_if_fast()
  fake_chain <- function() list(
    draws = matrix(stats::rnorm(20), 10, 2, dimnames = list(NULL, c("a", "b"))),
    accept_prob = stats::runif(10), divergent = rep(FALSE, 10),
    treedepth = rep(3L, 10), epsilon = 0.1)
  chains <- list(fake_chain(), fake_chain())
  sc <- list(n.iter = 10L, n.warmup = 5L, n.chains = 2L, seed = 7L)

  f1 <- tulpaObs:::.tobs_nuts_attach_convergence(
    list(draws = NULL), chains, par_names = c("a", "b"), sampler_control = sc)
  expect_identical(f1$sampler_control, sc)

  rc <- list(draws = do.call(rbind, lapply(chains, `[[`, "draws")),
             accept = unlist(lapply(chains, `[[`, "accept_prob")),
             divergent = unlist(lapply(chains, `[[`, "divergent")),
             treedepth = unlist(lapply(chains, `[[`, "treedepth")),
             epsilon = 0.1, chains = lapply(chains, `[[`, "draws"))
  f2 <- tulpaObs:::.ms_ocs_finalize_nuts_fit(
    list(draws = NULL), rc, lay = list(mu = 1:2), n_chains = 2L,
    sampler_control = sc)
  expect_identical(f2$sampler_control, sc)
})
