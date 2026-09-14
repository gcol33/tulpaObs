# =============================================================================
# test-sampled-fit-reporting.R - what a sampled fit reports about itself.
#
#   - a pg_gibbs fit carries the per-parameter convergence record (Rhat, bulk
#     and tail ESS) that summary() / print() read, decided on the split-Rhat <
#     1.01 rule the NUTS paths use (#319)
#   - logLik() / AIC() / BIC() are finite on a pg_gibbs fit, and a field-free
#     Laplace or NUTS fit keeps the value it reports (#320)
#   - glance()$converged agrees with converged() on every engine (#321)
# =============================================================================

.sfr_fits <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    sim <- simulate_occu(N = 120, J = 4, n_occ_covs = 1, n_det_covs = 1,
                         seed = 1)
    mk <- function(m, ctl) suppressMessages(suppressWarnings(
      tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ 1,
           y = sim$y, method = m, control = ctl)))
    cache <<- list(
      laplace  = mk("laplace", list(progress = FALSE)),
      pg_gibbs = mk("pg_gibbs", list(n.iter = 600L, n.warmup = 200L,
                                     n.chains = 2L, seed = 1L,
                                     progress = FALSE)),
      nuts     = mk("nuts", list(n.iter = 150L, n.warmup = 150L,
                                 n.chains = 2L, seed = 1L, progress = FALSE,
                                 verbose = FALSE)))
    cache
  }
})

test_that("a pg_gibbs fit carries the chain convergence record", {
  skip_if_fast()
  g <- .sfr_fits()$pg_gibbs
  expect_true(all(c("parameter", "rhat", "ess_bulk", "ess_tail") %in%
                    names(g$convergence)))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% colnames(summary(g))))
  expect_true(is.finite(g$max_rhat))
  expect_true(is.finite(g$min_ess))
  expect_identical(g$convergence$converged, all(g$convergence$rhat < 1.01))
  expect_true(any(grepl("Convergence:", capture.output(print(g)))))
})

test_that("logLik() is finite on a pg_gibbs fit, and AIC() / BIC() refuse it", {
  skip_if_fast()
  g <- .sfr_fits()$pg_gibbs
  ll <- logLik(g)
  expect_true(is.finite(as.numeric(ll)))
  # A sampled fit reports a log posterior mean, which is not a maximised
  # log-likelihood, so tulpa's information criteria decline it by name.
  expect_identical(attr(ll, "quantity"), "log_posterior_mean")
  expect_error(AIC(g), "maximised log-likelihood")
  expect_error(BIC(g), "maximised log-likelihood")
})

test_that("a fit that already reports a logLik() keeps it", {
  skip_if_fast()
  # The tail fill runs only where logLik() has no finite value. A NUTS fit
  # reports mean(log_prob) over its draws with `log_lik` left NA, so a gate on
  # the slot rather than on logLik() would overwrite it with the posterior-mean
  # marginal. Asserted as a no-op on the fit itself: a hard-coded value would
  # differ across BLAS / platform and could not tell the two failures apart.
  f <- .sfr_fits()
  fill <- tulpaObs:::.tobs_attach_sampled_loglik
  for (nm in c("laplace", "nuts", "pg_gibbs")) {
    expect_identical(fill(f[[nm]]), f[[nm]], info = nm)
  }
  expect_equal(as.numeric(logLik(f$nuts)), mean(f$nuts$log_prob))
  expect_gt(stats::sd(f$nuts$log_prob), 0)
  expect_null(f$laplace$model$field_eta_offset)
})

test_that("glance()$converged agrees with converged() on every engine", {
  skip_if_fast()
  for (ft in .sfr_fits()) {
    expect_identical(glance(ft)$converged, converged(ft))
  }
})
