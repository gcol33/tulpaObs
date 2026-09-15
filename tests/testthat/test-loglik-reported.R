# =============================================================================
# test-loglik-reported.R -- logLik() is a finite value, or NA with a reason.
#
# logLik.tulpa_fit() takes mean(log_prob) whenever log_prob is non-NULL, so a
# sampler that stores an all-NA placeholder reports NaN, and a fit that records
# its log marginal only in `log_lik` reports NA. The tobs() tail fills a missing
# value through the family's pointwise kernel at the posterior mean, and when
# there is no kernel to fill it with, drops the placeholder and records why.
# =============================================================================

.ll_finite <- function(fit, info) {
  ll <- stats::logLik(fit)
  expect_true(is.finite(as.numeric(ll)), info = info)
  expect_true(is.finite(attr(ll, "nobs")), info = info)
}

test_that("a fit with no value to report declines with a reason, not NaN", {
  fit <- structure(list(log_prob = rep(NA_real_, 5L), draws = NULL,
                        model = list(model_type = "none")),
                   class = c("tobs_fit", "tulpa_fit"))
  out <- .tobs_attach_sampled_loglik(fit)
  expect_null(out$log_prob)
  expect_identical(out$log_evidence_declined, "no_posterior_draws")
  ll <- stats::logLik(out)
  expect_true(is.na(as.numeric(ll)) && !is.nan(as.numeric(ll)))
  expect_identical(attr(ll, "declined"), "no_posterior_draws")
})

test_that("t_occu() PG fits report logLik and WAIC from the per-(site, season) marginal", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_t_occu(N = 60L, T_seasons = 5L, J = 3L, beta_occ = c(0.2, 0.6),
                         p = 0.4, rho = 0.6, sigma = 0.7, seed = 31L)
  fit <- tobs(~ x, data = sim$data, family = t_occu(), detection = ~ 1,
              y = sim$y, method = "pg_gibbs",
              control = list(verbose = FALSE, progress = FALSE, n.iter = 600L,
                             n.warmup = 300L, n.chains = 2L))
  .ll_finite(fit, "t_occu")

  pl <- .tobs_ploglik_t_occu(fit, n.draws = 50L)
  expect_identical(dim(pl), c(50L, sum(fit$model$nvis > 0L)))
  expect_true(all(is.finite(pl)) && all(pl <= 0))

  # One cell against the two-state marginal written out by hand.
  d <- 1L; cell <- which(fit$model$nvis > 0L)[1L]
  i <- (cell - 1L) %% fit$model$n_sites + 1L
  t <- (cell - 1L) %/% fit$model$n_sites + 1L
  th <- fit$draws[d, ]
  psi <- stats::plogis(sum(fit$model$X_occ[i, ] * th[1:2]) + fit$temporal_field[t])
  p <- stats::plogis(th[[3L]])
  n <- fit$model$nvis[cell]; k <- fit$model$kdet[cell]
  ref <- log(psi * p^k * (1 - p)^(n - k) + (1 - psi) * (k == 0L))
  expect_equal(.tobs_ploglik_t_occu(fit, n.draws = 1L)[1L, 1L], ref,
               tolerance = 1e-10)

  expect_true(is.finite(waic(fit)$estimates["waic", "Estimate"]))
})

test_that("ms_abun() Laplace fits surface their log marginal and N", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_ms_abun(n_species = 6L, N = 80L, J = 3L, seed = 0L)
  fit <- tobs(~ 1, data = sim$data, family = ms_abun(mixture = "poisson"),
              detection = ~ 1, y = sim$y, species = paste0("sp", 1:6),
              method = "laplace", control = list(verbose = FALSE))
  .ll_finite(fit, "ms_abun")
  expect_identical(as.numeric(stats::logLik(fit)), fit$log_lik)
  expect_identical(fit$N, sum(!is.na(sim$y)))
})

test_that("jsdm() / ms_count() sampler fits report a finite logLik", {
  skip_if_fast()
  skip_on_cran()
  smp <- list(verbose = FALSE, progress = FALSE, n.iter = 200L, n.warmup = 200L,
              n.chains = 2L, seed = 1L)
  sj <- simulate_jsdm(N = 80L, n_species = 6L, seed = 0L)
  fj <- tobs(~ x, data = sj$data, family = jsdm(), y = sj$y,
             species = paste0("sp", 1:6), method = "nuts", control = smp)
  .ll_finite(fj, "jsdm nuts")
  fg <- tobs(~ x, data = sj$data, family = jsdm(), y = sj$y,
             species = paste0("sp", 1:6), method = "pg_gibbs",
             control = list(verbose = FALSE, progress = FALSE, n.iter = 800L,
                            n.warmup = 300L, n.chains = 2L))
  .ll_finite(fg, "jsdm pg_gibbs")
  sc <- simulate_ms_count(N = 80L, n_species = 6L, response = "poisson",
                          seed = 0L)
  fc <- tobs(~ x, data = sc$data, family = ms_count("poisson"), y = sc$y,
             species = paste0("sp", 1:6), method = "nuts", control = smp)
  .ll_finite(fc, "ms_count nuts")
  expect_equal(as.numeric(stats::logLik(fc)),
               sum(.tobs_loglik_at_mean(fc, n.draws = 1L)))
})
