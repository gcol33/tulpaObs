# NUTS sampler controls (dotted names) + cross-chain convergence diagnostics.
# Verifies that n.chains pools draws, n.thin reduces retained draws, the
# renamed controls reach the sampler, and split-Rhat / ESS are surfaced on the
# fit and in summary(). Iteration counts are kept small for speed; the model
# is well-identified so a short run still mixes.

simulate_occu_small <- function(seed = 1, N = 60L, J = 4L,
                                psi = 0.6, p = 0.45) {
  set.seed(seed)
  z <- rbinom(N, 1, psi)
  y <- matrix(rbinom(N * J, 1, rep(z, times = J) * p), N, J)
  list(y = y, data = data.frame(x = rnorm(N)))
}

test_that("n.chains pools chains and reports convergence diagnostics", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_small()
  fit <- tobs(~ 1, data = sim$data, y = sim$y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.chains = 2L, n.iter = 400L, n.warmup = 200L,
                             adapt.delta = 0.9, max.treedepth = 8L,
                             seed = 7L, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$n_chains, 2L)
  expect_equal(nrow(fit$draws), 2L * (400L - 200L))
  expect_false(is.null(fit$convergence))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(fit$convergence)))
  expect_true(all(is.finite(fit$convergence$rhat)))
  expect_lt(max(fit$convergence$rhat), 1.2)
  expect_true(all(fit$convergence$ess_bulk > 0))
  # chain_id maps every pooled draw row to a chain
  expect_length(fit$chain_id, nrow(fit$draws))
  expect_setequal(unique(fit$chain_id), c(1L, 2L))
})

test_that("n.thin reduces retained draws per chain", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_small(seed = 2)
  fit <- tobs(~ 1, data = sim$data, y = sim$y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.chains = 1L, n.iter = 400L, n.warmup = 200L,
                             n.thin = 2L, seed = 9L, verbose = FALSE))
  expect_identical(fit$n_thin, 2L)
  expect_equal(nrow(fit$draws), length(seq.int(1L, 400L - 200L, by = 2L)))
})

test_that("summary() carries Rhat / ESS columns for a NUTS fit", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_small(seed = 3)
  fit <- tobs(~ 1, data = sim$data, y = sim$y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.chains = 2L, n.iter = 300L, n.warmup = 150L,
                             seed = 11L, verbose = FALSE))
  s <- summary(fit)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(s)))
  expect_true(any(is.finite(s$rhat)))
})

test_that("single-chain NUTS still works and reports a 1-chain convergence table", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_small(seed = 4)
  fit <- tobs(~ 1, data = sim$data, y = sim$y, detection = ~ 1,
              family = occu(), method = "nuts",
              control = list(n.iter = 300L, n.warmup = 150L,
                             seed = 13L, verbose = FALSE))
  expect_identical(fit$n_chains, 1L)
  expect_equal(nrow(fit$draws), 150L)
  # split-Rhat is defined for one chain (posterior splits it in two)
  expect_false(is.null(fit$convergence))
  expect_true(any(is.finite(fit$convergence$rhat)))
})
