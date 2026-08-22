# =============================================================================
# test-nuts-sampler-knobs.R -- `n.thin`, `n.threads` and `n.threads.grad` reach
# every NUTS route, and capping the gradient threads does not move the chain.
#
# `n.thin` used to be read by the occupancy route alone: `?tobs` documented it,
# the control validator admitted it and the engine table carried it, so a user
# thinning a long chain got the unthinned chain and no warning. It now runs
# through one pair of helpers (.tobs_nuts_run_parallel / .tobs_nuts_thin_chain)
# that every family's chain loop calls, and it thins the per-iteration
# diagnostics by the same stride so `divergent` still lines up with `draws`.
#
# `n.threads.grad` is the OpenMP thread count inside ONE gradient evaluation of
# a community target, whose per-species loop is parallel. The per-species
# reduction after that loop is serial and order-fixed, so the chain is the same
# at any count -- asserted here to the bit, since a thread count that changed
# the posterior would be a reduction bug, not a tuning knob.
# =============================================================================

ctl <- function(...) c(list(n.iter = 120L, n.warmup = 120L, verbose = FALSE,
                            progress = FALSE, seed = 7L), list(...))

test_that("the sampler table answers for both thinning and thread knobs", {
  d <- .tobs_engine_defaults("nuts")
  expect_identical(d$n.thin, 1L)
  expect_identical(d$n.threads, 1L)
  expect_identical(d$n.threads.grad, 0L)
  # `.tobs_fit_model()` forwards explicit values, so a literal there would be
  # the live answer and would not track the table.
  expect_null(eval(formals(.tobs_fit_model)$n.threads))
  expect_null(eval(formals(.tobs_fit_model)$n.thin))
  expect_null(eval(formals(.tobs_fit_model)$sigma.logr))
})

test_that("n.thin reaches a single-species count NUTS fit", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_abun(N = 40, J = 3, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(5), 0.3), beta_p = c(0.4, -0.2),
                       seed = 2)
  f <- function(...) tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1,
                          y = sim$y, family = abun(), method = "nuts",
                          control = ctl(...))
  f1 <- f(); f2 <- f(n.thin = 2L)
  expect_identical(nrow(f1$draws), 120L)
  expect_identical(nrow(f2$draws), 60L)
  # The kept draws are the stride-2 subset of the unthinned chain, not a
  # separately seeded run.
  expect_equal(unname(f2$draws),
               unname(f1$draws[seq.int(1L, 120L, by = 2L), , drop = FALSE]))
  # A diagnostic left unthinned would be silently misaligned with the draws.
  expect_identical(length(f2$nuts$divergent), 60L)
  expect_identical(length(f2$nuts$accept), 60L)
})

test_that("control$sigma.logr reaches the negbin NUTS prior", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_abun(N = 40, J = 3, n_abund_covs = 1, n_det_covs = 1,
                       beta_lambda = c(log(5), 0.3), beta_p = c(0.4, -0.2),
                       mixture = "negbin", size = 3, seed = 2)
  f <- function(...) tobs(~ abund_cov1, data = sim$data, detection = ~ det_cov1,
                          y = sim$y, family = abun(mixture = "negbin"),
                          method = "nuts", control = ctl(...))
  wide <- f()                      # table default 1.5
  tight <- f(sigma.logr = 0.05)
  # The knob used to vanish into `...` while a hardcoded 1.5 won, so the two
  # fits were identical; a prior 30x tighter must pull the log_r posterior in.
  expect_lt(stats::sd(tight$draws[, "log_r"]),
            0.5 * stats::sd(wide$draws[, "log_r"]))
})

test_that("n.thin and n.threads.grad reach a community NUTS fit", {
  skip_if_fast()
  skip_on_cran()
  S   <- 4L
  sim <- simulate_ms_occu(N = 40, J = 3, n_species = S,
                          beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                          alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                          seed = 3)
  f <- function(...) tobs(~ x, data = sim$data, detection = ~ 1, y = sim$y,
                          family = ms_occu(), species = paste0("sp", seq_len(S)),
                          method = "nuts", control = ctl(...))
  base <- f()
  expect_identical(base$nuts$draws, f(n.threads.grad = 1L)$nuts$draws)
  expect_identical(nrow(f(n.thin = 2L)$nuts$draws),
                   as.integer(nrow(base$nuts$draws) / 2L))
})
