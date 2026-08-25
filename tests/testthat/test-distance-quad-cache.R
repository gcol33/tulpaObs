# Per-fit quadrature caching.
#
# cpp_distance_site_sweep / cpp_distance_total_log_lik / cpp_distance_laplace_fixed
# / cpp_distance_grouped_oracle / cpp_distance_ploglik_batch / cpp_distance_nuts
# used to Newton-Raphson root-find the Gauss-Legendre quadrature fresh on every
# .Call() even though (cutpoints, transect, quad_order) are fixed for a fit's
# lifetime. They now take a quad_xptr built once by cpp_distance_build_quad()
# and reused across every repeated call. dist_build_quad() is a pure function of
# its three arguments, so a quad built once and reused must give byte-identical
# results to a quad rebuilt fresh for every call with the same arguments -- what
# this file asserts (the caching introduces no numerical difference).

test_that("a reused quad_xptr matches a freshly-rebuilt one, site sweep", {
  set.seed(201)
  cuts <- seq(0, 1, length.out = 6)
  n <- 40
  y <- matrix(rpois(n * 5, 3), n, 5)
  eta_lam <- rnorm(n, log(20), 0.3)
  eta_sig <- rnorm(n, log(0.4), 0.2)

  reused <- tulpaObs:::cpp_distance_build_quad(cuts, 0L, 64L)
  out_reused <- tulpaObs:::cpp_distance_site_sweep(
    y, eta_lam, eta_sig, reused, K_max = 200L, nb = FALSE, r = Inf)

  # A SEPARATE build from the identical (cutpoints, transect, quad_order):
  # mirrors what "rebuild on every call" used to do, since dist_build_quad()
  # is a pure function of its three arguments.
  fresh <- tulpaObs:::cpp_distance_build_quad(cuts, 0L, 64L)
  out_fresh <- tulpaObs:::cpp_distance_site_sweep(
    y, eta_lam, eta_sig, fresh, K_max = 200L, nb = FALSE, r = Inf)

  expect_identical(out_reused$log_lik, out_fresh$log_lik)
  expect_identical(out_reused$grad_lam, out_fresh$grad_lam)
  expect_identical(out_reused$grad_sig, out_fresh$grad_sig)
  expect_identical(out_reused$info_lam, out_fresh$info_lam)

  # And calling the SAME reused xptr twice in a row (the actual caching
  # pattern: one build, many calls) reproduces itself exactly.
  out_reused2 <- tulpaObs:::cpp_distance_site_sweep(
    y, eta_lam, eta_sig, reused, K_max = 200L, nb = FALSE, r = Inf)
  expect_identical(out_reused$log_lik, out_reused2$log_lik)
})

test_that("a reused quad_xptr matches a freshly-rebuilt one, total log-lik + grad", {
  set.seed(202)
  cuts <- seq(0, 1.5, length.out = 7)
  n <- 25
  y <- matrix(rpois(n * 6, 2), n, 6)
  eta_lam <- rnorm(n, log(15), 0.3)
  eta_sig <- rnorm(n, log(0.5), 0.2)

  reused <- tulpaObs:::cpp_distance_build_quad(cuts, 1L, 48L)
  o1 <- tulpaObs:::cpp_distance_total_log_lik(
    y, eta_lam, eta_sig, 0, reused, key = 0L, K_max = 200L, r = Inf)
  fresh <- tulpaObs:::cpp_distance_build_quad(cuts, 1L, 48L)
  o2 <- tulpaObs:::cpp_distance_total_log_lik(
    y, eta_lam, eta_sig, 0, fresh, key = 0L, K_max = 200L, r = Inf)

  expect_identical(o1$log_lik, o2$log_lik)
  expect_identical(o1$grad_eta_lambda, o2$grad_eta_lambda)
  expect_identical(o1$grad_eta_sigma, o2$grad_eta_sigma)
})

test_that("distance_laplace()'s internal quad build matches an externally-supplied one", {
  sim <- simulate_distance(N = 60, key = "halfnorm", transect = "line",
                           beta_lambda = c(log(30), 0.3),
                           beta_sigma  = c(log(0.45), 0.2), seed = 9)
  Xl <- model.matrix(~ abund_cov1, sim$data)
  Xs <- model.matrix(~ sigma_cov1, sim$data)

  fit_internal <- tulpaObs:::distance_laplace(
    sim$y, Xl, Xs, sim$cutpoints, key = "halfnorm", transect = "line",
    mixture = "P", verbose = FALSE)

  qptr <- tulpaObs:::cpp_distance_build_quad(as.numeric(sim$cutpoints), 0L, 64L)
  fit_external <- tulpaObs:::distance_laplace(
    sim$y, Xl, Xs, sim$cutpoints, key = "halfnorm", transect = "line",
    mixture = "P", verbose = FALSE, quad_xptr = qptr)

  expect_equal(fit_internal$beta_lambda, fit_external$beta_lambda)
  expect_equal(fit_internal$beta_sigma, fit_external$beta_sigma)
  expect_equal(fit_internal$log_lik, fit_external$log_lik)
})

test_that("ms_distance()'s shared per-species engine quad matches a standalone build", {
  set.seed(203)
  sim <- simulate_ms_distance(n_species = 3, N = 30, key = "halfnorm",
                              n_abund_covs = 1, seed = 203)
  model <- tulpaObs:::.tobs_build_ms_distance(
    ~ abund_cov1, ~ 1, sim$data, sim$y, sim$species, sim$cutpoints,
    key = "halfnorm", transect = "line")
  eng <- tulpaObs:::.tobs_ms_distance_engine(model)

  n <- model$n_sites
  eta_lam <- rnorm(n, log(20), 0.3)
  eta_sig <- rnorm(n, log(40), 0.2)
  sw_engine <- eng$sweep(1, eta_lam, eta_sig, 0)

  standalone_quad <- tulpaObs:::cpp_distance_build_quad(
    as.numeric(model$cutpoints), 0L, as.integer(model$quad_order))
  sw_standalone <- tulpaObs:::cpp_distance_site_sweep(
    eng$y_s[[1]], eta_lam, eta_sig, standalone_quad, eng$K_max,
    nb = FALSE, r = Inf)

  expect_identical(sw_engine$log_lik, sw_standalone$log_lik)
  expect_identical(sw_engine$grad_lam, sw_standalone$grad_lam)
})


# ---------------------------------------------------------------------------
# The SIMULATOR builds the rule at the fit's own order.
# cpp_simulate_distance had no quad_order parameter and hardcoded 64, while
# quad.order is a user-facing knob every likelihood path reads from
# model$quad_order -- so a fit at any other order was simulated from per-bin
# detection probabilities it was never fit against, and an SBC run at a
# non-default order compared ranks against a different generative pi.
# dist_build_quad() is a pure function of (cutpoints, transect, quad_order), so
# passing the order IS passing the rule.
# ---------------------------------------------------------------------------

# A detection scale well below the bin width, so g() has real curvature INSIDE
# a bin. That is what a coarse rule cannot integrate: 2-node Gauss-Legendre is
# exact for cubics, so a gently curved integrand (sigma at or above the bin
# width) gives the same answer at every order from 2 up and separates nothing.
.qsim_cuts <- c(0, 25, 50, 75, 100)
.qsim_draws <- matrix(c(log(200), log(15)), nrow = 1L)

test_that("the simulator's quadrature order is an argument, and it bites", {
  n <- 40L; B <- length(.qsim_cuts) - 1L
  X1 <- matrix(1, n, 1L)
  sim <- function(order, seed = 5L) {
    set.seed(seed)
    tulpaObs:::cpp_simulate_distance(X1, X1, .qsim_draws, as.numeric(.qsim_cuts),
                                     tulpaObs:::.dist_key_code("halfnorm"), 0L,
                                     as.integer(order), 0,
                                     n, B, 1L, 1L, FALSE, NA_real_, 1L)[[1L]]
  }
  # Same order, same rule, same stream: byte-identical.
  expect_identical(sim(64L), sim(64L))
  # A two-node rule does not integrate the same pi, so the allocation moves.
  expect_false(identical(sim(2L), sim(64L)))
  # And it converges by three nodes and stays there, which is what says the
  # difference above is the coarse rule rather than the plumbing.
  expect_identical(sim(3L), sim(64L))
  expect_identical(sim(512L), sim(64L))
})

test_that("the simulator refuses a quadrature order below one", {
  cuts <- c(0, 25, 50, 75, 100)
  X1 <- matrix(1, 4L, 1L)
  dr <- matrix(c(log(10), log(40)), nrow = 1L)
  expect_error(
    tulpaObs:::cpp_simulate_distance(X1, X1, dr, as.numeric(cuts),
                                     tulpaObs:::.dist_key_code("halfnorm"), 0L,
                                     0L, 0, 4L, 4L, 1L, 1L, FALSE, NA_real_, 1L),
    "quad_order")
})

test_that("simulate() on a distance fit uses that fit's quad_order", {
  cuts <- .qsim_cuts
  n <- 40L; B <- length(cuts) - 1L
  X1 <- matrix(1, n, 1L)
  # The handler reads model$quad_order; everything else here is the smallest
  # object it touches.
  fake <- function(order) structure(list(
    model = list(n_sites = n, n_bins = B, cutpoints = cuts, key = "halfnorm",
                 transect = "line", quad_order = as.integer(order),
                 X_processes = list(X1, X1),
                 process_info = list(list(p = 1L), list(p = 1L))),
    draws = .qsim_draws,
    nmix_dispersion = list(r = NULL), distance_shape = list(shape = NULL)),
    class = c("tobs_fit", "tulpa_fit"))
  direct <- function(order, seed = 5L) {
    set.seed(seed)
    tulpaObs:::cpp_simulate_distance(X1, X1, .qsim_draws, as.numeric(cuts),
                                     tulpaObs:::.dist_key_code("halfnorm"), 0L,
                                     as.integer(order), 0,
                                     n, B, 1L, 1L, FALSE, NA_real_, 1L)[[1L]]
  }
  for (order in c(2L, 64L)) {
    set.seed(5L)
    got <- tulpaObs:::.tobs_simulate_distance(fake(order), nsim = 1L)
    expect_identical(got, direct(order),
                     info = paste("quad_order", order))
  }
  # The two orders are genuinely different fits, so the loop above is not
  # comparing one answer with itself.
  expect_false(identical(direct(2L), direct(64L)))
})
