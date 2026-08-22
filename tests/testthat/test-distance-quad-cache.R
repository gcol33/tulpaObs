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
