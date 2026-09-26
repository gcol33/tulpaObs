# `control$max.iter` / `control$tol` reach the optimiser on every Laplace route.
# The BFGS routes over an exact marginal (dyn_abun(), fp_occu()) keep their own
# budget when the caller sets none, and take the caller's when one is set.

cap_fit_dyn_abun <- function(ctrl) {
  sim <- simulate_dyn_abun(N = 80L, T = 4L, J = 3L,
                           beta_lambda = c(log(5), 0.3), p = 0.5,
                           omega = 0.6, gamma = 1, seed = 23L)
  suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = dyn_abun(),
                        detection = ~ 1, omega = ~ 1, gamma = ~ 1, y = sim$y,
                        method = "laplace",
                        control = c(list(verbose = FALSE, progress = FALSE), ctrl)))
}

cap_fit_fp_occu <- function(ctrl) {
  sim <- simulate_fp_occu(N = 120L, J = 5L, seed = 16L)
  suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                        detection = ~ occ_cov1, y = sim$y, method = "laplace",
                        control = c(list(verbose = FALSE, progress = FALSE), ctrl)))
}

test_that("max.iter reaches the dyn_abun() BFGS", {
  skip_on_cran()
  capped <- cap_fit_dyn_abun(list(max.iter = 2L))
  expect_false(capped$convergence$converged)
  # Unset, the route keeps its own budget and converges. `n_iter` counts
  # optim's evaluations, which run past its iteration cap, so it is compared
  # against the unset fit rather than against the cap itself.
  full <- cap_fit_dyn_abun(list())
  expect_true(full$convergence$converged)
  expect_lt(capped$convergence$n_iter, full$convergence$n_iter)
})

test_that("max.iter reaches the plain fp_occu() BFGS", {
  skip_on_cran()
  capped <- cap_fit_fp_occu(list(max.iter = 2L))
  expect_false(capped$convergence$converged)
  full <- cap_fit_fp_occu(list())
  expect_true(full$convergence$converged)
  expect_lt(capped$convergence$n_iter, full$convergence$n_iter)
})
