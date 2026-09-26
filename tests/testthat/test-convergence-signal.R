# One convergence record and one signal per fit, on every family
# (R/convergence_signal.R).
#
# An optimiser fit that stops without meeting its criterion warns once, names
# the control to raise, and records the verdict at `fit$convergence`; the flat
# top-level `fit$converged` some fitters wrote is not left beside it.

conv_warnings <- function(expr) {
  w <- character()
  fit <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(fit = fit, w = w[grepl("convergence criterion", w)])
}

conv_ctl <- function(...) c(list(verbose = FALSE, progress = FALSE), list(...))

test_that("a capped occu() Laplace fit warns once and records the verdict", {
  sim <- simulate_occu(N = 60L, J = 4L, seed = 11L)
  r <- conv_warnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                          detection = ~ det_cov1, y = sim$y, method = "laplace",
                          control = conv_ctl(max.iter = 2L)))
  expect_length(r$w, 1L)
  expect_match(r$w, "occu\\(\\).*control\\$max.iter")
  expect_false(converged(r$fit))
  expect_false(r$fit$convergence$converged)
  expect_null(r$fit$converged)
})

test_that("a converged fit is silent", {
  sim <- simulate_occu(N = 60L, J = 4L, seed = 11L)
  r <- conv_warnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                          detection = ~ det_cov1, y = sim$y, method = "laplace",
                          control = conv_ctl()))
  expect_length(r$w, 0L)
  expect_true(converged(r$fit))
})

test_that("abun() takes the same single signal its fitter used to raise itself", {
  sim <- simulate_abun(N = 60L, J = 4L, seed = 13L)
  r <- conv_warnings(tobs(~ abund_cov1, data = sim$data, family = abun(),
                          detection = ~ det_cov1, y = sim$y, method = "laplace",
                          control = conv_ctl(max.iter = 2L)))
  expect_length(r$w, 1L)
  expect_false(converged(r$fit))
})

test_that("fp_occu() signals without verbose", {
  skip_on_cran()
  sim <- simulate_fp_occu(N = 120L, J = 5L, seed = 16L)
  r <- conv_warnings(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                          detection = ~ occ_cov1, y = sim$y, method = "laplace",
                          control = conv_ctl(max.iter = 2L)))
  expect_length(r$w, 1L)
  expect_false(converged(r$fit))
})

test_that("a flat-only record is completed and the flat slot removed", {
  fit <- list(converged = FALSE, n_iter = 7L, method = "laplace")
  expect_warning(out <- .tobs_finalize_convergence(fit, list(name = "cover")),
                 "cover\\(\\).*after 7 iterations")
  expect_identical(out$convergence$converged, FALSE)
  expect_identical(out$convergence$n_iter, 7L)
  expect_null(out$converged)
})

test_that("a sampler's Rhat verdict is recorded but not signalled", {
  fit <- list(convergence = list(converged = FALSE, n_iter = 1000L),
              method = "nuts")
  expect_no_warning(out <- .tobs_finalize_convergence(fit, list(name = "occu")))
  expect_false(out$convergence$converged)
})
