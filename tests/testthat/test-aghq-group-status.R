# The per-group AGHQ solve status on a fit.
#
# tulpa::tulpa_re_aghq() reports, per group, whether that group's posterior mode
# search and precision factorization succeeded. A group it marks failed comes
# back with NA BLUPs, so every per-group quantity read off it is NA and every
# community-level quantity summing over groups is missing that group's
# information. What is tested here is that the status reaches the fit: a caller
# filters on `convergence(fit)$groups_failed` rather than on the text of the
# warning the engine raises, and a fit carrying one is not reported as
# converged.

test_that("the status helper reads the engine's group_ok", {
  ok  <- .tobs_aghq_group_status(list(group_ok = c(TRUE, TRUE, TRUE)))
  expect_true(ok$all_ok)
  expect_identical(ok$n_groups, 3L)
  expect_length(ok$failed, 0L)
  expect_null(ok$failed_names)

  bad <- .tobs_aghq_group_status(
    list(group_ok = c(TRUE, FALSE, TRUE, FALSE)),
    group_names = c("sp1", "sp2", "sp3", "sp4"))
  expect_false(bad$all_ok)
  expect_identical(bad$failed, c(2L, 4L))
  expect_identical(bad$failed_names, c("sp2", "sp4"))

  # Names that do not describe the groups are dropped rather than mis-aligned.
  expect_null(.tobs_aghq_group_status(list(group_ok = c(TRUE, FALSE)),
                                      group_names = c("a", "b", "c"))$failed_names)
})

test_that("a failed group takes the fit out of converged, whatever optim said", {
  # The optimizer stopping cleanly is not the question: a fit missing a group's
  # information is not an estimate of what that group contributes to.
  expect_true(.tobs_aghq_converged(list(converged = TRUE,
                                        group_ok = c(TRUE, TRUE))))
  expect_false(.tobs_aghq_converged(list(converged = TRUE,
                                         group_ok = c(TRUE, FALSE))))
  expect_false(.tobs_aghq_converged(list(converged = FALSE,
                                         group_ok = c(TRUE, TRUE))))
})

test_that("the reported iteration count is the optimizer's own", {
  # The joint AGHQ driver is one stats::optim call; BFGS counts evaluations, and
  # they are reported in the units the package's other optim-driven fitters use.
  expect_identical(.tobs_aghq_n_iter(list(counts = c(`function` = 31, gradient = 9))),
                   31L)
  # An unnamed count vector still reads the function count.
  expect_identical(.tobs_aghq_n_iter(list(counts = c(31, 9))), 31L)
  # An engine that reports none gives NA rather than a fabricated number.
  expect_identical(.tobs_aghq_n_iter(list()), NA_integer_)
})

test_that("the convergence record carries the status through to the fit", {
  rec <- .tobs_aghq_convergence_record(
    list(converged = FALSE, n_iter = 44L, group_ok = c(TRUE, FALSE, TRUE),
         groups_failed = 2L),
    group_names = c("sp1", "sp2", "sp3"))
  expect_false(rec$converged)
  expect_identical(rec$n_iter, 44L)
  expect_identical(rec$groups_failed, 2L)
  expect_identical(rec$groups_failed_names, "sp2")

  # A fit from a path with no per-group solve carries the two fields it always
  # had, and nothing that would read as an all-succeeded status.
  plain <- .tobs_aghq_convergence_record(list(converged = TRUE, n_iter = 12L))
  expect_true(plain$converged)
  expect_null(plain$group_ok)
  expect_null(plain$groups_failed)
})

test_that("a community N-mixture fit records the per-species solve status", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_abun(n_species = 3, N = 12, J = 2,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(3), 0.2), mu_p = c(0.5, -0.2),
                          sd_lambda = 0.3, sd_p = 0.3,
                          mixture = "negbin", size = 5, sigma_logr = 0.4,
                          seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1,
              species = sim$species, method = "laplace",
              control = list(verbose = FALSE, n.quad = 3L))

  cv <- convergence(fit)
  expect_length(cv$group_ok, 3L)
  expect_true(all(cv$group_ok))
  expect_length(cv$groups_failed, 0L)
  # The pair the fit reports: a converged fit says how much work it took rather
  # than reporting NA iterations alongside TRUE.
  expect_true(cv$converged)
  expect_false(is.na(cv$n_iter))
  expect_gt(cv$n_iter, 0L)
})

test_that("a failed species is on the fit, not only in a warning", {
  skip_on_cran()
  skip_if_fast()
  # A species whose solve fails takes the AGHQ objective to its failure
  # sentinel at the same parameter, and the engine refuses such an optimum
  # outright -- so a fit REPORTING one is not reachable from
  # data. The status is injected at the extractor instead, which is the boundary
  # it crosses on its way to this package; what the engine does with a genuinely
  # unsolvable group is pinned upstream.
  bad  <- 2L
  # (fetched from the namespace rather than with `:::`, which R CMD check
  # flags on a call into another package.)
  real <- get("cpp_aghq_blups", envir = asNamespace("tulpa"))
  fail_one <- function(par, oracle, nc, full) {
    bl <- real(par, oracle, nc, full)
    bl$group_ok[bad] <- FALSE
    bl$bhat[bad, ]   <- NA_real_
    bl$bvar[bad, ]   <- NA_real_
    bl$bcov[bad, , ] <- NA_real_
    bl
  }
  testthat::local_mocked_bindings(cpp_aghq_blups = fail_one, .package = "tulpa")

  sim <- simulate_ms_abun(n_species = 3, N = 12, J = 2,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(3), 0.2), mu_p = c(0.5, -0.2),
                          sd_lambda = 0.3, sd_p = 0.3,
                          mixture = "negbin", size = 5, sigma_logr = 0.4,
                          seed = 1)
  fit <- suppressWarnings(
    tobs(~ abund_cov1, data = sim$data, y = sim$y,
         family = ms_abun(mixture = "negbin"),
         detection = ~ det_cov1,
         species = sim$species, method = "laplace",
         control = list(verbose = FALSE, n.quad = 3L)))

  # The same fixture, uninjected, converges (the test above), so the FALSE below
  # is the injected species and not an exhausted iteration budget.
  cv <- convergence(fit)
  expect_identical(cv$groups_failed, bad)
  expect_identical(cv$groups_failed_names, fit$model$species_names[bad])
  expect_false(cv$group_ok[bad])
  expect_true(all(cv$group_ok[-bad]))
  # The whole point: a fit carrying a failed species does not report itself as
  # converged, so a caller filtering on `converged()` drops it without knowing
  # anything about AGHQ.
  expect_false(cv$converged)
  expect_false(converged(fit))
})
