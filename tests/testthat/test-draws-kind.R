# =============================================================================
# test-draws-kind.R -- every fit states what its `$draws` are.
#
# tulpa's chain-vs-iid diagnostic gate reads `draws_kind` and treats an unstamped
# fit as a chain. A Laplace-family fit reports i.i.d. draws from its
# approximation, on which split-Rhat sits at ~1 and ESS at n_draws by
# construction, so without the stamp `check_diagnostics()` prints a convergence
# pass that says nothing. The stamp comes from the route table, at the one tail
# every tobs() fit passes.
# =============================================================================

ctl <- list(verbose = FALSE, progress = FALSE)

test_that("every public method declares the provenance of its draws", {
  kinds <- vapply(.tobs_method_table, function(r) r$draws %||% NA_character_,
                  character(1))
  expect_false(anyNA(kinds))
  expect_true(all(kinds %in% c("chain", "iid")))
  expect_setequal(names(kinds)[kinds == "chain"], c("nuts", "pg_gibbs"))
})

test_that("a Laplace fit is stamped iid and chain diagnostics are withheld", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_occu(N = 60L, J = 4L, seed = 11L)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, method = "laplace",
              control = ctl)
  expect_identical(fit$draws_kind, "iid")
  expect_true(is.na(tulpa::check_diagnostics(fit, quiet = TRUE)))
  d <- suppressMessages(tulpa::diagnostics(fit))
  expect_false(any(c("rhat", "ess_bulk") %in% names(d)))
})

test_that("a sampler fit is stamped chain and keeps its Rhat / ESS", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_occu(N = 60L, J = 4L, seed = 11L)
  fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
              detection = ~ det_cov1, y = sim$y, method = "nuts",
              control = list(verbose = FALSE, progress = FALSE, n.iter = 200L,
                             n.warmup = 200L, n.chains = 2L, seed = 1L))
  expect_identical(fit$draws_kind, "chain")
  d <- tulpa::diagnostics(fit)
  expect_true(all(c("rhat", "ess_bulk") %in% names(d)))
  expect_true(all(is.finite(d$rhat)))
})
