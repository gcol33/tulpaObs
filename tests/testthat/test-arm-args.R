# Extra-arm formula arguments of tobs() and the community simulators use bare
# process / symbol names, matching detection / positive / availability:
#   dyn_occu   colonization =   extinction =
#   dyn_abun   omega =          gamma =
#   fp_occu    p10 =            certainty =   (the `b` arm; `b` would partial-
#                                              match tobs()'s `by`)
# Community simulators match their family constructors:
#   simulate_ms_dyn_occu()   simulate_ms_int_occu()

test_that("community simulators carry their family-constructor names", {
  d <- simulate_ms_dyn_occu(N = 8, J = 2, n_species = 2, n_seasons = 3, seed = 1)
  expect_type(d, "list")
  expect_true(all(c("y", "data", "truth") %in% names(d)))

  i <- simulate_ms_int_occu(N = 8, J = c(2, 2), n_species = 2, seed = 1)
  expect_type(i, "list")
  expect_true("y" %in% names(i))
})

test_that("arm-formula arguments fit end to end under their bare names", {
  skip_if_fast()
  skip_on_cran()

  # fp_occu: the p10 (false-positive) and certainty (b) arms.
  sim <- simulate_fp_occu(N = 120, J = 4, seed = 7)
  fp <- tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
             detection = ~ 1, y = sim$y, p10 = ~ 1, certainty = ~ 1,
             control = list(verbose = FALSE))
  expect_s3_class(fp, "tobs_fit")
  # coef() is a per-arm list (psi / p11 / p10 / b); the p10 and b arms are the
  # ones set by the p10 = and certainty = arguments.
  expect_true(all(c("psi", "p11", "p10", "b") %in% names(coef(fp))))
  expect_true(all(is.finite(unlist(coef(fp)))))

  # dyn_occu: the colonization / extinction arms.
  sd <- simulate_dyn_occu(N = 60, J = 3, n_seasons = 3, seed = 8)
  dy <- tobs(~ 1, data = sd$data, family = dyn_occu(), detection = ~ 1,
             y = sd$y, colonization = ~ 1, extinction = ~ 1,
             control = list(verbose = FALSE))
  expect_s3_class(dy, "tobs_fit")
  # coef() is a per-arm list; gamma / epsilon are the colonization / extinction
  # arms set by the colonization = and extinction = arguments.
  expect_true(all(c("psi1", "p", "gamma", "epsilon") %in% names(coef(dy))))
  expect_true(all(is.finite(unlist(coef(dy)))))
})

test_that("a missing required arm argument errors with a pointer", {
  sd <- simulate_dyn_occu(N = 30, J = 2, n_seasons = 2, seed = 3)
  expect_error(
    tobs(~ 1, data = sd$data, family = dyn_occu(), detection = ~ 1,
         y = sd$y, extinction = ~ 1, control = list(verbose = FALSE)),
    "colonization")
})
