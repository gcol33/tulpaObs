# Backward-compatibility for the 0.0.109 argument / simulator renames.
#
# Two rename families, both with deprecated pass-throughs:
#   * community simulators   simulate_dyn_ms_occu -> simulate_ms_dyn_occu
#                            simulate_int_ms_occu -> simulate_ms_int_occu
#   * extra-arm formula args col_formula  -> colonization   ext_formula -> extinction
#                            omega_formula -> omega          gamma_formula -> gamma
#                            fp_formula    -> p10            b_formula    -> certainty
# The new names are canonical; the old names still work but warn once.

# Evaluate `expr`, returning its value and any warning messages (muffled).
catch_w <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      ws[[length(ws) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = ws)
}
has_dep <- function(ws) any(grepl("deprecated", ws))

test_that("renamed community simulators keep the old names as deprecated aliases", {
  new_d <- simulate_ms_dyn_occu(N = 8, J = 2, n_species = 2, n_seasons = 3, seed = 1)
  expect_type(new_d, "list")
  expect_true(all(c("y", "data", "truth") %in% names(new_d)))

  old_d <- catch_w(simulate_dyn_ms_occu(N = 8, J = 2, n_species = 2,
                                        n_seasons = 3, seed = 1))
  expect_true(has_dep(old_d$warnings))
  expect_identical(dim(old_d$value$y), dim(new_d$y))

  new_i <- simulate_ms_int_occu(N = 8, J = c(2, 2), n_species = 2, seed = 1)
  old_i <- catch_w(simulate_int_ms_occu(N = 8, J = c(2, 2), n_species = 2, seed = 1))
  expect_true(has_dep(old_i$warnings))
  expect_equal(length(old_i$value$y), length(new_i$y))
})

test_that(".tobs_arm_formula resolves new names, warns on the deprecated ones", {
  f <- ~ x
  # formulas carry their environment, so compare the deparsed form
  eq <- function(a, b) expect_identical(deparse(a), deparse(b))

  # new name: taken silently
  new <- catch_w(.tobs_arm_formula(list(omega = f), "omega", "omega_formula"))
  expect_false(has_dep(new$warnings)); eq(new$value, f)

  # deprecated name: honoured, but warns
  old <- catch_w(.tobs_arm_formula(list(omega_formula = f), "omega", "omega_formula"))
  expect_true(has_dep(old$warnings)); eq(old$value, f)

  # neither supplied: default
  eq(.tobs_arm_formula(list(), "omega", "omega_formula"), ~ 1)

  # new name wins over deprecated, no warning
  win <- catch_w(.tobs_arm_formula(list(p10 = ~ z, fp_formula = ~ w),
                                   "p10", "fp_formula"))
  expect_false(has_dep(win$warnings)); eq(win$value, ~ z)
})

test_that("deprecated arm-formula spellings still fit end to end", {
  skip_if_fast()
  skip_on_cran()

  # Reset the RNG before each fit: the pseudo-draws advance base R's RNG, so the
  # deprecated-name fit must start from the same stream to isolate the rename as
  # the only difference.

  # fp_occu: new p10 = vs deprecated fp_formula = must give the same fit.
  sim <- simulate_fp_occu(N = 120, J = 4, seed = 7)
  set.seed(100)
  new_fit <- tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                  detection = ~ 1, y = sim$y, p10 = ~ 1,
                  control = list(verbose = FALSE))
  set.seed(100)
  old_fit <- catch_w(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                          detection = ~ 1, y = sim$y, fp_formula = ~ 1,
                          control = list(verbose = FALSE)))
  expect_true(has_dep(old_fit$warnings))
  expect_equal(coef(new_fit), coef(old_fit$value), tolerance = 1e-6)

  # dyn_occu: new colonization/extinction vs deprecated col_formula/ext_formula.
  sd <- simulate_dyn_occu(N = 60, J = 3, n_seasons = 3, seed = 8)
  set.seed(100)
  new_d <- tobs(~ 1, data = sd$data, family = dyn_occu(), detection = ~ 1,
                y = sd$y, colonization = ~ 1, extinction = ~ 1,
                control = list(verbose = FALSE))
  set.seed(100)
  old_d <- catch_w(tobs(~ 1, data = sd$data, family = dyn_occu(), detection = ~ 1,
                        y = sd$y, col_formula = ~ 1, ext_formula = ~ 1,
                        control = list(verbose = FALSE)))
  expect_true(has_dep(old_d$warnings))
  expect_equal(coef(new_d), coef(old_d$value), tolerance = 1e-6)
})
