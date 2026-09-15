# =============================================================================
# test-sbc-gates.R -- one rejection gate, called by every SBC spec.
#
# The two checks a registered spec needs were two separate calls, and the
# preamble had already drifted: twelve specs called both, six community ones
# called the structured-term gate and not the visit-design one, and
# occu_categorical called neither.
#
# The omission bit on ms_abun(), whose binder stores `X_det_visit` and appends
# those columns to the detection coefficient names while `formulas$det` keeps
# only the site-level arm. SBC on such a fit drew theta from a posterior
# carrying visit-level coefficients and refit a model that had none -- a
# different model on each side, with nothing in the ranks saying so.
# =============================================================================

test_that("every SBC spec routes through the one rejection gate", {
  ns <- asNamespace("tulpaObs")
  specs <- grep("^[.]tobs_sbc_spec_", ls(ns, all.names = TRUE), value = TRUE)
  expect_gt(length(specs), 15L)
  for (nm in specs) {
    src <- paste(deparse(body(get(nm, envir = ns))), collapse = " ")
    expect_true(grepl(".tobs_sbc_reject_unsupported", src, fixed = TRUE),
                info = nm)
    # No spec keeps a second, separate call that could drift again.
    expect_false(grepl(".tobs_sbc_reject_visit_design(fit)", src, fixed = TRUE),
                 info = nm)
  }
  src <- paste(deparse(body(.tobs_sbc_reject_unsupported)), collapse = " ")
  expect_true(grepl(".tobs_sbc_reject_visit_design", src, fixed = TRUE))
})

test_that("a visit-level observation design is refused on every family", {
  mk <- function(fam, ncol_visit = 2L) {
    f <- structure(
      list(model = list(structured_terms = NULL,
                        X_det_visit = matrix(0, 4L, ncol_visit))),
      class = "tobs_fit")
    attr(f, "tobs_family") <- list(name = fam)
    f
  }
  for (fam in c("ms_abun", "ms_count", "jsdm", "ms_distance", "ms_dyn_occu",
                "ms_int_occu", "occu_categorical", "occu")) {
    expect_error(.tobs_sbc_reject_unsupported(mk(fam)), "visit-level",
                 info = fam)
  }
  # An EMPTY visit matrix is not a visit-level design, so the gate is a check
  # on the columns rather than on the slot's presence.
  expect_null(.tobs_sbc_reject_unsupported(mk("ms_abun", 0L)))
})

test_that("the structured-term half of the gate still fires", {
  f <- structure(list(model = list(structured_terms = list(spatial = TRUE))),
                 class = "tobs_fit")
  attr(f, "tobs_family") <- list(name = "occu")
  expect_error(.tobs_sbc_reject_unsupported(f), "structured term")
})

test_that("the RE half of the gate fires from `fit$re` alone (#340)", {
  # `.tobs_em_nested_laplace()` stamps `fit$re` with the term spec (the same
  # rule `fit$spatial`/`fit$temporal` are caught by) but never populates
  # `fit$re_effects` -- checking only `re_effects` let a nested_laplace RE fit
  # pass this gate and reach the nested route's own "requires at least one
  # latent block" error instead of the intended refusal.
  f_re <- structure(list(model = list(structured_terms = NULL), re = list(1)),
                    class = "tobs_fit")
  attr(f_re, "tobs_family") <- list(name = "occu")
  expect_error(.tobs_sbc_reject_unsupported(f_re), "structured term \\(re\\)")

  f_reeff <- structure(
    list(model = list(structured_terms = NULL), re_effects = list(g = 1)),
    class = "tobs_fit")
  attr(f_reeff, "tobs_family") <- list(name = "occu")
  expect_error(.tobs_sbc_reject_unsupported(f_reeff), "structured term \\(re\\)")
})

test_that("sbc() on a nested_laplace fit with (1 | g) refuses like laplace (#340)", {
  skip_on_cran()
  sim <- simulate_occu(N = 60, J = 4, n_occ_covs = 1, seed = 11)
  d <- sim$data
  set.seed(6)
  d$g <- factor(sample(letters[1:6], nrow(d), TRUE))
  fr <- tobs(~ occ_cov1 + (1 | g), data = d, family = occu(),
             detection = ~ det_cov1, y = sim$y, method = "nested_laplace",
             control = list(verbose = FALSE))
  expect_false(is.null(fr$re))
  err <- tryCatch(sbc(fr, n.sim = 1L, n.draws = 10L, n.ref = 5L),
                  error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "structured term \\(re\\)")
  expect_no_match(conditionMessage(err), "requires at least one latent block")
})
