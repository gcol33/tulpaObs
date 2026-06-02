# Per-family backend coverage: `.tobs_family_methods` is the single source of
# truth for which method each working family supports. tobs() validates the
# resolved method against it and errors with a pointer to the supported set,
# instead of silently downgrading the engine (the old nested_laplace ->
# single-Laplace fall-back, which also mislabelled `fit$method`) or relying on
# per-dispatcher rejections. These tests pin the contract; they fire before any
# fitting, so the data is deliberately minimal.

# ---- validator unit level -------------------------------------------------

test_that(".tobs_validate_family_method accepts every supported method", {
  for (fam_name in names(.tobs_family_methods)) {
    fam <- get(fam_name, envir = asNamespace("tulpaObs"))()
    for (m in .tobs_family_methods[[fam_name]]) {
      expect_silent(.tobs_validate_family_method(m, fam))
    }
  }
})

test_that(".tobs_validate_family_method rejects an unsupported method", {
  # nested_laplace is wired for occu / int_occu / dyn_occu (the multi-block
  # driver) and the cover joint path -- but NOT jsdm (no latent state to
  # integrate; the response is observed directly).
  expect_error(.tobs_validate_family_method("nested_laplace", jsdm()),
               "not available for jsdm")
  # nested_laplace_sla (skew on the nested path) is occu + cover only.
  expect_error(.tobs_validate_family_method("nested_laplace_sla", dyn_occu()),
               "not available for dyn_occu")
  # cover has no HMC likelihood and no EM correction engine.
  expect_error(.tobs_validate_family_method("nuts", cover()),
               "not available for cover")
  expect_error(.tobs_validate_family_method("laplace_gibbs", cover()),
               "not available for cover")
})

test_that(".tobs_validate_family_method now accepts nested_laplace for the multi-block families", {
  # gcol33/tulpaObs: nested-Laplace generalised beyond single-season occupancy.
  expect_silent(.tobs_validate_family_method("nested_laplace", int_occu()))
  expect_silent(.tobs_validate_family_method("nested_laplace", dyn_occu()))
})

test_that("community occupancy families are laplace-only (gcol33/tulpaObs#30)", {
  # ms_occu / ms_dyn_occu / ms_int_occu fit by the dedicated community Laplace-EM
  # (independent per-arm REs); the generic-engine nested_laplace / nuts community
  # paths were removed with the mis-specified legacy RE path.
  for (fam in list(ms_occu(), ms_dyn_occu(), ms_int_occu())) {
    expect_silent(.tobs_validate_family_method("laplace", fam))
    expect_error(.tobs_validate_family_method("nested_laplace", fam),
                 "not available")
    expect_error(.tobs_validate_family_method("nuts", fam), "not available")
  }
})

test_that("the rejection lists the family's supported methods", {
  err <- tryCatch(.tobs_validate_family_method("nuts", cover()),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Supported:")
  expect_match(err, "nested_laplace", fixed = TRUE)
  expect_false(grepl("\"nuts\"[^.]*Supported", err))  # nuts is not in the set
})

test_that("planned families are a no-op (they error earlier via dispatch)", {
  expect_silent(.tobs_validate_family_method("laplace", abun()))
  expect_silent(.tobs_validate_family_method("nuts", distance()))
})

# ---- public tobs() path ---------------------------------------------------

test_that("tobs() passes the method gate for dyn_occu + nested_laplace (now supported)", {
  # dyn_occu now supports nested_laplace, so the central method check passes and
  # dispatch proceeds; with no latent term the nested driver reports the
  # missing-block error -- proving the method gate let it through rather than
  # rejecting it as "not available".
  set.seed(1)
  y  <- array(rbinom(4 * 2 * 2, 1, 0.3), dim = c(4L, 2L, 2L))
  df <- data.frame(x = rnorm(4))
  expect_error(
    tobs(~ x, data = df, family = dyn_occu(), detection = ~ 1, y = y,
         col_formula = ~ 1, ext_formula = ~ 1, method = "nested_laplace",
         control = list(verbose = FALSE)),
    "at least one latent block"
  )
})

test_that("tobs() rejects nuts and gibbs/mi for the cover hurdle", {
  df <- data.frame(x = c(0, 1, 0, 1))
  yc <- c(0, 0.3, 0, 0.7)
  expect_error(
    tobs(~ x, data = df, family = cover(), y = yc, method = "nuts"),
    "not available for cover"
  )
  expect_error(
    tobs(~ x, data = df, family = cover(), y = yc, method = "laplace_mi"),
    "not available for cover"
  )
})
