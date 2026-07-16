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
  # nested_laplace is wired across the areal families but NOT the non-spatial
  # community joint cover hurdle (the per-species RE layered on the shared
  # coupled field needs upstream tulpa support; gcol33/tulpaObs#47).
  expect_error(.tobs_validate_family_method("nested_laplace", ms_occu_cover()),
               "not available for ms_occu_cover")
  # nested_laplace_sla (skew on the nested path) is occu + cover only.
  expect_error(.tobs_validate_family_method("nested_laplace_sla", dyn_occu()),
               "not available for dyn_occu")
  # cover has no EM correction engine (no gibbs/mi); it does have a NUTS path.
  expect_error(.tobs_validate_family_method("laplace_gibbs", cover()),
               "not available for cover")
})

test_that(".tobs_validate_family_method now accepts nested_laplace for the multi-block families", {
  # gcol33/tulpaObs: nested-Laplace generalised beyond single-season occupancy.
  expect_silent(.tobs_validate_family_method("nested_laplace", int_occu()))
  expect_silent(.tobs_validate_family_method("nested_laplace", dyn_occu()))
  # The shared areal field on the JSDM latent occupancy (gcol33/tulpaObs#76).
  expect_silent(.tobs_validate_family_method("nested_laplace", jsdm()))
})

test_that("ms_occu supports laplace/nuts/nested_laplace; ms_int_occu laplace-only", {
  # ms_occu gained a community NUTS sampler (gcol33/tulpaObs#69) and a shared
  # areal field (gcol33/tulpaObs#75); ms_int_occu remains on the dedicated
  # community Laplace-EM only.
  for (m in c("laplace", "nuts", "nested_laplace"))
    expect_silent(.tobs_validate_family_method(m, ms_occu()))
  expect_silent(.tobs_validate_family_method("laplace", ms_int_occu()))
  expect_error(.tobs_validate_family_method("nested_laplace", ms_int_occu()),
               "not available")
  expect_error(.tobs_validate_family_method("nuts", ms_int_occu()),
               "not available")
})

test_that("ms_dyn_occu supports laplace + nested_laplace (stMsPGOcc, #123)", {
  # A shared areal field on the first-season occupancy arm added the
  # nested_laplace route (gcol33/tulpaObs#123); NUTS remains a follow-up.
  expect_silent(.tobs_validate_family_method("laplace", ms_dyn_occu()))
  expect_silent(.tobs_validate_family_method("nested_laplace", ms_dyn_occu()))
  expect_error(.tobs_validate_family_method("nuts", ms_dyn_occu()),
               "not available")
})

test_that("the rejection lists the family's supported methods", {
  err <- tryCatch(.tobs_validate_family_method("nuts", ms_dyn_occu()),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Supported:")
  expect_match(err, "laplace", fixed = TRUE)
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
         colonization = ~ 1, extinction = ~ 1, method = "nested_laplace",
         control = list(verbose = FALSE)),
    "at least one latent block"
  )
})

test_that("tobs() rejects gibbs/mi for the cover hurdle (no EM correction engine)", {
  df <- data.frame(x = c(0, 1, 0, 1))
  yc <- c(0, 0.3, 0, 0.7)
  expect_error(
    tobs(~ x, data = df, family = cover(), y = yc, method = "laplace_gibbs"),
    "not available for cover"
  )
  expect_error(
    tobs(~ x, data = df, family = cover(), y = yc, method = "laplace_mi"),
    "not available for cover"
  )
})
