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
  # coupled field needs upstream tulpa support).
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
  # Nested-Laplace generalised beyond single-season occupancy.
  expect_silent(.tobs_validate_family_method("nested_laplace", int_occu()))
  expect_silent(.tobs_validate_family_method("nested_laplace", dyn_occu()))
  # The shared areal field on the JSDM latent occupancy.
  expect_silent(.tobs_validate_family_method("nested_laplace", jsdm()))
})

test_that("ms_occu / ms_int_occu support laplace + pg_gibbs + nuts, not areal", {
  # ms_occu gained a community NUTS sampler and a shared areal field;
  # ms_int_occu gained the multi-source community NUTS sampler but has no
  # areal path.
  for (m in c("laplace", "nuts", "nested_laplace"))
    expect_silent(.tobs_validate_family_method(m, ms_occu()))
  for (m in c("laplace", "pg_gibbs", "nuts"))
    expect_silent(.tobs_validate_family_method(m, ms_int_occu()))
  expect_error(.tobs_validate_family_method("nested_laplace", ms_int_occu()),
               "not available")
})

test_that("ms_dyn_occu supports laplace + nested_laplace + nuts", {
  # A shared areal field on the first-season occupancy arm added the
  # nested_laplace route; the community HMM-forward NUTS sampler added the
  # nuts route (0.0.158).
  for (m in c("laplace", "nested_laplace", "nuts"))
    expect_silent(.tobs_validate_family_method(m, ms_dyn_occu()))
})

test_that("the rejection lists the family's supported methods", {
  # occu_categorical is Laplace-only, so nuts is genuinely rejected.
  err <- tryCatch(.tobs_validate_family_method("nuts", occu_categorical()),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Supported:")
  expect_match(err, "laplace", fixed = TRUE)
  expect_false(grepl("\"nuts\"[^.]*Supported", err))  # nuts is not in the set
})

test_that("a family with no roster entry is a no-op for the validator", {
  # Every family constructor the package exports has a `.tobs_family_methods`
  # entry, so the no-entry branch is reached only by a hand-built
  # `obs_family()`. The validator passes it through; tobs() rejects it at
  # dispatch.
  fam <- obs_family(name        = "not_a_family",
                    class_long  = "hand-built family object",
                    latent      = "bernoulli",
                    observation = "binomial_detection")
  expect_null(.tobs_family_methods[[fam$name]])
  expect_silent(.tobs_validate_family_method("laplace", fam))
})

test_that("the #116 laplace-only families give the friendly registry error", {
  # royle_nichols / occu_ttd / occu_multi / double_observer / gdistremoval /
  # distsamp_open are laplace-only. Each is in `.tobs_family_methods`, so an
  # unsupported method= is rejected by the central registry with a "Supported:"
  # pointer BEFORE dispatch, rather than falling through to a per-dispatcher
  # `.map_engine` internal-error.
  laplace_only <- c("royle_nichols", "occu_ttd", "occu_multi",
                    "double_observer", "gdistremoval", "distsamp_open")
  for (fam_name in laplace_only) {
    fam <- get(fam_name, envir = asNamespace("tulpaObs"))()
    expect_identical(.tobs_family_methods[[fam_name]], "laplace",
                     info = fam_name)
    err <- tryCatch(
      .tobs_validate_family_method("nested_laplace", fam),
      error = function(e) conditionMessage(e)
    )
    expect_match(err, sprintf("not available for %s", fam_name), fixed = TRUE)
    expect_match(err, "Supported:", fixed = TRUE)
    expect_match(err, "\"laplace\"", fixed = TRUE)
    # The internal .map_engine fall-through message must NOT surface.
    expect_false(grepl("Internal error", err), info = fam_name)
  }
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

# ---- method = "auto" ------------------------------------------------------

test_that("every family's default_engine resolves under method = 'auto'", {
  # The "auto" branch had no `pg_gibbs` case, so `t_occu()` -- the one family
  # declaring it -- could not be fitted at the documented default at all.
  ctors <- list(occu = occu, dyn_occu = dyn_occu, int_occu = int_occu,
                abun = abun, cover = cover, occu_cover = occu_cover,
                t_occu = t_occu, count = count, jsdm = jsdm)
  for (nm in names(ctors)) {
    fam <- ctors[[nm]]()
    route <- .tobs_resolve_method("auto", fam)
    expect_equal(route$engine, fam$default_engine,
                 info = sprintf("family %s", nm))
    # The resolved public name, not the literal "auto" the caller passed.
    expect_false(identical(route$method, "auto"), info = sprintf("family %s", nm))
  }
})

test_that("a fit records the resolved route, not the literal 'auto'", {
  skip_on_cran()
  # `fit$method` is branched on downstream (`identical(object$method, "nuts")`),
  # so recording "auto" for a default fit clobbers the concrete label the fitter
  # already set.
  sim <- simulate_occu(N = 80, J = 3, seed = 1)
  auto <- tobs(~ occ_cov1, data = sim$data, detection = ~ det_cov1, y = sim$y,
               family = occu(), control = list(verbose = FALSE, progress = FALSE))
  expect_equal(auto$method, "laplace")
  explicit <- tobs(~ occ_cov1, data = sim$data, detection = ~ det_cov1,
                   y = sim$y, family = occu(), method = "laplace",
                   control = list(verbose = FALSE, progress = FALSE))
  expect_equal(explicit$method, "laplace")
  expect_equal(auto$means, explicit$means)
})
