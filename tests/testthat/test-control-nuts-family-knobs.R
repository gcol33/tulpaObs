# NUTS knobs that only some families read are admitted only on those families
# (#384). Each lives in a family-opted control group hosted on the nuts route;
# every other NUTS family rejects it as a key it does not use.

# The NUTS fitter each opting family reaches. A family that opts into a
# family-scoped nuts group must be listed here, so a new opt-in cannot land
# without its fitter being checked for the formal.
nuts_group_fitters <- c(
  abun          = ".tobs_fit_abun_nuts",
  removal       = ".tobs_fit_removal_nuts",
  distance      = ".tobs_fit_distance_nuts",
  ms_count      = ".tobs_fit_ms_count_nuts",
  jsdm          = ".tobs_fit_ms_count_nuts",
  ms_abun       = ".tobs_fit_ms_abun_nuts",
  ms_occu       = ".tobs_fit_ms_occu_nuts",
  ms_dyn_occu   = ".tobs_fit_ms_dyn_occu_nuts",
  ms_occu_cover = ".tobs_fit_ms_occu_cover_nuts")

nuts_family_groups <- names(Filter(function(h) identical(h, "nuts"),
                                   tulpaObs:::.tobs_family_group_hosts))

family_ctor <- function(nm) {
  f <- get(nm, envir = asNamespace("tulpaObs"))
  tryCatch(f(), error = function(e) f("beta"))
}

nuts_families <- names(Filter(function(m) "nuts" %in% m,
                              tulpaObs:::.tobs_family_methods))

test_that("every key of an opted nuts group is a formal of the family's NUTS fitter", {
  for (nm in nuts_families) {
    fam <- family_ctor(nm)
    groups <- intersect(fam$control_groups, nuts_family_groups)
    if (!length(groups)) next
    expect_true(nm %in% names(nuts_group_fitters), info = nm)
    fitter <- get(nuts_group_fitters[[nm]], envir = asNamespace("tulpaObs"))
    keys <- unlist(tulpaObs:::.tobs_control_groups[groups], use.names = FALSE)
    expect_true(all(keys %in% names(formals(fitter))),
                info = paste(nm, ":", paste(setdiff(keys, names(formals(fitter))),
                                            collapse = ", ")))
  }
})

test_that("a family-scoped nuts key is rejected on a family that does not read it", {
  route <- tulpaObs:::.tobs_resolve_method("nuts", occu())
  for (k in c("dispersion.re", "sigma.ld.init", "n.threads.grad", "sigma.logr"))
    expect_error(tulpaObs:::.tobs_validate_control(setNames(list(1), k), route, occu()),
                 sprintf("'%s' is not used by occu\\(\\)", k))
  route_mi <- tulpaObs:::.tobs_resolve_method("nuts", ms_int_occu())
  expect_error(tulpaObs:::.tobs_validate_control(list(n.threads.grad = 2L), route_mi,
                                                 ms_int_occu()),
               "not used by ms_int_occu")
})

test_that("an opting family accepts its keys under nuts and not under laplace", {
  fam <- ms_abun()
  ok <- list(n.threads.grad = 2L, sigma.logr = 1)
  expect_silent(tulpaObs:::.tobs_validate_control(
    ok, tulpaObs:::.tobs_resolve_method("nuts", fam), fam))
  expect_error(tulpaObs:::.tobs_validate_control(
    list(sigma.logr = 1), tulpaObs:::.tobs_resolve_method("laplace", fam), fam),
    "sigma.logr")
  fam_oc <- ms_occu_cover("beta")
  expect_silent(tulpaObs:::.tobs_validate_control(
    list(dispersion.re = TRUE, sigma.ld.init = 0.3),
    tulpaObs:::.tobs_resolve_method("nuts", fam_oc), fam_oc))
})
