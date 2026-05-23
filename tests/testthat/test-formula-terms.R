# Tests for the formula-term registry + parser (R/formula_terms.R,
# R/formula_parse.R). These exercise the parsing/registry layer in isolation,
# before the builder/fitter wiring.

make_dat <- function(n = 20L) {
  set.seed(1)
  data.frame(
    elev     = rnorm(n),
    forest   = rnorm(n),
    lon      = runif(n),
    lat      = runif(n),
    observer = factor(sample(letters[1:4], n, replace = TRUE)),
    year     = sample(2001:2005, n, replace = TRUE)
  )
}

chain_adj <- function(n = 20L) {
  a <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) { a[i, i + 1L] <- 1L; a[i + 1L, i] <- 1L }
  a
}

test_that("plain formula has no structured terms and keeps fixed effects", {
  dat <- make_dat()
  p <- tulpaObs:::.tobs_parse_formula(~ elev + forest, data = dat)
  expect_length(p$terms, 0L)
  expect_setequal(attr(stats::terms(p$fe_formula), "term.labels"),
                  c("elev", "forest"))
})

test_that("icar() is stripped from fixed effects and yields a spatial spec", {
  dat <- make_dat()
  adj <- chain_adj(20L)
  p <- tulpaObs:::.tobs_parse_formula(~ elev + icar(graph = adj), data = dat)
  expect_setequal(attr(stats::terms(p$fe_formula), "term.labels"), "elev")
  expect_length(p$terms, 1L)
  expect_s3_class(p$terms[[1]], "tobs_spatial")
  expect_identical(p$terms[[1]]$type, "icar")
  expect_identical(p$terms[[1]]$n_units, 20L)
})

test_that("gp() resolves bare coordinate columns from data", {
  dat <- make_dat()
  p <- tulpaObs:::.tobs_parse_formula(~ gp(lon, lat), data = dat)
  expect_length(p$terms, 1L)
  expect_s3_class(p$terms[[1]], "tobs_spatial")
  expect_identical(p$terms[[1]]$type, "gp")
  expect_identical(p$terms[[1]]$n_obs, 20L)
})

test_that("re() and temporal() resolve grouping/time columns to codes", {
  dat <- make_dat()
  pr <- tulpaObs:::.tobs_parse_formula(~ re(observer), data = dat)
  expect_s3_class(pr$terms[[1]], "tobs_re")
  expect_identical(pr$terms[[1]]$n_groups, 4L)

  pt <- tulpaObs:::.tobs_parse_formula(~ temporal(year, type = "rw1"), data = dat)
  expect_s3_class(pt$terms[[1]], "tobs_temporal")
  expect_identical(pt$terms[[1]]$type, "rw1")
  expect_identical(pt$terms[[1]]$n_times, 5L)
})

test_that("non-registered calls stay in the fixed-effects design", {
  dat <- make_dat()
  p <- tulpaObs:::.tobs_parse_formula(~ log(elev + 5) + re(observer), data = dat)
  expect_setequal(attr(stats::terms(p$fe_formula), "term.labels"), "log(elev + 5)")
  expect_length(p$terms, 1L)
})

test_that("intercept removal is preserved through parsing", {
  dat <- make_dat()
  adj <- chain_adj(20L)
  p <- tulpaObs:::.tobs_parse_formula(~ elev - 1 + icar(graph = adj), data = dat)
  expect_identical(attr(stats::terms(p$fe_formula), "intercept"), 0L)
})

test_that("term in one process is tagged to that process only", {
  dat <- make_dat()
  adj <- chain_adj(20L)
  parsed <- tulpaObs:::.tobs_parse_processes(
    list(psi = ~ elev + icar(graph = adj), p = ~ forest),
    data = dat, env = environment()
  )
  resolved <- tulpaObs:::.tobs_resolve_terms(parsed$terms)
  expect_length(resolved, 1L)
  expect_identical(resolved[[1]]$processes, 1L)
})

test_that("copy() shares one realization across processes", {
  dat <- make_dat()
  parsed <- tulpaObs:::.tobs_parse_processes(
    list(psi = ~ gp(lon, lat, id = "u"), p = ~ forest + copy("u")),
    data = dat, env = environment()
  )
  resolved <- tulpaObs:::.tobs_resolve_terms(parsed$terms)
  expect_length(resolved, 1L)
  expect_identical(resolved[[1]]$spec$id, "u")
  expect_identical(resolved[[1]]$processes, c(1L, 2L))
})

test_that("copy() to a missing id is an error", {
  dat <- make_dat()
  parsed <- tulpaObs:::.tobs_parse_processes(
    list(psi = ~ elev, p = ~ copy("nope")),
    data = dat, env = environment()
  )
  expect_error(tulpaObs:::.tobs_resolve_terms(parsed$terms), "no term with id")
})
