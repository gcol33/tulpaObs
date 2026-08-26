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
  p <- tulpaObs:::.tobs_parse_formula(
    ~ gp(lon, lat, prior_range = c(0.1, 0.05)), data = dat)
  expect_length(p$terms, 1L)
  expect_s3_class(p$terms[[1]], "tobs_spatial")
  expect_identical(p$terms[[1]]$type, "gp")
  expect_identical(p$terms[[1]]$n_obs, 20L)
})

# --- spatial() umbrella: dispatches to the specific constructors -----------

# Compare two specs on their substantive fields (the recorded call text and
# label differ between `spatial(..., model = "x")` and `x(...)`).
spec_fields <- function(spec) {
  spec[setdiff(names(spec), c("term_call", "label"))]
}

test_that("spatial(model = 'bym2') is identical to bym2() (areal)", {
  dat <- make_dat()
  adj <- chain_adj(20L)
  u <- tulpaObs:::.tobs_parse_formula(
    ~ elev + spatial(graph = adj, model = "bym2"), data = dat)
  r <- tulpaObs:::.tobs_parse_formula(~ elev + bym2(graph = adj), data = dat)
  expect_length(u$terms, 1L)
  expect_s3_class(u$terms[[1]], "tobs_spatial")
  expect_identical(u$terms[[1]]$type, "bym2")
  expect_identical(spec_fields(u$terms[[1]]), spec_fields(r$terms[[1]]))
  # the umbrella is stripped from the fixed-effects design like any term
  expect_setequal(attr(stats::terms(u$fe_formula), "term.labels"), "elev")
})

test_that("spatial(lon, lat, model = 'gp') is identical to gp() (continuous)", {
  dat <- make_dat()
  u <- tulpaObs:::.tobs_parse_formula(
    ~ spatial(lon, lat, model = "gp", prior_range = c(0.1, 0.05)), data = dat)
  r <- tulpaObs:::.tobs_parse_formula(
    ~ gp(lon, lat, prior_range = c(0.1, 0.05)), data = dat)
  expect_identical(u$terms[[1]]$type, "gp")
  expect_identical(spec_fields(u$terms[[1]]), spec_fields(r$terms[[1]]))
})

test_that("spatial() forwards per-model arguments and id", {
  dat <- make_dat()
  u <- tulpaObs:::.tobs_parse_formula(
    ~ spatial(lon, lat, model = "gp", nn = 5, id = "f",
              prior_range = c(0.1, 0.05)), data = dat)
  r <- tulpaObs:::.tobs_parse_formula(
    ~ gp(lon, lat, nn = 5, id = "f", prior_range = c(0.1, 0.05)), data = dat)
  expect_equal(u$terms[[1]]$nn, 5)
  expect_identical(u$terms[[1]]$id, "f")
  expect_identical(spec_fields(u$terms[[1]]), spec_fields(r$terms[[1]]))
})

test_that("spatial() rejects an unknown / non-spatial model", {
  dat <- make_dat()
  expect_error(
    tulpaObs:::.tobs_parse_formula(~ spatial(lon, lat, model = "re"), data = dat))
  expect_error(
    tulpaObs:::.tobs_parse_formula(~ spatial(lon, lat, model = "nope"), data = dat))
})

test_that("spatial() rejects a typo'd or wrong-model argument", {
  dat <- make_dat()
  adj <- chain_adj(20L)
  # continuous term: `...` would otherwise swallow the typo as a coordinate
  expect_error(
    tulpaObs:::.tobs_parse_formula(
      ~ spatial(lon, lat, model = "gp", nnn = 5), data = dat),
    "unknown argument")
  # areal argument passed to a continuous model
  expect_error(
    tulpaObs:::.tobs_parse_formula(
      ~ spatial(lon, lat, model = "gp", graph = adj), data = dat),
    "unknown argument")
  # continuous argument passed to an areal model (no `...`, named in message)
  expect_error(
    tulpaObs:::.tobs_parse_formula(
      ~ spatial(graph = adj, model = "bym2", nn = 5), data = dat),
    "unknown argument")
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

test_that("print methods on structured-term specs emit their term summaries", {
  dat <- make_dat()
  adj <- chain_adj(20L)

  p_icar <- tulpaObs:::.tobs_parse_formula(~ elev + icar(graph = adj), data = dat)
  expect_output(print(p_icar$terms[[1]]), "tobs spatial term: icar (20 units)",
                fixed = TRUE)

  p_re <- tulpaObs:::.tobs_parse_formula(~ re(observer), data = dat)
  expect_output(print(p_re$terms[[1]]), "tobs re term: intercept (iid, 4 groups)",
                fixed = TRUE)

  p_temp <- tulpaObs:::.tobs_parse_formula(~ temporal(year, type = "rw1"), data = dat)
  expect_output(print(p_temp$terms[[1]]), "tobs temporal term: rw1 (5 times)",
                fixed = TRUE)

  p_copy <- tulpaObs:::.tobs_parse_formula(~ elev + copy("u"), data = dat)
  expect_output(print(p_copy$terms[[1]]), "tobs copy term: -> u")
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
    list(psi = ~ gp(lon, lat, id = "u", prior_range = c(0.1, 0.05)),
         p = ~ forest + copy("u")),
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

# --- lme4-style bar syntax: sugar over re() --------------------------------

# Compare two tobs_re specs on their substantive fields (ignore the recorded
# call text, which differs between `(1|g)` and `re(g)`).
re_fields <- function(spec) {
  spec[setdiff(names(spec), "term_call")]
}

test_that("(1 | g) desugars to re(g, type = 'intercept')", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ elev + (1 | observer), data = dat)
  ref <- tulpaObs:::.tobs_parse_formula(~ elev + re(observer), data = dat)
  expect_length(bar$terms, 1L)
  expect_s3_class(bar$terms[[1]], "tobs_re")
  expect_identical(bar$terms[[1]]$type, "intercept")
  expect_identical(re_fields(bar$terms[[1]]), re_fields(ref$terms[[1]]))
  # the bar leaves elev in the fixed-effects design
  expect_setequal(attr(stats::terms(bar$fe_formula), "term.labels"), "elev")
})

test_that("(x | g) is a correlated random intercept + slope", {
  dat <- make_dat()
  bar  <- tulpaObs:::.tobs_parse_formula(~ (elev | observer), data = dat)
  ref  <- tulpaObs:::.tobs_parse_formula(
    ~ re(observer, type = "slope", covariate = cbind(elev)), data = dat)
  expect_identical(bar$terms[[1]]$type, "slope")
  expect_true(bar$terms[[1]]$correlated)
  expect_true(bar$terms[[1]]$intercept)
  # one slope column, named after the covariate
  expect_identical(colnames(bar$terms[[1]]$covariate), "elev")
  expect_identical(re_fields(bar$terms[[1]]), re_fields(ref$terms[[1]]))
})

test_that("(1 + x | g) equals (x | g)", {
  dat <- make_dat()
  a <- tulpaObs:::.tobs_parse_formula(~ (1 + elev | observer), data = dat)
  b <- tulpaObs:::.tobs_parse_formula(~ (elev | observer), data = dat)
  expect_identical(re_fields(a$terms[[1]]), re_fields(b$terms[[1]]))
})

test_that("(x || g) drops the intercept/slope correlation", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (elev || observer), data = dat)
  expect_identical(bar$terms[[1]]$type, "slope")
  expect_false(bar$terms[[1]]$correlated)
})

test_that("multiple bars yield multiple re terms alongside fixed effects", {
  dat <- make_dat()
  p <- tulpaObs:::.tobs_parse_formula(
    ~ forest + (1 | observer) + (elev || year), data = dat)
  expect_length(p$terms, 2L)
  expect_true(all(vapply(p$terms, inherits, logical(1), "tobs_re")))
  expect_setequal(attr(stats::terms(p$fe_formula), "term.labels"), "forest")
})

test_that("bars work in a process formula and tag to that process", {
  dat <- make_dat()
  parsed <- tulpaObs:::.tobs_parse_processes(
    list(psi = ~ elev, p = ~ effort + (1 | observer)),
    data = cbind(dat, effort = rnorm(nrow(dat))), env = environment()
  )
  resolved <- tulpaObs:::.tobs_resolve_terms(parsed$terms)
  expect_length(resolved, 1L)
  expect_s3_class(resolved[[1]]$spec, "tobs_re")
  expect_identical(resolved[[1]]$processes, 2L)
})

test_that("(0 + x | g) is a slope-only block (no intercept)", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (0 + elev | observer), data = dat)
  expect_length(bar$terms, 1L)
  expect_identical(bar$terms[[1]]$type, "slope")
  expect_false(bar$terms[[1]]$intercept)
  expect_identical(colnames(bar$terms[[1]]$covariate), "elev")
})

test_that("(1 + x + z | g) stacks multiple slopes in one block", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (1 + elev + forest | observer), data = dat)
  expect_length(bar$terms, 1L)
  expect_identical(bar$terms[[1]]$type, "slope")
  expect_true(bar$terms[[1]]$intercept)
  expect_identical(colnames(bar$terms[[1]]$covariate), c("elev", "forest"))
})

test_that("(1 | g:h) groups over the interaction factor", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (1 | observer:year), data = dat)
  expect_length(bar$terms, 1L)
  expect_identical(bar$terms[[1]]$type, "intercept")
  # interaction(observer, year, drop = TRUE) has one group per observed combo
  combos <- length(unique(interaction(dat$observer, dat$year, drop = TRUE)))
  expect_identical(bar$terms[[1]]$n_groups, combos)
})

test_that("(1 | g/h) expands to nested grouping (g and g:h)", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (1 | observer/year), data = dat)
  expect_length(bar$terms, 2L)
  expect_true(all(vapply(bar$terms, inherits, logical(1), "tobs_re")))
  expect_identical(bar$terms[[1]]$n_groups, nlevels(dat$observer))
  combos <- length(unique(interaction(dat$observer, dat$year, drop = TRUE)))
  expect_identical(bar$terms[[2]]$n_groups, combos)
})

test_that("the LHS slopes distribute across nested grouping factors", {
  dat <- make_dat()
  bar <- tulpaObs:::.tobs_parse_formula(~ (1 + elev | observer/year), data = dat)
  expect_length(bar$terms, 2L)
  expect_true(all(vapply(bar$terms, function(s) s$type == "slope", logical(1))))
  expect_true(all(vapply(bar$terms,
                         function(s) identical(colnames(s$covariate), "elev"),
                         logical(1))))
})
