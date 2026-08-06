# The engine-defaults table (#183). These assertions pin the resolved default
# for each (engine, family), so an edit to the table shows up as a failing test
# naming the knob that moved rather than as one family quietly sampling a
# different chain from its siblings. Values here are the ones that were written
# out at the call sites before the move.

test_that("the NUTS profile resolves to the shipped defaults", {
  d <- .tobs_engine_defaults("nuts")
  expect_identical(d$n.iter, 1000L)
  expect_identical(d$n.warmup, 1000L)
  expect_identical(d$n.chains, 1L)
  expect_identical(d$max.treedepth, 10L)
  expect_identical(d$adapt.delta, 0.9)
  expect_identical(d$seed, 1L)
  expect_identical(d$sigma.beta, 5)
  expect_identical(d$sigma.logr, 1.5)
})

test_that("the pg_gibbs profile resolves to the shipped defaults", {
  d <- .tobs_engine_defaults("pg_gibbs")
  expect_identical(d$n.iter, 3000L)
  expect_identical(d$n.warmup, 1500L)
  expect_identical(d$n.chains, 2L)
  expect_identical(d$n.thin, 1L)
  expect_identical(d$seed, 1L)
  expect_identical(d$sigma.beta, 2.5)
  # adapt.delta / max.treedepth are HMC knobs; a Gibbs sweep has neither, and
  # that absence is the answer to "what is the default adapt.delta here".
  expect_false("adapt.delta" %in% names(d))
  expect_false("max.treedepth" %in% names(d))
})

test_that("the log-link community families carry the wider NUTS coefficient prior", {
  # Count / abundance sample on a log link and use sigma.beta = 10; the
  # logit-link occupancy families use 5. This mirrors the C++ model defaults.
  for (fam in c("ms_count", "jsdm", "ms_abun")) {
    expect_identical(.tobs_default("nuts", "sigma.beta", fam), 10,
                     info = fam)
  }
  for (fam in c("ms_occu", "ms_dyn_occu", "ms_int_occu", "ms_occu_cover")) {
    expect_identical(.tobs_default("nuts", "sigma.beta", fam), 5, info = fam)
  }
  # Everything else stays on the engine profile.
  expect_identical(.tobs_engine_defaults("nuts", "ms_count")$n.iter, 1000L)
})

test_that("the spatial-factor community sampler keeps its own warmup + adaptation", {
  d <- .tobs_engine_defaults("nuts", "ms_occu_cover_spatial")
  expect_identical(d$n.warmup, 500L)
  expect_identical(d$adapt.delta, 0.95)
  expect_identical(d$n.iter, 1000L)          # the rest is the shared profile
  expect_identical(d$sigma.beta, 5)
})

test_that("the Laplace-EM ridge is one value", {
  expect_identical(.tobs_default("laplace", "sigma.beta"), 5)
  # max.iter / tol are per-route, not profile-shaped, so the table does not
  # claim them; asking must error rather than return a plausible NULL.
  expect_error(.tobs_default("laplace", "max.iter"), "no default")
  expect_error(.tobs_default("nuts", "tol"), "no default")
})

test_that("a user control value always wins over the profile", {
  ctl <- .tobs_control_defaults(list(n.iter = 42L, adapt.delta = 0.99),
                                "nuts", "ms_occu")
  expect_identical(ctl$n.iter, 42L)
  expect_identical(ctl$adapt.delta, 0.99)
  expect_identical(ctl$n.warmup, 1000L)       # unset knobs still fill
  expect_identical(ctl$sigma.beta, 5)
})

test_that("an unset or NULL knob counts as absent, matching the %||% it replaced", {
  expect_identical(.tobs_control_defaults(list(seed = NULL), "nuts")$seed, 1L)
  expect_identical(.tobs_control_defaults(NULL, "pg_gibbs")$n.chains, 2L)
  # A knob the caller set to a falsy-but-real value is kept.
  expect_identical(.tobs_control_defaults(list(seed = 0L), "nuts")$seed, 0L)
})

test_that("an engine with no profile leaves control untouched", {
  ctl <- list(a = 1)
  expect_identical(.tobs_control_defaults(ctl, "nested_laplace"), ctl)
  expect_identical(.tobs_engine_defaults("no_such_engine"), list())
})
