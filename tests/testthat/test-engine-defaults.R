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

test_that("the Laplace-EM prior / regularization scales are one value each", {
  expect_identical(.tobs_default("laplace", "sigma.beta"), 5)
  # A loading prior width and an LKJ shape are engine-level values of the same
  # kind as the ridge (gcol33/tulpaObs#189). sd.load in particular had three
  # copies that must move together: the auto-K ladder selects a rank by marginal
  # evidence under it, so a value drifting between the selection fit and the
  # final fit selects a rank the fit does not use.
  expect_identical(.tobs_default("laplace", "sd.load"), 1.0)
  expect_identical(.tobs_default("laplace", "re.lkj"), 1.5)
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

# =========================================================================== #
# Past the table: what a FITTER runs (gcol33/tulpaObs#188)                     #
#                                                                             #
# Everything above asserts what the table RESOLVES. None of it could see a     #
# fitter that ignores the table, and 21 of them did -- each restating the      #
# knobs in its own formals, three with values that disagreed. The assertions   #
# below are on the fitters.                                                   #
# =========================================================================== #

# Every fitter that samples. A new one added without routing through the table
# fails the structural test below rather than shipping a silently different
# chain, which is the failure mode the whole table exists to prevent.
.ed_sampler_fitters <- c(
  # single-species entry (.tobs_fit_model forwards explicit values to these)
  ".tobs_fit_abun_nuts", ".tobs_fit_abun_nuts_spatial",
  ".tobs_fit_distance_nuts", ".tobs_fit_distance_nuts_spatial",
  ".tobs_fit_dyn_abun_nuts", ".tobs_fit_dyn_abun_nuts_spatial",
  ".tobs_fit_dyn_abun_nuts_temporal",
  ".tobs_fit_fp_occu_nuts", ".tobs_fit_fp_occu_nuts_spatial",
  ".tobs_fit_removal_nuts", ".tobs_fit_removal_nuts_spatial",
  ".tobs_fit_occu_pg_gibbs", ".tobs_fit_occu_pg_gibbs_spatial",
  # own dispatch
  ".tobs_fit_cover_nuts", ".tobs_fit_occu_cover_nuts",
  ".tobs_fit_occu_cover_nuts_spatial",
  ".tobs_fit_occu_multiscale_cover_nuts",
  # community
  ".tobs_fit_ms_abun_nuts", ".tobs_fit_ms_abun_nuts_spatial",
  ".tobs_fit_ms_count_nuts", ".tobs_fit_ms_occu_nuts",
  ".tobs_fit_ms_dyn_occu_nuts", ".tobs_fit_ms_int_occu_nuts",
  ".tobs_fit_ms_occu_cover_nuts",
  ".tobs_fit_ms_occu_pg_gibbs", ".tobs_fit_ms_count_pg_gibbs",
  ".tobs_fit_ms_dyn_occu_pg_gibbs", ".tobs_fit_ms_int_occu_pg_gibbs",
  ".tobs_fit_t_occu_pg_gibbs"
)

.ed_sampler_knobs <- c("n.iter", "n.warmup", "n.chains", "n.thin",
                       "max.treedepth", "adapt.delta", "seed",
                       "sigma.beta", "sigma.logr")

test_that("no sampler fitter states its own default for a sampler knob", {
  # A literal in a fitter's formals is a SECOND answer to "what is the default
  # n.iter". For the fitters below `.tobs_fit_model()` reaches, it was an
  # unreachable second answer (that entry forwards explicit values, so the
  # formal never applied); for the rest it was the live one and disagreed with
  # the table. Either way the fix is the same: the formal is a NULL sentinel and
  # the table is the only answer.
  offenders <- character(0)
  for (nm in .ed_sampler_fitters) {
    f <- get(nm, envir = asNamespace("tulpaObs"))
    fm <- formals(f)
    for (k in intersect(.ed_sampler_knobs, names(fm))) {
      if (!is.null(fm[[k]]) && !identical(fm[[k]], quote(expr = ))) {
        offenders <- c(offenders, sprintf("%s(%s = %s)", nm, k,
                                          deparse(fm[[k]])))
      }
    }
  }
  expect_identical(offenders, character(0))
})

test_that("every sampler fitter resolves the table as its first act", {
  # The NULL sentinels above are only correct because something fills them.
  for (nm in .ed_sampler_fitters) {
    src <- paste(deparse(body(get(nm, envir = asNamespace("tulpaObs")))),
                 collapse = " ")
    expect_true(grepl(".tobs_fill_sampler", src, fixed = TRUE),
                info = nm)
  }
})

test_that(".tobs_fill_sampler fills only the knobs a fitter declares", {
  # One call serves fitters with different knob sets: a pg_gibbs fitter has no
  # adapt.delta, and filling one into its frame would invent a knob it never
  # reads.
  e <- new.env()
  assign("n.iter", NULL, envir = e)
  assign("seed", 99L, envir = e)                       # caller-supplied: wins
  .tobs_fill_sampler(e, "nuts")
  expect_identical(get("n.iter", envir = e), 1000L)
  expect_identical(get("seed", envir = e), 99L)
  expect_false(exists("adapt.delta", envir = e, inherits = FALSE))
})

test_that("the single-species entry's departures are on the record", {
  # `.tobs_fit_model()` serves occupancy and every observation family through
  # one entry, and its departures belong to the ENTRY, not to the nine families
  # passing through it. A wider ridge, a looser adaptation target and a
  # different stream seed are kept because each is the value the single-species
  # recovery / coverage tests were calibrated against.
  d <- .tobs_single_species_defaults("nuts")
  expect_identical(d$sigma.beta, 10)
  expect_identical(d$adapt.delta, 0.8)
  expect_identical(d$seed, 42)
  # n.iter is NOT among them. It read 2000 only because commit 8975470 changed
  # `n.iter` from meaning the TOTAL run to meaning kept post-warmup draws on
  # exactly these paths and left the literal alone -- a default that had always
  # kept 2000 - 1000 = 1000 draws silently began keeping 2000. The table's 1000
  # is the count these paths were calibrated at.
  expect_identical(d$n.iter, 1000L)
  expect_identical(d$n.warmup, 1000L)
  # The Laplace / nested-Laplace routes share the entry's wider ridge and carry
  # no chain knobs at all.
  expect_identical(.tobs_single_species_defaults("laplace")$sigma.beta, 10)
  expect_identical(.tobs_single_species_defaults("nested_laplace")$sigma.beta, 10)
  expect_false("n.iter" %in% names(.tobs_single_species_defaults("laplace")))
})

test_that("occu(pg_gibbs) keeps the same chain as its ms_* siblings", {
  # It kept 1000 draws where every community sibling kept 1500, with no reason
  # recorded -- verbatim the failure this table exists to prevent. The entry now
  # carries NO pg_gibbs departure, so the two resolve identically.
  ss <- .tobs_single_species_defaults("pg_gibbs")
  ms <- .tobs_engine_defaults("pg_gibbs")
  expect_identical(ss, ms)
  expect_identical(ss$n.iter, 3000L)
  expect_identical(ss$n.warmup, 1500L)
})

test_that("the pg_gibbs n.iter convention is the opposite of the NUTS one", {
  # Under "nuts", n.iter is the KEPT count and the run is n.iter + n.warmup.
  # Under "pg_gibbs", n.iter is the TOTAL and warmup comes out of it. Every
  # pg_gibbs fitter computes its kept count as
  # length(seq.int(n.warmup + 1L, n.iter, by = n.thin)); this pins the profile
  # against that formula so the two cannot drift.
  d <- .tobs_engine_defaults("pg_gibbs")
  kept <- length(seq.int(d$n.warmup + 1L, d$n.iter, by = d$n.thin))
  expect_identical(kept, 1500L)
  expect_lt(d$n.warmup, d$n.iter)   # a total below its own warmup keeps nothing
})

test_that("n.quad names one control across several marginals, each enumerated", {
  # Not one number, and deliberately so: each route integrates a different
  # marginal over a different latent dimension. What was wrong was that no
  # reader could find out which number applied (gcol33/tulpaObs#189).
  expect_identical(.tobs_n_quad("re_aghq"), 9L)
  expect_identical(.tobs_n_quad("ms_nmix"), 1L)
  expect_identical(.tobs_n_quad("ms_nmix_scalar"), 2L)
  expect_identical(.tobs_n_quad("ms_occu_cover"), 5L)
  expect_identical(.tobs_n_quad("community_latent"), 5L)
  # The lognormal per-unit cover marginal is closed form, so it needs no
  # quadrature; the beta one is not.
  expect_identical(.tobs_n_quad("cover_latent_lognormal"), 1L)
  expect_identical(.tobs_n_quad("cover_latent_beta"), 15L)
  # An unknown route errors rather than returning a plausible NULL.
  expect_error(.tobs_n_quad("no_such_route"), "no n.quad default")
})
