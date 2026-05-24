# Control-option validation. `tobs()` rejects control names that do not apply
# to the resolved method (misapplied or misspelled) instead of silently
# swallowing them into `.tobs_fit_model()`'s `...`. The validator is exercised
# directly via the internal helper (no fitting needed) and once end-to-end
# through tobs() to confirm it is wired ahead of dispatch.

route <- function(method) tulpaObs:::.tobs_resolve_method(method, occu())
check <- function(control, method) {
  tulpaObs:::.tobs_validate_control(control, route(method))
}

test_that("Laplace rejects sampler controls", {
  for (key in c("n.chains", "n.iter", "n.warmup", "n.thin", "n.threads",
                "adapt.delta", "max.treedepth", "sigma.re.scale")) {
    ctrl <- setNames(list(1L), key)
    expect_error(check(ctrl, "laplace"), key, fixed = TRUE)
  }
})

test_that("Laplace rejects seed (deterministic route), points at stochastic methods", {
  err <- expect_error(check(list(seed = 1L), "laplace"))
  expect_match(conditionMessage(err), "seed", fixed = TRUE)
  # the hint should name the routes that actually consume a seed
  expect_match(conditionMessage(err), "laplace_gibbs")
  expect_match(conditionMessage(err), "nuts")
})

test_that("NUTS rejects Laplace-only controls", {
  for (key in c("max.iter", "tol", "damping", "n.gibbs", "n.imputations")) {
    ctrl <- setNames(list(1L), key)
    expect_error(check(ctrl, "nuts"), key, fixed = TRUE)
  }
})

test_that("stochastic-correction routes accept seed + correction draw counts", {
  expect_silent(check(list(seed = 1L, n.gibbs = 5L), "laplace_gibbs"))
  expect_silent(check(list(seed = 1L, n.imputations = 5L), "laplace_mi"))
  # but still reject sampler controls
  expect_error(check(list(n.chains = 2L), "laplace_gibbs"), "n.chains", fixed = TRUE)
})

test_that("n.seeds is allowed only on stochastic routes", {
  expect_silent(check(list(n.seeds = 3L), "nuts"))
  expect_silent(check(list(n.seeds = 3L), "laplace_gibbs"))
  expect_silent(check(list(n.seeds = 3L), "laplace_mi"))
  # deterministic routes reject it (seed-variants would be identical)
  err <- expect_error(check(list(n.seeds = 3L), "laplace"))
  expect_match(conditionMessage(err), "n.seeds", fixed = TRUE)
  expect_match(conditionMessage(err), "nuts")
  expect_error(check(list(n.seeds = 3L), "nested_laplace"), "n.seeds", fixed = TRUE)
})

test_that("valid controls pass for each engine family", {
  expect_silent(check(list(max.iter = 50L, tol = 1e-5, damping = 0.5,
                           sigma.beta = 5, verbose = FALSE), "laplace"))
  expect_silent(check(list(max.iter = 50L), "nested_laplace"))
  expect_silent(check(list(n.chains = 4L, n.iter = 2000L, n.warmup = 1000L,
                           n.thin = 2L, n.threads = 2L, adapt.delta = 0.9,
                           max.treedepth = 12L, seed = 7L, sigma.beta = 5,
                           sigma.re.scale = 1, verbose = FALSE), "nuts"))
  expect_silent(check(list(), "laplace"))
})

test_that("unknown control names error with a fuzzy suggestion", {
  err <- expect_error(check(list(niter = 5000L), "nuts"))
  expect_match(conditionMessage(err), "not a known control option")
  expect_match(conditionMessage(err), "n.iter")

  expect_error(check(list(nonsense = 1L), "nuts"), "not a known control option")
})

test_that("unnamed control list is rejected", {
  expect_error(check(list(2000L), "nuts"), "fully named list")
})

test_that("all offending names are reported together", {
  # both max.iter and tol are Laplace-only -> both invalid under nuts
  err <- expect_error(check(list(max.iter = 10L, tol = 1e-5), "nuts"))
  expect_match(conditionMessage(err), "max.iter")
  expect_match(conditionMessage(err), "tol")
})

test_that("validation fires through tobs() ahead of dispatch", {
  # data/formula are immaterial: the control error precedes the detection/y
  # checks in the occu dispatcher.
  expect_error(
    tobs(~ 1, data = data.frame(x = 1), family = occu(),
         method = "laplace", control = list(n.chains = 2L)),
    "n.chains", fixed = TRUE
  )
})
