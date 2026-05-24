# LOO-stacked ensembles: tobs_stack() weights, input validation, the combined
# fitted/predict methods, and the n.seeds sugar in tobs(). Stacking weights and
# PSIS-LOO come from the `loo` package; tests are skipped when it is absent.
# Laplace fits are deterministic and fast, so the core tobs_stack() assertions
# run off-CRAN; the n.seeds path uses a stochastic route and is gated.

sim_occu <- function(seed = 1, N = 80L, J = 4L, psi = 0.6, p = 0.45) {
  set.seed(seed)
  z <- rbinom(N, 1, psi)
  y <- matrix(rbinom(N * J, 1, rep(z, times = J) * p), N, J)
  list(y = y, data = data.frame(x = rnorm(N), w = rnorm(N)))
}

fit_lap <- function(form, sim) {
  tobs(form, data = sim$data, y = sim$y, detection = ~ 1,
       family = occu(), method = "laplace", control = list(verbose = FALSE))
}

test_that(".tobs_pointwise_loglik returns a finite [n_draws x n_obs] matrix", {
  sim <- sim_occu()
  f <- fit_lap(~ x, sim)
  ll <- tulpaObs:::.tobs_pointwise_loglik(f)
  expect_true(is.matrix(ll))
  expect_equal(nrow(ll), nrow(f$draws))
  expect_equal(ncol(ll), nrow(sim$y))
  expect_true(all(is.finite(ll)))
})

test_that(".tobs_pointwise_loglik errors on an unknown model_type", {
  fake <- list(draws = matrix(0, 10, 2), model = list(model_type = "bogus"))
  expect_error(tulpaObs:::.tobs_pointwise_loglik(fake), "not implemented")
})

test_that(".tobs_pointwise_loglik covers every family (shape + finite)", {
  # community: site x species replicated occupancy
  set.seed(2); ns <- 18L; nsp <- 3L; mv <- 3L; yl <- list()
  for (s in seq_len(nsp)) {
    zz <- rbinom(ns, 1, 0.5); ys <- matrix(0L, ns, mv)
    for (i in seq_len(ns)) if (zz[i]) ys[i, ] <- rbinom(mv, 1, 0.4)
    yl[[paste0("sp", s)]] <- ys
  }
  fc <- tobs(~ 1, data = data.frame(x = rnorm(ns)), family = ms_occu(),
             detection = ~ 1, y = yl, species = TRUE, method = "laplace",
             control = list(verbose = FALSE))
  ll <- tulpaObs:::.tobs_pointwise_loglik(fc)
  expect_equal(ncol(ll), ns * nsp)
  expect_true(all(is.finite(ll)))

  # dynamic: per-site HMM marginal
  set.seed(3); ns <- 24L; Tn <- 3L; ya <- array(0L, c(ns, mv, Tn))
  z <- matrix(0, ns, Tn); z[, 1] <- rbinom(ns, 1, 0.6)
  for (t in 2:Tn) z[, t] <- z[, t-1]*(1-rbinom(ns,1,0.1)) + (1-z[, t-1])*rbinom(ns,1,0.2)
  for (i in seq_len(ns)) for (t in seq_len(Tn)) if (z[i, t]) ya[i, , t] <- rbinom(mv, 1, 0.5)
  fd <- tobs(~ 1, data = data.frame(x = rnorm(ns)), family = dyn_occu(),
             detection = ~ 1, y = ya, col_formula = ~ 1, ext_formula = ~ 1,
             method = "laplace", control = list(verbose = FALSE))
  lld <- tulpaObs:::.tobs_pointwise_loglik(fd)
  expect_equal(ncol(lld), ns)
  expect_true(all(is.finite(lld)))

  # cover hurdle (beta): occurrence Bernoulli + beta positive density
  set.seed(7); N <- 70L; xx <- rnorm(N)
  occ <- rbinom(N, 1, plogis(-0.3 + 0.7*xx)); mu <- plogis(-0.5 + 0.4*xx)
  yb <- numeric(N); yb[occ==1] <- rbeta(sum(occ), mu[occ==1]*8, (1-mu[occ==1])*8)
  fb <- tobs(~ x, data = data.frame(x = xx), family = cover("beta"), y = yb,
             method = "laplace", control = list(verbose = FALSE))
  llb <- tulpaObs:::.tobs_pointwise_loglik(fb)
  expect_equal(ncol(llb), N)
  expect_true(all(is.finite(llb)))
})

test_that("cover stacking errors on the nested-joint path", {
  # the nested-joint cover fit lacks per-arm mode/Hessian
  fake <- structure(list(encoding = list(), occ = list(mode = 1)),
                    class = c("cover_fit", "tobs_fit"))
  expect_error(tulpaObs:::.tobs_ploglik_cover(fake), "nested-joint")
})

test_that("tobs_stack() returns weights that sum to 1 over named members", {
  skip_if_not_installed("loo")
  sim <- sim_occu()
  f1 <- fit_lap(~ x,     sim)
  f2 <- fit_lap(~ x + w, sim)
  ens <- tobs_stack(simple = f1, full = f2)

  expect_s3_class(ens, "tobs_stack")
  expect_named(ens$weights, c("simple", "full"))
  expect_equal(sum(ens$weights), 1, tolerance = 1e-6)
  expect_true(all(ens$weights >= 0))
  expect_setequal(ens$comparison$model, c("simple", "full"))
  expect_length(ens$fits, 2L)
})

test_that("tobs_stack() accepts a single list argument", {
  skip_if_not_installed("loo")
  sim <- sim_occu()
  ens <- tobs_stack(list(fit_lap(~ x, sim), fit_lap(~ x + w, sim)))
  expect_s3_class(ens, "tobs_stack")
  expect_named(ens$weights, c("model1", "model2"))
})

test_that("tobs_stack() rejects bad inputs", {
  skip_if_not_installed("loo")
  sim <- sim_occu()
  f1 <- fit_lap(~ x, sim)
  expect_error(tobs_stack(f1), "at least two")
  expect_error(tobs_stack(f1, "not a fit"), "tobs_fit")

  # members fit to a different number of observations
  sim2 <- sim_occu(seed = 2, N = 50L)
  f2 <- fit_lap(~ x, sim2)
  expect_error(tobs_stack(f1, f2), "same observations")
})

test_that("fitted.tobs_stack is the weight-combined member fitted values", {
  skip_if_not_installed("loo")
  sim <- sim_occu()
  f1 <- fit_lap(~ x,     sim)
  f2 <- fit_lap(~ x + w, sim)
  ens <- tobs_stack(f1, f2)
  w <- ens$weights

  fv <- fitted(ens)
  expect_named(fv, c("psi", "p", "z"))
  expect_length(fv$psi, nrow(sim$y))
  expect_true(all(fv$psi >= 0 & fv$psi <= 1))

  manual <- w[[1]] * fitted(f1)$psi + w[[2]] * fitted(f2)$psi
  expect_equal(unname(fv$psi), unname(manual), tolerance = 1e-8)

  # predict() with no X.0 is the in-sample fitted mixture
  expect_equal(predict(ens), fv)
})

test_that("predict.tobs_stack with X.0 needs matching member designs", {
  skip_if_not_installed("loo")
  sim <- sim_occu()
  f1 <- fit_lap(~ x,     sim)   # 2 occ coefs
  f2 <- fit_lap(~ x + w, sim)   # 3 occ coefs
  ens <- tobs_stack(f1, f2)
  X0 <- model.matrix(~ x, sim$data)[1:5, , drop = FALSE]
  expect_error(predict(ens, X.0 = X0), "same occupancy design")
})

test_that("n.seeds builds a tobs_stack of K members and predicts out-of-sample", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("loo")
  sim <- sim_occu()
  ens <- tobs(~ x, data = sim$data, y = sim$y, detection = ~ 1,
              family = occu(), method = "laplace_gibbs",
              control = list(n.seeds = 3L, n.gibbs = 5L, seed = 42L,
                             verbose = FALSE))
  expect_s3_class(ens, "tobs_stack")
  expect_length(ens$fits, 3L)
  expect_named(ens$weights, c("seed42", "seed43", "seed44"))
  expect_equal(sum(ens$weights), 1, tolerance = 1e-6)

  # members share the occupancy design, so X.0 prediction works
  X0 <- model.matrix(~ x, sim$data)[1:5, , drop = FALSE]
  pr <- predict(ens, X.0 = X0)
  expect_equal(nrow(pr), 5L)
  expect_true(all(c("mean", "sd", "q2.5", "q50", "q97.5") %in% names(pr)))
  expect_true(all(pr$mean >= 0 & pr$mean <= 1))
})

test_that("n.seeds is rejected on the deterministic Laplace route", {
  sim <- sim_occu()
  expect_error(
    tobs(~ x, data = sim$data, y = sim$y, detection = ~ 1, family = occu(),
         method = "laplace", control = list(n.seeds = 3L)),
    "n.seeds", fixed = TRUE
  )
})
