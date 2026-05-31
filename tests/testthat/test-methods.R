.fit_simple <- function(formula = ~ elev, det = ~ 1, n = 50, seed = 42,
                       method = "laplace") {
  set.seed(seed)
  d <- data.frame(elev = rnorm(n))
  psi <- plogis(0.5 + 0.5 * d$elev)
  z <- rbinom(n, 1, psi)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- tobs(formula = formula, data = d, family = occu(),
              detection = det, y = y, method = method,
              control = list(verbose = FALSE))
  list(fit = fit, y = y, d = d, n = n)
}

test_that("S3 methods work on single-season fit", {
  res <- .fit_simple()
  fit <- res$fit; y <- res$y; n <- res$n

  cf <- coef(fit)
  expect_type(cf, "list")
  expect_length(cf$psi, 2)
  expect_length(cf$p, 1)

  ci <- confint(fit)
  expect_true(nrow(ci) >= 3)

  V <- vcov(fit)
  expect_true(all(diag(V) > 0))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")

  expect_equal(nobs(fit), sum(y >= 0, na.rm = TRUE))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "p", "z"))
  expect_length(fv$psi, n)
  expect_true(all(fv$psi >= 0 & fv$psi <= 1))

  r <- residuals(fit)
  expect_named(r, c("occ", "det"))
  expect_length(r$occ, n)
  expect_equal(dim(r$det), c(n, 3))

  y_sim <- simulate(fit, nsim = 1, seed = 1)
  expect_equal(dim(y_sim), dim(y))

  pred <- predict(fit)
  expect_named(pred, c("psi", "p", "z"))

  X0 <- model.matrix(~ elev, data.frame(elev = c(-1, 0, 1)))
  pred_dm <- predict(fit, X.0 = X0)
  expect_equal(nrow(pred_dm), 3)
  expect_true(all(pred_dm$mean >= 0 & pred_dm$mean <= 1))

  td <- tulpa::tidy(fit)
  expect_s3_class(td, "data.frame")

  gl <- tulpa::glance(fit)
  expect_s3_class(gl, "data.frame")

  re <- tulpa::ranef(fit)
  expect_true(is.data.frame(re) || is.list(re))
})

test_that("non-NUTS fits report NA sampler diagnostics, NUTS reports numeric", {
  # tulpaObs#17: a Laplace / nested-Laplace fit ran no HMC trajectory, so the
  # NUTS-only sampler-health fields (acceptance, divergence, tree depth, step
  # size) must be NA rather than the constants 1 / 0 / 0 / 0 -- otherwise a user
  # checking sampler health reads "no sampler ran" as "sampler ran cleanly".
  fit_lap <- .fit_simple(method = "laplace")$fit
  expect_identical(fit_lap$method, "laplace")
  expect_true(all(is.na(fit_lap$accept_prob)))
  expect_true(all(is.na(fit_lap$divergent)))
  expect_true(all(is.na(fit_lap$treedepth)))
  expect_true(is.na(fit_lap$epsilon))

  # Areal nested-Laplace N-mixture build path (abun.R) -- the same rule.
  set.seed(11)
  side <- 5L; ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) adj[i, j] <- 1L
  phi <- as.numeric(scale(rnorm(ng))) * 0.5; phi <- phi - mean(phi)
  x_ab <- rnorm(ng)
  N <- rpois(ng, exp(log(5) + 0.5 * x_ab + phi))
  yA <- matrix(NA_integer_, ng, 5L)
  for (i in seq_len(ng)) yA[i, ] <- rbinom(5L, N[i], plogis(0.4))
  fit_nl <- tobs(~ abund_cov1 + icar(graph = adj),
                 data = data.frame(abund_cov1 = x_ab), family = abun(),
                 detection = ~ 1, y = yA, method = "nested_laplace",
                 control = list(verbose = FALSE))
  expect_identical(fit_nl$method, "nested_laplace")
  expect_true(all(is.na(fit_nl$accept_prob)))
  expect_true(all(is.na(fit_nl$divergent)))
  expect_true(all(is.na(fit_nl$treedepth)))
  expect_true(is.na(fit_nl$epsilon))

  # A NUTS fit, by contrast, carries real numeric diagnostics.
  fit_nuts <- .fit_simple(method = "nuts")$fit
  expect_identical(fit_nuts$method, "nuts")
  expect_true(any(is.finite(fit_nuts$accept_prob)))
  expect_false(all(is.na(fit_nuts$divergent)))
})

test_that("WAIC works on single-season fit", {
  res <- .fit_simple(formula = ~ elev, n = 30, seed = 42)
  w <- tobs_waic(res$fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(w$elpd))
  expect_true(w$p_waic >= 0)
})

test_that("PPC works on single-season fit", {
  res <- .fit_simple(formula = ~ 1, n = 30, seed = 42)
  ppc <- tobs_ppc(res$fit, n.samples = 50)
  expect_length(ppc$fit.y, 50)
  expect_length(ppc$fit.y.rep, 50)
  expect_true(ppc$bayesian.p >= 0 && ppc$bayesian.p <= 1)
})

test_that("compare_models works", {
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)

  fit1 <- tobs(~ 1, d, family = occu(), detection = ~ 1, y = y,
               method = "nuts",
               control = list(n.iter = 200, n.warmup = 100, seed = 42, verbose = FALSE))
  fit2 <- tobs(~ x, d, family = occu(), detection = ~ 1, y = y,
               method = "nuts",
               control = list(n.iter = 200, n.warmup = 100, seed = 42, verbose = FALSE))

  comp <- tulpa::compare_models(null = fit1, elev = fit2)
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 2)
})

test_that("simulation functions work", {
  sim <- simulate_occu(N = 20, J = 3, seed = 42)
  expect_equal(dim(sim$y), c(20, 3))
  expect_equal(nrow(sim$data), 20)

  sim_ms <- simulate_ms_occu(N = 10, J = 3, n_species = 3, seed = 42)
  expect_equal(dim(sim_ms$y), c(10, 3, 3))

  sim_t <- simulate_dyn_occu(N = 10, J = 3, n_seasons = 4, seed = 42)
  expect_equal(dim(sim_t$y), c(10, 3, 4))
})

test_that("tobs_data long format conversion works", {
  df <- expand.grid(site = 1:5, visit = 1:3)
  df$detected <- rbinom(15, 1, 0.3)
  df$effort <- rnorm(15)
  df$habitat <- rep(c("forest", "grass", "forest", "grass", "forest"), each = 3)

  od <- tobs_data(df, y = "detected", site = "site", visit = "visit",
                  occ.covs = "habitat", det.covs = "effort")
  expect_s3_class(od, "tobs_data")
  expect_equal(nrow(od$y), 5)
  expect_equal(ncol(od$y), 3)
  expect_equal(nrow(od$occ.covs), 5)
})
