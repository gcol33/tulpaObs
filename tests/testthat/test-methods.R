.fit_simple <- function(formula = ~ elev, det = ~ 1, n = 50, seed = 42,
                       engine = "laplace") {
  set.seed(seed)
  d <- data.frame(elev = rnorm(n))
  psi <- plogis(0.5 + 0.5 * d$elev)
  z <- rbinom(n, 1, psi)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- tobs(formula = formula, data = d, family = occu(),
              detection = det, y = y, engine = engine,
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

test_that("WAIC works on single-season fit", {
  res <- .fit_simple(formula = ~ elev, n = 30, seed = 42)
  w <- waicOccu(res$fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(w$elpd))
  expect_true(w$p_waic >= 0)
})

test_that("PPC works on single-season fit", {
  res <- .fit_simple(formula = ~ 1, n = 30, seed = 42)
  ppc <- ppcOccu(res$fit, n.samples = 50)
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
               engine = "nuts",
               control = list(iter = 200, warmup = 100, seed = 42, verbose = FALSE))
  fit2 <- tobs(~ x, d, family = occu(), detection = ~ 1, y = y,
               engine = "nuts",
               control = list(iter = 200, warmup = 100, seed = 42, verbose = FALSE))

  comp <- tulpa::compare_models(null = fit1, elev = fit2)
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 2)
})

test_that("simulation functions work", {
  sim <- simulate_occu(N = 20, J = 3, seed = 42)
  expect_equal(dim(sim$y), c(20, 3))
  expect_equal(nrow(sim$data), 20)

  sim_ms <- simMsOcc(N = 10, J = 3, n_species = 3, seed = 42)
  expect_equal(dim(sim_ms$y), c(10, 3, 3))

  sim_t <- simTOcc(N = 10, J = 3, n_seasons = 4, seed = 42)
  expect_equal(dim(sim_t$y), c(10, 3, 4))
})

test_that("occu_data long format conversion works", {
  df <- expand.grid(site = 1:5, visit = 1:3)
  df$detected <- rbinom(15, 1, 0.3)
  df$effort <- rnorm(15)
  df$habitat <- rep(c("forest", "grass", "forest", "grass", "forest"), each = 3)

  od <- occu_data(df, y = "detected", site = "site", visit = "visit",
                  occ.covs = "habitat", det.covs = "effort")
  expect_s3_class(od, "occu_data")
  expect_equal(nrow(od$y), 5)
  expect_equal(ncol(od$y), 3)
  expect_equal(nrow(od$occ.covs), 5)
})
