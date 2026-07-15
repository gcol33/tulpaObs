# test-community-methods.R - predict() / residuals() for the community-occupancy
# families and fitted() / predict() / residuals() for jsdm(), which previously
# errored (a whole family cluster with no post-fit surface).

test_that("predict() / residuals() work on ms_occu", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 50, J = 3, n_species = 5, seed = 1)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:5), method = "laplace",
              control = list(verbose = FALSE))
  psi <- predict(fit)
  expect_equal(dim(psi), c(50L, 5L))
  expect_true(all(psi >= 0 & psi <= 1))
  expect_equal(colnames(psi), fit$model$species_names)
  # newdata recomputes from the per-species coefficients.
  pn <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
  expect_equal(dim(pn), c(3L, 5L))
  # detection prediction.
  pd <- predict(fit, type = "detection")
  expect_equal(dim(pd), c(50L, 5L))
  rr <- residuals(fit)
  expect_equal(dim(rr$occ), c(50L, 5L))
  expect_true(all(is.finite(rr$occ)))
})

test_that("predict() / residuals() work on ms_dyn_occu and ms_int_occu", {
  skip_on_cran()
  skip_if_fast()
  sd_ <- simulate_ms_dyn_occu(N = 40, J = 3, n_species = 5, n_seasons = 3,
                              gamma = 0.2, epsilon = 0.1, seed = 2)
  fd <- tobs(~ 1, data = sd_$data, family = ms_dyn_occu(), detection = ~ 1,
             y = sd_$y, species = paste0("sp", 1:5), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(dim(predict(fd)), c(40L, 5L))
  expect_true(all(is.finite(residuals(fd)$occ)))

  si <- simulate_ms_int_occu(N = 60, J = c(3, 4), n_species = 5, seed = 3)
  fi <- tobs(~ 1, data = si$data, family = ms_int_occu(), detection = ~ 1,
             y = si$y, species = paste0("sp", 1:5), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(dim(predict(fi)), c(60L, 5L))
  expect_true(all(is.finite(residuals(fi)$occ)))
})

test_that("fitted() / predict() / residuals() work on jsdm", {
  skip_on_cran()
  # Since gcol33/tulpaObs#121 jsdm() IS the community GLMM with a logit link, so
  # it shares the ms_count() post-fit surface: fitted()/predict() return the
  # per-(site, species) mean on the response scale (a probability here) under
  # `$mu`, and residuals() returns `$mu`.
  sim <- simulate_ms_occu(N = 40, J = 1, n_species = 5, seed = 4)
  yj  <- apply(sim$y, c(1, 3),
               function(v) as.integer(any(v[!is.na(v)] == 1)))
  fit <- tobs(~ x, data = sim$data, family = jsdm(), y = yj,
              species = paste0("sp", 1:5), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  f <- fitted(fit)
  expect_equal(dim(f$mu), c(40L, 5L))
  expect_true(all(f$mu >= 0 & f$mu <= 1))
  expect_equal(dim(predict(fit)), c(40L, 5L))
  expect_equal(dim(predict(fit, newdata = data.frame(x = c(-1, 0, 1)))),
               c(3L, 5L))
  rr <- residuals(fit)
  expect_equal(dim(rr$mu), c(40L, 5L))
  expect_true(all(is.finite(rr$mu)))
})
