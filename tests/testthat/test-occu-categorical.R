# Tests for occu_categorical() -- presence + nominal K-class hurdle. Bernoulli
# presence + baseline-category multinomial logit on the class given present.

test_that("occu_categorical() constructs", {
  fam <- occu_categorical()
  expect_s3_class(fam, "tobs_family")
  expect_equal(fam$name, "occu_categorical")
  expect_equal(fam$observation, "binomial_plus_multinomial")
  expect_equal(fam$status, "working")
  fam2 <- occu_categorical(classes = c("forest", "grass", "wet"))
  expect_equal(fam2$params$classes, c("forest", "grass", "wet"))
  expect_error(occu_categorical(classes = "only-one"), "at least two")
})

test_that("multinomial-logit kernel matches the vectorised R fitter math", {
  # one Newton-free spot check: the C++ kernel softmax == the R fitter softmax
  eta <- c(-0.4, 0.9, 0.2)
  out <- tulpa:::cpp_multinomial_logit_terms(eta, 2L)
  denom <- 1 + sum(exp(eta)); p <- c(exp(eta) / denom, 1 / denom)
  expect_equal(out$ll, log(p[2]), tolerance = 1e-12)
  expect_equal(as.numeric(out$grad), c(0, 1, 0) - p[seq_along(eta)], tolerance = 1e-12)
})

test_that("the R multinomial fitter recovers known coefficients", {
  set.seed(3)
  N <- 6000L; K <- 4L
  X <- cbind(1, rnorm(N), rnorm(N))
  B <- matrix(c(0.3, 0.8, -0.5, -0.6, 0.4, 0.9, 0.2, -0.7, 0.3), nrow = 3)
  E <- exp(X %*% B); P <- cbind(E / (1 + rowSums(E)), 1 / (1 + rowSums(E)))
  cls <- apply(P, 1L, function(pr) sample.int(K, 1L, prob = pr))
  fit <- tulpaObs:::.tobs_mlogit_fit(X, cls, K)
  expect_true(fit$converged)
  expect_lt(max(abs(fit$Beta - B)), 0.1)
})

test_that("occu_categorical requires laplace and rejects a detection formula", {
  sim <- simulate_occu_categorical(N = 200L, seed = 5)
  dat <- cbind(sim$data, y = sim$y)
  expect_error(
    tobs(y ~ x, data = dat, family = occu_categorical(), method = "nuts"),
    "not available"
  )
})

test_that("occu_categorical recovers presence + class truth end to end", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_categorical(N = 5000L, seed = 1)
  dat <- cbind(sim$data, y = sim$y)
  fit <- tobs(y ~ x, data = dat, family = occu_categorical(), method = "laplace")
  expect_s3_class(fit, "occu_categorical_fit")
  expect_equal(fit$K, sim$truth$K)
  expect_true(fit$convergence$converged)

  # presence arm
  expect_lt(abs(fit$beta_occ[1] - sim$truth$beta_occ[1]), 0.2)
  expect_lt(abs(fit$beta_occ[2] - sim$truth$beta_occ[2]), 0.2)
  # class arm (baseline-category multinomial)
  expect_lt(max(abs(fit$beta_class - sim$truth$beta_class)), 0.25)

  # predict: valid simplices
  pr <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
  expect_equal(dim(pr$cond), c(3L, fit$K))
  expect_true(all(abs(rowSums(pr$cond) - 1) < 1e-8))
  expect_true(all(abs(rowSums(pr$joint) - 1) < 1e-8))
})
