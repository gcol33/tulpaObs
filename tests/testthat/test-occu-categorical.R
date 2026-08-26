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

test_that("occu_categorical() recovers presence + class truth (multi-seed) and CIs cover (#276)", {
  # Mirrors the multi-seed recovery + CI-coverage pattern every other roster
  # family carries (test-occu-ttd.R, test-double-observer.R, ...); a single
  # seed cannot separate a working estimator from one that lands close by luck.
  skip_on_cran()
  skip_if_fast()
  n_seed <- 15L
  beta_occ_truth <- c(0.2, 0.8)
  bo1 <- bo2 <- bc <- rep(NA_real_, n_seed)
  hit_o <- hit_c <- tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu_categorical(N = 1500L, beta_occ = beta_occ_truth,
                                     seed = 900L + s)
    dat <- cbind(sim$data, y = sim$y)
    fit <- tryCatch(
      tobs(y ~ x, data = dat, family = occu_categorical(), method = "laplace"),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    tot <- tot + 1L
    bo1[tot] <- fit$beta_occ[1]; bo2[tot] <- fit$beta_occ[2]
    bc[tot]  <- fit$beta_class[2, 1]   # slope on x, first non-baseline class
    if (abs(fit$beta_occ[2] - beta_occ_truth[2]) <= 1.96 * fit$se_occ[2])
      hit_o <- hit_o + 1L
    if (abs(fit$beta_class[2, 1] - sim$truth$beta_class[2, 1]) <=
          1.96 * fit$se_class[2, 1])
      hit_c <- hit_c + 1L
  }
  expect_true(tot >= floor(0.8 * n_seed))
  expect_lt(abs(mean(bo1[seq_len(tot)]) - beta_occ_truth[1]), 0.15)
  expect_lt(abs(mean(bo2[seq_len(tot)]) - beta_occ_truth[2]), 0.15)
  expect_lt(abs(mean(bc[seq_len(tot)])  - sim$truth$beta_class[2, 1]), 0.2)
  expect_gte(hit_o / tot, 0.8)
  expect_gte(hit_c / tot, 0.8)
})

test_that("occu_categorical_fit S3 surface: coef/vcov/confint/logLik/glance/tidy/summary/print/nobs (#276)", {
  sim <- simulate_occu_categorical(N = 300L, seed = 9)
  dat <- cbind(sim$data, y = sim$y)
  fit <- tobs(y ~ x, data = dat, family = occu_categorical(), method = "laplace")

  expect_identical(nobs(fit), 300L)

  # coef.occu_categorical_fit() overrides the generic multiarm coef() and
  # returns the raw (p x (K-1)) matrix, not the flattened arm-blocks vector.
  cf <- coef(fit)
  expect_named(cf, c("presence", "class"))
  expect_equal(unname(cf$presence), unname(fit$beta_occ))
  expect_equal(cf$class, fit$beta_class)

  # vcov()/confint()/tidy() have no family override, so they DO go through
  # `.tobs_arm_blocks_categorical()`, which reshapes the (p x (K-1)) beta_class
  # / se_class matrices column-major into a flat vector; check the flatten
  # actually lines each name up with the matrix cell it came from (#276 item 4).
  cls_labels  <- colnames(fit$beta_class); coef_labels <- rownames(fit$beta_class)
  V <- vcov(fit)
  se_flat <- sqrt(diag(V))
  for (cl in cls_labels) for (co in coef_labels) {
    key <- paste0("class:", cl, ":", co)
    expect_equal(unname(V[key, key]), unname(fit$se_class[co, cl])^2, tolerance = 1e-8)
    expect_equal(unname(se_flat[[key]]), unname(fit$se_class[co, cl]), tolerance = 1e-8)
  }

  ci <- confint(fit)
  expect_true(all(ci[, 1] <= ci[, 2]))

  ll <- logLik(fit)
  expect_true(is.finite(as.numeric(ll)))
  expect_identical(attr(ll, "nobs"), 300L)

  g <- glance(fit)
  expect_true(all(c("n", "logLik", "df", "converged") %in% names(g)))

  td <- tidy(fit)
  expect_true(all(c("arm", "term", "estimate", "std.error",
                    "conf.low", "conf.high") %in% names(td)))
  expect_setequal(unique(td$arm), c("presence", "class"))

  expect_output(print(fit), "occu_categorical")
  expect_output(summary(fit), "presence arm")
})
