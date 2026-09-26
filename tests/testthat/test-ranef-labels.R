test_that("ranef() reports the grouping factor's own labels (state arm)", {
  skip_if_fast()
  skip_on_cran()
  set.seed(5)
  ng <- 12L; per <- 15L; J <- 5L; N <- ng * per
  lab <- sprintf("grp%02d", sample.int(ng))
  gi <- rep(seq_len(ng), each = per)
  b <- rnorm(ng, sd = 0.9)
  x <- rnorm(N)
  psi <- plogis(0.2 + 0.5 * x + b[gi])
  z <- rbinom(N, 1, psi)
  y <- matrix(rbinom(N * J, 1, z * 0.6), N, J)
  d <- data.frame(x = x, g = factor(lab[gi]))

  fit <- tobs(~ x + (1 | g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  rf <- ranef(fit)

  expect_setequal(as.character(rf$level), lab)
  expect_equal(nrow(rf), ng)
  est <- rf$estimate[match(lab, rf$level)]
  expect_gt(stats::cor(est, b), 0.3)
})

test_that("ranef() reports the grouping factor's own labels (detection arm)", {
  skip_if_fast()
  skip_on_cran()
  set.seed(1)
  ng <- 10L; N <- 200L; J <- 5L
  lab <- sprintf("obs%02d", sample.int(ng))
  gi <- rep_len(seq_len(ng), N)
  b <- rnorm(ng, sd = 0.8)
  x <- rnorm(N)
  z <- rbinom(N, 1, plogis(0.3 + 0.5 * x))
  y <- matrix(rbinom(N * J, 1, z * plogis(0.2 + b[gi])), N, J)
  d <- data.frame(x = x, observer = factor(lab[gi]))

  fit <- tobs(~ x, detection = ~ (1 | observer), data = d, y = y,
              family = occu(), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  rf <- ranef(fit)

  expect_setequal(as.character(rf$level), lab)
  est <- rf$estimate[match(lab, rf$level)]
  expect_gt(stats::cor(est, b), 0.3)
})

test_that("random-effect design carries the factor levels in code order", {
  d <- data.frame(g = factor(c("b", "a", "c", "a", "b", "c")))
  codes <- tulpaObs:::.tobs_index_codes(d$g, "re", "group")
  expect_identical(codes$levels, c("a", "b", "c"))
  expect_identical(codes$levels[codes$idx], as.character(d$g))
  codes_num <- tulpaObs:::.tobs_index_codes(c(10, 2, 10, 2), "re", "group")
  expect_identical(codes_num$levels, c("2", "10"))
})
