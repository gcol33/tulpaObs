# predict(newdata = ) on a random-effect fit (#372): a row whose grouping
# column names a fitted level takes that group's effect, an unseen level the
# population mean (0), and a newdata without the grouping column -- or a bare
# X.0 -- is population-level.

re_occu_fixture <- function() {
  set.seed(3); G <- 10; per <- 20; N <- G * per
  g <- factor(sprintf("grp%02d", rep(seq_len(G), each = per)))
  x <- rnorm(N); b <- rnorm(G, 0, 1)
  z <- rbinom(N, 1, plogis(0.2 + 0.5 * x + b[as.integer(g)]))
  y <- matrix(rbinom(N * 4, 1, z * 0.6), N, 4)
  list(d = data.frame(x = x, g = g), y = y)
}

test_that("occu() (1|g): a seen level adds its effect, an unseen one adds 0", {
  skip_if_fast()
  fx <- re_occu_fixture()
  fit <- suppressWarnings(tobs(~ x + (1 | g), data = fx$d, y = fx$y,
                               detection = ~ 1, family = occu(),
                               method = "laplace",
                               control = list(progress = FALSE, verbose = FALSE)))
  re <- ranef(fit)
  lev <- re$level[which.max(re$estimate)]
  nd  <- data.frame(x = c(0, 0), g = c(lev, "never_seen"))
  pg  <- predict(fit, newdata = nd)
  pop <- predict(fit, newdata = data.frame(x = c(0, 0)))
  pX  <- predict(fit, X.0 = cbind(1, c(0, 0)))
  # Population-level without the grouping column, and the same through X.0.
  expect_equal(pop[[1]], pX[[1]])
  # The unseen level is the population prediction.
  expect_equal(pg[[1]][2], pop[[1]][2])
  # The seen level moves by its own effect, on the logit scale.
  b <- re$estimate[re$level == lev]
  expect_gt(b, 0)
  expect_equal(qlogis(pg[[1]][1]) - qlogis(pop[[1]][1]), b, tolerance = 0.05)
})

test_that("abun() (1|g) on the abundance arm takes the group effect", {
  skip_if_fast()
  sa <- simulate_abun(N = 60L, J = 4L, seed = 13L)
  set.seed(4)
  g <- factor(rep(sprintf("s%d", 1:6), 10))
  bg <- rnorm(6, 0, 0.6)
  sa$data$g <- g
  fit <- suppressWarnings(tobs(~ abund_cov1 + (1 | g), data = sa$data,
                               family = abun(), detection = ~ det_cov1,
                               y = sa$y, method = "laplace",
                               control = list(progress = FALSE, verbose = FALSE)))
  re <- ranef(fit)
  lev <- re$level[1]
  pg  <- predict(fit, newdata = data.frame(abund_cov1 = 0, g = lev))
  pop <- predict(fit, newdata = data.frame(abund_cov1 = 0))
  b <- re$estimate[re$level == lev]
  expect_equal(log(pg$mean) - log(pop$mean), b, tolerance = 0.05)
})
