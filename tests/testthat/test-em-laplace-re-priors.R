test_that("RE Laplace path warns on a user prior and has no dead damping knob (#427)", {
  expect_false("damping" %in% names(formals(tulpaObs:::.tobs_em_laplace_re)))
  skip_on_cran()
  set.seed(4)
  N <- 120L; J <- 4L; G <- 12L
  g <- factor(rep(seq_len(G), length.out = N))
  x <- stats::rnorm(N)
  z <- stats::rbinom(N, 1, stats::plogis(0.2 + 0.5 * x + stats::rnorm(G, 0, 0.6)[g]))
  y <- matrix(stats::rbinom(N * J, 1, 0.5 * z), N, J)
  d <- data.frame(x = x, g = g)
  ctl <- list(verbose = FALSE, progress = FALSE)
  expect_warning(
    tobs(~ x + (1 | g), data = d, family = occu(), detection = ~ 1, y = y,
         method = "laplace", priors = occu_priors(), control = ctl),
    "not applied on the random-effect Laplace path")
  expect_no_warning(
    tobs(~ x + (1 | g), data = d, family = occu(), detection = ~ 1, y = y,
         method = "laplace", control = ctl))
})
