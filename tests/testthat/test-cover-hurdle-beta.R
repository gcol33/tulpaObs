## Tests for cover(response = "beta") via tulpa_laplace_beta.

simulate_beta_cover <- function(N = 500, beta_occ = c(-0.4, 0.8),
                                beta_pos = c(0.5, -1.2), phi = 30,
                                seed = 1) {
  set.seed(seed)
  x <- runif(N, -2, 2)
  eta_occ <- beta_occ[1] + beta_occ[2] * x
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x
  mu      <- plogis(eta_pos)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos]  <- rbeta(sum(is_pos), mu[is_pos] * phi, (1 - mu[is_pos]) * phi)
  y[!is_pos] <- 0
  y <- pmin(pmax(y, 0), 1 - 1e-6)
  list(
    data  = data.frame(x = x),
    y     = y,
    truth = list(beta_occ = beta_occ, beta_pos = beta_pos, phi = phi,
                 occur = occur)
  )
}

test_that("cover(response='beta') recovers betas and phi on simulated data", {
  sim <- simulate_beta_cover(N = 800, seed = 11)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )
  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "beta")
  expect_true(fit$converged)

  expect_lt(abs(fit$beta_occ[1] - sim$truth$beta_occ[1]), 0.3)
  expect_lt(abs(fit$beta_occ[2] - sim$truth$beta_occ[2]), 0.3)
  expect_lt(abs(fit$beta_pos[1] - sim$truth$beta_pos[1]), 0.2)
  expect_lt(abs(fit$beta_pos[2] - sim$truth$beta_pos[2]), 0.2)
  expect_lt(abs(fit$phi_pos / sim$truth$phi - 1), 0.25)
})

test_that("predict() respects the beta back-transform", {
  sim <- simulate_beta_cover(N = 400, seed = 5)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )
  newdata <- data.frame(x = seq(-2, 2, length.out = 21))

  X <- model.matrix(~ x, newdata)
  expected_mu <- plogis(as.numeric(X %*% fit$beta_pos))
  expected_p  <- plogis(as.numeric(X %*% fit$beta_occ))

  expect_equal(predict(fit, newdata, type = "occupancy"),  expected_p)
  expect_equal(predict(fit, newdata, type = "conditional"), expected_mu)
  expect_equal(predict(fit, newdata, type = "expected"),
               expected_p * expected_mu)
})

test_that("cover(response='beta') matches a separate-fit pipeline mirroring the reference INLA hurdle", {
  skip_if_not_installed("betareg")

  sim <- simulate_beta_cover(N = 800, seed = 19)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )

  # Replicate Michael's MOT_abund_data.Rmd separate-hurdle pattern (binomial
  # on occurrence, beta on positive subset) using glm + betareg as a check.
  dat <- sim$data
  dat$occur <- as.integer(sim$y > 0)
  m_occ <- glm(occur ~ x, family = binomial(), data = dat)
  pos_idx <- which(dat$occur == 1)
  m_pos <- betareg::betareg(y ~ x, data = data.frame(y = sim$y[pos_idx],
                                                     x = dat$x[pos_idx]))

  expect_lt(abs(fit$beta_occ[1] - coef(m_occ)[1]), 0.05)
  expect_lt(abs(fit$beta_occ[2] - coef(m_occ)[2]), 0.05)
  expect_lt(abs(fit$beta_pos[1] - coef(m_pos)[1]), 0.05)
  expect_lt(abs(fit$beta_pos[2] - coef(m_pos)[2]), 0.05)
  expect_lt(abs(fit$phi_pos / coef(m_pos)[3] - 1), 0.05)
})
