# One-inflated Beta cover family.
#
# cover(response = "beta_oi") models plots recorded at exactly full cover
# (y = 1) as a genuine point mass rather than clamping them to 1 - 1e-6. With a
# constant inflation probability the likelihood factorizes: pi is the share of
# positive plots at the ceiling (a binomial proportion), and the interior Beta
# is fit on the (0, 1) plots. These tests prove (a) the ceiling clamp the plain
# Beta suffers biases its precision while beta_oi does not, (b) pi and the
# interior Beta coefficients recover with calibrated intervals, (c) the spatial
# nested-Laplace path carries it.

# Truth: presence ~ Bernoulli; among positives a share pi sit at cover = 1, the
# rest follow a Beta(mu, phi) on (0, 1).
sim_beta_oi <- function(N = 600, pi_one = 0.25, phi = 25,
                        beta_occ = c(0.3, 0.6), beta_pos = c(0.2, -0.5), seed = 1) {
  set.seed(seed)
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(beta_occ[1] + beta_occ[2] * x))
  mu <- plogis(beta_pos[1] + beta_pos[2] * x)
  y <- numeric(N); is_pos <- occur == 1L; np <- sum(is_pos)
  at_ceiling <- rbinom(np, 1, pi_one) == 1L
  yp <- numeric(np); yp[at_ceiling] <- 1
  yp[!at_ceiling] <- rbeta(sum(!at_ceiling), mu[is_pos][!at_ceiling] * phi,
                           (1 - mu[is_pos][!at_ceiling]) * phi)
  y[is_pos] <- yp
  list(data = data.frame(x = x), y = y,
       truth = list(pi_one = pi_one, phi = phi, beta_pos = beta_pos))
}

ci_cover <- function(est, se, truth, z = 1.96) {
  truth >= est - z * se && truth <= est + z * se
}

test_that("beta_oi reports pi_one and the structure", {
  sim <- sim_beta_oi(N = 600, pi_one = 0.3, seed = 11)
  fit <- tobs(~ x, data = sim$data, family = cover("beta_oi"), y = sim$y,
              method = "laplace")
  expect_s3_class(fit, "cover_fit")
  expect_identical(fit$positive, "beta_oi")
  expect_true(is.finite(fit$pi_one) && fit$pi_one > 0 && fit$pi_one < 1)
  expect_true(is.finite(fit$pi_one_sd))
  expect_gt(fit$n_ceiling, 0)
  # The conditional cover mixes the ceiling mass with the interior Beta mean.
  cond <- predict(fit, data.frame(x = 0), type = "conditional")
  mu0  <- plogis(fit$beta_pos[[1]])
  expect_equal(unname(cond), unname(fit$pi_one + (1 - fit$pi_one) * mu0),
               tolerance = 1e-8)
})

test_that("the ceiling clamp biases plain beta's precision; beta_oi does not", {
  # ~25% of positives at cover = 1. The plain Beta clamps them to 1 - 1e-6,
  # which forces a tiny precision; beta_oi sets them aside and recovers phi.
  sim <- sim_beta_oi(N = 1200, pi_one = 0.25, phi = 25, seed = 3)
  fit_oi <- tobs(~ x, data = sim$data, family = cover("beta_oi"), y = sim$y,
                 method = "laplace")
  fit_b  <- tobs(~ x, data = sim$data, family = cover("beta"), y = sim$y,
                 method = "laplace")
  # beta_oi recovers phi within 25%; the clamped plain Beta is far too low.
  expect_lt(abs(fit_oi$phi_pos - 25) / 25, 0.25)
  expect_lt(fit_b$phi_pos, 5)
  expect_gt(fit_oi$phi_pos, 3 * fit_b$phi_pos)
})

test_that("beta_oi recovers pi and the interior coefficients with calibrated CIs", {
  skip_if_fast(); skip_on_cran()
  truth_pi <- 0.25; truth_slope <- -0.5; truth_phi <- 25
  n_seeds <- 20L
  hit_pi <- logical(n_seeds); hit_slope <- logical(n_seeds)
  err_pi <- numeric(n_seeds); err_phi <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    sim <- sim_beta_oi(N = 700, pi_one = truth_pi, phi = truth_phi,
                       beta_pos = c(0.2, truth_slope), seed = 100 + s)
    fit <- tobs(~ x, data = sim$data, family = cover("beta_oi"), y = sim$y,
                method = "laplace")
    hit_pi[s]    <- ci_cover(fit$pi_one, fit$pi_one_sd, truth_pi)
    hit_slope[s] <- ci_cover(fit$beta_pos[["x"]], fit$se_pos[["x"]], truth_slope)
    err_pi[s]    <- abs(fit$pi_one - truth_pi)
    err_phi[s]   <- abs(fit$phi_pos - truth_phi) / truth_phi
  }
  expect_lt(median(err_pi), 0.04)              # pi recovery
  expect_lt(median(err_phi), 0.2)              # interior precision recovery
  expect_gte(mean(c(hit_pi, hit_slope)), 0.85) # pooled CI coverage
})

test_that("beta_oi carries through the spatial nested-Laplace path", {
  skip_if_fast(); skip_on_cran()
  set.seed(9)
  n_s <- 25L
  nbr <- lapply(seq_len(n_s), function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L
  N <- 700L
  w <- 0.6 * (sqrt(0.7) * rnorm(n_s) + sqrt(0.3) * rnorm(n_s))
  region <- sample.int(n_s, N, replace = TRUE); x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(0.2 + 0.6 * x + w[region]))
  mu <- plogis(0.3 - 0.5 * x + w[region])
  y <- numeric(N); is_pos <- occur == 1L; np <- sum(is_pos)
  ceil <- rbinom(np, 1, 0.25) == 1L; yp <- numeric(np); yp[ceil] <- 1
  yp[!ceil] <- rbeta(sum(!ceil), mu[is_pos][!ceil] * 25, (1 - mu[is_pos][!ceil]) * 25)
  y[is_pos] <- yp
  dat <- data.frame(x = x, region = factor(region))

  # The simulator puts the SAME field w in both means -- the occurrence logit
  # and the interior cover mean -- so the fit has to be told the arms share it.
  # A field without share() loads on occurrence alone, which leaves w in the
  # cover arm's residual and reads back as over-dispersion: phi comes out near
  # 10 against a truth of 25, the same value a fit with no field at all gives.
  # Truth is one field at unit loading, so the alpha axis straddles 1.
  fit <- tobs(formula = ~ x + bym2(graph = adj, group_var = "region") +
                share(spatial(), alpha = grid(c(0.5, 1.0, 1.5))),
              data = dat, family = cover("beta_oi"), y = y,
              method = "nested_laplace",
              control = list(sigma.grid = c(0.4, 0.8), rho.grid = c(0.5, 0.9)))
  expect_s3_class(fit, "cover_fit")
  expect_identical(fit$positive, "beta_oi")
  expect_true(fit$converged)
  expect_lt(abs(fit$pi_one - 0.25), 0.08)
  expect_lt(abs(fit$phi_pos - 25) / 25, 0.3)
  # Spatial predict projects the one-inflated conditional cover.
  nd <- data.frame(x = 0, region = factor(1, levels = levels(dat$region)))
  pr <- predict(fit, nd, type = "cover_cond")
  expect_s3_class(pr, "tobs_prediction")
})
