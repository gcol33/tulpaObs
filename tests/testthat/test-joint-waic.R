# WAIC / pointwise log-likelihood for the cover-family joint fits:
# occu_cover() (non-spatial Laplace + spatial nested Laplace) and the
# nested-joint cover() shared-field fit.

.jw_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  adj
}

.jw_occu_cover_long <- function(sim, N, J) {
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  tobs_data(long, y = "y", site = "site_id", visit = "visit",
            det.covs = c("det_cov1", "pos_cov1"))
}

test_that("occu_cover() non-spatial: WAIC + pointwise log-lik (#26)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(11)
  N <- 120L; J <- 5L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", seed = 21L)
  od <- .jw_occu_cover_long(sim, N, J)
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- tobs(formula = ~ occ_cov1, data = cell_dat,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = od$y, y_pos = y_pos, visits = od$det.covs,
              method = "laplace", control = list(verbose = FALSE))

  ll <- tulpaObs:::.tobs_pointwise_loglik(fit)
  expect_equal(ncol(ll), N)                       # per-site pointwise
  expect_true(all(is.finite(ll)))

  w <- waic(fit)
  expect_true(is.finite(w$waic) && is.finite(w$elpd))
  expect_gte(w$p_waic, 0)
  # lppd is a sane per-observation magnitude for a hurdle (not absurd).
  expect_lt(abs(w$lppd / N), 5)
})

test_that("occu_cover() beta arm: pointwise log-lik is finite (#26)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(12)
  N <- 120L; J <- 5L
  sim <- simulate_occu_cover(N = N, J = J, positive = "beta", phi = 30, seed = 22L)
  od <- .jw_occu_cover_long(sim, N, J)
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  y_pos <- pmin(pmax(y_pos, 1e-6), 1 - 1e-6)

  fit <- tobs(formula = ~ occ_cov1, data = cell_dat, family = occu_cover("beta"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = od$y, y_pos = y_pos, visits = od$det.covs,
              method = "laplace", control = list(verbose = FALSE))
  w <- waic(fit)
  expect_true(is.finite(w$waic))
  expect_equal(ncol(tulpaObs:::.tobs_pointwise_loglik(fit)), N)
})

test_that("occu_cover() spatial joint: WAIC + pointwise log-lik (#26)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(13)
  N <- 30L; J <- 5L
  adj <- .jw_chain_adj(N)
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 0.8, alpha = 0.6, seed = 23L)
  od <- .jw_occu_cover_long(sim, N, J)
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 200L, engine = "joint",
                   sigma.grid = c(0.5, 1.0), alpha.grid = c(0, 0.5))
  ))

  ll <- tulpaObs:::.tobs_pointwise_loglik(fit)
  expect_equal(ncol(ll), N)
  expect_true(all(is.finite(ll)))
  w <- waic(fit)
  expect_true(is.finite(w$waic))
  expect_gte(w$p_waic, 0)
})

test_that("cover() nested-joint: WAIC + pointwise log-lik (#26)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(14)
  N <- 200L; n_s <- 25L
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  adj <- .jw_chain_adj(n_s)
  phi <- rnorm(n_s); theta <- rnorm(n_s)
  w_s <- 0.6 * (sqrt(0.7) * phi + sqrt(0.3) * theta)
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.3 + 0.7 * x + w_s[spatial_idx]))
  y <- ifelse(occur == 1L,
              pmin(exp(rnorm(N, 0.4 - 0.5 * x + w_s[spatial_idx], 0.4)),
                   1 - 1e-6), 0)
  dat <- data.frame(x = x, region = factor(spatial_idx))

  fit <- tobs(formula = ~ x + bym2(graph = adj, group_var = "region"),
              data = dat, family = cover("lognormal"), y = y,
              method = "nested_laplace",
              control = list(sigma.grid = c(0.4, 0.8), rho.grid = c(0.5, 0.9)))

  expect_false(is.null(fit$spi_full))
  ll <- tulpaObs:::.tobs_pointwise_loglik(fit)
  expect_equal(ncol(ll), N)                       # per-observation pointwise
  expect_true(all(is.finite(ll)))
  w <- waic(fit)
  expect_true(is.finite(w$waic) && is.finite(w$elpd))
  expect_lt(abs(w$lppd / N), 5)

  # Without the stored spatial-unit index the joint log-lik errors clearly.
  fit2 <- fit; fit2$spi_full <- NULL
  expect_error(waic(fit2), "spatial-unit index")
})

test_that("cover() separate-Laplace WAIC is unaffected (#26)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(15)
  N <- 150L
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.2 + 0.6 * x))
  y <- ifelse(occur == 1L, pmin(exp(rnorm(N, 0.3 - 0.4 * x, 0.4)), 1 - 1e-6), 0)
  fit <- tobs(formula = ~ x, data = data.frame(x = x),
              family = cover("lognormal"), y = y, method = "laplace")
  w <- waic(fit)
  expect_true(is.finite(w$waic))
  expect_equal(ncol(tulpaObs:::.tobs_pointwise_loglik(fit)), N)
})
