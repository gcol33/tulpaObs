# DIC / CPO / LPML via the engine criteria layer, and family-specific PIT + PPC
# for the cover hurdle and occu_cover (gcol33/tulpa#47, gcol33/tulpaObs#27).

.cc_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  adj
}

.cc_occu_cover_long <- function(sim, N, J) {
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  tobs_data(long, y = "y", site = "site_id", visit = "visit",
            det.covs = c("det_cov1", "pos_cov1"))
}

test_that("cover() separate-Laplace: WAIC/DIC/CPO + PIT + PPC", {
  skip_on_cran()
  skip_if_fast()
  set.seed(101)
  N <- 200L
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.2 + 0.6 * x))
  y <- ifelse(occur == 1L, pmin(exp(rnorm(N, 0.3 - 0.4 * x, 0.4)), 1 - 1e-6), 0)
  fit <- tobs(formula = ~ x, data = data.frame(x = x),
              family = cover("lognormal"), y = y, method = "laplace")

  w <- waic(fit)
  expect_true(is.finite(w$waic) && is.finite(w$elpd))

  d <- dic(fit, n.draws = 300L)
  expect_true(is.finite(d$dic) && is.finite(d$p_dic))
  expect_true(d$p_dic >= -1)               # effective parameters ~ small positive

  cp <- cpo(fit, n.draws = 300L)
  expect_true(is.finite(cp$lpml))
  expect_equal(cp$lpml, cp$elpd_loo, tolerance = 1e-8)
  expect_length(cp$pointwise$cpo, N)
  expect_true(all(cp$pointwise$cpo > 0))

  pit <- pit_residuals(fit, n.samples = 300L)
  expect_length(pit, N)
  expect_true(all(pit >= 0 & pit <= 1))

  ppc <- ppc(fit, n.samples = 200L)
  expect_true(ppc$bayesian.p >= 0 && ppc$bayesian.p <= 1)
  # Correct model: not in the extreme tails.
  expect_gt(ppc$bayesian.p, 0.001)
  expect_lt(ppc$bayesian.p, 0.999)
})

test_that("cover() beta arm: criteria + PIT + PPC run", {
  skip_on_cran()
  skip_if_fast()
  set.seed(102)
  N <- 200L
  x <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.1 + 0.5 * x))
  mu <- plogis(0.2 - 0.3 * x); phi <- 25
  y <- ifelse(occur == 1L,
              pmin(pmax(rbeta(N, mu * phi, (1 - mu) * phi), 1e-6), 1 - 1e-6), 0)
  fit <- tobs(formula = ~ x, data = data.frame(x = x),
              family = cover("beta"), y = y, method = "laplace")
  expect_true(is.finite(dic(fit, n.draws = 200L)$dic))
  pit <- pit_residuals(fit, n.samples = 200L)
  expect_true(all(pit >= 0 & pit <= 1))
  expect_true(ppc(fit, n.samples = 150L)$bayesian.p <= 1)
})

test_that("cover() nested-joint: PIT + PPC project the shared field", {
  skip_on_cran()
  skip_if_fast()
  set.seed(103)
  N <- 200L; n_s <- 25L
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  adj <- .cc_chain_adj(n_s)
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
  pit <- pit_residuals(fit, n.samples = 200L)
  expect_length(pit, N)
  expect_true(all(is.finite(pit) & pit >= 0 & pit <= 1))
  expect_true(is.finite(dic(fit, n.draws = 200L)$dic))
  expect_true(ppc(fit, n.samples = 150L)$bayesian.p <= 1)
})

test_that("occu_cover(): DIC/CPO + PIT + PPC", {
  skip_on_cran()
  skip_if_fast()
  set.seed(104)
  N <- 120L; J <- 5L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", seed = 41L)
  od <- .cc_occu_cover_long(sim, N, J)
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  fit <- tobs(formula = ~ occ_cov1, data = cell_dat,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = od$y, y_pos = y_pos, visits = od$det.covs,
              method = "laplace", control = list(verbose = FALSE))

  d <- dic(fit, n.draws = 300L)
  expect_true(is.finite(d$dic) && is.finite(d$p_dic))
  cp <- cpo(fit, n.draws = 300L)
  expect_equal(cp$lpml, cp$elpd_loo, tolerance = 1e-8)
  expect_length(cp$pointwise$cpo, N)

  pit <- pit_residuals(fit, n.samples = 250L)
  expect_length(pit, N)
  expect_true(all(pit >= 0 & pit <= 1))

  ppc <- ppc(fit, n.samples = 200L)
  expect_true(ppc$bayesian.p >= 0 && ppc$bayesian.p <= 1)
})
