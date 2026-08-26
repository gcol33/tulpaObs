# Cover-hurdle spatially-varying trend declared as a weighted areal formula
# term. A spatially varying trend is model structure, so it lives in the
# formula as a second weighted icar() term, not in control$trend (removed,
# breaking).

# Build a grid ICAR adjacency.
# Simulate a two-field cover hurdle: a shared intercept ICAR field plus a
# time-weighted spatially-varying trend field, both copied onto the cover arm.
simulate_cover_trend <- function(g = 5L, N = 2500L,
                                  sigma1 = 0.8, alpha1 = 1.0,
                                  sigma2 = 0.7, alpha2 = 1.4,
                                  beta_occ = c(0.2, 0.3),
                                  beta_pos = c(-1.6, 0.2),
                                  sd_pos = 0.4, seed = 7) {
  set.seed(seed)
  n_cells <- g * g
  adj <- rook_adj(g)
  smooth_field <- function() {
    z <- rnorm(n_cells)
    for (it in 1:40) {
      zn <- z
      for (i in seq_len(n_cells)) {
        nb <- which(adj[i, ] == 1L); zn[i] <- 0.4 * z[i] + 0.6 * mean(z[nb])
      }
      z <- zn
    }
    z <- z - mean(z); z / stats::sd(z)
  }
  z1 <- smooth_field(); z2 <- smooth_field()

  cell <- sample.int(n_cells, N, replace = TRUE)
  time <- as.numeric(scale(rnorm(N)))
  eta_occ <- beta_occ[1] + beta_occ[2] * time + sigma1 * z1[cell] +
             time * (sigma2 * z2[cell])
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * time +
             (alpha1 * sigma1) * z1[cell] +
             time * ((alpha2 * sigma2) * z2[cell])
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, sd_pos)), 1 - 1e-6), 0)

  list(adj = adj,
       data = data.frame(cell = cell, time = time, cover = cover),
       y = cover,
       truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                    sigma2 = sigma2, alpha2 = alpha2))
}

trend_control <- list(verbose = FALSE, n.threads = 1L, adaptive.grid = TRUE,
                      sigma.grid = exp(seq(log(0.2), log(2), length.out = 5)),
                      alpha.grid = c(0, exp(seq(log(0.2), log(3), length.out = 4))))

# ---- Recovery: the formula-driven trend path fits and recovers the trend ----

test_that("cover() recovers a trend declared as a weighted areal formula term", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_cover_trend(seed = 7)

  fit <- tobs(
    formula = ~ time + icar(graph = sim$adj, group_var = "cell") +
                icar(graph = sim$adj, weight = time, group_var = "cell"),
    data = sim$data, family = cover(response = "lognormal"), y = sim$y,
    method = "nested_laplace", control = trend_control)

  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)
  expect_identical(fit$n_fields, 2L)
  expect_identical(fit$trend_weight, "time")

  # Fixed effects recovered on the right side of zero.
  expect_gt(fit$beta_occ[2], 0)   # truth 0.3
  expect_lt(fit$beta_pos[1], 0)   # truth -1.6

  # The coupled trend block is identified: a finite positive amplitude
  # (b2.sigma * b2.alpha), with alpha_trend on the right side of zero.
  expect_true(is.finite(fit$sigma_trend) && fit$sigma_trend > 0)
  expect_true(is.finite(fit$alpha_trend) && fit$alpha_trend > 0.3)
})

# ---- Equivalence: bare icar(weight=) == umbrella spatial(model, weight=) ----

test_that("spatial(model='icar', weight=) resolves identically to icar(weight=)", {
  skip_if_fast()
  skip_on_cran()
  sim <- simulate_cover_trend(seed = 7)

  fit_bare <- tobs(
    formula = ~ time + icar(graph = sim$adj, group_var = "cell") +
                icar(graph = sim$adj, weight = time, group_var = "cell"),
    data = sim$data, family = cover(response = "lognormal"), y = sim$y,
    method = "nested_laplace", control = trend_control)

  fit_umb <- tobs(
    formula = ~ time +
                spatial(graph = sim$adj, model = "icar", group_var = "cell") +
                spatial(graph = sim$adj, model = "icar", weight = time,
                        group_var = "cell"),
    data = sim$data, family = cover(response = "lognormal"), y = sim$y,
    method = "nested_laplace", control = trend_control)

  expect_equal(fit_bare$beta_occ,    fit_umb$beta_occ)
  expect_equal(fit_bare$beta_pos,    fit_umb$beta_pos)
  expect_equal(fit_bare$sigma_trend, fit_umb$sigma_trend)
  expect_equal(fit_bare$alpha_trend, fit_umb$alpha_trend)
})

# ---- control$trend is removed (breaking) ------------------------------------

test_that("control$trend errors with a migration pointer", {
  adj <- rook_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)

  expect_error(
    tobs(formula = ~ time + icar(graph = adj, group_var = "cell"),
         data = df, family = cover(response = "lognormal"), y = y,
         method = "nested_laplace",
         control = list(trend = list(weight = "time"))),
    "control\\$trend is no longer supported")
})

# ---- New parser guards on the weighted areal term ---------------------------

test_that("a weighted areal term without an intercept field errors", {
  adj <- rook_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)

  expect_error(
    tobs(formula = ~ time + icar(graph = adj, weight = time, group_var = "cell"),
         data = df, family = cover(response = "lognormal"), y = y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "unweighted intercept field")
})

test_that("a weighted trend term requires method = 'nested_laplace'", {
  adj <- rook_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)

  expect_error(
    tobs(formula = ~ time + icar(graph = adj, group_var = "cell") +
                     icar(graph = adj, weight = time, group_var = "cell"),
         data = df, family = cover(response = "lognormal"), y = y,
         method = "laplace", control = list(verbose = FALSE)),
    "requires method = 'nested_laplace'")
})

test_that("two unweighted areal terms error (one intercept field only)", {
  adj <- rook_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)

  expect_error(
    tobs(formula = ~ time + icar(graph = adj, group_var = "cell") +
                     bym2(graph = adj, group_var = "cell"),
         data = df, family = cover(response = "lognormal"), y = y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "exactly one unweighted intercept field")
})
