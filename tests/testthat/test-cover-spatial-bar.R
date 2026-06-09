# Single-term varying-coefficient spatial bar in cover() (gcol33/tulpaObs#61).
#
# The shared `spatial(~ 1 + time || cell, graph = adj, to = c("presence",
# "positive"))` form is sugar: it desugars to the existing two weighted-areal-
# term coupled trend path (#59) -- the intercept column -> the unweighted
# intercept field, the covariate column -> the weight-scaled trend field, both on
# the bar node index, presence-anchored and copied to the positive arm with an
# estimated coupling alpha. These tests prove the sugar produces the SAME fit as
# the two-term form (byte-identical), that `to =` is order-free and defaults to
# both arms, and that the #61 scope gates fire.

# Rook-adjacency on a g x g grid (self-contained so the file runs in isolation).
.bar_grid_adj <- function(g) {
  n <- g * g
  co <- expand.grid(r = seq_len(g), c = seq_len(g))
  adj <- matrix(0L, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i < j && abs(co$r[i] - co$r[j]) + abs(co$c[i] - co$c[j]) == 1L) {
      adj[i, j] <- 1L; adj[j, i] <- 1L
    }
  }
  adj
}

# Two-field cover hurdle: a shared intercept ICAR field plus a time-weighted
# spatially-varying trend field, both copied onto the positive arm.
.bar_sim_cover_trend <- function(g = 5L, N = 2500L,
                                 sigma1 = 0.8, alpha1 = 1.0,
                                 sigma2 = 0.7, alpha2 = 1.4,
                                 beta_occ = c(0.2, 0.3),
                                 beta_pos = c(-1.6, 0.2),
                                 sd_pos = 0.4, seed = 7) {
  set.seed(seed)
  n_cells <- g * g
  adj <- .bar_grid_adj(g)
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
       y = cover)
}

.bar_trend_control <- list(
  verbose = FALSE, n.threads = 1L, adaptive.grid = TRUE,
  sigma.grid = exp(seq(log(0.2), log(2), length.out = 5)),
  alpha.grid = c(0, exp(seq(log(0.2), log(3), length.out = 4))))

# ---- Byte-identical: bar sugar == two-term coupled form --------------------

test_that("the shared spatial() bar is byte-identical to the two-term form", {
  skip_if_fast()
  skip_on_cran()
  sim <- .bar_sim_cover_trend(seed = 7)

  fit_two <- tobs(
    formula = ~ time + icar(graph = sim$adj, group_var = "cell") +
                icar(graph = sim$adj, weight = time, group_var = "cell"),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  fit_bar <- tobs(
    formula = ~ time +
                spatial(~ 1 + time || cell, graph = sim$adj,
                        to = c("presence", "positive")),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  # Core correctness proof: the desugaring must reproduce the existing
  # machinery's fitted parameter vector, log-marginal, and field exactly.
  expect_equal(fit_bar$beta_occ,      fit_two$beta_occ)
  expect_equal(fit_bar$beta_pos,      fit_two$beta_pos)
  expect_equal(fit_bar$se_occ,        fit_two$se_occ)
  expect_equal(fit_bar$se_pos,        fit_two$se_pos)
  expect_equal(fit_bar$sigma_pos,     fit_two$sigma_pos)
  expect_equal(fit_bar$sigma_trend,   fit_two$sigma_trend)
  expect_equal(fit_bar$alpha_trend,   fit_two$alpha_trend)
  expect_equal(fit_bar$log_marginal,  fit_two$log_marginal)
  expect_equal(fit_bar$joint$modes,   fit_two$joint$modes)
  expect_equal(fit_bar$joint$weights, fit_two$joint$weights)
  expect_identical(fit_bar$n_fields, 2L)
  expect_identical(fit_bar$trend_weight, "time")
})

# ---- to= is order-free (presence anchor regardless of order) ---------------

test_that("to= is unordered: c('positive','presence') == c('presence','positive')", {
  skip_if_fast()
  skip_on_cran()
  sim <- .bar_sim_cover_trend(seed = 7)

  fit_ab <- tobs(
    formula = ~ time +
                spatial(~ 1 + time || cell, graph = sim$adj,
                        to = c("presence", "positive")),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  fit_ba <- tobs(
    formula = ~ time +
                spatial(~ 1 + time || cell, graph = sim$adj,
                        to = c("positive", "presence")),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  expect_equal(fit_ab$beta_occ,    fit_ba$beta_occ)
  expect_equal(fit_ab$beta_pos,    fit_ba$beta_pos)
  expect_equal(fit_ab$sigma_trend, fit_ba$sigma_trend)
  expect_equal(fit_ab$alpha_trend, fit_ba$alpha_trend)
})

# ---- to= omitted defaults to both arms -------------------------------------

test_that("to= omitted defaults to both cover arms (same fit)", {
  skip_if_fast()
  skip_on_cran()
  sim <- .bar_sim_cover_trend(seed = 7)

  fit_explicit <- tobs(
    formula = ~ time +
                spatial(~ 1 + time || cell, graph = sim$adj,
                        to = c("presence", "positive")),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  fit_default <- tobs(
    formula = ~ time + spatial(~ 1 + time || cell, graph = sim$adj),
    data = sim$data, family = cover(positive = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  expect_equal(fit_default$beta_occ,    fit_explicit$beta_occ)
  expect_equal(fit_default$beta_pos,    fit_explicit$beta_pos)
  expect_equal(fit_default$sigma_trend, fit_explicit$sigma_trend)
  expect_equal(fit_default$alpha_trend, fit_explicit$alpha_trend)
})

# ---- Validation / scope gates (no fit, always run) -------------------------

.bar_small_data <- function(n_cells = 16L) {
  adj <- .bar_grid_adj(as.integer(round(sqrt(n_cells))))
  df  <- data.frame(cell = rep(seq_len(n_cells), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  list(adj = adj, df = df, y = y)
}

test_that("an unknown to= arm label errors listing the valid arms", {
  d <- .bar_small_data()
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = d$adj,
                             to = c("occ", "cover")),
         data = d$df, family = cover(positive = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "unknown arm label")
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = d$adj, to = "occurrence"),
         data = d$df, family = cover(positive = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "presence.*positive|positive.*presence")
})

test_that("a correlated `|` bar fits an MCAR field (gcol33/tulpaObs#64)", {
  d <- .bar_small_data()
  # Tiny smoke data: the outer CCD over Sigma is weakly identified and declines
  # to the tensor grid (a benign grid-size note); the assertion is plumbing
  # (structure + summary shape), with parameter recovery in
  # test-cover-spatial-bar-mcar.R.
  fit <- suppressWarnings(tobs(
              formula = ~ time +
                          spatial(~ 1 + time | cell, graph = d$adj,
                                  to = c("presence", "positive")),
              data = d$df, family = cover(positive = "lognormal"), y = d$y,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE, max.iter = 40L,
                             integration = "grid")))
  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$mcar))
  # The cross-covariance summary carries one SD per field and the cross-corr.
  expect_length(fit$sigma_mcar, 2L)
  expect_length(fit$rho_mcar, 1L)
  expect_true(is.finite(fit$alpha_mcar))
})

test_that("a correlated `|` bar requires both cover arms on `to`", {
  d <- .bar_small_data()
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time | cell, graph = d$adj, to = "presence"),
         data = d$df, family = cover(positive = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "both cover arms|arm-specific|separate latent")
})

test_that("a correlated `|` bar cannot co-exist with another areal term", {
  d <- .bar_small_data()
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time | cell, graph = d$adj,
                             to = c("presence", "positive")) +
                     icar(graph = d$adj, group_var = "cell"),
         data = d$df, family = cover(positive = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "whole spatial structure|other areal terms")
})

test_that("a single-arm to= errors pointing to the shared form", {
  d <- .bar_small_data()
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = d$adj, to = "presence"),
         data = d$df, family = cover(positive = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "both cover arms|arm-specific|separate latent")
})

test_that("a node index not matching the graph dimension errors", {
  adj <- .bar_grid_adj(4L)   # 16 nodes
  # cell index runs to 25 -> exceeds the 16-node graph.
  df  <- data.frame(cell = rep(seq_len(25L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = adj,
                             to = c("presence", "positive")),
         data = df, family = cover(positive = "lognormal"), y = y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "node")
})
