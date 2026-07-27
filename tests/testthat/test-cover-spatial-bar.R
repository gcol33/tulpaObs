# Single-term varying-coefficient spatial bar in cover() (gcol33/tulpaObs#61).
#
# A `spatial(~ 1 + time || cell, graph = adj)` bar in a shared cover() formula is
# sugar: it desugars to the existing two weighted-areal-term coupled trend path
# (#59) -- the intercept column -> the unweighted intercept field, the covariate
# column -> the weight-scaled trend field, both on the bar node index,
# presence-anchored and copied to the positive arm with an estimated coupling
# alpha. A bar in a shared formula reaches both cover arms; an arm-specific field
# is written in that arm's per-arm formula (placement). These tests prove the
# sugar produces the SAME fit as the two-term form (byte-identical) and that the
# #61 scope gates fire.

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
    data = sim$data, family = cover(response = "lognormal"), y = sim$y,
    method = "nested_laplace", control = .bar_trend_control)

  fit_bar <- tobs(
    formula = ~ time +
                spatial(~ 1 + time || cell, graph = sim$adj),
    data = sim$data, family = cover(response = "lognormal"), y = sim$y,
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

test_that("the coupled intercept+trend || bar recovers both field amplitudes + occurrence slope (#140)", {
  skip_if_fast()
  skip_on_cran()
  # End-to-end recovery for the flagship coupled intercept+trend bar on the cover
  # hurdle (the byte-identical test above only proves it reproduces the two-term
  # machinery). The cover-arm FE time slope beta_pos[2] is CONFOUNDED with the
  # time-weighted trend field -- both multiply `time` -- so it attenuates (measured
  # mean ~0.12 vs 0.2, coverage ~0.2); that is a model property (a linear trend FE
  # competes with a trend field), not a fitter bug, and is deliberately not
  # asserted here.
  #
  # The field amplitudes are asserted as a RATIO, not against their simulated
  # values. `sigma` / `sigma_trend` are ICAR CONDITIONAL-scale hyperparameters,
  # and the map from those onto a field's marginal SD is an engine
  # parameterisation, not a property of the data. Re-measured over seeds 701-712
  # against tulpa 0.0.95, BOTH hyperparameters sit a common factor below their
  # simulated marginal SDs -- sigma 0.538 vs 0.8 and sigma_trend 0.478 vs 0.7,
  # -33% and -32%, sd 0.031 with no seed anywhere near 0.7 -- while the FE slope
  # (0.275), its coverage (0.92) and alpha_trend (1.14) are where they were. A
  # uniform factor on two independently simulated fields with different truths is
  # the signature of a scale convention, and it cancels in the ratio: measured
  # sigma/sigma_trend 1.129 against a simulated 0.8/0.7 = 1.143, 1.2% over 12
  # seeds and 4.1% over the 6 run here. So the ratio is the part of the two
  # amplitudes that survives a reparameterisation, and that is what is asserted.
  # The cover joint fit surfaces no per-cell field vector, so the
  # correlation-with-truth idiom the occu() SVC recovery suites use
  # (test-occu-spatial-svc-recovery.R) is not available on this path.
  n_seeds <- 6L
  bo2 <- st <- sg <- numeric(n_seeds); co <- logical(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- .bar_sim_cover_trend(seed = 700L + r)
    fit <- tobs(~ time + spatial(~ 1 + time || cell, graph = sim$adj),
                data = sim$data, family = cover(response = "lognormal"), y = sim$y,
                method = "nested_laplace", control = .bar_trend_control)
    bo2[r] <- fit$beta_occ[2]
    st[r]  <- fit$sigma_trend
    sg[r]  <- fit$hyperpar$spatial[["b1.sigma"]]
    co[r]  <- abs(fit$beta_occ[2] - 0.3) <= 1.96 * fit$se_occ[2]
  }
  expect_lt(abs(mean(bo2) - 0.3), 0.12)          # occurrence slope recovers
  expect_gte(mean(co), 0.6)                       # occurrence-slope CI coverage

  # Relative amplitude of the two fields: convention-free, so this is the
  # assertion that survives an engine reparameterisation.
  expect_lt(abs(mean(sg / st) / (0.8 / 0.7) - 1), 0.15)

  # Both amplitudes stay positive and finite, and inside the bracket measured on
  # the current parameterisation (sigma_trend mean 0.478, sd 0.031 over 12 seeds).
  # A gross-regression guard, NOT a recovery claim -- the recovery claim is the
  # ratio above.
  expect_true(all(is.finite(c(sg, st))) && all(c(sg, st) > 0))
  expect_lt(abs(mean(st) - 0.48), 0.15)
})

# ---- Validation / scope gates (no fit, always run) -------------------------

.bar_small_data <- function(n_cells = 16L) {
  adj <- .bar_grid_adj(as.integer(round(sqrt(n_cells))))
  df  <- data.frame(cell = rep(seq_len(n_cells), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  list(adj = adj, df = df, y = y)
}

test_that("`to =` is not a bar argument (retired: arm chosen by placement)", {
  d <- .bar_small_data()
  # `to =` is gone from the call surface: the bar form takes only `graph`/`by`,
  # so an explicit arm argument is an unknown argument. The arm is chosen by
  # placement (write the field in that arm's formula) and shared with copy().
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = d$adj,
                             to = c("presence", "positive")),
         data = d$df, family = cover(response = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "unexpected argument `to`")
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = d$adj, to = "positive"),
         data = d$df, family = cover(response = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "unexpected argument `to`")
})

test_that("a correlated `|` bar fits an MCAR field (gcol33/tulpaObs#64)", {
  d <- .bar_small_data()
  # Tiny smoke data: the outer CCD over Sigma is weakly identified and declines
  # to the tensor grid (a benign grid-size note); the assertion is plumbing
  # (structure + summary shape), with parameter recovery in
  # test-cover-spatial-bar-mcar.R.
  fit <- suppressWarnings(tobs(
              formula = ~ time +
                          spatial(~ 1 + time | cell, graph = d$adj),
              data = d$df, family = cover(response = "lognormal"), y = d$y,
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

test_that("a single-arm correlated `|` bar fits on the occurrence arm alone (#109)", {
  d <- .bar_small_data()
  # Placement: the correlated bar written in the presence formula only is the
  # free-Sigma field on that arm, with no cross-arm copy.
  fit <- suppressWarnings(tobs(
    presence = ~ time + spatial(~ 1 + time | cell, graph = d$adj),
    positive = ~ time,
    data = d$df, family = cover(response = "lognormal"), y = d$y,
    method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE)))
  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$mcar))
  expect_length(fit$sigma_mcar, 2L)
  expect_length(fit$rho_mcar, 1L)
  expect_true(is.na(fit$alpha_mcar))   # single arm: no cross-arm copy
})

test_that("a correlated `|` bar cannot co-exist with another areal term", {
  d <- .bar_small_data()
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time | cell, graph = d$adj) +
                     icar(graph = d$adj, group_var = "cell"),
         data = d$df, family = cover(response = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "whole spatial structure|other areal terms")
})

test_that("a single-arm || placement is wired as an arm-specific separate latent", {
  # gcol33/tulpaObs#65: an INDEPENDENT (`||`) bar placed in one arm's formula
  # fits an arm-specific separate field on that arm only, with its own precision
  # and no cross-arm copy. Recovery lives in
  # test-cover-spatial-bar-armspecific.R; here the assertion is plumbing.
  d <- .bar_small_data()
  fit <- suppressWarnings(tobs(
    presence = ~ time +
                spatial(~ 1 + time || cell, graph = d$adj),
    positive = ~ time,
    data = d$df, family = cover(response = "lognormal"), y = d$y,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE, integration = "grid")))
  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$armspecific))
  expect_identical(fit$armspec_blocks[[1L]]$arm, "presence")
})

test_that("a node index not matching the graph dimension errors", {
  adj <- .bar_grid_adj(4L)   # 16 nodes
  # cell index runs to 25 -> exceeds the 16-node graph.
  df  <- data.frame(cell = rep(seq_len(25L), length.out = 64L),
                    time = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  expect_error(
    tobs(formula = ~ time +
                     spatial(~ 1 + time || cell, graph = adj),
         data = df, family = cover(response = "lognormal"), y = y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "node")
})
