# Control threading on the community latent routes (gcol33/tulpaObs#158).
#
# `.tobs_validate_control()` admits `max.outer` / `factor.starts` (the
# `block_coordinate` group) and `n.quad` (`laplace_em`) on both Laplace routes for
# every family, but whether a value is then handed to
# `.tobs_community_latent_ascent()` is decided separately in each dispatcher. It
# was not: `factor.starts` reached only ms_abun and `n.quad` reached nothing, so a
# passed value was accepted and silently replaced by the driver default.
#
# The knobs leave no trace on a fitted object either, which is why the drift
# needed no mocking to go unnoticed. So the driver now returns the settings it
# RESOLVED and the fit carries them as `fit$latent_control`, and these tests read
# that: it is the same quantity a user has to read to know a control took effect.
#
# Every fit here is deliberately run at `max.outer = 1L` -- these test the
# plumbing, not recovery, and the factor block does not converge early (its
# per-cell offset change decays about 2% a pass), so a default-budget fit would
# cost 150 passes to assert an integer.

.clc_count_sim <- function(N = 40L, S = 5L, Q = 1L, seed = 1L) {
  set.seed(seed)
  d   <- data.frame(x = stats::rnorm(N))
  X   <- stats::model.matrix(~ x, d)
  bs  <- cbind(stats::rnorm(S, 0.8, 0.3), stats::rnorm(S, 0.4, 0.2))
  lam <- matrix(stats::rnorm(S * Q, 0, 0.5), S, Q)
  eta <- matrix(stats::rnorm(N * Q), N, Q)
  mu  <- X %*% t(bs) + eta %*% t(lam)
  list(data = d,
       y_pois = matrix(stats::rpois(N * S, exp(mu)), N, S,
                       dimnames = list(NULL, paste0("sp", seq_len(S)))),
       y_bern = matrix(stats::rbinom(N * S, 1, stats::plogis(mu)), N, S,
                       dimnames = list(NULL, paste0("sp", seq_len(S)))))
}

.clc_occu_sim <- function(N = 40L, S = 5L, J = 3L, Q = 1L, seed = 2L) {
  set.seed(seed)
  d   <- data.frame(x = stats::rnorm(N))
  X   <- stats::model.matrix(~ x, d)
  bp  <- cbind(stats::rnorm(S, 0.2, 0.3), stats::rnorm(S, 0.6, 0.2))
  lam <- matrix(stats::rnorm(S * Q, 0, 0.8), S, Q)
  zt  <- matrix(stats::rnorm(N * Q), N, Q)
  psi <- stats::plogis(X %*% t(bp) + zt %*% t(lam))
  p   <- stats::plogis(matrix(stats::rnorm(S, 0.5, 0.2), N, S, byrow = TRUE))
  y <- array(0L, c(N, J, S), dimnames = list(NULL, NULL, paste0("sp", 1:S)))
  for (s in seq_len(S)) {
    z <- stats::rbinom(N, 1, psi[, s])
    for (j in seq_len(J)) y[, j, s] <- stats::rbinom(N, 1, z * p[, s])
  }
  list(data = d, y = y)
}

.clc_abun_sim <- function(N = 30L, S = 4L, J = 3L, Q = 1L, seed = 3L) {
  set.seed(seed)
  d   <- data.frame(x = stats::rnorm(N))
  X   <- stats::model.matrix(~ x, d)
  bl  <- cbind(stats::rnorm(S, log(3), 0.3), stats::rnorm(S, 0.3, 0.2))
  bp  <- stats::rnorm(S, 0.5, 0.2)
  lam <- matrix(stats::rnorm(S * Q, 0, 0.4), S, Q)
  zt  <- matrix(stats::rnorm(N * Q), N, Q)
  lambda <- exp(X %*% t(bl) + zt %*% t(lam))
  y <- array(NA_integer_, c(N, J, S),
             dimnames = list(NULL, NULL, paste0("sp", seq_len(S))))
  for (s in seq_len(S)) {
    Nn <- stats::rpois(N, lambda[, s])
    for (j in seq_len(J)) y[, j, s] <- stats::rbinom(N, Nn, stats::plogis(bp[s]))
  }
  list(data = d, y = y)
}

# The three settings a caller can set, on the families whose oracle is a
# closed-form density and so cheap enough to fit three times here.
test_that("max.outer / factor.starts / n.quad reach the latent driver", {
  ctl <- list(max.outer = 1L, factor.starts = 3L, n.quad = 7L,
              verbose = FALSE, progress = FALSE)

  d <- .clc_count_sim()

  # ms_count() + latent() -- lfMsAbund
  f_count <- tobs(~ x + latent(1), data = d$data, family = ms_count(),
                  y = d$y_pois, species = colnames(d$y_pois),
                  method = "laplace", control = ctl)
  expect_identical(f_count$latent_control$max.outer, 1L)
  expect_identical(f_count$latent_control$factor.starts, 3L)
  expect_identical(f_count$latent_control$n.quad, 7L)

  # jsdm() -- lfJSDM; shares the ms_count binder and fitter, so this covers the
  # second dispatch site into the same fitter rather than a second code path.
  f_jsdm <- tobs(~ x + latent(1), data = d$data, family = jsdm(),
                 y = d$y_bern, species = colnames(d$y_bern),
                 method = "laplace", control = ctl)
  expect_identical(f_jsdm$latent_control$max.outer, 1L)
  expect_identical(f_jsdm$latent_control$factor.starts, 3L)
  expect_identical(f_jsdm$latent_control$n.quad, 7L)

  # ms_occu() + latent() -- lfMsPGOcc
  o <- .clc_occu_sim()
  f_occu <- tobs(~ x + latent(1), detection = ~ 1, data = o$data,
                 family = ms_occu(), y = o$y, species = dimnames(o$y)[[3L]],
                 method = "laplace", control = ctl)
  expect_identical(f_occu$latent_control$max.outer, 1L)
  expect_identical(f_occu$latent_control$factor.starts, 3L)
  expect_identical(f_occu$latent_control$n.quad, 7L)
})

# The per-family defaults are the point of the report: `max.outer` resolves
# against the family's own measured `factor.outer` when a factor block is present,
# so what a fit ran under is not readable off the call.
test_that("the reported settings carry each family's own resolved defaults", {
  d <- .clc_count_sim()

  # ms_count sets factor.outer = 150 from its measured intercept-bias curve
  # (#156). Reported even when max.outer is overridden, which is what makes the
  # default auditable without paying for it.
  f <- tobs(~ x + latent(1), data = d$data, family = ms_count(), y = d$y_pois,
            species = colnames(d$y_pois), method = "laplace",
            control = list(max.outer = 1L, verbose = FALSE, progress = FALSE))
  expect_identical(f$latent_control$factor.outer, 150L)
  expect_identical(f$latent_control$max.outer, 1L)
  # Driver defaults, unset by this call.
  expect_identical(f$latent_control$factor.starts, 8L)
  expect_identical(f$latent_control$n.quad, 5L)

  # A field-only fit has no factor block, so the factor knobs read NA rather than
  # the number they would have taken -- reporting one would claim an effect the
  # fit does not carry. A field block DOES reach `tol` and break early, so its
  # budget is a cap and the default 25 is cheap.
  A <- matrix(0L, 16L, 16L)
  for (i in 1:15) { A[i, i + 1L] <- 1L; A[i + 1L, i] <- 1L }
  cell <- rep(seq_len(16L), length.out = nrow(d$data))
  dd <- cbind(d$data, cell = cell)
  ff <- tobs(~ x + icar(graph = A, group_var = "cell"), data = dd,
             family = ms_count(), y = d$y_pois, species = colnames(d$y_pois),
             method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  expect_identical(ff$latent_control$max.outer, 25L)
  expect_true(is.na(ff$latent_control$factor.starts))
  expect_true(is.na(ff$latent_control$factor.outer))
  expect_true(is.na(ff$latent_control$n.quad))
})

# The two families whose oracle marginalises a latent state. They are the reason
# `factor.starts` is a per-family setting at all -- each candidate direction costs
# a full loading EM against that oracle -- so they are also the ones where a dead
# knob costs the most.
test_that("the latent-state families thread the same controls", {
  skip_if_fast()
  skip_on_cran()
  ctl <- list(max.outer = 1L, factor.starts = 2L, n.quad = 7L,
              verbose = FALSE, progress = FALSE)

  a <- .clc_abun_sim()
  f_abun <- tobs(~ x + latent(1), detection = ~ 1, data = a$data,
                 family = ms_abun(), y = a$y,
                 species = dimnames(a$y)[[3L]], method = "laplace",
                 control = ctl)
  expect_identical(f_abun$latent_control$max.outer, 1L)
  expect_identical(f_abun$latent_control$factor.starts, 2L)
  expect_identical(f_abun$latent_control$n.quad, 7L)

  cut <- c(0, 25, 50, 75, 100)
  dd <- simulate_ms_distance(n_species = 5, N = 50, cutpoints = cut,
                             n_factors = 1, load_sd = 0.5, seed = 11)
  f_dist <- tobs(~ abund_cov1 + latent(1), detection = ~ 1, data = dd$data,
                 family = ms_distance(cutpoints = cut), y = dd$y,
                 species = dd$species, method = "laplace", control = ctl)
  expect_identical(f_dist$latent_control$max.outer, 1L)
  expect_identical(f_dist$latent_control$factor.starts, 2L)
  expect_identical(f_dist$latent_control$n.quad, 7L)

  # ms_abun sets factor.starts = 1 as its family default: 16 seeds at N=80, S=8,
  # Q=2 produced no magnitude failure at one start, and re-fitting the three worst
  # by residual correlation at eight moved nothing (largest change 0.0027) for a
  # consistent 2.0-2.3x. So the default differs from the driver's 8, and which one
  # a fit used is exactly what is otherwise unreadable.
  f_def <- tobs(~ x + latent(1), detection = ~ 1, data = a$data,
                family = ms_abun(), y = a$y, species = dimnames(a$y)[[3L]],
                method = "laplace",
                control = list(max.outer = 1L, verbose = FALSE,
                               progress = FALSE))
  expect_identical(f_def$latent_control$factor.starts, 1L)
})

# Widening the candidate set cannot lower the objective it is selected on: the
# one-start set (the cosine direction) is a subset of any wider one, and the
# winner is whichever converged loading EM lands highest on the joint marginal.
test_that("a wider factor.starts cannot lower the selected marginal", {
  skip_if_fast()
  skip_on_cran()
  d <- .clc_count_sim(N = 60L, S = 6L, seed = 7L)
  fit1 <- tobs(~ x + latent(1), data = d$data, family = ms_count(),
               y = d$y_pois, species = colnames(d$y_pois), method = "laplace",
               control = list(max.outer = 1L, factor.starts = 1L,
                              verbose = FALSE, progress = FALSE))
  fit4 <- tobs(~ x + latent(1), data = d$data, family = ms_count(),
               y = d$y_pois, species = colnames(d$y_pois), method = "laplace",
               control = list(max.outer = 1L, factor.starts = 4L,
                              verbose = FALSE, progress = FALSE))
  m1 <- fit1$ms_factor$marginal_loglik
  m4 <- fit4$ms_factor$marginal_loglik
  expect_true(is.finite(m1) && is.finite(m4))
  expect_gte(m4, m1 - 1e-6)
})
