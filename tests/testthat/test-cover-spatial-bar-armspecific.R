# Arm-specific separate spatial latents in cover() (gcol33/tulpaObs#65).
#
# An INDEPENDENT (`||`) spatial bar placed in ONE arm's formula puts an areal
# field on that arm only, with its own precision and NO cross-arm copy. A field in
# each arm's formula = independent per-arm fields, each its own precision, no
# coupling. This is the FREE counterpart to #61's shared (copied) field (a field
# in the shared formula, or copy()): there the engine anchors one field on
# presence and copies it to positive with an estimated alpha; here each arm
# carries its own field.
#
# Engine mechanism (no engine change needed): the joint multi-block driver reads a
# per-arm `spatial_idx` of 0 as "this arm's rows do not see this block" (the
# `l_b > 0` scatter guard in tulpa's nested_laplace_joint_multi.h). An arm-
# specific block is a NON-copied areal block whose other-arm spatial_idx is the
# all-zero sentinel; its precision integrates on the outer grid.
#
# These tests prove parameter recovery (a statistical fitter is not validated by
# smoke alone): (a) a POSITIVE(cover)-only field is recovered and the presence arm
# shows no spurious field; (b) two INDEPENDENT per-arm fields are both recovered
# with no cross-arm coupling; (c) on genuinely independent per-arm fields the arm-
# specific fit recovers the arm the shared (copied) fit must miss.

# Rook-adjacency on a g x g grid (self-contained so the file runs in isolation).
.as_grid_adj <- function(g) {
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

# A smooth, mean-zero, unit-SD ICAR-like field over the grid graph.
.as_smooth_field <- function(adj) {
  n <- nrow(adj); z <- rnorm(n)
  for (it in 1:60) {
    zn <- z
    for (i in seq_len(n)) {
      nb <- which(adj[i, ] == 1L); zn[i] <- 0.35 * z[i] + 0.65 * mean(z[nb])
    }
    z <- zn
  }
  z <- z - mean(z); z / stats::sd(z)
}

.as_control <- list(
  verbose = FALSE, progress = FALSE, n.threads = 1L, adaptive.grid = TRUE,
  integration = "grid",
  sigma.grid = exp(seq(log(0.3), log(2.5), length.out = 6)))

# Posterior-mean per-cell field (mean-centred) and per-arm amplitudes for the
# b-th block of a fitted arm-specific cover_fit, via the joint-draw projection.
.as_block_field <- function(fit, b, n = 400L) {
  bundle <- .tobs_joint_draws(fit, n = n)
  blk <- bundle$blocks[[b]]
  z <- colMeans(blk$z); z <- z - mean(z)
  list(z = z, amp_occ = mean(blk$amp_occ), amp_pos = mean(blk$amp_pos))
}


# ---- (a) positive(cover)-only field; presence has no spatial structure ------

test_that("a positive-only single-arm field recovers, presence shows none", {
  skip_if_fast()
  skip_on_cran()
  set.seed(3)
  g <- 6L; n_cells <- g * g
  adj <- .as_grid_adj(g)
  z_pos <- .as_smooth_field(adj)
  N <- 4000L
  cell <- sample.int(n_cells, N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sigma_field <- 0.9; beta_occ <- c(0.3, 0.4); beta_pos <- c(-1.0, 0.25)

  eta_occ <- beta_occ[1] + beta_occ[2] * x                 # presence: NO field
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sigma_field * z_pos[cell]
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.4)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, x = x, cover = cover)

  fit <- tobs(presence = ~ x, positive = ~ x + spatial(~ 1 || cell, graph = adj),
              data = dat, family = cover(response = "lognormal"), y = dat$cover,
              method = "nested_laplace", control = .as_control)

  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$armspecific))
  expect_identical(fit$n_fields, 1L)
  expect_identical(fit$armspec_blocks[[1L]]$arm, "positive")

  # Fixed effects recover.
  expect_equal(unname(fit$beta_occ), beta_occ, tolerance = 0.1)
  expect_equal(unname(fit$beta_pos), beta_pos, tolerance = 0.15)

  # The positive-arm field is recovered (shape); the presence arm carries NONE.
  fld <- .as_block_field(fit, 1L)
  truth <- sigma_field * (z_pos - mean(z_pos))
  expect_gt(cor(fld$z, truth), 0.85)
  expect_gt(fld$amp_pos, 0.1)          # positive arm carries the field
  expect_equal(fld$amp_occ, 0)         # presence arm carries no field (no copy)
  # Field amplitude (precision) integrated on the grid, surfaced and positive.
  expect_true(is.finite(fit$sigma_armspecific[["positive.Intercept"]]))
  expect_gt(fit$sigma_armspecific[["positive.Intercept"]], 0)
})


# ---- (b) two independent per-arm fields, no cross-arm coupling --------------

test_that("two independent single-arm fields both recover, no cross-arm copy", {
  skip_if_fast()
  skip_on_cran()
  set.seed(21)
  g <- 8L; n_cells <- g * g
  adj <- .as_grid_adj(g)
  z_pre <- .as_smooth_field(adj)
  # Orthogonalize the positive field against presence -> genuinely independent
  # realizations (cor ~ 0), so a copied field could not fit both arms.
  z_raw <- .as_smooth_field(adj)
  z_pos <- z_raw - sum(z_raw * z_pre) / sum(z_pre * z_pre) * z_pre
  z_pos <- z_pos - mean(z_pos); z_pos <- z_pos / stats::sd(z_pos)
  expect_lt(abs(cor(z_pre, z_pos)), 0.05)

  N <- 5000L
  cell <- sample.int(n_cells, N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sig_pre <- 1.0; sig_pos <- 0.9
  beta_occ <- c(0.2, 0.3); beta_pos <- c(-0.8, 0.2)

  eta_occ <- beta_occ[1] + beta_occ[2] * x + sig_pre * z_pre[cell]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sig_pos * z_pos[cell]
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.35)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, x = x, cover = cover)

  fit <- tobs(presence = ~ x + spatial(~ 1 || cell, graph = adj),
              positive = ~ x + spatial(~ 1 || cell, graph = adj),
              data = dat, family = cover(response = "lognormal"), y = dat$cover,
              method = "nested_laplace", control = .as_control)

  expect_true(isTRUE(fit$armspecific))
  expect_identical(fit$n_fields, 2L)
  # Block order follows formula order: presence block first, positive second.
  arms <- vapply(fit$armspec_blocks, function(m) m$arm, character(1))
  expect_identical(arms, c("presence", "positive"))

  pre <- .as_block_field(fit, 1L)
  pos <- .as_block_field(fit, 2L)

  # Each per-arm field is recovered on its own arm.
  expect_gt(cor(pre$z, sig_pre * (z_pre - mean(z_pre))), 0.85)
  expect_gt(cor(pos$z, sig_pos * (z_pos - mean(z_pos))), 0.85)

  # No cross-arm copy: the presence block contributes ZERO to the positive arm
  # and vice versa (the 0-sentinel spatial_idx), and the recovered presence field
  # is NOT the positive field (independent realizations stay separate).
  expect_equal(pre$amp_pos, 0)
  expect_equal(pos$amp_occ, 0)
  expect_lt(abs(cor(pre$z, sig_pos * (z_pos - mean(z_pos)))), 0.3)
  expect_lt(abs(cor(pos$z, sig_pre * (z_pre - mean(z_pre)))), 0.3)
})


# ---- (c) arm-specific recovers where the shared (copied) fit cannot ---------

test_that("arm-specific fit recovers an arm the shared copied fit must miss", {
  skip_if_fast()
  skip_on_cran()
  set.seed(21)
  g <- 8L; n_cells <- g * g
  adj <- .as_grid_adj(g)
  z_pre <- .as_smooth_field(adj)
  z_raw <- .as_smooth_field(adj)
  z_pos <- z_raw - sum(z_raw * z_pre) / sum(z_pre * z_pre) * z_pre
  z_pos <- z_pos - mean(z_pos); z_pos <- z_pos / stats::sd(z_pos)

  N <- 5000L
  cell <- sample.int(n_cells, N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sig_pre <- 1.0; sig_pos <- 0.9
  beta_occ <- c(0.2, 0.3); beta_pos <- c(-0.8, 0.2)

  eta_occ <- beta_occ[1] + beta_occ[2] * x + sig_pre * z_pre[cell]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sig_pos * z_pos[cell]
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.35)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, x = x, cover = cover)

  fit_arm <- tobs(presence = ~ x + spatial(~ 1 || cell, graph = adj),
                  positive = ~ x + spatial(~ 1 || cell, graph = adj),
                  data = dat, family = cover(response = "lognormal"),
                  y = dat$cover, method = "nested_laplace", control = .as_control)

  fit_shared <- tobs(formula = ~ x +
                       spatial(~ 1 || cell, graph = adj),
                     data = dat, family = cover(response = "lognormal"),
                     y = dat$cover, method = "nested_laplace",
                     control = .as_control)

  # Arm-specific recovers BOTH per-arm fields.
  pre <- .as_block_field(fit_arm, 1L)
  pos <- .as_block_field(fit_arm, 2L)
  cor_arm_pre <- cor(pre$z, z_pre - mean(z_pre))
  cor_arm_pos <- cor(pos$z, z_pos - mean(z_pos))
  expect_gt(cor_arm_pre, 0.85)
  expect_gt(cor_arm_pos, 0.85)

  # The shared fit forces ONE field copied across both arms. On genuinely
  # independent per-arm truths it can latch onto at most one arm's structure;
  # the OTHER arm's field is then badly fit. The arm-specific fit recovers both,
  # so it strictly beats the shared fit on whichever arm the shared fit misses.
  bsh <- .tobs_joint_draws(fit_shared, n = 400L)
  zsh <- colMeans(bsh$blocks[[1L]]$z); zsh <- zsh - mean(zsh)
  cor_sh_pre <- abs(cor(zsh, z_pre - mean(z_pre)))
  cor_sh_pos <- abs(cor(zsh, z_pos - mean(z_pos)))
  shared_missed <- min(cor_sh_pre, cor_sh_pos)
  armspecific_on_missed <- if (cor_sh_pre <= cor_sh_pos) cor_arm_pre else cor_arm_pos

  # On the arm the shared fit misses, the arm-specific fit recovers it much better.
  expect_lt(shared_missed, 0.5)
  expect_gt(armspecific_on_missed, 0.85)
  expect_gt(armspecific_on_missed - shared_missed, 0.3)
})


# ---- (d) predict projects the per-arm fields, not flat maps (#95) -----------

test_that("arm-specific predict() projects each per-arm field (not flat)", {
  skip_if_fast()
  skip_on_cran()
  set.seed(7)
  g <- 6L; n_cells <- g * g
  adj <- .as_grid_adj(g)
  z_pre <- .as_smooth_field(adj)
  z_raw <- .as_smooth_field(adj)
  z_pos <- z_raw - sum(z_raw * z_pre) / sum(z_pre * z_pre) * z_pre
  z_pos <- z_pos - mean(z_pos); z_pos <- z_pos / stats::sd(z_pos)

  N <- 4000L
  cell <- sample.int(n_cells, N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sig_pre <- 1.0; sig_pos <- 0.9
  beta_occ <- c(0.2, 0.3); beta_pos <- c(-0.8, 0.2)

  eta_occ <- beta_occ[1] + beta_occ[2] * x + sig_pre * z_pre[cell]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sig_pos * z_pos[cell]
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.35)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, x = x, cover = cover)

  fit <- tobs(presence = ~ x + spatial(~ 1 || cell, graph = adj),
              positive = ~ x + spatial(~ 1 || cell, graph = adj),
              data = dat, family = cover(response = "lognormal"),
              y = dat$cover, method = "nested_laplace", control = .as_control)

  nd <- data.frame(cell = seq_len(n_cells), x = 0)
  occ  <- as.data.frame(predict(fit, newdata = nd, type = "occurrence",
                                nsim = 300L, draws = FALSE))$mean
  cond <- as.data.frame(predict(fit, newdata = nd, type = "cover_cond",
                                nsim = 300L, draws = FALSE))$mean

  # Per-cell predictions vary (the bug projected a flat map: n_distinct == 1).
  expect_gt(length(unique(round(occ, 6))), 1L)
  expect_gt(length(unique(round(cond, 6))), 1L)
  expect_gt(stats::sd(occ), 1e-4)
  expect_gt(stats::sd(cond), 1e-4)

  # The projected map tracks the arm's own field, not the other arm's: occurrence
  # follows the presence field, conditional cover follows the positive field.
  expect_gt(cor(qlogis(pmin(pmax(occ, 1e-6), 1 - 1e-6)), z_pre - mean(z_pre)), 0.8)
  expect_gt(cor(log(cond), z_pos - mean(z_pos)), 0.7)
})

test_that("intercept-only arm-specific predict() needs no time_col (#95)", {
  skip_if_fast()
  skip_on_cran()
  set.seed(11)
  adj <- .as_grid_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L), x = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  fit <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 || cell, graph = adj),
    positive = ~ x + spatial(~ 1 || cell, graph = adj),
    data = df, family = cover(response = "lognormal"), y = y,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE, integration = "grid")))
  nd <- data.frame(cell = seq_len(16L), x = 0)
  # No trend field is present, so predict must not demand a `time_col`.
  expect_no_error(
    occ <- as.data.frame(predict(fit, newdata = nd, type = "occurrence",
                                 nsim = 100L, draws = FALSE))$mean)
  expect_gt(length(unique(round(occ, 6))), 1L)
})


# ---- Scope gates (no fit, always run) --------------------------------------

.as_small <- function() {
  adj <- .as_grid_adj(4L)
  df  <- data.frame(cell = rep(seq_len(16L), length.out = 64L), x = rnorm(64L))
  y   <- ifelse(rbinom(64L, 1, 0.5) == 1L, runif(64L, 0.01, 0.9), 0)
  list(adj = adj, df = df, y = y)
}

test_that("a single-arm || bar requires the nested-Laplace method", {
  d <- .as_small()
  expect_error(
    tobs(presence = ~ x, positive = ~ x + spatial(~ 1 || cell, graph = d$adj),
         data = d$df, family = cover(response = "lognormal"), y = d$y,
         method = "laplace", control = list(verbose = FALSE)),
    "arm-specific spatial bar|nested_laplace")
})

test_that("a single-arm correlated `|` bar fits on that arm alone (no copy, #109)", {
  d <- .as_small()
  fit <- suppressWarnings(tobs(
    presence = ~ x, positive = ~ x + spatial(~ 1 + x | cell, graph = d$adj),
    data = d$df, family = cover(response = "lognormal"), y = d$y,
    method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE)))
  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$mcar))
  expect_length(fit$sigma_mcar, 2L)
  expect_length(fit$rho_mcar, 1L)
  expect_true(is.na(fit$alpha_mcar))   # single arm: no cross-arm copy
})

test_that("an arm-specific bar cannot be mixed with a shared field", {
  d <- .as_small()
  # A shared field (presence field copied onto positive) plus an arm-specific
  # positive field in the same fit is the disallowed mix.
  expect_error(
    tobs(presence = ~ x + spatial(~ 1 || cell, graph = d$adj),
         positive = ~ x + copy(spatial()) +
                    spatial(~ 0 + x || cell, graph = d$adj),
         data = d$df, family = cover(response = "lognormal"), y = d$y,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "separate per-arm structure|cannot be combined")
})

test_that("two arm-specific fields on the SAME arm error", {
  d <- .as_small()
  expect_error(
    suppressWarnings(tobs(
      presence = ~ 1,
      positive = ~ spatial(~ 1 || cell, graph = d$adj) +
                   spatial(~ 0 + x || cell, graph = d$adj),
      data = d$df, family = cover(response = "lognormal"), y = d$y,
      method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))),
    "distinct arm")
})

test_that("a presence-only single-arm bar is accepted and fits", {
  d <- .as_small()
  fit <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 || cell, graph = d$adj), positive = ~ x,
    data = d$df, family = cover(response = "lognormal"), y = d$y,
    method = "nested_laplace",
    control = list(verbose = FALSE, progress = FALSE, integration = "grid")))
  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$armspecific))
  expect_identical(fit$armspec_blocks[[1L]]$arm, "presence")
})


# ---- (e) BYM2 arm-specific field on the cover arm (gcol33/tulpaObs#107) ------

test_that("a positive-only BYM2 single-arm field recovers the rho-mixed field", {
  skip_if_fast()
  skip_on_cran()
  set.seed(7)
  g <- 6L; n_cells <- g * g
  adj <- .as_grid_adj(g)
  z_pos <- .as_smooth_field(adj)          # mostly structured -> high rho
  N <- 4000L
  cell <- sample.int(n_cells, N, replace = TRUE)
  x    <- as.numeric(scale(rnorm(N)))
  sigma_field <- 0.9; beta_occ <- c(0.3, 0.4); beta_pos <- c(-1.0, 0.25)

  eta_occ <- beta_occ[1] + beta_occ[2] * x                 # presence: NO field
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + sigma_field * z_pos[cell]
  cover   <- ifelse(occur == 1L,
                    pmin(exp(eta_pos + rnorm(N, 0, 0.4)), 1 - 1e-6), 0)
  dat <- data.frame(cell = cell, x = x, cover = cover)

  fit <- tobs(
    presence = ~ x,
    positive = ~ x + spatial(~ 1 || cell, graph = adj, model = "bym2"),
    data = dat, family = cover(response = "lognormal"), y = dat$cover,
    method = "nested_laplace",
    control = utils::modifyList(.as_control,
                                list(rho.grid = c(0.2, 0.5, 0.8, 0.95))))

  expect_s3_class(fit, "cover_fit")
  expect_true(isTRUE(fit$armspecific))
  expect_identical(fit$n_fields, 1L)
  expect_identical(fit$armspec_blocks[[1L]]$arm, "positive")
  expect_identical(fit$armspec_blocks[[1L]]$type, "bym2")

  expect_equal(unname(fit$beta_occ), beta_occ, tolerance = 0.1)
  expect_equal(unname(fit$beta_pos), beta_pos, tolerance = 0.15)

  # The reconstructed rho-mixed unit field (phi + theta) recovers the truth, and
  # the presence arm carries none. This exercises the BYM2 projection.
  fld <- .as_block_field(fit, 1L)
  truth <- sigma_field * (z_pos - mean(z_pos))
  expect_gt(cor(fld$z, truth), 0.85)
  expect_gt(fld$amp_pos, 0.1)
  expect_equal(fld$amp_occ, 0)
  expect_true(is.finite(fit$sigma_armspecific[["positive.Intercept"]]))
  expect_gt(fit$sigma_armspecific[["positive.Intercept"]], 0)

  # Structured fraction: a smooth field is mostly structured, so the grid-weighted
  # posterior rho sits in the upper half of the grid.
  tg <- fit$joint$theta_grid; w <- fit$joint$weights
  fin <- is.finite(w) & w > 0; w <- w[fin] / sum(w[fin]); tg <- tg[fin, , drop = FALSE]
  rho_post <- sum(w * as.numeric(tg[, "b1.rho"]))
  expect_gt(rho_post, 0.5)
})
