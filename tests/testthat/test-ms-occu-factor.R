# Community latent-factor occupancy -- ms_occu() + latent() (the spOccupancy
# lfMsPGOcc analogue; gcol33/tulpaObs#119). Residual species co-occurrence on the
# occupancy arm via Q per-site latent factors + per-species loadings, fit by block
# coordinate ascent (community occupancy EM with the factor offset <-> a two-state
# marginal factor update; R/ms_occu_field.R over the shared driver in
# R/community_latent.R). The loadings / factors are identified only up to
# rotation, so recovery is judged on the residual species correlation matrix
# (Sigma_res = lambda lambda'), which IS identified.
#
# Thresholds are set from a measured 10-seed run at N=250, S=16, J=5
# (dev_notes/probe_ms_occu_factor_seeds.R): median residual correlation 0.885,
# min 0.740. A detection history carries far less information per (site, species)
# than a count, so these sit below the count family's (0.85 median) bar.

# logit psi_{s,i} = X_i (mu + b_s) + sum_q lambda_{s,q} zeta_{q,i}
.msof_sim <- function(N = 250L, S = 16L, Q = 2L, J = 5L, load_sd = 0.8,
                      seed = 1L) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, d)
  b_psi <- vapply(1:2, function(j) stats::rnorm(S, c(0, 0.8)[j], c(0.4, 0.3)[j]),
                  numeric(S))
  b_p   <- stats::rnorm(S, 0.4, 0.3)
  lam   <- matrix(stats::rnorm(S * Q, 0, load_sd), S, Q)
  zeta  <- matrix(stats::rnorm(N * Q), N, Q)
  psi <- stats::plogis(X %*% t(b_psi) + zeta %*% t(lam))
  p   <- stats::plogis(matrix(b_p, N, S, byrow = TRUE))
  y <- array(0L, c(N, J, S))
  for (s in seq_len(S)) {
    z <- stats::rbinom(N, 1, psi[, s])
    for (j in seq_len(J)) y[, j, s] <- stats::rbinom(N, 1, z * p[, s])
  }
  dimnames(y) <- list(NULL, NULL, paste0("sp", seq_len(S)))
  list(y = y, data = d, S = S, beta_psi = c(0, 0.8),
       cor_res = stats::cov2cor(tcrossprod(lam) + diag(1e-8, S)))
}

.msof_fit <- function(d, ...) {
  tobs(~ x + latent(2), data = d$data, family = ms_occu(), detection = ~ 1,
       y = d$y, species = paste0("sp", seq_len(d$S)), method = "laplace",
       control = list(verbose = FALSE, progress = FALSE), ...)
}


test_that("ms_occu() + latent() gates unsupported combinations", {
  d <- .msof_sim(N = 60L, S = 8L, J = 3L, seed = 3L)
  sp <- paste0("sp", seq_len(8))

  # A factor-only model is the block-coordinate Laplace-EM, not nested_laplace
  # (which needs a field) and not the non-spatial NUTS sampler.
  expect_error(
    tobs(~ x + latent(2), data = d$data, family = ms_occu(), detection = ~ 1,
         y = d$y, species = sp, method = "nested_laplace"),
    "block-coordinate|laplace")
  expect_error(
    tobs(~ x + latent(2), data = d$data, family = ms_occu(), detection = ~ 1,
         y = d$y, species = sp, method = "nuts"),
    "block-coordinate|laplace")
  # n_factors must be < n_species
  expect_error(
    tobs(~ x + latent(8), data = d$data, family = ms_occu(), detection = ~ 1,
         y = d$y, species = sp, method = "laplace"),
    "n_factors")
})

test_that("lfMsPGOcc recovers residual co-occurrence and wires S3", {
  skip_on_cran()
  d   <- .msof_sim(seed = 4L)
  fit <- .msof_fit(d)

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "laplace")
  expect_identical(fit$ms_factor$n_factors, 2L)
  expect_equal(dim(fit$ms_factor$loadings), c(16L, 2L))
  expect_equal(dim(fit$ms_factor$residual_cov), c(16L, 16L))
  expect_identical(rownames(fit$ms_factor$loadings), paste0("sp", 1:16))

  # residual species-correlation recovery (identified up to rotation)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.8)
  # community occupancy means still recovered alongside the factors
  expect_equal(unname(fit$means[1:2]), d$beta_psi, tolerance = 0.25)
  # fitted() is factor-aware: psi picks up the per-(site, species) offset
  expect_false(is.null(fit$model$occu_factor_offset))
  expect_equal(dim(fitted(fit)$psi), c(250L, 16L))
  # WAIC scores the factor structure (community_ploglik.R adds the offset)
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("lfMsPGOcc recovers the residual correlation over seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 10L
  rc <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d <- .msof_sim(seed = 300L + s)
    fit <- .msof_fit(d)
    off <- upper.tri(d$cor_res)
    rc[s] <- stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off])
  }
  expect_gt(stats::median(rc), 0.80)
})


# --- spatial-factor composition: a shared ICAR field AND latent factors --------
#
# sfMsPGOcc. The centred loadings (sum_s lambda_sq = 0) make the shared field
# (loading == 1) and the factors orthogonal in species space, so both are
# identified -- but they compete for the same per-site signal, and occupancy is
# information-poor. A measured sweep (dev_notes/probe_sf_occu_diagnose.R) isolates
# the cost: it is NOT a convergence failure (max.outer 20 vs 60 changes the
# residual correlation by 0.002) and NOT the field's presence (a field in the
# model but absent from the data costs ~nothing). It scales with FIELD STRENGTH
# (median rc 0.827 / 0.715 / 0.649 at field_sd 0.5 / 1.0 / 2.0): a strong field
# saturates psi toward 0/1, where the binary curvature psi(1 - psi) collapses and
# a detection history carries almost no information about species-specific
# deviations. More species / visits restores it (J=10 S=30 -> median 0.908,
# min 0.774), which is where these thresholds are set. The field itself recovers
# throughout (>= 0.95).

.msosf_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# logit psi_{s,i} = X_i (mu + b_s) + f_i + sum_q lambda_{s,q} zeta_{q,i}
.msosf_sim <- function(side = 16L, S = 30L, Q = 2L, J = 10L, seed = 1L) {
  set.seed(seed)
  A <- .msosf_grid_graph(side); Ns <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f  <- f - mean(f)
  d  <- data.frame(x = stats::rnorm(Ns))
  X  <- stats::model.matrix(~ x, d)
  b_psi <- vapply(1:2, function(j) stats::rnorm(S, c(0, 0.8)[j], c(0.4, 0.3)[j]),
                  numeric(S))
  b_p   <- stats::rnorm(S, 0.4, 0.3)
  lam  <- scale(matrix(stats::rnorm(S * Q, 0, 0.7), S, Q), scale = FALSE)
  zeta <- matrix(stats::rnorm(Ns * Q), Ns, Q)
  psi <- stats::plogis(X %*% t(b_psi) + matrix(f, Ns, S) + zeta %*% t(lam))
  p   <- stats::plogis(matrix(b_p, Ns, S, byrow = TRUE))
  y <- array(0L, c(Ns, J, S))
  for (s in seq_len(S)) {
    z <- stats::rbinom(Ns, 1, psi[, s])
    for (j in seq_len(J)) y[, j, s] <- stats::rbinom(Ns, 1, z * p[, s])
  }
  dimnames(y) <- list(NULL, NULL, paste0("sp", seq_len(S)))
  list(y = y, data = d, graph = A, f = f, S = S,
       cor_res = stats::cov2cor(tcrossprod(lam) + diag(1e-8, S)))
}

.msosf_fit <- function(d) {
  tobs(~ x + icar(graph = d$graph) + latent(2), data = d$data,
       family = ms_occu(), detection = ~ 1, y = d$y,
       species = paste0("sp", seq_len(d$S)), method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

test_that("sfMsPGOcc recovers BOTH the shared field and the factors", {
  skip_if_fast()
  skip_on_cran()
  d   <- .msosf_sim(seed = 41L)
  fit <- .msosf_fit(d)

  expect_identical(fit$method, "nested_laplace")
  # both latents present: the centred loadings separate the shared spatial mean
  # (the field) from the between-species residual structure (the factors)
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$ms_factor))
  expect_length(fit$spatial_field, 256L)
  expect_gt(stats::cor(fit$spatial_field, d$f), 0.9)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.7)
  # fitted() adds BOTH offsets (field per-site, factors per-(site, species))
  expect_false(is.null(fit$model$occu_field_offset))
  expect_false(is.null(fit$model$occu_factor_offset))
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("sfMsPGOcc recovers both structures over seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 4L
  rc <- numeric(n_seed); fc <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d   <- .msosf_sim(seed = 40L + s)
    fit <- .msosf_fit(d)
    off <- upper.tri(d$cor_res)
    rc[s] <- stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off])
    fc[s] <- stats::cor(fit$spatial_field, d$f)
  }
  expect_gt(stats::median(fc), 0.90)
  expect_gt(stats::median(rc), 0.75)
})
