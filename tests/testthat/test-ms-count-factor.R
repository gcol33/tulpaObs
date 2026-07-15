# Community latent-factor count -- ms_count() + latent() (the spAbundance
# lfMsAbund analogue, Poisson; gcol33/tulpaObs#117). Residual species
# co-occurrence via Q per-site latent factors + per-species loadings, fit by
# block coordinate ascent (community EM with the factor offset <-> a Poisson
# factor update; R/ms_count_factor.R). The loadings / factors are identified
# only up to rotation, so recovery is judged on the residual species correlation
# matrix (Sigma_res = lambda lambda'), which IS identified.

# log mu_{s,i} = X_i (mu + b_s) + sum_q lambda_{s,q} eta_{q,i}.
.mscf_sim <- function(N = 160L, S = 14L, Q = 2L, load_sd = 0.6, seed = 1L) {
  set.seed(seed)
  d  <- data.frame(x = stats::rnorm(N))
  X  <- stats::model.matrix(~ x, d)
  bs <- vapply(1:2, function(j) stats::rnorm(S, c(1, 0.5)[j], c(0.4, 0.3)[j]),
               numeric(S))
  lam <- matrix(stats::rnorm(S * Q, 0, load_sd), S, Q)
  eta <- matrix(stats::rnorm(N * Q), N, Q)
  y <- matrix(stats::rpois(N * S, exp(pmin(X %*% t(bs) + eta %*% t(lam), 700))),
              N, S, dimnames = list(NULL, paste0("sp", seq_len(S))))
  cor_res <- stats::cov2cor(tcrossprod(lam) + diag(1e-8, S))
  list(y = y, data = d, cor_res = cor_res, beta = c(1, 0.5), S = S, N = N)
}


test_that("ms_count() + latent() gates unsupported combinations", {
  d <- .mscf_sim(N = 60L, S = 8L, seed = 3L)

  # nested_laplace / nuts are not the factor engine
  expect_error(
    tobs(~ x + latent(2), data = d$data, family = ms_count(), y = d$y,
         species = colnames(d$y), method = "nested_laplace"),
    "block-coordinate|laplace")
  # A per-site latent structure supports the responses with no dispersion
  # parameter -- Poisson (ms_count) and Bernoulli (jsdm). A negbin size /
  # Gaussian residual variance is not identified against it.
  expect_error(
    tobs(~ x + latent(2), data = d$data, family = ms_count("negbin"), y = d$y,
         species = colnames(d$y), method = "laplace"),
    "not identified|Poisson")
})

test_that("a latent-factor count fit recovers residual co-occurrence + S3", {
  skip_on_cran()
  d <- .mscf_sim(N = 160L, S = 14L, Q = 2L, seed = 4L)
  fit <- tobs(~ x + latent(2), data = d$data, family = ms_count(), y = d$y,
              species = colnames(d$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$ms_factor$n_factors, 2L)
  expect_equal(dim(fit$ms_factor$loadings), c(14L, 2L))
  expect_equal(dim(fit$ms_factor$residual_cov), c(14L, 14L))
  # residual species-correlation recovery (identified up to rotation)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.8)
  # community means recovered
  expect_equal(unname(unlist(coef(fit))), d$beta, tolerance = 0.2)
  # fitted() is factor-aware
  expect_gt(stats::cor(as.numeric(fitted(fit)$mu), as.numeric(d$y)), 0.7)
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("latent-factor count recovers the residual correlation over seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 12L
  rc <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d <- .mscf_sim(N = 160L, S = 14L, Q = 2L, seed = 200 + s)
    fit <- tobs(~ x + latent(2), data = d$data, family = ms_count(), y = d$y,
                species = colnames(d$y), method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    off  <- upper.tri(d$cor_res)
    rc[s] <- stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off])
  }
  expect_gt(stats::median(rc), 0.85)
})


# --- spatial-factor composition: a shared field AND latent factors together ----

.mscsf_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# log mu_{s,i} = X_i (mu + b_s) + f_i + sum_q lambda_{s,q} eta_{q,i}, with a shared
# ICAR field f AND Q centred-loading factors (residual co-occurrence).
.mscsf_sim <- function(side = 11L, S = 16L, Q = 2L, seed = 1L) {
  set.seed(seed)
  A <- .mscsf_grid_graph(side); Ns <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- 0.6 * scale(sin(co$r/side*pi) + cos(co$c/side*pi))[, 1]; f <- f - mean(f)
  d <- data.frame(x = stats::rnorm(Ns))
  X <- stats::model.matrix(~ x, d)
  bs  <- vapply(1:2, function(j) stats::rnorm(S, c(1, 0.5)[j], c(0.4, 0.3)[j]),
                numeric(S))
  lam <- scale(matrix(stats::rnorm(S * Q, 0, 0.5), S, Q), scale = FALSE)
  eta <- matrix(stats::rnorm(Ns * Q), Ns, Q)
  lin <- X %*% t(bs) + matrix(f, Ns, S) + eta %*% t(lam)
  y <- matrix(stats::rpois(Ns * S, exp(pmin(lin, 700))), Ns, S,
              dimnames = list(NULL, paste0("sp", seq_len(S))))
  cor_res <- stats::cov2cor(tcrossprod(lam) + diag(1e-8, S))
  list(y = y, data = d, graph = A, f = f, cor_res = cor_res, Ns = Ns)
}

test_that("spatial-factor count recovers BOTH the shared field and the factors", {
  skip_on_cran()
  d <- .mscsf_sim(side = 11L, S = 16L, Q = 2L, seed = 6L)
  fit <- tobs(~ x + icar(graph = d$graph) + latent(2), data = d$data,
              family = ms_count(), y = d$y, species = colnames(d$y),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  # both latents present + recovered (the centred loadings separate them)
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$ms_factor))
  expect_gt(stats::cor(fit$spatial_field, d$f), 0.8)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.8)
  # fitted() adds both offsets
  expect_gt(stats::cor(as.numeric(fitted(fit)$mu), as.numeric(d$y)), 0.7)
  expect_true(is.finite(tobs_waic(fit)$waic))
})
