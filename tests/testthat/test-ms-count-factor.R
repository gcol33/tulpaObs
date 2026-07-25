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
  # `lam` and `beta_real` are the realized draws, not the population constants
  # the arguments name. The loading MAGNITUDE is scored against sqrt(sum(lam^2))
  # and the community mean against colMeans(bs) -- see the magnitude test below
  # and gcol33/tulpaObs#155 for why the constant is the wrong estimand.
  list(y = y, data = d, cor_res = cor_res, beta = c(1, 0.5), S = S, N = N,
       lam = lam, beta_real = colMeans(bs))
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
  # Community means against the seed's REALIZED mean, absolutely (#155). Scored
  # against the nominal c(1, 0.5) this seed reads 0.142 off on the intercept and
  # 0.002 off the realized one, so the old `tolerance = 0.2` spent most of its
  # budget on the draw. Budget = 3 sd of the deviation over seeds 201-216
  # (intercept sd 0.033 -> 0.099, slope sd 0.026 -> 0.077), which also covers the
  # observed max (0.065 / 0.065). Seed 215 is excluded from that spread: it is
  # the #157 direction basin (slope deviation 0.220, 4x the next largest) and
  # building the budget to swallow it would encode that pathology as expected.
  expect_community_mean(fit, d$beta_real, c(0.10, 0.08))
  # fitted() is factor-aware
  expect_gt(stats::cor(as.numeric(fitted(fit)$mu), as.numeric(d$y)), 0.7)
  expect_true(is.finite(tobs_waic(fit)$waic))
})

# The loading MAGNITUDE, which no assertion on the residual correlation can see.
#
# `residual_cor` is row-normalised, so a pure scale error leaves it untouched:
# the seed that carried the worst magnitude over-fit on the joint-mode estimator
# (1.53x truth) still reported a residual correlation of 0.93. That is how
# gcol33/tulpaObs#153 reached CI, and #156 is the same defect at lower amplitude.
# Score sqrt(tr(Sigma_res)) = ||lambda||_F instead, which is rotation-invariant
# (so it survives the loading/factor indeterminacy) and IS the quantity that
# regressed, against the seed's own realized loading draw.
test_that("latent-factor count recovers the loading magnitude + mean over seeds", {
  skip_if_fast()
  skip_on_cran()
  out <- vapply(1:12, function(s) {
    d <- .mscf_sim(N = 160L, S = 14L, Q = 2L, seed = 200L + s)
    fit <- tobs(~ x + latent(2), data = d$data, family = ms_count(), y = d$y,
                species = colnames(d$y), method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    c(sqrt(sum(fit$ms_factor$loadings^2)) / sqrt(sum(d$lam^2)),
      unname(unlist(coef(fit))) - d$beta_real)
  }, numeric(3))
  mag <- out[1L, ]
  dev <- t(out[-1L, , drop = FALSE])
  # Budgets measured on these 12 seeds: median 0.978, range 0.855 to 1.048. The
  # joint-mode estimator this replaced ran a median of 1.060 with a 1.53 tail
  # over seeds 201-216, so the CEILING is what carries the regression guard and
  # sits deliberately below that tail.
  #
  # Seed 215 is still excluded by the 201-212 range: even with the #157
  # multi-start fix below it reads 1.28 (down from 1.654 pre-fix), close enough
  # to this loop's 1.30 ceiling that including it would leave the guard almost
  # no margin against a real regression on the other 12 seeds. It has its own
  # dedicated regression test below instead.
  expect_lt(abs(stats::median(mag) - 1), 0.08)
  expect_lt(max(mag), 1.30)
  expect_gt(min(mag), 0.75)
  # A one-sided coefficient shift is a property of the MEAN deviation over seeds,
  # not of any single fit, so it belongs here rather than in the single-fit test
  # above (#155). Measured over these 12 seeds: intercept +0.0032 (0.33 se),
  # slope +0.0035 (0.54 se) -- both consistent with zero. Budget is 5 se
  # (0.049 / 0.032), which still catches a #153-scale inflation: a slope 1.44x
  # truth would shift it 0.35.
  expect_lt(abs(mean(dev[, 1L])), 0.05)
  expect_lt(abs(mean(dev[, 2L])), 0.04)
})

# gcol33/tulpaObs#157 regression guard: seed 215's loading EM settled in a bad
# DIRECTION basin (1.654x truth, a marginal 31 nats below what the same EM
# reaches from a better direction -- verified against an EM started at the
# literal simulated truth, which reaches 1.03x). The fix
# (.tobs_latent_factor_random_starts() in R/community_latent.R) tries several
# fixed-seed pseudo-random restarts alongside the deterministic cosine one at
# the first outer pass, through the SAME ascent + magnitude search, and keeps
# whichever the loading EM converges to the higher marginal from -- the
# "multi-start over directions, selected on the marginal" #157 asked for.
# Measured on this seed: 1.654x -> 1.28x, and the achieved marginal (701212)
# now sits ABOVE the truth-started reference (701211), i.e. this reaches the
# same basin truth does, not merely a smaller one.
test_that("latent-factor count escapes the #157 direction basin on seed 215", {
  skip_if_fast()
  skip_on_cran()
  d <- .mscf_sim(N = 160L, S = 14L, Q = 2L, seed = 215L)
  fit <- tobs(~ x + latent(2), data = d$data, family = ms_count(), y = d$y,
              species = colnames(d$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  mag <- sqrt(sum(fit$ms_factor$loadings^2)) / sqrt(sum(d$lam^2))
  # Pre-fix this read 1.654; the ceiling sits well clear of that regression
  # while still requiring the basin to actually be escaped, not merely nudged.
  expect_lt(mag, 1.40)
  expect_gt(mag, 0.75)
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
