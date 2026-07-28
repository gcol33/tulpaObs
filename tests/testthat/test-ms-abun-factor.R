# Community N-mixture + latent() factors -- the spAbundance lfMsNMix analogue,
# and the spatial-factor case (a shared field alongside the factors).
# gcol33/tulpaObs#117. Poisson.
#
#   N_{s,i} ~ Poisson(lambda_{s,i}),  y_{s,i,j} | N ~ Binomial(N_{s,i}, p_{s,i,j})
#   log lambda_{s,i} = X_i (mu + b_s) [+ f_i] + sum_q lambda_{s,q} zeta_{q,i}
#
# The latent N still integrates out in closed form per species-site, so the whole
# latent structure sits on eta_lambda and the family reduces to one working
# oracle over the Royle marginal (R/ms_abun_latent.R):
#   score = grad_eta_lambda,  curv = info_eta_lambda - var_N * score_wt_lambda^2
# (the Louis (1982) abundance curvature with the detection arm profiled out),
# which drives the shared block-coordinate engine in R/community_latent.R.
#
# The loadings / factors are identified only up to rotation, so recovery is
# judged on the residual species correlation (Sigma_res = lambda lambda'), which
# IS identified -- as for every other community factor family.

.msaf_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

.msaf_sim <- function(N = 120L, S = 12L, J = 4L, Q = 2L, load_sd = 0.5,
                      field = NULL, seed = 1L) {
  set.seed(seed)
  d  <- data.frame(x = stats::rnorm(N))
  X  <- stats::model.matrix(~ x, d)
  bl <- cbind(stats::rnorm(S, log(4), 0.4), stats::rnorm(S, 0.4, 0.3))
  bp <- matrix(stats::rnorm(S, 0.5, 0.3), S, 1L)
  lam <- matrix(stats::rnorm(S * Q, 0, load_sd), S, Q)
  # With a shared field the loadings are centred across species: the field owns
  # the shared spatial mean, the factors the between-species residual.
  if (!is.null(field)) lam <- scale(lam, scale = FALSE)
  zeta <- matrix(stats::rnorm(N * Q), N, Q)
  f <- if (is.null(field)) numeric(N) else field
  y <- array(NA_integer_, dim = c(N, J, S),
             dimnames = list(NULL, NULL, paste0("sp", seq_len(S))))
  for (s in seq_len(S)) {
    l  <- exp(as.numeric(X %*% bl[s, ]) + f + as.numeric(zeta %*% lam[s, ]))
    Ns <- stats::rpois(N, l)
    p  <- stats::plogis(bp[s, 1L])
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, Ns[i], p)
  }
  list(y = y, data = d, species = dimnames(y)[[3]], f = f,
       cor_res = stats::cov2cor(tcrossprod(lam) + diag(1e-8, S)),
       mu_lambda = colMeans(bl), mu_p = mean(bp), S = S, N = N)
}


test_that("abun() / ms_abun() reject a non-numeric K_max", {
  # `K_max` is the first formal, so abun("negbin") -- reaching for the mixing
  # distribution -- binds the string to K_max. Unchecked it coerces to NA and
  # resurfaces as an unrelated comparison error deep inside a kernel.
  expect_error(ms_abun("negbin"), "K_max.*mixture|mixture = ")
  expect_error(abun("negbin"),    "K_max.*mixture|mixture = ")
  expect_error(abun(K_max = -5),  "K_max")
  # the honest spellings still work
  expect_s3_class(ms_abun(mixture = "negbin"), "tobs_family")
  expect_s3_class(abun(K_max = 50), "tobs_family")
})

test_that("ms_abun() + latent() gates unsupported combinations", {
  d <- .msaf_sim(N = 40L, S = 6L, J = 3L, seed = 3L)

  # A negbin size is a second per-site dispersion: not identified against a
  # per-site latent structure.
  expect_error(
    tobs(~ x + latent(2), detection = ~ 1, data = d$data,
         family = ms_abun(mixture = "negbin"), y = d$y, species = d$species,
         method = "laplace"),
    "not identified|Poisson")
  # The sampler carries no latent term.
  expect_error(
    tobs(~ x + latent(2), detection = ~ 1, data = d$data,
         family = ms_abun(), y = d$y, species = d$species, method = "nuts"),
    "not yet wired|laplace")
  # A factor-only model is the plain block-coordinate Laplace-EM.
  expect_error(
    tobs(~ x + latent(2), detection = ~ 1, data = d$data,
         family = ms_abun(), y = d$y, species = d$species,
         method = "nested_laplace"),
    "laplace")
})

# Sized for runtime, not for the ceiling: N = 120 / S = 12 recovers the residual
# correlation at 0.994, but a fit is several minutes (each block-coordinate pass
# re-enters the community EM, and every likelihood evaluation sums over the
# latent N). N = 80 / S = 8 keeps the recovery well clear of the threshold.
test_that("a latent-factor community N-mixture recovers residual co-occurrence", {
  skip_if_fast()
  skip_on_cran()
  d <- .msaf_sim(N = 80L, S = 8L, J = 4L, Q = 2L, seed = 4L)
  fit <- tobs(~ x + latent(2), detection = ~ 1, data = d$data,
              family = ms_abun(), y = d$y, species = d$species,
              method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$ms_factor$n_factors, 2L)
  expect_equal(dim(fit$ms_factor$loadings), c(8L, 2L))
  expect_equal(dim(fit$ms_factor$residual_cov), c(8L, 8L))
  # residual species-correlation recovery (identified up to rotation)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.7)
  # community means recovered on both arms
  expect_equal(unname(coef(fit)$lambda), d$mu_lambda, tolerance = 0.3)
  expect_equal(unname(coef(fit)$p), d$mu_p, tolerance = 0.35)
  # per-species structure + S3
  expect_equal(dim(fit$ms_community$coef_lambda), c(8L, 2L))
  expect_true(all(is.finite(unlist(vcov(fit)))))
  rf <- ranef(fit)
  expect_true(all(c("species", "arm", "term", "estimate") %in% names(rf)))
  # fitted() is factor-aware: the factor offset enters lambda
  ft <- fitted(fit)
  expect_equal(dim(ft$lambda), c(80L, 8L))
  expect_gt(stats::cor(as.numeric(ft$lambda),
                       as.numeric(apply(d$y, c(1, 3), mean))), 0.7)
})

test_that("spatial-factor community N-mixture recovers the field and the factors", {
  skip_if_fast()
  skip_on_cran()
  side <- 10L
  A  <- .msaf_grid_graph(side)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- 0.6 * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f  <- f - mean(f)
  d  <- .msaf_sim(N = side * side, S = 10L, J = 4L, Q = 2L, field = f, seed = 6L)
  fit <- tobs(~ x + icar(graph = A) + latent(2), detection = ~ 1, data = d$data,
              family = ms_abun(), y = d$y, species = d$species,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  # both latents present and separated by the centred loadings
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$ms_factor))
  expect_gt(stats::cor(fit$spatial_field, d$f), 0.7)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off]), 0.7)
  # fitted() adds both offsets
  expect_gt(stats::cor(as.numeric(fitted(fit)$lambda),
                       as.numeric(apply(d$y, c(1, 3), mean))), 0.7)
})

test_that("latent-factor community N-mixture recovers over seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 6L
  rc <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d <- .msaf_sim(N = 80L, S = 8L, J = 4L, Q = 2L, seed = 300 + s)
    fit <- tobs(~ x + latent(2), detection = ~ 1, data = d$data,
                family = ms_abun(), y = d$y, species = d$species,
                method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    off <- upper.tri(d$cor_res)
    rc[s] <- stats::cor(fit$ms_factor$residual_cor[off], d$cor_res[off])
  }
  expect_gt(stats::median(rc), 0.8)
})
