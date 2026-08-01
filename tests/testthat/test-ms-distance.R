# Community binned distance sampling -- ms_distance() (the spAbundance msDS
# analogue), with latent() factors (lfMsDS) and a shared field (sfMsDS).
# gcol33/tulpaObs#117. Poisson.
#
#   N_{s,i} ~ Poisson(lambda_{s,i})
#   y_{s,i,.} | N ~ Multinomial(N_{s,i}; pi_{s,i,1..B}, 1 - p_det)
#   log lambda_{s,i} = X_i (mu_lambda + b_lambda_s) [+ f_i] [+ sum_q l_{s,q} z_{q,i}]
#   log sigma_{s,i}  = X_sig_i (mu_sigma + b_sigma_s)
#
# The latent N integrates out per species-site in closed form, so the family adds
# no C++: every fit is driven by the existing cpp_distance_site_sweep kernel
# (R/ms_distance.R). The community EM reads its per-species score from that
# sweep, and the latent driver's working oracle is the same Louis (1982) formula
# the community N-mixture uses:
#   score = grad_lam,  curv = info_lam - var_N * swl^2.
# Factors are identified only up to rotation, so factor recovery is judged on the
# residual species correlation (Sigma_res = lambda lambda'), which IS identified.

.msds_cut <- c(0, 25, 50, 75, 100)

.msds_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}


test_that("simulate_ms_distance() returns a well-formed community design", {
  d <- simulate_ms_distance(n_species = 5, N = 30, cutpoints = .msds_cut,
                            seed = 1)
  expect_equal(dim(d$y), c(30L, 4L, 5L))
  expect_length(d$species, 5L)
  expect_equal(nrow(d$data), 30L)
  expect_true(all(d$y >= 0L))
  expect_equal(dim(d$truth$beta_lambda), c(5L, 2L))
  # factors are optional truth
  df <- simulate_ms_distance(n_species = 6, N = 30, cutpoints = .msds_cut,
                             n_factors = 2, seed = 1)
  expect_equal(dim(df$truth$loadings), c(6L, 2L))
  expect_equal(dim(df$truth$cor_res), c(6L, 6L))
})

test_that("ms_distance() gates unsupported combinations", {
  d <- simulate_ms_distance(n_species = 4, N = 25, cutpoints = .msds_cut,
                            seed = 3)
  # cutpoints are required on the family
  expect_error(
    tobs(~ abund_cov1, detection = ~ 1, data = d$data, family = ms_distance(),
         y = d$y, species = d$species, method = "laplace"),
    "cutpoints")
  # negbin is not yet a per-species dispersion RE
  expect_error(
    tobs(~ abund_cov1, detection = ~ 1, data = d$data,
         family = ms_distance(cutpoints = .msds_cut, mixture = "negbin"),
         y = d$y, species = d$species, method = "laplace"),
    "Poisson")
  # a field needs nested_laplace; nested_laplace needs a field
  expect_error(
    tobs(~ abund_cov1, detection = ~ 1, data = d$data,
         family = ms_distance(cutpoints = .msds_cut), y = d$y,
         species = d$species, method = "nested_laplace"),
    "needs a shared field|laplace")
  # NUTS is not wired for this family
  expect_error(
    tobs(~ abund_cov1, detection = ~ 1, data = d$data,
         family = ms_distance(cutpoints = .msds_cut), y = d$y,
         species = d$species, method = "nuts"),
    "nuts|method")
})

# Cheap plumbing companion (gcol33/tulpaObs#159): a small design with no truth
# thresholds to re-calibrate, so it stays ungated and keeps this path exercised
# on every push while the recovery block below moves to the recovery tier.
test_that("a small msDS fit wires the community S3 surface", {
  d <- simulate_ms_distance(n_species = 3, N = 25, cutpoints = .msds_cut,
                            seed = 4)
  fit <- tobs(~ abund_cov1, detection = ~ 1, data = d$data,
              family = ms_distance(cutpoints = .msds_cut), y = d$y,
              species = d$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "laplace")
  expect_identical(fit$model$model_type, "ms_distance")
  expect_equal(dim(fit$ms_community$coef_lambda), c(3L, 2L))
  expect_equal(dim(fit$ms_community$Sigma_lambda), c(2L, 2L))
  expect_identical(rownames(fit$ms_community$coef_lambda), d$species)
  ft <- fitted(fit)
  expect_equal(dim(ft$lambda), c(25L, 3L))
  expect_true(is.finite(nobs(fit)))
})

test_that("msDS recovers the community means and per-species structure", {
  skip_if_fast()
  skip_on_cran()
  d <- simulate_ms_distance(n_species = 10, N = 100, cutpoints = .msds_cut,
                            seed = 4)
  fit <- tobs(~ abund_cov1, detection = ~ 1, data = d$data,
              family = ms_distance(cutpoints = .msds_cut), y = d$y,
              species = d$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$model$model_type, "ms_distance")
  # Community means against the seed's REALIZED mean, absolutely (#155). The old
  # `tolerance = 0.25` compared to mu_lambda = c(log(30), 0.4) = c(3.40, 0.4) and
  # mu_sigma = log(40) = 3.69: both intercepts exceed 0.25 in magnitude, so
  # `all.equal.numeric` treated the tolerance as RELATIVE there, silently
  # allowing +-0.85 / +-0.92 absolute -- effectively no check on either
  # intercept. Budget below is 3-4x the sd of the deviation from
  # colMeans(truth$beta_lambda) / colMeans(truth$beta_sigma) over seeds 1-16 of
  # this fixture (lambda intercept sd 0.029, slope sd 0.008, sigma intercept sd
  # 0.054), widened where needed to also cover this seed. Seed 4 itself sits at
  # the sample MAX on both intercepts (-0.103 lambda / +0.197 sigma) -- the
  # classic lambda/sigma trade-off direction the multi-seed test below already
  # documents (one realisation can move along that ridge without either
  # estimator arm being wrong), so the budget is set to cover it with margin
  # rather than only the typical seed.
  expect_community_mean(fit, c(colMeans(d$truth$beta_lambda),
                               colMeans(d$truth$beta_sigma)),
                        c(0.13, 0.025, 0.23))
  # per-species coefficients + community covariances
  expect_equal(dim(fit$ms_community$coef_lambda), c(10L, 2L))
  expect_equal(dim(fit$ms_community$Sigma_lambda), c(2L, 2L))
  expect_true(all(is.finite(unlist(vcov(fit)))))
  expect_true(all(is.finite(unlist(confint(fit)))))
  # per-species abundance surface tracks the truth
  ft <- fitted(fit)
  expect_equal(dim(ft$lambda), c(100L, 10L))
  expect_gt(stats::cor(as.numeric(ft$lambda), as.numeric(d$truth$lambda)), 0.8)
  # S3
  rf <- ranef(fit)
  expect_true(all(c("species", "arm", "term", "estimate") %in% names(rf)))
  expect_true(all(c("lambda", "sigma") %in% unique(rf$arm)))
  expect_true(is.finite(nobs(fit)))
})

# A single fit cannot separate bias from a draw: seed 4 alone lands mu_lambda
# 0.18 low and mu_sigma 0.17 high (the classic lambda / sigma trade-off
# direction), which looks like bias but is one realisation. Averaged over seeds
# both community means sit on the truth. Note the estimand: each species' beta is
# drawn around the community mean, so even a perfect fit scatters around
# mu_lambda with SD ~ sd_lambda / sqrt(n_species).
test_that("msDS community means are unbiased over seeds with nominal coverage", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 10L
  est <- matrix(NA_real_, n_seed, 3L)
  hit <- matrix(NA, n_seed, 3L)
  truth <- NULL
  for (s in seq_len(n_seed)) {
    d <- simulate_ms_distance(n_species = 10, N = 100, cutpoints = .msds_cut,
                              seed = 100 + s)
    truth <- c(d$truth$mu_lambda, d$truth$mu_sigma)
    fit <- tobs(~ abund_cov1, detection = ~ 1, data = d$data,
                family = ms_distance(cutpoints = .msds_cut), y = d$y,
                species = d$species, method = "laplace",
                control = list(verbose = FALSE, progress = FALSE))
    e <- c(coef(fit)$lambda, coef(fit)$sigma)
    est[s, ] <- e
    hit[s, ] <- abs(e - truth) < 1.96 * fit$sds
  }
  # Bias of the community means, judged against the Monte Carlo error of the
  # seed average rather than an absolute tolerance.
  bias <- colMeans(est) - truth
  mcse <- apply(est, 2, stats::sd) / sqrt(n_seed)
  expect_lt(max(abs(bias / mcse)), 4)
  # 95% Wald coverage, on the rubric's pooled floor (not a literal 0.95, which
  # is flaky at 10 seeds).
  expect_gt(mean(hit), 0.8)
})

# Sized for runtime: each block-coordinate pass re-enters the community EM, whose
# per-species Hessian is finite-differenced here (the distance kernel exposes the
# per-site Louis pieces but no assembled block, so unlike ms_abun this family
# does not yet pass `sp_info`).
test_that("lfMsDS recovers residual species co-occurrence", {
  skip_if_fast()
  skip_on_cran()
  d <- simulate_ms_distance(n_species = 8, N = 80, cutpoints = .msds_cut,
                            n_factors = 2, load_sd = 0.5, seed = 5)
  fit <- tobs(~ abund_cov1 + latent(2), detection = ~ 1, data = d$data,
              family = ms_distance(cutpoints = .msds_cut), y = d$y,
              species = d$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$ms_factor$n_factors, 2L)
  expect_equal(dim(fit$ms_factor$loadings), c(8L, 2L))
  off <- upper.tri(d$truth$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off],
                       d$truth$cor_res[off]), 0.7)
  expect_equal(unname(coef(fit)$lambda), d$truth$mu_lambda, tolerance = 0.3)
})

test_that("sfMsDS recovers the shared field alongside the factors", {
  skip_if_fast()
  skip_on_cran()
  side <- 9L
  A  <- .msds_grid_graph(side)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- 0.5 * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f  <- f - mean(f)
  d  <- simulate_ms_distance(n_species = 8, cutpoints = .msds_cut,
                             n_factors = 2, field = f, seed = 6)
  fit <- tobs(~ abund_cov1 + icar(graph = A) + latent(2), detection = ~ 1,
              data = d$data, family = ms_distance(cutpoints = .msds_cut),
              y = d$y, species = d$species, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$ms_factor))
  expect_gt(stats::cor(fit$spatial_field, d$truth$field), 0.7)
  off <- upper.tri(d$truth$cor_res)
  expect_gt(stats::cor(fit$ms_factor$residual_cor[off],
                       d$truth$cor_res[off]), 0.6)
})
