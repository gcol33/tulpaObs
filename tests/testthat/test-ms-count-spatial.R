# Community-spatial count -- ms_count() + a shared areal field icar() (the
# spAbundance sfMsAbund analogue, Poisson; gcol33/tulpaObs#117). Fit by block
# coordinate ascent: the community Laplace-EM (shared field as a per-site offset)
# alternated with a self-contained Poisson-ICAR field update (R/ms_count_spatial.R),
# no C++. The shared field is well identified (informed by every species at each
# site), so it recovers cleanly alongside the community means.

.msc_grid_graph <- function(side) {
  N <- side * side
  A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# Community-spatial Poisson counts: log mu_{s,i} = X_i (mu + b_s) + f_i.
.msc_sim <- function(side = 10L, S = 12L, beta = c(1, 0.5),
                     beta_sd = c(0.4, 0.3), field_sd = 0.7, seed = 1L) {
  set.seed(seed)
  A <- .msc_grid_graph(side); Ns <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f <- f - mean(f)
  d <- data.frame(x = stats::rnorm(Ns))
  X <- stats::model.matrix(~ x, d)
  bs <- vapply(seq_along(beta),
               function(j) stats::rnorm(S, beta[j], beta_sd[j]), numeric(S))
  y <- matrix(NA_real_, Ns, S, dimnames = list(NULL, paste0("sp", seq_len(S))))
  for (s in seq_len(S)) y[, s] <- stats::rpois(Ns, exp(as.numeric(X %*% bs[s, ]) + f))
  list(y = y, data = d, graph = A, field = f, beta = beta, Ns = Ns, S = S)
}


test_that("ms_count() reports nested_laplace and gates the areal field", {
  expect_true("nested_laplace" %in% tulpaObs:::.tobs_family_methods$ms_count)
  d <- .msc_sim(side = 6L, S = 6L, seed = 3L)

  # non-spatial + nested_laplace, and areal + laplace, both error
  expect_error(
    tobs(~ x, data = d$data, family = ms_count(), y = d$y,
         species = colnames(d$y), method = "nested_laplace"),
    "needs a shared areal field")
  expect_error(
    tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_count(),
         y = d$y, species = colnames(d$y), method = "laplace"),
    "needs method = .nested_laplace.")
  # negbin + areal is gated (Poisson-only: field vs dispersion not identified)
  expect_error(
    tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_count("negbin"),
         y = d$y, species = colnames(d$y), method = "nested_laplace"),
    "Poisson-only|not identified")
})

test_that("a community-spatial count fit recovers the field and wires S3", {
  skip_on_cran()
  d   <- .msc_sim(side = 10L, S = 12L, seed = 5L)
  fit <- tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_count(),
              y = d$y, species = colnames(d$y), method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_equal(length(fit$spatial_field), d$Ns)
  expect_gt(stats::cor(fit$spatial_field, d$field), 0.8)
  expect_true(is.finite(fit$spatial_hyper$sigma) && fit$spatial_hyper$sigma > 0)

  expect_length(unlist(coef(fit)), 2L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  # fitted() is field-aware -> tracks the counts
  ft <- fitted(fit)$mu
  expect_equal(dim(ft), c(d$Ns, d$S))
  expect_gt(stats::cor(as.numeric(ft), as.numeric(d$y)), 0.6)
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("community-spatial count recovers community means + field over seeds", {
  skip_if_fast()
  skip_on_cran()
  beta <- c(1, 0.5)
  n_seed <- 20L
  cover <- matrix(FALSE, n_seed, length(beta))
  est   <- matrix(NA_real_, n_seed, length(beta))
  fcor  <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d   <- .msc_sim(side = 10L, S = 12L, beta = beta, seed = s)
    fit <- tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_count(),
                y = d$y, species = colnames(d$y), method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
    fcor[s]    <- stats::cor(fit$spatial_field, d$field)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.06)
  expect_gte(mean(cover), 0.85)
  expect_gt(stats::median(fcor), 0.9)
})


# --- community SVC (svcMsAbund): intercept + varying-coefficient field ---------

# log mu_{s,i} = X_i (mu + b_s) + f0_i + w_i * f1_i, two independent ICAR fields.
.msc_svc_sim <- function(side = 10L, S = 14L, seed = 1L) {
  set.seed(seed)
  A <- .msc_grid_graph(side); Ns <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f0 <- 0.6 * scale(sin(co$r/side*pi) + cos(co$c/side*pi))[, 1]; f0 <- f0 - mean(f0)
  f1 <- 0.6 * scale(cos(co$r/side*pi*1.3) + sin(co$c/side*pi*0.7))[, 1]; f1 <- f1 - mean(f1)
  d <- data.frame(x = stats::rnorm(Ns), w = stats::rnorm(Ns), cell = seq_len(Ns))
  X <- stats::model.matrix(~ x, d)
  bs <- vapply(1:2, function(j) stats::rnorm(S, c(1, 0.5)[j], c(0.4, 0.3)[j]),
               numeric(S))
  y <- matrix(NA_real_, Ns, S, dimnames = list(NULL, paste0("sp", seq_len(S))))
  for (s in seq_len(S))
    y[, s] <- stats::rpois(Ns, exp(as.numeric(X %*% bs[s, ]) + f0 + d$w * f1))
  list(y = y, data = d, graph = A, f0 = f0, f1 = f1, Ns = Ns, S = S)
}

test_that("community field recovers under bym2 (scaled structured + iid)", {
  skip_on_cran()
  d <- .msc_sim(side = 9L, S = 12L, seed = 5L)
  fit <- tobs(~ x + bym2(graph = d$graph), data = d$data,
              family = ms_count(), y = d$y, species = colnames(d$y),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$spatial_hyper$type, "bym2")
  expect_true(fit$spatial_hyper$rho > 0 && fit$spatial_hyper$rho <= 1)
  expect_true(fit$spatial_hyper$sigma > 0)
  # the combined field recovers a smooth structured truth
  expect_gt(stats::cor(fit$spatial_field, d$field), 0.8)
  expect_equal(unname(unlist(coef(fit))), c(1, 0.5), tolerance = 0.15)
  # bym2 is the single shared intercept field only: an SVC bar errors
  expect_error(
    tobs(~ x + bym2(graph = d$graph, weight = x, group_var = "cell"),
         data = cbind(d$data, cell = seq_len(d$Ns)), family = ms_count(),
         y = d$y, species = colnames(d$y), method = "nested_laplace",
         control = list(progress = FALSE)),
    "single shared intercept|bym2")
})

test_that("community field recovers under group_var (sites > cells)", {
  skip_on_cran()
  set.seed(8)
  side <- 8L; Acell <- .msc_grid_graph(side); Ncell <- nrow(Acell)
  R <- 2L; Ns <- Ncell * R; S <- 12L
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  fcell <- 0.7 * scale(sin(co$r/side*pi) + cos(co$c/side*pi))[, 1]; fcell <- fcell - mean(fcell)
  cell_of_site <- rep(seq_len(Ncell), each = R)
  d <- data.frame(x = stats::rnorm(Ns), cell = cell_of_site)
  X <- stats::model.matrix(~ x, d)
  bs <- vapply(1:2, function(j) stats::rnorm(S, c(1, 0.5)[j], c(0.4, 0.3)[j]), numeric(S))
  y <- matrix(NA_real_, Ns, S, dimnames = list(NULL, paste0("sp", seq_len(S))))
  for (s in seq_len(S))
    y[, s] <- stats::rpois(Ns, exp(as.numeric(X %*% bs[s, ]) + fcell[cell_of_site]))
  fit <- tobs(~ x + icar(graph = Acell, group_var = "cell"), data = d,
              family = ms_count(), y = y, species = colnames(y),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  # the field has one node per CELL (fewer than sites) and recovers the cell field
  expect_equal(length(fit$spatial_field), Ncell)
  expect_gt(stats::cor(fit$spatial_field, fcell), 0.85)
  expect_equal(unname(unlist(coef(fit))), c(1, 0.5), tolerance = 0.15)
})

test_that("community field recovers under car_proper (proper CAR)", {
  skip_on_cran()
  d <- .msc_sim(side = 10L, S = 12L, seed = 5L)
  fit <- tobs(~ x + car_proper(graph = d$graph), data = d$data,
              family = ms_count(), y = d$y, species = colnames(d$y),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_identical(fit$spatial_hyper$type, "car_proper")
  expect_true(fit$spatial_hyper$rho > 0 && fit$spatial_hyper$rho < 1)
  expect_gt(stats::cor(fit$spatial_field, d$field), 0.8)
  # the intercept stays identified (the field is sum-to-zero deviations)
  expect_equal(unname(unlist(coef(fit))), c(1, 0.5), tolerance = 0.15)
})

test_that("community SVC (svcMsAbund) recovers the intercept + trend fields", {
  skip_on_cran()
  d <- .msc_svc_sim(side = 10L, S = 14L, seed = 7L)
  fit <- tobs(~ x + spatial(~ 1 + w || cell, graph = d$graph), data = d$data,
              family = ms_count(), y = d$y, species = colnames(d$y),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  # both fields recover
  expect_gt(stats::cor(fit$spatial_field, d$f0), 0.8)    # intercept field
  expect_gt(stats::cor(fit$trend_field,   d$f1), 0.8)    # varying-coefficient field
  expect_equal(fit$spatial_hyper$field_labels, c("intercept", "w"))
  expect_true(is.finite(tobs_waic(fit)$waic))
})

test_that("community SVC recovers both fields across seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 12L
  c0 <- numeric(n_seed); c1 <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d <- .msc_svc_sim(side = 10L, S = 14L, seed = 100 + s)
    fit <- tobs(~ x + spatial(~ 1 + w || cell, graph = d$graph), data = d$data,
                family = ms_count(), y = d$y, species = colnames(d$y),
                method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    c0[s] <- stats::cor(fit$spatial_field, d$f0)
    c1[s] <- stats::cor(fit$trend_field,   d$f1)
  }
  expect_gt(stats::median(c0), 0.85)
  expect_gt(stats::median(c1), 0.85)
})
