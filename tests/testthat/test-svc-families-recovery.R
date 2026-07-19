# Continuous NNGP varying coefficients (svc()) on the observation families'
# Laplace backends (gcol33/tulpaObs#144).
#
# svc() surfaces are latent field blocks on the state arm, so every family whose
# marginal hands the shared areal-BFGS driver a per-site eta gradient carries them
# with no family-specific code: `.tobs_build_field_spec()` appends one NNGP block
# per `indices` entry and the family's existing eval() is untouched. These tests
# fit a KNOWN smooth surface and assert it is recovered.
#
# The recovery target is the SURFACE (correlation against truth) plus the
# coefficients, not the hyperparameters: the marginal SD and range sit on the GP
# ridge and trade off against each other at these site counts (the same reason
# gcol33/tulpaObs#119 asserts divergences / phi / sigma on the NUTS route rather
# than surface correlation -- each backend asserts what its own bottleneck leaves
# identified). Count families are far more informative per site than binary
# occupancy ones, so their tolerances differ by design.

# One smooth coordinate surface, shared by every fit below.
.svc_fam_truth <- function(n, seed) {
  set.seed(seed)
  co <- cbind(stats::runif(n), stats::runif(n))
  z <- 1.2 * sin(3 * co[, 1]) + 1.2 * cos(3 * co[, 2])
  list(co = co, z = z - mean(z),
       df = data.frame(lon = co[, 1], lat = co[, 2], x = stats::rnorm(n)))
}

.svc_fam_formula <- function(covariate = TRUE) {
  if (covariate)
    ~ 1 + x + svc(lon, lat, indices = 1, nn = 10, prior_range = c(0.3, 0.5))
  else
    ~ 1 + svc(lon, lat, indices = 1, nn = 10, prior_range = c(0.3, 0.5))
}


test_that("removal() recovers a known svc() surface on the abundance arm", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(120L, 1L)
  N  <- nrow(tr$df)
  lambda <- exp(1.6 + 0.5 * tr$df$x + tr$z)
  Ntrue  <- stats::rpois(N, lambda)
  K <- 3L; p_cap <- 0.45
  y <- matrix(0L, N, K); remaining <- Ntrue
  for (k in seq_len(K)) {
    y[, k] <- stats::rbinom(N, remaining, p_cap)
    remaining <- remaining - y[, k]
  }
  fit <- tobs(.svc_fam_formula(), data = tr$df, y = y, family = removal(),
              detection = ~ 1, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))

  surf <- as.numeric(fit$svc_field)
  expect_length(surf, N)
  expect_true(all(is.finite(surf)))
  # Counts identify the surface sharply.
  expect_gt(stats::cor(surf, tr$z), 0.85)
  cf <- coef(fit)
  expect_lt(abs(cf$lambda[["(Intercept)"]] - 1.6), 0.5)
  expect_lt(abs(cf$lambda[["x"]] - 0.5), 0.25)
  expect_lt(abs(stats::plogis(cf$p[["(Intercept)"]]) - p_cap), 0.12)

  # The surface is reported the way the occu() Laplace and NUTS routes report it.
  expect_s3_class(fit$svc, "tobs_svc")
  expect_identical(fit$svc_field_arm, "abundance")
  expect_true(all(c("sigma", "phi") %in% names(fit$svc_hyper[[1L]])))
  # No areal term was supplied, so no areal field is invented alongside it.
  expect_null(fit$spatial_field)
})


test_that("removal() recovers a covariate-weighted svc() surface", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(120L, 4L)
  N <- nrow(tr$df)
  # The surface multiplies x: a spatially-varying slope, not a varying intercept.
  lambda <- exp(1.7 + tr$z * tr$df$x)
  Ntrue  <- stats::rpois(N, lambda)
  K <- 3L; p_cap <- 0.45
  y <- matrix(0L, N, K); remaining <- Ntrue
  for (k in seq_len(K)) {
    y[, k] <- stats::rbinom(N, remaining, p_cap)
    remaining <- remaining - y[, k]
  }
  fit <- tobs(~ 1 + x + svc(lon, lat, indices = 2, nn = 10,
                            prior_range = c(0.3, 0.5)),
              data = tr$df, y = y, family = removal(), detection = ~ 1,
              method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  surf <- as.numeric(fit$svc_field)
  expect_length(surf, N)
  # A weighted surface is identified only where |x| carries signal, so it recovers
  # less sharply than the intercept surface above.
  expect_gt(stats::cor(surf, tr$z), 0.5)
  expect_identical(fit$svc_indices, 2L)
})


test_that("fp_occu() recovers a known svc() surface on the psi arm", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(120L, 1L)
  N <- nrow(tr$df); J <- 6L
  psi <- stats::plogis(-0.2 + 0.8 * tr$df$x + tr$z)
  z   <- stats::rbinom(N, 1, psi)
  y   <- matrix(0L, N, J)
  for (i in seq_len(N)) for (j in seq_len(J)) {
    y[i, j] <- if (z[i] == 1L) sample(0:2, 1L, prob = c(0.4, 0.3, 0.3))
               else sample(0:1, 1L, prob = c(0.95, 0.05))
  }
  fit <- tobs(.svc_fam_formula(), data = tr$df, y = y, family = fp_occu(),
              detection = ~ 1, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))

  surf <- as.numeric(fit$svc_field)
  expect_length(surf, N)
  expect_true(all(is.finite(surf)))
  # One binary site per surface node: information-bounded, not backend-bounded --
  # the same regime the single-season occu() Laplace route measures (#143).
  expect_gt(stats::cor(surf, tr$z), 0.55)
  cf <- coef(fit)
  expect_lt(abs(stats::plogis(cf$p11[["(Intercept)"]]) - 0.6), 0.15)
  expect_lt(stats::plogis(cf$p10[["(Intercept)"]]), 0.2)
  expect_identical(fit$svc_field_arm, "occupancy")
})


test_that("distance() recovers a known svc() surface on the abundance arm", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(60L, 2L)
  N <- nrow(tr$df)
  cutpoints <- c(0, 25, 50, 75, 100)
  sigma <- 50
  lambda <- exp(1.8 + tr$z)
  Ntrue  <- stats::rpois(N, lambda)
  # Half-normal detection: per-bin masses over the covered strip.
  pb <- vapply(seq_len(4L), function(b)
    stats::integrate(function(u) exp(-u^2 / (2 * sigma^2)),
                     cutpoints[b], cutpoints[b + 1L])$value / 100, 0)
  y <- t(vapply(seq_len(N), function(i)
    stats::rmultinom(1L, Ntrue[i], c(pb, 1 - sum(pb)))[1:4, 1], numeric(4L)))

  fit <- tobs(.svc_fam_formula(covariate = FALSE), data = tr$df, y = y,
              family = distance(cutpoints = cutpoints, transect = "line"),
              detection = ~ 1, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  surf <- as.numeric(fit$svc_field)
  expect_length(surf, N)
  expect_gt(stats::cor(surf, tr$z), 0.8)
  cf <- coef(fit)
  expect_lt(abs(cf$lambda[["(Intercept)"]] - 1.8), 0.5)
  expect_lt(abs(exp(cf$sigma[["(Intercept)"]]) - sigma) / sigma, 0.3)
  expect_identical(fit$svc_field_arm, "abundance")
})


test_that("dyn_abun() recovers a known svc() surface on the initial-abundance arm", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(45L, 2L)
  N <- nrow(tr$df); n_seasons <- 3L; J <- 3L
  omega <- 0.7; gamma <- 1.2; p_det <- 0.5
  Nmat <- matrix(0L, N, n_seasons)
  Nmat[, 1L] <- stats::rpois(N, exp(1.2 + tr$z))
  for (t in 2:n_seasons)
    Nmat[, t] <- stats::rbinom(N, Nmat[, t - 1L], omega) + stats::rpois(N, gamma)
  y <- array(0L, dim = c(N, J, n_seasons))
  for (t in seq_len(n_seasons)) for (j in seq_len(J))
    y[, j, t] <- stats::rbinom(N, Nmat[, t], p_det)

  fit <- tobs(.svc_fam_formula(covariate = FALSE), data = tr$df, y = y,
              family = dyn_abun(K_max = 25L), detection = ~ 1,
              omega = ~ 1, gamma = ~ 1, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  surf <- as.numeric(fit$svc_field)
  expect_length(surf, N)
  expect_gt(stats::cor(surf, tr$z), 0.8)
  cf <- coef(fit)
  expect_lt(abs(cf$lambda[["(Intercept)"]] - 1.2), 0.5)
  expect_lt(abs(stats::plogis(cf$p[["(Intercept)"]]) - p_det), 0.15)
  expect_lt(abs(stats::plogis(cf$omega[["(Intercept)"]]) - omega), 0.2)
  expect_identical(fit$svc_field_arm, "abundance")
})


test_that("an svc() surface composes with an areal field on the same arm", {
  skip_on_cran(); skip_if_fast()
  tr <- .svc_fam_truth(100L, 7L)
  N <- nrow(tr$df)
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  tr$df$cell <- seq_len(N)
  lambda <- exp(1.6 + tr$z)
  Ntrue  <- stats::rpois(N, lambda)
  K <- 3L; y <- matrix(0L, N, K); remaining <- Ntrue
  for (k in seq_len(K)) {
    y[, k] <- stats::rbinom(N, remaining, 0.45)
    remaining <- remaining - y[, k]
  }
  fit <- tobs(~ 1 + icar(graph = adj, group_var = "cell") +
                svc(lon, lat, indices = 1, nn = 10, prior_range = c(0.3, 0.5)),
              data = tr$df, y = y, family = removal(), detection = ~ 1,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  # Both blocks are fit and reported under their own slots; the areal field and
  # the continuous surface share the abundance arm and split its structure, so
  # neither is asserted against the truth here -- only that both are present and
  # finite, and that the coefficients survive the composition.
  expect_length(as.numeric(fit$svc_field), N)
  expect_length(fit$spatial_field, N)
  expect_true(all(is.finite(as.numeric(fit$svc_field))))
  expect_true(all(is.finite(fit$spatial_field)))
  expect_lt(abs(coef(fit)$lambda[["(Intercept)"]] - 1.6), 0.6)
})


test_that("svc() on the detection formula errors on the areal-BFGS families", {
  set.seed(11)
  n <- 30L
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n))
  y  <- matrix(stats::rpois(n * 3L, 2), n, 3L)
  expect_error(
    tobs(~ 1, data = df, y = y, family = removal(),
         detection = ~ svc(lon, lat, indices = 1, prior_range = c(0.3, 0.5)),
         method = "laplace"),
    "detection|abundance")
})


test_that("svc() under NUTS still errors on the areal-BFGS families", {
  set.seed(12)
  n <- 30L
  df <- data.frame(lon = stats::runif(n), lat = stats::runif(n))
  y  <- matrix(stats::rpois(n * 3L, 2), n, 3L)
  # The in-tree NUTS targets carry no NNGP block; the guard fires before any fit.
  expect_error(
    tobs(~ svc(lon, lat, indices = 1, prior_range = c(0.3, 0.5)), data = df,
         y = y, family = removal(), detection = ~ 1, method = "nuts"),
    "svc|areal|spatially-varying")
})
