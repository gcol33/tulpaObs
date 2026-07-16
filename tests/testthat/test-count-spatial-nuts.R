# NUTS + areal car_proper field on the count / occupancy observation families
# (gcol33/tulpaObs#72). Each family carries a FIXED-HYPER non-centered proper-CAR
# field on its abundance / occupancy arm: the field precision (tau, rho) is fixed
# at the nested-Laplace areal posterior mean and the whitened raw ~ N(0, I) is
# sampled jointly with the coefficients via the shared field block
# (src/nuts_field_block.h). The recovery invariant is that NUTS reproduces the
# integrated nested-Laplace field (cor high), recovers the abundance / occupancy
# slope, calibrates the coefficient SDs to the nested-Laplace SEs (fixed-hyper,
# same field precision), and samples without divergences.
#
# The dispatcher gate for the four observation families (removal / distance /
# fp_occu / dyn_abun) lives in R/occu_fit.R (NUTS + spatial -> stop), so these
# tests drive the family spatial-NUTS fitters directly; the one-line dispatch
# rewiring is the occu_fit.R follow-up.

.csn_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

# Smoothed, demeaned ICAR-like field on the grid.
.csn_field <- function(adj, sd_phi = 0.6, seed = 1L) {
  set.seed(seed); ng <- nrow(adj); phi <- as.numeric(scale(rnorm(ng)))
  for (r in 1:3) {
    pn <- phi
    for (i in seq_len(ng)) { nb <- which(adj[i, ] == 1L); pn[i] <- 0.5 * phi[i] + 0.5 * mean(phi[nb]) }
    phi <- pn
  }
  phi <- sd_phi * as.numeric(scale(phi)); phi - mean(phi)
}

# Resolve the fully-built spatial term from a fitted model (the spatial-NUTS
# fitters take the tobs_spatial object the dispatcher would have passed them).
.csn_spatial <- function(fit) tulpaObs:::.tobs_structures_from_model(fit$model)$spatial

# SE of an arm coefficient from the nested-Laplace grid-integrated vcov.
.csn_nl_se <- function(fit, nm) sqrt(diag(fit$vcov))[[nm]]

# A covariate orthogonal to the latent field, so the abundance / occupancy slope
# is identified separately from the field (a covariate that happens to correlate
# with the smoothed field would confound the slope, defeating slope recovery).
.csn_orth_cov <- function(phi) {
  x <- rnorm(length(phi))
  as.numeric(scale(residuals(lm(x ~ phi))))
}

# Hazard-rate detection g(x) = 1 - exp(-(x/sigma)^(-shape)); per-bin detection
# probability via the bin-midpoint x line-transect uniform density (a rough
# generator, enough to identify the abundance field + the scalar shape).
.csn_haz_pi <- function(cut, sigma, shape, W) {
  mid <- (head(cut, -1) + tail(cut, -1)) / 2
  (1 - exp(-(mid / sigma)^(-shape))) * diff(cut) / W
}


test_that("removal() NUTS + car_proper field recovers + calibrates to nested-Laplace", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 11L); ng <- nrow(adj)
  set.seed(11); xab <- .csn_orth_cov(phi)
  lambda <- exp(log(8) + 0.6 * xab + phi); p <- 0.45; K <- 4L
  N <- rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- N
  for (k in 1:K) { y[, k] <- rbinom(ng, rem, p); rem <- rem - y[, k] }
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), data = data.frame(abund_cov1 = xab),
             family = removal(), detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_removal_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 3L)
  expect_identical(nu$method, "nuts")
  expect_false(is.null(nu$spatial_field))
  expect_equal(mean(nu$divergent), 0)                                   # clean geometry
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.95)             # reproduces NL field
  expect_gt(cor(nu$spatial_field, phi), 0.6)                           # tracks truth
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.6) /
              nu$sds[["lambda_abund_cov1"]], 3)                        # slope recovered
  rr <- nu$sds[["lambda_abund_cov1"]] / .csn_nl_se(nl, "lambda_abund_cov1")
  expect_gt(rr, 0.6); expect_lt(rr, 1.6)                               # SD calibrated to NL
})


test_that("distance() NUTS + car_proper field recovers + calibrates to nested-Laplace", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 22L); ng <- nrow(adj)
  set.seed(22); xab <- .csn_orth_cov(phi); cut <- c(0, 10, 20, 30, 40); B <- 4L
  sig <- exp(3.0); lamd <- exp(log(40) + 0.5 * xab + phi)
  g_mid <- exp(-((head(cut, -1) + tail(cut, -1)) / 2)^2 / (2 * sig^2))
  pi_b <- g_mid * diff(cut) / 40 * 0.6
  Nd <- rpois(ng, lamd); y <- matrix(0L, ng, B)
  for (i in seq_len(ng)) y[i, ] <- rbinom(B, Nd[i], pi_b)
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), data = data.frame(abund_cov1 = xab),
             family = distance(key = "halfnorm", transect = "line", cutpoints = cut),
             detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_distance_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 4L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.95)
  expect_gt(cor(nu$spatial_field, phi), 0.6)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
              nu$sds[["lambda_abund_cov1"]], 3)
  rr <- nu$sds[["lambda_abund_cov1"]] / .csn_nl_se(nl, "lambda_abund_cov1")
  expect_gt(rr, 0.6); expect_lt(rr, 1.6)
})


test_that("fp_occu() NUTS + car_proper field reproduces the nested-Laplace field, 0 divergences", {
  skip_on_cran()
  skip_if_fast()
  # The occupancy field is weakly identified (one binary site per node), so the
  # invariant is that NUTS reproduces the integrated nested-Laplace field and
  # samples cleanly -- not that the field tracks the latent truth tightly (the
  # nested-Laplace field itself only weakly does at this replication).
  adj <- .csn_grid_adj(7L); phi <- .csn_field(adj, sd_phi = 0.8, seed = 33L); ng <- nrow(adj)
  set.seed(33); xab <- .csn_orth_cov(phi); psi <- plogis(0.3 + 0.7 * xab + phi)
  z <- rbinom(ng, 1, psi); J <- 6L; yf <- integer(ng * J); idx <- 1L
  for (i in seq_len(ng)) for (j in 1:J) {
    yf[idx] <- if (z[i] == 1) sample(0:2, 1, prob = c(0.45, 0.25, 0.30))
               else sample(0:1, 1, prob = c(0.92, 0.08))
    idx <- idx + 1L
  }
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), data = data.frame(abund_cov1 = xab),
             family = fp_occu(), detection = ~ 1, y = matrix(yf, ng, J, byrow = TRUE),
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_fp_occu_nuts_spatial(
    nl$model, .csn_spatial(nl), n.iter = 600L, n.warmup = 600L,
    adapt.delta = 0.95, seed = 5L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)             # reproduces NL field
  # NUTS reproduces the nested-Laplace occupancy slope (the field is weakly
  # identified, so the invariant is NUTS == nested-Laplace, not truth recovery).
  expect_lt(abs(nu$means[["psi_abund_cov1"]] - nl$means[["psi_abund_cov1"]]) /
              nu$sds[["psi_abund_cov1"]], 1)
})


test_that("dyn_abun() NUTS + car_proper field reproduces the nested-Laplace field, 0 divergences", {
  skip_on_cran()
  skip_if_fast()
  # The forward-HMM NUTS is the slowest per leapfrog, so a small grid with a
  # bounded tree depth keeps the recovery affordable while still exercising the
  # fixed-hyper field block on the initial-abundance arm.
  adj <- .csn_grid_adj(4L); phi <- .csn_field(adj, seed = 44L); ng <- nrow(adj)
  set.seed(44); xab <- .csn_orth_cov(phi); Td <- 3L; Jd <- 3L
  N1 <- rpois(ng, exp(log(7) + 0.5 * xab + phi)); ya <- array(0L, c(ng, Jd, Td)); Ncur <- N1
  for (t in 1:Td) {
    for (j in 1:Jd) ya[, j, t] <- rbinom(ng, Ncur, 0.5)
    if (t < Td) Ncur <- rbinom(ng, Ncur, 0.6) + rpois(ng, 2)
  }
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), data = data.frame(abund_cov1 = xab),
             family = dyn_abun(K_max = 20L), detection = ~ 1, y = ya,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_dyn_abun_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 200L, n.warmup = 200L, max.treedepth = 8L, seed = 6L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)
  # NUTS reproduces the nested-Laplace initial-abundance slope (the field-arm
  # confounding leaves the truth weakly identified at this grid / replication).
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - nl$means[["lambda_abund_cov1"]]) /
              nu$sds[["lambda_abund_cov1"]], 1)
})


# ---------------------------------------------------------------------------
# Intrinsic ICAR fields under NUTS (gcol33/tulpaObs#113). The #71 sum-to-zero
# reparameterisation drops the constant precision null direction, so the icar
# field's whitened raw ~ N(0, I_{n-1}) samples with the same well-conditioned
# geometry as the full-rank proper-CAR field (0 divergences, reproduces the
# nested-Laplace icar field). Same invariants as the car_proper tests above.
# ---------------------------------------------------------------------------

test_that("removal() NUTS + icar field recovers + reproduces the nested-Laplace field (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 11L); ng <- nrow(adj)
  set.seed(11); xab <- .csn_orth_cov(phi)
  lambda <- exp(log(8) + 0.6 * xab + phi); p <- 0.45; K <- 4L
  N <- rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- N
  for (k in 1:K) { y[, k] <- rbinom(ng, rem, p); rem <- rem - y[, k] }
  nl <- tobs(~ abund_cov1 + icar(graph = adj), data = data.frame(abund_cov1 = xab),
             family = removal(), detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_removal_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 3L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)                                  # sum-to-zero geometry
  expect_length(nu$spatial_field, ng)                                  # centred field, length n
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)                         # sum z = 0
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.9)             # reproduces NL field
  expect_gt(cor(nu$spatial_field, phi), 0.6)                          # tracks truth
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.6) /
              nu$sds[["lambda_abund_cov1"]], 3)                       # slope recovered
})


test_that("distance() NUTS + icar field recovers + reproduces the nested-Laplace field (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 22L); ng <- nrow(adj)
  set.seed(22); xab <- .csn_orth_cov(phi); cut <- c(0, 10, 20, 30, 40); B <- 4L
  sig <- exp(3.0); lamd <- exp(log(40) + 0.5 * xab + phi)
  g_mid <- exp(-((head(cut, -1) + tail(cut, -1)) / 2)^2 / (2 * sig^2))
  pi_b <- g_mid * diff(cut) / 40 * 0.6
  Nd <- rpois(ng, lamd); y <- matrix(0L, ng, B)
  for (i in seq_len(ng)) y[i, ] <- rbinom(B, Nd[i], pi_b)
  nl <- tobs(~ abund_cov1 + icar(graph = adj), data = data.frame(abund_cov1 = xab),
             family = distance(key = "halfnorm", transect = "line", cutpoints = cut),
             detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_distance_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 4L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.9)
  expect_gt(cor(nu$spatial_field, phi), 0.6)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
              nu$sds[["lambda_abund_cov1"]], 3)
})


test_that("fp_occu() NUTS + icar field reproduces the nested-Laplace field, 0 divergences (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(7L); phi <- .csn_field(adj, sd_phi = 0.8, seed = 33L); ng <- nrow(adj)
  set.seed(33); xab <- .csn_orth_cov(phi); psi <- plogis(0.3 + 0.7 * xab + phi)
  z <- rbinom(ng, 1, psi); J <- 6L; yf <- integer(ng * J); idx <- 1L
  for (i in seq_len(ng)) for (j in 1:J) {
    yf[idx] <- if (z[i] == 1) sample(0:2, 1, prob = c(0.45, 0.25, 0.30))
               else sample(0:1, 1, prob = c(0.92, 0.08))
    idx <- idx + 1L
  }
  nl <- tobs(~ abund_cov1 + icar(graph = adj), data = data.frame(abund_cov1 = xab),
             family = fp_occu(), detection = ~ 1, y = matrix(yf, ng, J, byrow = TRUE),
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_fp_occu_nuts_spatial(
    nl$model, .csn_spatial(nl), n.iter = 600L, n.warmup = 600L,
    adapt.delta = 0.95, seed = 5L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)
  expect_lt(abs(nu$means[["psi_abund_cov1"]] - nl$means[["psi_abund_cov1"]]) /
              nu$sds[["psi_abund_cov1"]], 1)
})


test_that("dyn_abun() NUTS + icar field reproduces the nested-Laplace field, 0 divergences (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(4L); phi <- .csn_field(adj, seed = 44L); ng <- nrow(adj)
  set.seed(44); xab <- .csn_orth_cov(phi); Td <- 3L; Jd <- 3L
  N1 <- rpois(ng, exp(log(7) + 0.5 * xab + phi)); ya <- array(0L, c(ng, Jd, Td)); Ncur <- N1
  for (t in 1:Td) {
    for (j in 1:Jd) ya[, j, t] <- rbinom(ng, Ncur, 0.5)
    if (t < Td) Ncur <- rbinom(ng, Ncur, 0.6) + rpois(ng, 2)
  }
  nl <- tobs(~ abund_cov1 + icar(graph = adj), data = data.frame(abund_cov1 = xab),
             family = dyn_abun(K_max = 20L), detection = ~ 1, y = ya,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_dyn_abun_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 200L, n.warmup = 200L, max.treedepth = 8L, seed = 6L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - nl$means[["lambda_abund_cov1"]]) /
              nu$sds[["lambda_abund_cov1"]], 1)
})


# ---------------------------------------------------------------------------
# Intrinsic BYM2 fields under NUTS (gcol33/tulpaObs#113). BYM2 = a structured
# ICAR component (sum-to-zero eigen-loading) PLUS an iid component; the whitened
# raw ~ N(0, I_{2n-1}) stacks both, reconstructed as z = sqrt(rho) sf phi +
# sqrt(1-rho) theta. The unit field is NOT sum-to-zero (the iid part), so these
# assert on field-shape recovery + 0 divergences, not on centring.
# ---------------------------------------------------------------------------

test_that("removal() NUTS + bym2 field recovers + reproduces the nested-Laplace field (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 11L); ng <- nrow(adj)
  set.seed(11); xab <- .csn_orth_cov(phi)
  lambda <- exp(log(8) + 0.6 * xab + phi); p <- 0.45; K <- 4L
  N <- rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- N
  for (k in 1:K) { y[, k] <- rbinom(ng, rem, p); rem <- rem - y[, k] }
  nl <- tobs(~ abund_cov1 + bym2(graph = adj), data = data.frame(abund_cov1 = xab),
             family = removal(), detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_removal_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 3L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_length(nu$spatial_field, ng)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.8)
  expect_gt(cor(nu$spatial_field, phi), 0.55)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.6) /
              nu$sds[["lambda_abund_cov1"]], 3)
})


test_that("distance() NUTS + bym2 field recovers + reproduces the nested-Laplace field (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 22L); ng <- nrow(adj)
  set.seed(22); xab <- .csn_orth_cov(phi); cut <- c(0, 10, 20, 30, 40); B <- 4L
  sig <- exp(3.0); lamd <- exp(log(40) + 0.5 * xab + phi)
  g_mid <- exp(-((head(cut, -1) + tail(cut, -1)) / 2)^2 / (2 * sig^2))
  pi_b <- g_mid * diff(cut) / 40 * 0.6
  Nd <- rpois(ng, lamd); y <- matrix(0L, ng, B)
  for (i in seq_len(ng)) y[i, ] <- rbinom(B, Nd[i], pi_b)
  nl <- tobs(~ abund_cov1 + bym2(graph = adj), data = data.frame(abund_cov1 = xab),
             family = distance(key = "halfnorm", transect = "line", cutpoints = cut),
             detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_distance_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 4L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.8)
  expect_gt(cor(nu$spatial_field, phi), 0.55)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
              nu$sds[["lambda_abund_cov1"]], 3)
})


test_that("fp_occu() NUTS + bym2 field reproduces the nested-Laplace field, 0 divergences (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(7L); phi <- .csn_field(adj, sd_phi = 0.8, seed = 33L); ng <- nrow(adj)
  set.seed(33); xab <- .csn_orth_cov(phi); psi <- plogis(0.3 + 0.7 * xab + phi)
  z <- rbinom(ng, 1, psi); J <- 6L; yf <- integer(ng * J); idx <- 1L
  for (i in seq_len(ng)) for (j in 1:J) {
    yf[idx] <- if (z[i] == 1) sample(0:2, 1, prob = c(0.45, 0.25, 0.30))
               else sample(0:1, 1, prob = c(0.92, 0.08))
    idx <- idx + 1L
  }
  nl <- tobs(~ abund_cov1 + bym2(graph = adj), data = data.frame(abund_cov1 = xab),
             family = fp_occu(), detection = ~ 1, y = matrix(yf, ng, J, byrow = TRUE),
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_fp_occu_nuts_spatial(
    nl$model, .csn_spatial(nl), n.iter = 600L, n.warmup = 600L,
    adapt.delta = 0.95, seed = 5L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)
  expect_lt(abs(nu$means[["psi_abund_cov1"]] - nl$means[["psi_abund_cov1"]]) /
              nu$sds[["psi_abund_cov1"]], 1)
})


test_that("dyn_abun() NUTS + bym2 field reproduces the nested-Laplace field, 0 divergences (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(4L); phi <- .csn_field(adj, seed = 44L); ng <- nrow(adj)
  set.seed(44); xab <- .csn_orth_cov(phi); Td <- 3L; Jd <- 3L
  N1 <- rpois(ng, exp(log(7) + 0.5 * xab + phi)); ya <- array(0L, c(ng, Jd, Td)); Ncur <- N1
  for (t in 1:Td) {
    for (j in 1:Jd) ya[, j, t] <- rbinom(ng, Ncur, 0.5)
    if (t < Td) Ncur <- rbinom(ng, Ncur, 0.6) + rpois(ng, 2)
  }
  nl <- tobs(~ abund_cov1 + bym2(graph = adj), data = data.frame(abund_cov1 = xab),
             family = dyn_abun(K_max = 20L), detection = ~ 1, y = ya,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_dyn_abun_nuts_spatial(
    nl$model, .csn_spatial(nl), mixture = "poisson",
    n.iter = 200L, n.warmup = 200L, max.treedepth = 8L, seed = 6L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.7)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - nl$means[["lambda_abund_cov1"]]) /
              nu$sds[["lambda_abund_cov1"]], 1)
})


# ---------------------------------------------------------------------------
# >= 20-seed fixed-effect 95%-interval coverage for the intrinsic icar field
# under NUTS (gcol33/tulpaObs#113 Definition of Done). removal() is the cheapest
# count-marginal family, so it carries the coverage evidence for the sum-to-zero
# reparam; the field is fixed-hyper (warmed from nested-Laplace), so the
# abundance-slope interval is the calibrated quantity.
# ---------------------------------------------------------------------------

test_that("removal() NUTS + icar field: abundance-slope 95% coverage over 20 seeds (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); ng <- nrow(adj)
  n_seeds <- 20L; beta_truth <- 0.6
  covered <- logical(n_seeds); diverg <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    phi <- .csn_field(adj, seed = 100L + s)
    set.seed(100L + s); xab <- .csn_orth_cov(phi)
    lambda <- exp(log(8) + beta_truth * xab + phi); p <- 0.45; K <- 4L
    N <- rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- N
    for (k in 1:K) { y[, k] <- rbinom(ng, rem, p); rem <- rem - y[, k] }
    nl <- tryCatch(tobs(~ abund_cov1 + icar(graph = adj),
                        data = data.frame(abund_cov1 = xab), family = removal(),
                        detection = ~ 1, y = y, method = "nested_laplace",
                        control = list(verbose = FALSE, progress = FALSE)),
                   error = function(e) NULL)
    if (is.null(nl)) { covered[s] <- NA; next }
    nu <- tryCatch(tulpaObs:::.tobs_fit_removal_nuts_spatial(
      nl$model, .csn_spatial(nl), mixture = "poisson",
      n.iter = 400L, n.warmup = 400L, seed = 1L), error = function(e) NULL)
    if (is.null(nu)) { covered[s] <- NA; next }
    m <- nu$means[["lambda_abund_cov1"]]; sd <- nu$sds[["lambda_abund_cov1"]]
    covered[s] <- abs(m - beta_truth) < 1.96 * sd
    diverg[s]  <- mean(nu$divergent)
  }
  ok <- !is.na(covered)
  expect_gte(mean(ok), 0.9)                       # nearly every seed fits
  expect_equal(mean(diverg[ok]), 0)               # clean sum-to-zero geometry
  expect_gte(mean(covered[ok]), 0.85)             # 95% Wald coverage at the floor
})


# ---------------------------------------------------------------------------
# Hazard-rate key distance() NUTS + areal field (gcol33/tulpaObs#114). The
# hazard key adds a single global log-shape coordinate eta_b, orthogonal to the
# fixed-hyper field block (the C++ target places it before the whitened field
# raw). The fit recovers the abundance field AND the scalar shape.
# ---------------------------------------------------------------------------

test_that("distance() hazard-key NUTS + car_proper recovers field + shape (#114)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 7L); ng <- nrow(adj)
  set.seed(7); xab <- .csn_orth_cov(phi); cut <- c(0, 10, 20, 30, 40); W <- 40
  sigma <- 18; shape <- 3.5; B <- 4L
  lam <- exp(log(45) + 0.5 * xab + phi); pib <- .csn_haz_pi(cut, sigma, shape, W)
  Nd <- rpois(ng, lam); y <- matrix(0L, ng, B)
  for (i in seq_len(ng)) y[i, ] <- rbinom(B, Nd[i], pib)
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), data = data.frame(abund_cov1 = xab),
             family = distance(key = "hazard", transect = "line", cutpoints = cut),
             detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_distance_nuts_spatial(
    nl$model, spatial = .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 4L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_gt(cor(nu$spatial_field, phi), 0.7)
  expect_lt(abs(nu$means[["lambda_abund_cov1"]] - 0.5) /
              nu$sds[["lambda_abund_cov1"]], 3)
  expect_lt(abs(nu$means[["log_shape"]] - log(shape)), 0.4)   # scalar shape recovered
})


test_that("distance() hazard-key NUTS + icar recovers field + shape (#114)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .csn_grid_adj(6L); phi <- .csn_field(adj, seed = 7L); ng <- nrow(adj)
  set.seed(7); xab <- .csn_orth_cov(phi); cut <- c(0, 10, 20, 30, 40); W <- 40
  sigma <- 18; shape <- 3.5; B <- 4L
  lam <- exp(log(45) + 0.5 * xab + phi); pib <- .csn_haz_pi(cut, sigma, shape, W)
  Nd <- rpois(ng, lam); y <- matrix(0L, ng, B)
  for (i in seq_len(ng)) y[i, ] <- rbinom(B, Nd[i], pib)
  nl <- tobs(~ abund_cov1 + icar(graph = adj), data = data.frame(abund_cov1 = xab),
             family = distance(key = "hazard", transect = "line", cutpoints = cut),
             detection = ~ 1, y = y, method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  nu <- tulpaObs:::.tobs_fit_distance_nuts_spatial(
    nl$model, spatial = .csn_spatial(nl), mixture = "poisson",
    n.iter = 500L, n.warmup = 500L, seed = 4L)
  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$divergent), 0)
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)                # sum-to-zero centred
  expect_gt(cor(nu$spatial_field, phi), 0.7)
  expect_lt(abs(nu$means[["log_shape"]] - log(shape)), 0.4)
})
