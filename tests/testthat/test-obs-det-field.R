# DETECTION-arm areal field on the observation families: removal / fp_occu /
# dyn_abun carry a spatially-varying field on their detection arm (capture logit /
# true-positive p11 / detection p) via the shared areal-BFGS nested-Laplace
# driver, the field routed to the detection arm instead of the abundance /
# occupancy arm. removal's detection design is per-pass, so its per-observation
# eta-gradient is summed to a per-site field gradient; fp_occu and dyn_abun carry
# per-site detection designs, so the field loads on eta directly (no aggregation).
# Detection fields are more weakly identified than abundance fields (a per-site
# binary-ish detection channel), so the recovery bar is a correlation with the
# truth surface, the arm slope, and a `spatial_field_arm` label of "detection". A
# detection-arm field under NUTS stays gated (the C++ field block loads on the
# abundance / occupancy arm); that gate is asserted too.

# Smoothed ICAR-like field on a side x side grid, demeaned (sum-to-zero).
.odf_field <- function(adj, sd_phi = 0.7, seed = 1L) {
  set.seed(seed); ng <- nrow(adj); phi <- as.numeric(scale(stats::rnorm(ng)))
  for (r in 1:3) { pn <- phi
    for (i in seq_len(ng)) { nb <- which(adj[i, ] == 1L); pn[i] <- 0.5*phi[i] + 0.5*mean(phi[nb]) }
    phi <- pn }
  phi <- sd_phi * as.numeric(scale(phi)); phi - mean(phi)
}


test_that("removal() detection-arm ICAR field recovers the capture surface (#114)", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(7L); ng <- nrow(adj); n_seeds <- 6L
  fcor <- slope <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    phi <- .odf_field(adj, sd_phi = 0.7, seed = 40L + s)
    set.seed(40L + s); x <- as.numeric(scale(stats::rnorm(ng)))
    lambda <- exp(log(9) + 0.5 * x)                 # abundance: no field
    p <- stats::plogis(0.2 + phi)                    # detection: spatial field
    K <- 4L; Nn <- stats::rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- Nn
    for (k in 1:K) { y[, k] <- stats::rbinom(ng, rem, p); rem <- rem - y[, k] }
    fit <- tryCatch(tobs(~ x, data = data.frame(x = x), family = removal(),
                         detection = ~ icar(graph = adj), y = y,
                         method = "nested_laplace",
                         control = list(verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) expect_identical(fit$spatial_field_arm, "detection")
    slope[s] <- fit$means[["lambda_x"]]
    fcor[s]  <- abs(stats::cor(fit$spatial_field, phi))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_lt(abs(mean(slope[ok]) - 0.5), 0.15)         # abundance slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.6)        # capture field recovered
})


test_that("fp_occu() detection-arm (p11) ICAR field recovers the surface (#114)", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(8L); ng <- nrow(adj); n_seeds <- 6L
  fcor <- slope <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    phi <- .odf_field(adj, sd_phi = 0.7, seed = 80L + s)
    set.seed(80L + s); x <- as.numeric(scale(stats::rnorm(ng)))
    psi <- stats::plogis(0.4 + 0.6 * x)              # occupancy: no field
    p11 <- stats::plogis(0.5 + phi)                   # true-positive detection: field
    J <- 10L; z <- stats::rbinom(ng, 1, psi); y <- matrix(0L, ng, J)
    for (i in seq_len(ng)) for (j in seq_len(J)) {
      y[i, j] <- if (z[i] == 1L)
        sample(0:2, 1, prob = c(1 - p11[i], 0.2 * p11[i], 0.8 * p11[i]))
      else sample(0:1, 1, prob = c(0.92, 0.08))
    }
    fit <- tryCatch(tobs(~ x, data = data.frame(x = x), family = fp_occu(),
                         detection = ~ icar(graph = adj), y = y,
                         method = "nested_laplace",
                         control = list(verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) expect_identical(fit$spatial_field_arm, "detection")
    slope[s] <- fit$means[["psi_x"]]
    fcor[s]  <- abs(stats::cor(fit$spatial_field, phi))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_lt(abs(mean(slope[ok]) - 0.6), 0.25)         # occupancy slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.55)       # p11 field recovered
})


test_that("dyn_abun() detection-arm ICAR field recovers the surface (#114)", {
  skip_on_cran()
  skip_if_fast()
  # dyn_abun's forward-HMM marginal over an areal grid is costly; a compact grid
  # and a few seeds validate the per-site detection-arm routing.
  adj <- rook_adj(5L); ng <- nrow(adj); n_seeds <- 4L
  fcor <- slope <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    phi <- .odf_field(adj, sd_phi = 0.7, seed = 90L + s)
    set.seed(90L + s); Tt <- 3L; J <- 4L; x <- as.numeric(scale(stats::rnorm(ng)))
    lambda <- exp(log(7) + 0.4 * x)                  # init abundance: no field
    p <- stats::plogis(0.2 + phi)                     # detection: spatial field
    omega <- 0.7; gamma <- 1.2
    ya <- array(0L, dim = c(ng, J, Tt)); Nt <- stats::rpois(ng, lambda)
    for (t in 1:Tt) { if (t > 1) Nt <- stats::rbinom(ng, Nt, omega) + stats::rpois(ng, gamma)
      for (j in 1:J) ya[, j, t] <- stats::rbinom(ng, Nt, p) }
    fit <- tryCatch(tobs(~ x, data = data.frame(x = x), family = dyn_abun(),
                         detection = ~ icar(graph = adj), y = ya,
                         method = "nested_laplace",
                         control = list(verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) expect_identical(fit$spatial_field_arm, "detection")
    slope[s] <- fit$means[["lambda_x"]]
    fcor[s]  <- abs(stats::cor(fit$spatial_field, phi))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_lt(abs(mean(slope[ok]) - 0.4), 0.25)         # abundance slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.55)       # detection field recovered
})


test_that("detection-arm areal field under NUTS stays gated (#114)", {
  adj <- rook_adj(4L); ng <- nrow(adj)
  set.seed(1); x <- stats::rnorm(ng); y <- matrix(stats::rpois(ng * 4L, 3), ng, 4L)
  expect_error(
    suppressWarnings(tobs(~ x, data = data.frame(x = x), family = removal(),
                          detection = ~ icar(graph = adj), y = y, method = "nuts",
                          control = list(verbose = FALSE, progress = FALSE))),
    "detection-arm|nested_laplace")
})
