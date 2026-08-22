# Temporal() field composed WITH the areal spatial field on the count observation
# families (removal / distance / fp_occu / dyn_abun) under method =
# "nested_laplace". The shared areal-BFGS driver (R/areal_bfgs.R) carries a
# second latent block (the AR1 temporal precision) alongside the spatial field,
# both grid-integrated. Each fit simulates a dataset with BOTH a smoothed
# ICAR-like spatial field on the abundance / occupancy arm and an AR1 temporal
# effect over seasons, then checks recovery of the fixed- effect slope, the
# spatial field shape, and -- where the field is identified -- the temporal field
# shape and autocorrelation.
#
# Identifiability differs by family. The count likelihoods (removal, distance)
# identify both fields well at a moderate grid, so they assert recovery of the
# slope + both field shapes + the AR1 rho. The occupancy field (fp_occu, one
# binary site per node) and the forward-HMM open N-mixture (dyn_abun) identify the
# temporal component only weakly at an affordable grid, so they assert the
# composition runs and recovers the fixed slope + the spatial field, with the
# temporal block estimated (length T) but not held to a tight truth correlation.

.cts_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

# Smoothed, demeaned ICAR-like field on the grid.
.cts_field <- function(adj, sd_phi = 0.7, seed = 1L) {
  set.seed(seed); ng <- nrow(adj); phi <- as.numeric(scale(rnorm(ng)))
  for (r in 1:3) {
    pn <- phi
    for (i in seq_len(ng)) { nb <- which(adj[i, ] == 1L); pn[i] <- 0.5 * phi[i] + 0.5 * mean(phi[nb]) }
    phi <- pn
  }
  phi <- sd_phi * as.numeric(scale(phi)); phi - mean(phi)
}

# A demeaned AR1 path over T seasons (the temporal truth).
.cts_ar1 <- function(Tt, rho, sd_u, seed) {
  set.seed(seed); u <- numeric(Tt); u[1] <- rnorm(1, 0, sd_u)
  for (t in 2:Tt) u[t] <- rho * u[t - 1] + rnorm(1, 0, sd_u * sqrt(1 - rho^2))
  u - mean(u)
}

# A covariate orthogonal to the spatial field, so the abundance / occupancy slope
# is identified separately from the field.
.cts_orth <- function(phi) {
  x <- rnorm(length(phi))
  as.numeric(scale(residuals(stats::lm(x ~ phi))))
}


test_that("removal() areal field + AR1 temporal recovers the slope + both fields", {
  skip_on_cran()
  skip_if_fast()
  side <- 6L; adj <- .cts_grid_adj(side); nsite <- nrow(adj)
  slope <- numeric(5); sp_cor <- numeric(5); tmp_cor <- numeric(5); rho_hat <- numeric(5)
  for (sd in 1:5) {
    phi <- .cts_field(adj, sd_phi = 0.7, seed = 10L + sd)
    Tt <- 4L; u <- .cts_ar1(Tt, 0.7, 0.9, 100L + sd)
    set.seed(sd); seas <- sample(rep(seq_len(Tt), length.out = nsite)); xab <- .cts_orth(phi)
    dat <- data.frame(abund_cov1 = xab, season = seas)
    lambda <- exp(log(10) + 0.6 * xab + phi + u[seas]); K <- 5L; p <- 0.5
    N <- rpois(nsite, lambda); y <- matrix(0L, nsite, K); rem <- N
    for (k in 1:K) { y[, k] <- rbinom(nsite, rem, p); rem <- rem - y[, k] }
    f <- tobs(~ abund_cov1 + icar(graph = adj) + temporal(season, type = "ar1"),
              data = dat, family = removal(), detection = ~ 1, y = y,
              method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
    expect_identical(f$method, "nested_laplace")
    expect_length(f$temporal_field, Tt)
    slope[sd]   <- unname(coef(f)$lambda["abund_cov1"])
    sp_cor[sd]  <- cor(f$spatial_field, phi)
    tmp_cor[sd] <- cor(f$temporal_field, u)
    rho_hat[sd] <- unname(f$temporal_hyper[["rho"]])
  }
  expect_lt(abs(median(slope) - 0.6), 0.12)      # abundance slope recovered
  expect_gt(median(sp_cor), 0.85)                # spatial field shape recovered
  expect_gt(median(tmp_cor), 0.6)                # temporal field shape recovered
  expect_gt(median(rho_hat), 0.3)                # positive AR1 autocorrelation
})


test_that("distance() areal field + AR1 temporal recovers the slope + both fields", {
  skip_on_cran()
  skip_if_fast()
  side <- 6L; adj <- .cts_grid_adj(side); nsite <- nrow(adj)
  cut <- c(0, 10, 20, 30, 40); B <- 4L; sig <- exp(3.0)
  g_mid <- exp(-((head(cut, -1) + tail(cut, -1)) / 2)^2 / (2 * sig^2))
  pi_b  <- g_mid * diff(cut) / 40 * 0.6
  slope <- numeric(5); sp_cor <- numeric(5); tmp_cor <- numeric(5); rho_hat <- numeric(5)
  for (sd in 1:5) {
    phi <- .cts_field(adj, sd_phi = 0.7, seed = 20L + sd)
    Tt <- 4L; u <- .cts_ar1(Tt, 0.7, 0.9, 200L + sd)
    set.seed(sd); seas <- sample(rep(seq_len(Tt), length.out = nsite)); xab <- .cts_orth(phi)
    dat <- data.frame(abund_cov1 = xab, season = seas)
    lamd <- exp(log(60) + 0.5 * xab + phi + u[seas])
    Nd <- rpois(nsite, lamd); y <- matrix(0L, nsite, B)
    for (i in seq_len(nsite)) y[i, ] <- rbinom(B, Nd[i], pi_b)
    f <- tobs(~ abund_cov1 + icar(graph = adj) + temporal(season, type = "ar1"),
              data = dat, family = distance(key = "halfnorm", transect = "line", cutpoints = cut),
              detection = ~ 1, y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
    expect_length(f$temporal_field, Tt)
    slope[sd]   <- unname(coef(f)$lambda["abund_cov1"])
    sp_cor[sd]  <- cor(f$spatial_field, phi)
    tmp_cor[sd] <- cor(f$temporal_field, u)
    rho_hat[sd] <- unname(f$temporal_hyper[["rho"]])
  }
  expect_lt(abs(median(slope) - 0.5), 0.12)
  expect_gt(median(sp_cor), 0.85)
  expect_gt(median(tmp_cor), 0.6)
  expect_gt(median(rho_hat), 0.3)
})


test_that("fp_occu() areal field + AR1 temporal composes and recovers the psi slope", {
  skip_on_cran()
  skip_if_fast()
  # The occupancy field is weakly identified (one binary site per node), so the
  # invariant is that the temporal block composes (length T) and the fixed psi
  # slope + spatial field recover -- not a tight temporal-truth correlation.
  side <- 7L; adj <- .cts_grid_adj(side); nsite <- nrow(adj)
  slope <- numeric(3); sp_cor <- numeric(3)
  for (sd in 1:3) {
    phi <- .cts_field(adj, sd_phi = 0.8, seed = 30L + sd)
    Tt <- 4L; u <- .cts_ar1(Tt, 0.7, 0.9, 300L + sd)
    set.seed(sd); seas <- sample(rep(seq_len(Tt), length.out = nsite)); xab <- .cts_orth(phi)
    dat <- data.frame(abund_cov1 = xab, season = seas)
    psi <- plogis(0.3 + 0.7 * xab + phi + u[seas]); z <- rbinom(nsite, 1, psi); J <- 6L
    yf <- matrix(0L, nsite, J)
    for (i in seq_len(nsite)) for (j in 1:J)
      yf[i, j] <- if (z[i] == 1L) sample(0:2, 1L, prob = c(0.45, 0.25, 0.30))
                  else sample(0:1, 1L, prob = c(0.92, 0.08))
    f <- tobs(~ abund_cov1 + icar(graph = adj) + temporal(season, type = "ar1"),
              data = dat, family = fp_occu(), detection = ~ 1, y = yf,
              method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
    expect_length(f$temporal_field, Tt)
    expect_false(is.null(f$temporal_hyper))
    slope[sd]  <- unname(coef(f)$psi["abund_cov1"])
    sp_cor[sd] <- cor(f$spatial_field, phi)
  }
  expect_lt(abs(median(slope) - 0.7), 0.4)       # psi slope in the right region
  expect_gt(median(sp_cor), 0.15)                # spatial field weakly recovered
})


test_that("dyn_abun() areal field + AR1 temporal composes and recovers the lambda slope", {
  skip_on_cran()
  skip_if_fast()
  # The open-population forward-HMM marginal identifies the temporal component
  # only weakly at an affordable grid, so the invariant is composition (length T)
  # + recovery of the fixed initial-abundance slope and the spatial field.
  side <- 4L; adj <- .cts_grid_adj(side); nsite <- nrow(adj)
  slope <- numeric(3); sp_cor <- numeric(3)
  for (sd in 1:3) {
    phi <- .cts_field(adj, sd_phi = 0.6, seed = 40L + sd)
    Tt <- 4L; u <- .cts_ar1(Tt, 0.7, 0.9, 400L + sd)
    set.seed(sd); seas <- sample(rep(seq_len(Tt), length.out = nsite)); xab <- .cts_orth(phi)
    dat <- data.frame(abund_cov1 = xab, season = seas); Td <- 3L; Jd <- 3L
    N1 <- rpois(nsite, exp(log(8) + 0.5 * xab + phi + u[seas]))
    ya <- array(0L, c(nsite, Jd, Td)); Ncur <- N1
    for (t in 1:Td) { for (j in 1:Jd) ya[, j, t] <- rbinom(nsite, Ncur, 0.5)
      if (t < Td) Ncur <- rbinom(nsite, Ncur, 0.6) + rpois(nsite, 2) }
    f <- tobs(~ abund_cov1 + icar(graph = adj) + temporal(season, type = "ar1"),
              data = dat, family = dyn_abun(K_max = 30L), detection = ~ 1, y = ya,
              method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
    expect_length(f$temporal_field, Tt)
    expect_false(is.null(f$temporal_hyper))
    slope[sd]  <- unname(coef(f)$lambda["abund_cov1"])
    sp_cor[sd] <- cor(f$spatial_field, phi)
  }
  expect_lt(abs(median(slope) - 0.5), 0.35)
  expect_gt(median(sp_cor), 0.6)
})


test_that("count() temporal() term errors with a pointer", {
  # The observation families (removal / distance / fp_occu / dyn_abun) all support
  # a temporal-only field now. The count() relative-abundance GLMM still gates
  # every non-areal structured term (temporal / re / svc / latent), so it carries
  # the not-yet-wired temporal pointer here.
  nsite <- 40L
  set.seed(1); seas <- sample(rep(1:3, length.out = nsite))
  dat <- data.frame(abund_cov1 = rnorm(nsite), season = seas)
  y <- rpois(nsite, 3)
  expect_error(
    tobs(~ abund_cov1 + temporal(season, type = "ar1"), data = dat,
         family = count(), y = y,
         control = list(verbose = FALSE, progress = FALSE)),
    "temporal"
  )
})
