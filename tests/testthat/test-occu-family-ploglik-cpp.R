# Parallel C++ pointwise log-likelihood kernels for the draw-matrix occupancy
# families (single-season, multi-season dynamic HMM, multi-source integrated).
# Each family's per-observation marginal was a pure-R loop; the ports mirror it
# and parallelise over the observation index. The oracles below are the original
# R loops, so the tests assert the C++ path reproduces them to libm rounding and
# is thread-count invariant (a divergence is a bug, not an improvement: both
# evaluate the same marginal).

.lp  <- function(e) stats::plogis(e, log.p = TRUE)
.l1  <- function(e) stats::plogis(-e, log.p = TRUE)
.lae <- function(a, b) {
  m <- pmax(a, b); out <- m + log1p(exp(pmin(a, b) - m))
  bi <- is.infinite(m) & (a == b); out[bi] <- m[bi]; out
}

test_that("single-season occupancy pointwise loglik: C++ == R oracle", {
  set.seed(31)
  for (rep in 1:3) {
    S <- 50L; N <- 120L; mv <- 5L
    eta_psi <- matrix(rnorm(S * N, 0, 1.2), S, N)
    eta_p   <- matrix(rnorm(S * N, 0, 1.2), S, N)
    y <- matrix(-1L, N, mv)
    for (i in seq_len(N)) { nv <- sample(0:mv, 1); if (nv > 0) y[i, seq_len(nv)] <- rbinom(nv, 1L, 0.4) }
    R <- matrix(0, S, N)
    for (i in seq_len(N)) {
      yi <- y[i, ]; valid <- yi >= 0; nv <- sum(valid); if (nv == 0L) next
      k1 <- sum(yi[valid] == 1); k0 <- nv - k1
      R[, i] <- if (k1 > 0L) .lp(eta_psi[, i]) + k1 * .lp(eta_p[, i]) + k0 * .l1(eta_p[, i])
                else .lae(.lp(eta_psi[, i]) + nv * .l1(eta_p[, i]), .l1(eta_psi[, i]))
    }
    C1 <- cpp_occu_single_ploglik(eta_psi, eta_p, y, 1L)
    C4 <- cpp_occu_single_ploglik(eta_psi, eta_p, y, 4L)
    expect_equal(C1, R, tolerance = 1e-9)
    expect_identical(C4, C1)
  }
})

test_that("dynamic (multi-season HMM) pointwise loglik: C++ == R oracle", {
  set.seed(41)
  n_sites <- 25L; mv <- 4L; Tn <- 3L; S <- 40L
  Xs <- replicate(4, cbind(1, rnorm(n_sites)), simplify = FALSE)
  y <- array(-1L, dim = c(n_sites, mv, Tn))
  nvis <- integer(n_sites * Tn); adet <- logical(n_sites * Tn)
  for (i in seq_len(n_sites)) for (t in seq_len(Tn)) {
    nv <- sample(0:mv, 1); idx <- (i - 1L) * Tn + t; nvis[idx] <- nv
    if (nv > 0) { yy <- rbinom(nv, 1L, 0.4); y[i, seq_len(nv), t] <- yy; adet[idx] <- any(yy == 1) }
  }
  model <- list(model_type = "dynamic", n_sites = n_sites, n_seasons = Tn,
                max_visits = mv, y = y, n_visits = nvis, any_detected = adet,
                process_info = lapply(Xs, function(X) list(p = ncol(X))),
                X_processes = Xs)
  draws <- matrix(rnorm(S * 8, 0, 0.8), S, 8)
  e1 <- .tobs_eta_draws(model, draws, 1L); ep <- .tobs_eta_draws(model, draws, 2L)
  eg <- .tobs_eta_draws(model, draws, 3L); ee <- .tobs_eta_draws(model, draws, 4L)
  NEG <- -1e10; R <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    lpp <- .lp(ep[, i]); l1mp <- .l1(ep[, i]); lgam <- .lp(eg[, i]); l1mgam <- .l1(eg[, i])
    leps <- .lp(ee[, i]); l1meps <- .l1(ee[, i]); a_occ <- .lp(e1[, i]); a_un <- .l1(e1[, i])
    sll <- numeric(S)
    for (t in seq_len(Tn)) {
      idx <- (i - 1L) * Tn + t; nv <- nvis[idx]
      if (nv > 0L) {
        yv <- y[i, , t]; valid <- yv >= 0; k1 <- sum(yv[valid] == 1); k0 <- sum(yv[valid] == 0)
        if (adet[idx]) { sll <- sll + a_occ + (k1 * lpp + k0 * l1mp); a_occ <- numeric(S); a_un <- rep(NEG, S) }
        else { t1 <- a_occ + nv * l1mp; ln <- .lae(t1, a_un); sll <- sll + ln; a_occ <- t1 - ln; a_un <- a_un - ln }
      }
      if (t < Tn) { no <- .lae(a_occ + l1meps, a_un + lgam); nu <- .lae(a_occ + leps, a_un + l1mgam); a_occ <- no; a_un <- nu }
    }
    R[, i] <- sll
  }
  C1 <- .tobs_ploglik_dynamic(model, draws, 1L)
  C4 <- .tobs_ploglik_dynamic(model, draws, 4L)
  expect_equal(C1, R, tolerance = 1e-9)
  expect_identical(C4, C1)
})

test_that("integrated (multi-source) pointwise loglik: C++ == R oracle", {
  set.seed(51)
  n_sites <- 30L; n_src <- 3L; mv <- 4L; S <- 40L
  Xs <- c(list(cbind(1, rnorm(n_sites))), replicate(n_src, cbind(1, rnorm(n_sites)), simplify = FALSE))
  site_maps <- vector("list", n_src); y_sources <- vector("list", n_src)
  for (s in seq_len(n_src)) {
    sites <- sort(sample(seq_len(n_sites), round(n_sites * 0.7)))
    site_maps[[s]] <- sites - 1L
    y_sources[[s]] <- t(vapply(sites, function(i) {
      nv <- sample(1:mv, 1); yy <- rep(-1L, mv); yy[seq_len(nv)] <- rbinom(nv, 1L, 0.35); yy
    }, integer(mv)))
  }
  model <- list(model_type = "integrated", n_sites = n_sites, n_sources = n_src,
                site_maps = site_maps, y_sources = y_sources,
                process_info = lapply(Xs, function(X) list(p = ncol(X))),
                X_processes = Xs)
  draws <- matrix(rnorm(S * (2 * (n_src + 1)), 0, 0.8), S, 2 * (n_src + 1))
  eta_psi <- .tobs_eta_draws(model, draws, 1L)
  lps <- lapply(seq_len(n_src), function(s) .lp(.tobs_eta_draws(model, draws, 1L + s)))
  l1ps <- lapply(seq_len(n_src), function(s) .l1(.tobs_eta_draws(model, draws, 1L + s)))
  R <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    ldo <- numeric(S); anyd <- FALSE
    for (s in seq_len(n_src)) {
      loc <- which(site_maps[[s]] + 1L == i); if (!length(loc)) next
      yv <- y_sources[[s]][loc[1L], ]; valid <- yv >= 0; nv <- sum(valid); if (nv == 0L) next
      k1 <- sum(yv[valid] == 1); k0 <- nv - k1
      ldo <- ldo + k1 * lps[[s]][, i] + k0 * l1ps[[s]][, i]; if (k1 > 0) anyd <- TRUE
    }
    R[, i] <- if (anyd) .lp(eta_psi[, i]) + ldo else .lae(.lp(eta_psi[, i]) + ldo, .l1(eta_psi[, i]))
  }
  C1 <- .tobs_ploglik_integrated(model, draws, 1L)
  C4 <- .tobs_ploglik_integrated(model, draws, 4L)
  expect_equal(C1, R, tolerance = 1e-9)
  expect_identical(C4, C1)
})
