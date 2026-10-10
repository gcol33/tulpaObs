# Parallel C++ pointwise log-likelihood for the cover() hurdle
# (cpp_cover_hurdle_ploglik, the WAIC / PSIS-LOO input). The R oracle is
# .tobs_cover_hurdle_ll; the kernel mirrors it draw for draw over all five
# positive families, so the two agree to libm rounding and the result is
# thread-count invariant (the draw loop has no shared writes).

.mk_cover_ploglik_case <- function(fam, N = 60L, S = 40L, seed = 11) {
  set.seed(seed)
  occur   <- rbinom(N, 1L, 0.5)
  present <- which(occur == 1L); Np <- length(present)
  eta_occ <- matrix(rnorm(S * N, 0, 1), S, N)
  eta_pos <- matrix(rnorm(S * Np, 0, 1), S, Np)
  disp    <- if (fam == "beta") runif(S, 3, 25) else runif(S, 0.3, 1.2)
  bounds  <- NULL
  if (fam == "beta") {
    y_pos <- runif(Np, 1e-3, 1 - 1e-3)
  } else if (fam == "ordinal") {
    y_pos <- rep(0, Np)
    lo <- runif(Np, -3, -0.5); bounds <- list(lower = lo, upper = lo + runif(Np, 0.3, 2))
  } else if (fam == "lognormal_trunc") {
    y_pos <- rnorm(Np, -1, 0.5); bounds <- list(trunc_upper = rep(0, Np))
  } else {
    y_pos <- rnorm(Np, -1, 0.5)
  }
  list(fam = fam, eta_occ = eta_occ, eta_pos = eta_pos, disp = disp,
       occur = occur, y_pos = y_pos, idx_pos = present, bounds = bounds)
}

test_that("C++ cover hurdle pointwise loglik matches the R oracle (all families)", {
  for (fam in c("lognormal", "lognormal_trunc", "ordinal", "beta", "gaussian")) {
    d <- .mk_cover_ploglik_case(fam)
    R <- .tobs_cover_hurdle_ll(d$eta_occ, d$eta_pos, d$disp, d$occur, d$y_pos,
                               d$idx_pos, fam, bounds = d$bounds)
    C <- .cover_hurdle_ploglik_core(d$eta_occ, d$eta_pos, d$disp, d$occur,
                                    d$y_pos, d$idx_pos, fam, bounds = d$bounds,
                                    n_threads = 1L)
    expect_equal(C, R, tolerance = 1e-8, info = fam)
  }
})

test_that("C++ cover hurdle pointwise loglik is thread-count invariant", {
  d  <- .mk_cover_ploglik_case("beta", N = 120L, S = 64L)
  C1 <- .cover_hurdle_ploglik_core(d$eta_occ, d$eta_pos, d$disp, d$occur,
                                   d$y_pos, d$idx_pos, "beta", bounds = d$bounds,
                                   n_threads = 1L)
  C8 <- .cover_hurdle_ploglik_core(d$eta_occ, d$eta_pos, d$disp, d$occur,
                                   d$y_pos, d$idx_pos, "beta", bounds = d$bounds,
                                   n_threads = 8L)
  expect_identical(C8, C1)
})

test_that("the positive arm scores the density tulpa's Laplace fit evaluates", {
  eta  <- c(-40, -34, -3, 0, 2.5, 34, 40)
  occur <- rep(1L, length(eta))
  cases <- list(
    beta      = list(y = rep(0.37, length(eta)), disp = 6.5, fam = "beta"),
    lognormal = list(y = rep(-0.8, length(eta)), disp = 0.6, fam = "lognormal"),
    gaussian  = list(y = rep(1.4, length(eta)),  disp = 0.6, fam = "gaussian")
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    # eta_occ = 50 makes log p ~ 0, leaving the positive-arm density.
    C <- .cover_hurdle_ploglik_core(matrix(50, 1, length(eta)), matrix(eta, 1),
                                    cs$disp, occur, cs$y, seq_along(eta), nm,
                                    n_threads = 1L)
    y_engine <- if (nm == "lognormal") exp(cs$y) else cs$y
    engine <- vapply(seq_along(eta), function(i)
      tulpa:::cpp_family_terms(y_engine[i], 1L, eta[i], cs$fam, cs$disp)[["log_lik"]],
      numeric(1))
    expect_equal(as.vector(C) - plogis(50, log.p = TRUE), engine,
                 tolerance = 1e-10, info = nm)
  }
})
