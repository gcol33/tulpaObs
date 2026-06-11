# Multi-seed recovery for the two "working" families the audit flagged as
# under-validated (gcol33/tulpaObs#84): dynamic occupancy (dyn_occu) had only a
# single-seed point check, and single-source integrated occupancy (int_occu) had
# FD / smoke / loglik-vs-EM checks but no |est - truth| recovery. These add
# multi-seed recovery so a calibration regression in either family is caught.

# --------------------------------------------------------------------------- #
# dyn_occu: multi-seed (psi1, gamma, epsilon, p) recovery.                      #
# (The dynamic-abundance family is the template that makes a dynamic family     #
# "solid"; dyn_occu now gets the same multi-seed treatment rather than one      #
# point check.)                                                                 #
# --------------------------------------------------------------------------- #

test_that("dyn_occu recovers (psi1, gamma, epsilon, p) across seeds", {
  skip_on_cran()
  skip_if_fast()

  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p_true <- 0.5
  truth <- c(psi1 = psi1, gam = gam, eps = eps, p = p_true)
  n_sites <- 200L; n_seasons <- 4L; n_visits <- 4L
  n_seed <- 12L

  rec <- matrix(NA_real_, n_seed, 4L,
                dimnames = list(NULL, c("psi1", "gam", "eps", "p")))
  for (s in seq_len(n_seed)) {
    set.seed(2000L + s)
    z <- matrix(NA_integer_, n_sites, n_seasons)
    z[, 1] <- rbinom(n_sites, 1, psi1)
    for (t in 2:n_seasons) {
      z[, t] <- ifelse(z[, t - 1] == 1, rbinom(n_sites, 1, 1 - eps),
                       rbinom(n_sites, 1, gam))
    }
    y <- array(0L, dim = c(n_sites, n_visits, n_seasons))
    for (i in seq_len(n_sites)) for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t]) rbinom(n_visits, 1, p_true) else 0L
    }
    fit <- tobs(~ 1, data = data.frame(idx = seq_len(n_sites)),
                family = dyn_occu(), detection = ~ 1, y = y,
                col_formula = ~ 1, ext_formula = ~ 1,
                control = list(verbose = FALSE))
    rec[s, ] <- c(
      plogis(fit$means[["psi1_(Intercept)"]]),
      plogis(fit$means[["gamma_(Intercept)"]]),
      plogis(fit$means[["epsilon_(Intercept)"]]),
      plogis(fit$means[["p_(Intercept)"]]))
  }

  # Across-seed MEDIAN recovers each parameter (robust to per-seed MC noise).
  med <- apply(rec, 2L, median)
  expect_lt(abs(med[["psi1"]] - truth[["psi1"]]), 0.08)
  expect_lt(abs(med[["gam"]]  - truth[["gam"]]),  0.08)
  expect_lt(abs(med[["eps"]]  - truth[["eps"]]),  0.08)
  expect_lt(abs(med[["p"]]    - truth[["p"]]),    0.08)
})

# --------------------------------------------------------------------------- #
# dyn_occu HMM marginal: independent R forward-recursion anchor.               #
# The fit integrates out the latent occupancy chain with a forward recursion   #
# (HMM filtering). Re-derive that marginal in plain R for one site and check    #
# the engine's per-site contribution matches -- catches a transition / emission #
# transcription error the recovery test alone could miss.                       #
# --------------------------------------------------------------------------- #

# Independent R forward recursion for the dynamic-occupancy marginal P(y | params)
# of one site (T x J detection history):
#   alpha_1(z) = pi(z) b_1(z); alpha_t(z) = b_t(z) sum_{z'} alpha_{t-1}(z') A(z',z)
#   marginal = sum_z alpha_T(z); b_t(1) = prod_j p^{y}(1-p)^{1-y}, b_t(0) = 1{all y=0}
# A(z_{t-1}, z_t): from unoccupied colonize w.p. gamma; from occupied go extinct
# w.p. epsilon -- the unmarked::colext / MacKenzie (2003) convention.
.dm_fwd_R <- function(y_TJ, psi1, gam, eps, p) {
  Tt <- nrow(y_TJ)
  A  <- matrix(c(1 - gam, gam, eps, 1 - eps), 2, 2, byrow = TRUE)  # rows z_{t-1}=0,1
  emit <- function(yv, occ) {
    if (occ) prod(p^yv * (1 - p)^(1 - yv)) else as.numeric(all(yv == 0))
  }
  a <- c(1 - psi1, psi1) * c(emit(y_TJ[1, ], FALSE), emit(y_TJ[1, ], TRUE))
  for (t in 2:Tt) {
    b <- c(emit(y_TJ[t, ], FALSE), emit(y_TJ[t, ], TRUE))
    a <- b * as.numeric(crossprod(A, a))     # a_t(z) = b_t(z) sum_z' a_{t-1}(z') A(z',z)
  }
  sum(a)
}

test_that("dyn_occu marginal recursion is internally consistent", {
  set.seed(7)
  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p <- 0.6
  yTJ <- rbind(c(1, 0), c(0, 0), c(1, 1))    # one site, T = 3, J = 2
  m <- .dm_fwd_R(yTJ, psi1, gam, eps, p)
  expect_true(is.finite(log(m)) && m > 0 && m <= 1)
  y0 <- matrix(0L, 3, 2)
  m0 <- .dm_fwd_R(y0, psi1, gam, eps, p)
  expect_gt(m0, m)                           # all-zero is the most likely history
  expect_lt(m0, 1)
})

# Strong correctness anchor: the independent R forward recursion reproduces the
# reference implementation (unmarked::colext) log-likelihood to numerical
# precision at colext's own MLE. This pins the dynamic-occupancy marginal *formula*
# against a trusted external reference (analogous to dm_fwd_R vs
# cpp_dyn_abun_total_log_lik in test-dyn_abun.R), independent of how tulpaObs's EM
# optimises it.
test_that("dyn_occu marginal recursion matches unmarked::colext at its MLE", {
  skip_if_fast()
  skip_if_not_installed("unmarked")
  suppressPackageStartupMessages(library(unmarked))
  set.seed(101)
  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p_true <- 0.5
  n_sites <- 120L; T <- 4L; J <- 3L
  z <- matrix(NA_integer_, n_sites, T); z[, 1] <- rbinom(n_sites, 1, psi1)
  for (t in 2:T)
    z[, t] <- ifelse(z[, t - 1] == 1, rbinom(n_sites, 1, 1 - eps),
                     rbinom(n_sites, 1, gam))
  y <- array(0L, c(n_sites, J, T))
  for (i in seq_len(n_sites)) for (t in seq_len(T))
    y[i, , t] <- if (z[i, t]) rbinom(J, 1, p_true) else 0L

  Y <- matrix(0L, n_sites, J * T)            # unmarked layout: sites x (J*T), season-major
  for (t in seq_len(T)) Y[, (t - 1) * J + seq_len(J)] <- y[, , t]
  fm <- unmarked::colext(~ 1, ~ 1, ~ 1, ~ 1,
                         data = unmarked::unmarkedMultFrame(y = Y, numPrimary = T),
                         se = FALSE)
  uc <- unmarked::coef(fm)                    # logit psi, col(gamma), ext(epsilon), p
  ps <- plogis(uc[1]); gm <- plogis(uc[2]); ep <- plogis(uc[3]); pp <- plogis(uc[4])
  r_marg <- sum(vapply(seq_len(n_sites),
    function(i) log(.dm_fwd_R(t(y[i, , ]), ps, gm, ep, pp)), numeric(1)))
  expect_equal(r_marg, as.numeric(logLik(fm)), tolerance = 1e-5)
})

# Head-to-head coefficient equivalence with unmarked::colext. The exact
# forward-backward E-step (gcol33/tulpaObs#86) brings tulpaObs's dynamic-occupancy
# EM to colext's MLE (before the fix the EM converged ~3.4 logLik below it, with
# biased colonisation / extinction). The residual gap is the EM pseudo-count
# discretisation, so the tolerance is looser than the N-mixture-vs-pcount gate.
test_that("dyn_occu coefficients match unmarked::colext (tulpaObs#86)", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("unmarked")
  suppressPackageStartupMessages(library(unmarked))
  set.seed(101)
  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p_true <- 0.5
  n_sites <- 120L; T <- 4L; J <- 3L
  z <- matrix(NA_integer_, n_sites, T); z[, 1] <- rbinom(n_sites, 1, psi1)
  for (t in 2:T)
    z[, t] <- ifelse(z[, t - 1] == 1, rbinom(n_sites, 1, 1 - eps),
                     rbinom(n_sites, 1, gam))
  y <- array(0L, c(n_sites, J, T))
  for (i in seq_len(n_sites)) for (t in seq_len(T))
    y[i, , t] <- if (z[i, t]) rbinom(J, 1, p_true) else 0L

  fit <- tobs(~ 1, data = data.frame(idx = seq_len(n_sites)), family = dyn_occu(),
              detection = ~ 1, y = y, col_formula = ~ 1, ext_formula = ~ 1,
              priors = FALSE, control = list(verbose = FALSE))
  Y <- matrix(0L, n_sites, J * T)
  for (t in seq_len(T)) Y[, (t - 1) * J + seq_len(J)] <- y[, , t]
  fm <- unmarked::colext(~ 1, ~ 1, ~ 1, ~ 1,
                         data = unmarked::unmarkedMultFrame(y = Y, numPrimary = T),
                         se = FALSE)
  uc <- unmarked::coef(fm)   # logit psi1, colonisation (gamma), extinction (eps), p
  expect_lt(abs(fit$means[["psi1_(Intercept)"]]    - uc[1]), 0.1)
  expect_lt(abs(fit$means[["gamma_(Intercept)"]]   - uc[2]), 0.1)
  expect_lt(abs(fit$means[["epsilon_(Intercept)"]] - uc[3]), 0.1)
  expect_lt(abs(fit$means[["p_(Intercept)"]]       - uc[4]), 0.1)
})

# --------------------------------------------------------------------------- #
# Single-source int_occu: multi-seed POINT recovery.                            #
#                                                                               #
# NOTE (a finding this test surfaces): single-source int_occu's deterministic   #
# Laplace standard errors are severely overconfident -- `fit$sds` for the        #
# occupancy intercept is ~0.005 here (an order of magnitude or more too small    #
# for ~200 sites), so a per-seed SD-coverage gate fails even though the point    #
# estimate is unbiased. The recovery target is therefore the across-seed MEDIAN  #
# point estimate (robust to the uncalibrated SE); the SE-calibration gap is a    #
# separate kernel issue worth its own fix.                                       #
# --------------------------------------------------------------------------- #

test_that("single-source int_occu recovers (psi intercept, psi slope) across seeds", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 10L
  psi_int <- psi_slope <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_int_occu(N_total = 250, n_data = 1L, J = 5L,
                             beta_occ = c(0, 0.4), beta_det = list(c(0, -0.3)),
                             seed = 30L + s)
    fit <- tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1,
                y = sim$y, method = "laplace", control = list(verbose = FALSE))
    psi_int[s]   <- fit$means[["psi_(Intercept)"]]
    psi_slope[s] <- fit$means[["psi_x"]]
  }
  # Across-seed median recovers the occupancy intercept (true 0) and slope (0.4).
  expect_lt(abs(median(psi_int)), 0.15)
  expect_lt(abs(median(psi_slope) - 0.4), 0.15)
})
