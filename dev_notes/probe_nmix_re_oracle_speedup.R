# Benchmark the native NMixGroupedOracle path in .tobs_nmix_re_aghq() against
# the documented R-closure baseline (~30-90 s per fit at this scale, per the
# pre-oracle test-abun-re.R header). One fit each at the suite's recovery
# scale (N = 100, J = 4, 10 groups), Poisson + NB, on both arms.

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
devtools::load_all(quiet = TRUE)

set.seed(1234)
N    <- 100L
J    <- 4L
ngrp <- 10L

# --- lambda arm, Poisson ---
sim_l <- local({
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = 0.6)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  eta_l  <- as.numeric(model.matrix(~ x1, data) %*% c(0.6, 0.4)) + b[grp]
  Nlat   <- stats::rpois(N, exp(eta_l))
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], plogis(0))
  list(y = y, data = data)
})

# --- p arm, Poisson ---
sim_p <- local({
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = 0.6)
  data <- data.frame(z1 = stats::rnorm(N), g = factor(grp))
  eta_p  <- as.numeric(model.matrix(~ z1, data) %*% c(0, 0.4)) + b[grp]
  Nlat   <- stats::rpois(N, exp(1.5))
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], plogis(eta_p[i]))
  list(y = y, data = data)
})

# --- lambda arm, NB ---
sim_nb <- local({
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = 0.5)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  eta_l  <- as.numeric(model.matrix(~ x1, data) %*% c(0.6, 0.3)) + b[grp]
  Nlat   <- stats::rnbinom(N, size = 4, mu = exp(eta_l))
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], plogis(0))
  list(y = y, data = data)
})

run <- function(label, expr) {
  t <- system.time(fit <- eval(expr, envir = parent.frame()))
  cat(sprintf("%-36s  user=%.2fs  elapsed=%.2fs  converged=%s\n",
              label, t[["user.self"]], t[["elapsed"]],
              isTRUE(fit$convergence$converged)))
  invisible(fit)
}

cat("---- native NMixGroupedOracle path, n.quad = 5 ----\n")
run("lambda  Poisson  n.quad=5", quote(
  tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
       data = sim_l$data, y = sim_l$y, family = abun(),
       method = "laplace", verbose = FALSE,
       control = list(n.quad = 5))))
run("p       Poisson  n.quad=5", quote(
  tobs(formula = ~ 1, detection = ~ z1 + (1 | g),
       data = sim_p$data, y = sim_p$y, family = abun(),
       method = "laplace", verbose = FALSE,
       control = list(n.quad = 5))))
run("lambda  NB       n.quad=3", quote(
  tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
       data = sim_nb$data, y = sim_nb$y,
       family = abun(mixture = "negbin"),
       method = "laplace", verbose = FALSE,
       control = list(n.quad = 3))))

cat("---- n.quad = 1 (joint Laplace, fastest setting) ----\n")
run("lambda  Poisson  n.quad=1", quote(
  tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
       data = sim_l$data, y = sim_l$y, family = abun(),
       method = "laplace", verbose = FALSE,
       control = list(n.quad = 1))))
