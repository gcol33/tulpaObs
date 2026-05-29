# Confirm n_quad is now threaded; small Poisson example.
suppressMessages({ devtools::load_all(".", quiet = TRUE) })

set.seed(2026)
N <- 100; J <- 4; ngrp <- 10
b <- stats::rnorm(ngrp, sd = 0.6); grp <- rep(seq_len(ngrp), length.out = N)
data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
eta_l <- as.numeric(model.matrix(~ x1, data) %*% c(0.5, 0.3)) + b[grp]
Nlat <- stats::rpois(N, exp(eta_l))
y <- matrix(NA_integer_, N, J)
for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], 0.5)

for (nq in c(1L, 5L, 9L)) {
  t0 <- Sys.time()
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
              data = data, y = y, family = abun(),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = nq))
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  cat(sprintf("n_quad = %d  | sigma_g = %.3f | beta_l = (%.3f, %.3f) | %.2fs\n",
              nq, fit$means["sigma_g1_(Intercept)"],
              fit$means["lambda_(Intercept)"], fit$means["lambda_x1"], el))
}
