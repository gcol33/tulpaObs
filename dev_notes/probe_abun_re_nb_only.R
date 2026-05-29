# Isolate the NB + RE path; print at each stage to find the hang.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

set.seed(7)
N <- 120; J <- 4; ngrp <- 12
b <- stats::rnorm(ngrp, sd = 0.6)
grp <- rep(seq_len(ngrp), length.out = N)
data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
eta_l <- as.numeric(model.matrix(~ x1, data) %*% c(0.5, 0.3)) + b[grp]
Nlat  <- stats::rnbinom(N, size = 4, mu = exp(eta_l))
p     <- 0.5
y <- matrix(NA_integer_, N, J)
for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p)

cat("[1] warm-start NB nmix_laplace (no RE)...\n"); flush.console()
t0 <- Sys.time()
warm <- nmix_laplace(y = as.integer(t(y)), site_idx = rep(seq_len(N), each = J),
                     X_lambda = model.matrix(~ x1, data),
                     X_p = matrix(1, N * J, 1),
                     mixture = "NB", K_max = max(y) + 100L, verbose = FALSE)
cat("    warm took", round(difftime(Sys.time(), t0, units = "secs"), 2),
    "s; beta_lambda =", round(warm$beta_lambda, 3),
    "; r =", round(warm$r, 3), "\n"); flush.console()

cat("[2] full NB + RE fit...\n"); flush.console()
t0 <- Sys.time()
fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
            data = data, y = y, family = abun(mixture = "negbin"),
            method = "laplace", verbose = TRUE,
            control = list(n.quad = 1))
cat("    NB+RE fit took", round(difftime(Sys.time(), t0, units = "secs"), 2),
    "s\n")
print(round(fit$means, 3))
cat("mixture:", fit$mixture, "  r:", fit$nmix_dispersion$r, "\n")
