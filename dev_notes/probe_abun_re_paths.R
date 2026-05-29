# Exercise all four cells of (arm in {lambda, p}) x (mixture in {P, NB}) +
# the n_quad > 1 sigma debias, for #13.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

sim_lambda_re <- function(N = 120, J = 4, ngrp = 12, sigma = 0.6, seed = 1) {
  set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b <- stats::rnorm(ngrp, sd = sigma)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  beta_lambda <- c(0.5, 0.3); beta_p <- 0
  eta_l <- as.numeric(model.matrix(~ x1, data) %*% beta_lambda) + b[grp]
  Nlat  <- stats::rpois(N, exp(eta_l))
  p <- plogis(beta_p)
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p)
  list(y = y, data = data, beta_lambda = beta_lambda, beta_p = beta_p,
       sigma = sigma, b = b)
}

sim_p_re <- function(N = 120, J = 4, ngrp = 12, sigma = 0.6, seed = 1) {
  set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b <- stats::rnorm(ngrp, sd = sigma)
  data <- data.frame(z1 = stats::rnorm(N), g = factor(grp))
  beta_lambda <- 1.5; beta_p <- c(0.0, 0.4)
  lambda <- exp(beta_lambda)
  Nlat   <- stats::rpois(N, lambda)
  eta_p  <- as.numeric(model.matrix(~ z1, data) %*% beta_p) + b[grp]
  p <- plogis(eta_p)
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p[i])
  list(y = y, data = data, beta_lambda = beta_lambda, beta_p = beta_p,
       sigma = sigma, b = b)
}

cat("\n--- Poisson, lambda-arm RE, n_quad = 1 ---\n")
s <- sim_lambda_re(seed = 7)
fit1 <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
             data = s$data, y = s$y, family = abun(),
             method = "laplace", verbose = FALSE,
             control = list(n.quad = 1))
print(round(fit1$means[1:4], 3))
cat("sigma truth =", s$sigma, "  estimated =",
    round(fit1$means["sigma_g1_(Intercept)"], 3), "\n")

cat("\n--- Poisson, lambda-arm RE, n_quad = 7 (AGHQ debias) ---\n")
fit2 <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
             data = s$data, y = s$y, family = abun(),
             method = "laplace", verbose = FALSE,
             control = list(n.quad = 7))
print(round(fit2$means[1:4], 3))
cat("sigma estimated =",
    round(fit2$means["sigma_g1_(Intercept)"], 3), "\n")

cat("\n--- Poisson, p-arm RE ---\n")
sp <- sim_p_re(seed = 13)
fitp <- tobs(formula = ~ 1, detection = ~ z1 + (1 | g),
             data = sp$data, y = sp$y, family = abun(),
             method = "laplace", verbose = FALSE,
             control = list(n.quad = 5))
print(round(fitp$means[1:5], 3))
cat("sigma_p truth =", sp$sigma, "  estimated =",
    round(fitp$means["sigma_p1_(Intercept)"], 3), "\n")

cat("\n--- NB, lambda-arm RE ---\n")
fitnb <- tryCatch(
  tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
       data = s$data, y = s$y, family = abun(mixture = "negbin"),
       method = "laplace", verbose = FALSE,
       control = list(n.quad = 5)),
  error = function(e) { cat("NB fit error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fitnb)) {
  print(round(fitnb$means[1:5], 3))
  cat("mixture:", fitnb$mixture, "\n")
  cat("dispersion r:", round(fitnb$nmix_dispersion$r %||% NA, 3), "\n")
}

cat("\n--- Error gates ---\n")
ge <- tryCatch(
  tobs(formula = ~ 1 + (1 | g) + icar(),
       detection = ~ 1, data = s$data, y = s$y, family = abun(),
       method = "nested_laplace", verbose = FALSE),
  error = function(e) conditionMessage(e))
cat("RE + spatial -> ", substr(ge, 1, 80), "\n")
