# verify_sla_1d.R
#
# Numerical sanity check for the simplified-Laplace 1D formula derived in
# dev_notes/simplified_laplace_derivation.md §2.4.
#
# Setup: intercept-only Bernoulli model.
#   x ~ N(0, tau^2)
#   y | x ~ Binomial(n, plogis(x))
#
# Predicted (from derivation):
#   gamma_SLA = -sigma_post^3 * n * p_star * (1 - p_star) * (1 - 2 * p_star)
#
# Compared against: direct numerical computation of the posterior skewness
# via fine-grid integration on the unnormalised posterior.

verify_one <- function(tau, n, y) {
  log_post <- function(x) {
    dnorm(x, 0, tau, log = TRUE) +
      dbinom(y, n, plogis(x), log = TRUE)
  }
  opt <- optimise(log_post, c(-20, 20), maximum = TRUE)
  x_star <- opt$maximum
  p_star <- plogis(x_star)
  Q <- 1 / tau^2 + n * p_star * (1 - p_star)
  sigma_post <- sqrt(1 / Q)

  l3 <- -n * p_star * (1 - p_star) * (1 - 2 * p_star)
  gamma_sla <- sigma_post^3 * l3

  # Direct posterior moments via fine grid
  grid <- seq(x_star - 10 * sigma_post, x_star + 10 * sigma_post, length.out = 4001)
  lp <- vapply(grid, log_post, numeric(1))
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  m1 <- sum(w * grid)
  m2 <- sum(w * (grid - m1)^2)
  m3 <- sum(w * (grid - m1)^3)
  gamma_exact <- m3 / m2^1.5

  list(
    x_star = x_star, p_star = p_star, sigma_post = sigma_post,
    gamma_sla = gamma_sla, gamma_exact = gamma_exact,
    rel_err = (gamma_sla - gamma_exact) / abs(gamma_exact)
  )
}

# Grid of test cases
cases <- expand.grid(tau = c(1, 2, 5), n = c(5, 20, 100), y = c(0, 1, 5, 10))
cases <- cases[cases$y <= cases$n, ]

res <- do.call(rbind, lapply(seq_len(nrow(cases)), function(k) {
  r <- verify_one(cases$tau[k], cases$n[k], cases$y[k])
  data.frame(
    tau = cases$tau[k], n = cases$n[k], y = cases$y[k],
    p_star = round(r$p_star, 3),
    sigma_post = round(r$sigma_post, 3),
    gamma_sla = round(r$gamma_sla, 4),
    gamma_exact = round(r$gamma_exact, 4),
    rel_err = round(r$rel_err, 3)
  )
}))

print(res)

cat("\nMax |rel_err| where |gamma_exact| > 0.05:",
    max(abs(res$rel_err[abs(res$gamma_exact) > 0.05])), "\n")
cat("Sign matches in", sum(sign(res$gamma_sla) == sign(res$gamma_exact)),
    "of", nrow(res), "cases\n")
