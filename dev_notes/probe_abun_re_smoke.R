# Smoke test for gcol33/tulpaObs#13: site-level RE on N-mixture abundance.
# Confirm the path returns a tobs_fit and the coefficients are in range.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

set.seed(2026)
N <- 100; J <- 4
ngrp <- 10
sigma_true <- 0.6
b_true <- stats::rnorm(ngrp, sd = sigma_true)
grp <- rep(seq_len(ngrp), length.out = N)
data <- data.frame(
  x1 = stats::rnorm(N),
  x2 = stats::rnorm(N),
  g  = factor(grp)
)

# True abundance arm: log lambda = 0.5 + 0.3*x1 + b[g]
beta_lambda_true <- c(0.5, 0.3)
eta_lambda <- as.numeric(model.matrix(~ x1, data) %*% beta_lambda_true) + b_true[grp]
lambda <- exp(eta_lambda)
Nlat <- stats::rpois(N, lambda)
# Detection arm: logit p = 0
p <- rep(0.5, N)
y <- matrix(NA_integer_, N, J)
for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p[i])

cat("simulated truth: beta_lambda =", beta_lambda_true,
    "  sigma_true =", sigma_true, "\n")
cat("range of lambda:", range(lambda), "\n")

t0 <- Sys.time()
fit <- tryCatch(
  tobs(formula   = ~ x1 + (1 | g),
       detection = ~ 1,
       data = data, y = y,
       family = abun(),
       method = "laplace", verbose = FALSE),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
t1 <- Sys.time()

if (!is.null(fit)) {
  cat("\n== fit ok in", round(as.numeric(t1 - t0), 2), "s ==\n")
  cat("class:", paste(class(fit), collapse = ", "), "\n")
  cat("means:\n"); print(round(fit$means, 3))
  cat("sds:\n"); print(round(fit$sds, 3))
  cat("converged:", fit$convergence$converged, "\n")
  if (!is.null(fit$nmix_re))   cat("nmix_re arm:", fit$nmix_re$arm, "\n")
  if (!is.null(fit$re_effects)) {
    cat("re_effects head:\n"); print(utils::head(fit$re_effects[[1]]))
  }
}
