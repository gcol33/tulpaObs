# Minimal repro: does tulpa::tulpa_laplace regress on a plain (no-RE) weighted
# binomial GLM -- the call tulpaObs's detection M-step makes -- after the tulpa
# ABI-24 upgrade? Compare against stats::glm (the exact MLE). If tulpa's
# coefficients diverge from glm, the no-RE detection bias seen in tulpaObs
# recovery tests is an upstream tulpa regression, not the RE work.
library(tulpa)
set.seed(1)
N <- 4000
x <- rnorm(N)
b <- c(-0.3, 0.8)               # intercept, slope
n <- rep(5L, N)                 # trials
p <- plogis(b[1] + b[2] * x)
y <- rbinom(N, n, p)
X <- cbind(1, x)

g <- glm(cbind(y, n - y) ~ x, family = binomial())
fo <- tulpa_laplace(y = y, n_trials = n, X = X, family = "binomial",
                    return_hessian = TRUE)

cat("truth          :", round(b, 4), "\n")
cat("glm   coef     :", round(unname(coef(g)), 4), "\n")
cat("tulpa coef     :", round(fo$mode[1:2], 4), "\n")
cat("abs diff t-glm :", round(abs(fo$mode[1:2] - unname(coef(g))), 5), "\n")

# Weighted version (detection M-step passes fractional weights).
w <- runif(N, 0.3, 1)
gw <- glm(cbind(y, n - y) ~ x, family = binomial(), weights = w)
fw <- tulpa_laplace(y = y, n_trials = n, X = X, weights = w,
                    family = "binomial", return_hessian = TRUE)
cat("\nweighted glm   :", round(unname(coef(gw)), 4), "\n")
cat("weighted tulpa :", round(fw$mode[1:2], 4), "\n")
cat("abs diff       :", round(abs(fw$mode[1:2] - unname(coef(gw))), 5), "\n")
