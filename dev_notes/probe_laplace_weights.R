# Localize: does tulpa_laplace() apply observation `weights` to the FIXED-effect
# block? The detection M-step relies on it. Compare against glm(weights=).
suppressMessages(library(tulpa))
set.seed(7)
cat("tulpa:", as.character(utils::packageVersion("tulpa")), "\n\n")

## A) Integer y, per-row weights, NO RE, NO duplication.
N <- 500L
X <- cbind(1, rnorm(N), rnorm(N))
beta <- c(0.5, -0.8, 0.6)
p <- plogis(X %*% beta)
y <- rbinom(N, 1, p)
wt <- runif(N, 0.2, 1)                       # nontrivial weights
fitA <- tulpa::tulpa_laplace(y = y, n_trials = rep(1, N), X = X,
                             weights = wt, family = "binomial")
glmA <- glm(y ~ X[, 2] + X[, 3], weights = wt, family = binomial)
cat("A) integer y + weights, no RE:\n")
cat("   tulpa:", round(fitA$mode[1:3], 4), "\n")
cat("   glm  :", round(unname(coef(glmA)), 4), "\n")
cat("   max abs diff:", round(max(abs(fitA$mode[1:3] - coef(glmA))), 5), "\n\n")

## B) Same data but NO weights (all 1) -> sanity that betas recover.
fitB <- tulpa::tulpa_laplace(y = y, n_trials = rep(1, N), X = X,
                             family = "binomial")
glmB <- glm(y ~ X[, 2] + X[, 3], family = binomial)
cat("B) integer y, no weights:\n")
cat("   tulpa:", round(fitB$mode[1:3], 4), "\n")
cat("   glm  :", round(unname(coef(glmB)), 4), "\n")
cat("   max abs diff:", round(max(abs(fitB$mode[1:3] - coef(glmB))), 5), "\n\n")

## C) Aggregated binomial: y = successes, n = trials (weights via counts).
## This is the exact weighted-Bernoulli WITHOUT a `weights` arg: encode the
## occupancy objective as two binomial rows per site with n_trials carrying
## an integer pseudo-count? No -- test plain integer binomial recovers beta.
ng <- 1L
agg_n <- rep(10L, N)
agg_y <- rbinom(N, agg_n, p)
fitC <- tulpa::tulpa_laplace(y = agg_y, n_trials = agg_n, X = X,
                             family = "binomial")
glmC <- glm(cbind(agg_y, agg_n - agg_y) ~ X[, 2] + X[, 3], family = binomial)
cat("C) integer binomial counts (n_trials=10):\n")
cat("   tulpa:", round(fitC$mode[1:3], 4), "\n")
cat("   glm  :", round(unname(coef(glmC)), 4), "\n")
cat("   max abs diff:", round(max(abs(fitC$mode[1:3] - coef(glmC))), 5), "\n")

cat("\n=== done ===\n")
