pkgbuild::clean_dll("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", recompile = TRUE))
set.seed(11)
N <- 200; J <- 6
g <- factor(sample(1:20, N, replace = TRUE))
x <- rnorm(N); z <- rnorm(N)
# Occupancy with a real correlated random intercept+slope by group.
ng <- 20
Sig <- matrix(c(0.8^2, 0.3, 0.2,
                0.3,   0.6^2, 0.1,
                0.2,   0.1,   0.5^2), 3, 3)
B <- MASS::mvrnorm(ng, mu = c(0,0,0), Sigma = Sig)
gi <- as.integer(g)
eta <- 0.2 + B[gi,1] + B[gi,2]*x + B[gi,3]*z
psi <- plogis(eta); zocc <- rbinom(N, 1, psi)
p <- 0.5
y <- matrix(0L, N, J)
for (i in 1:N) y[i,] <- rbinom(J, 1, zocc[i]*p)
d <- data.frame(g = g, x = x, z = z)

ctl <- list(iter = 120, warmup = 60, seed = 1, verbose = FALSE)
for (tag in c("(1+x+z|g)", "(0+x|g)")) {
  f <- if (tag == "(1+x+z|g)") ~ (1 + x + z | g) else ~ (0 + x | g)
  cat("\n==== ", tag, " (NUTS) ====\n", sep = "")
  fit <- tryCatch(
    tobs(formula = f, data = d, y = y, detection = ~ 1,
         family = occu(), engine = "nuts", control = ctl),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) next
  cat("ncol(draws):", ncol(fit$draws), " n_samples:", fit$n_samples, "\n")
  cat("names(means):\n"); print(names(fit$means))
  cat("colnames(draws):\n"); print(colnames(fit$draws))
}
cat("\n=== done ===\n")
