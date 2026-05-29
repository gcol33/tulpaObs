suppressMessages({ devtools::load_all(".", quiet = TRUE) })

# Walk through what the user would see at the public surface.
set.seed(42)
N <- 60; J <- 3; ngrp <- 6
b <- stats::rnorm(ngrp, sd = 0.5); grp <- rep(seq_len(ngrp), length.out = N)
data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
eta_l <- as.numeric(model.matrix(~ x1, data) %*% c(0.5, 0.3)) + b[grp]
y <- matrix(NA_integer_, N, J)
for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, stats::rpois(1, exp(eta_l[i])), 0.5)

fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
            data = data, y = y, family = abun(),
            method = "laplace", verbose = FALSE,
            control = list(n.quad = 1))

cat("class:", paste(class(fit), collapse = ", "), "\n")
cat("\ncoef:\n"); print(coef(fit))
cat("\nvcov dim:", dim(vcov(fit)), "\n")
cat("\nranef:\n"); print(ranef(fit))
cat("\nsummary names (first 10):\n"); print(head(names(fit$means), 10))
cat("\nfitted (head):\n"); print(utils::head(fitted(fit)$lambda))
cat("\nintercepts:\n"); print(fit$intercepts)
cat("\npredict (head):\n")
pred <- predict(fit, X.0 = cbind(1, 0:3))
print(pred)
