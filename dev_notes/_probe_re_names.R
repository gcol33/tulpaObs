suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs"))
set.seed(11)
N <- 300; J <- 5
g <- factor(sample(1:25, N, replace = TRUE))
x <- rnorm(N); z <- rnorm(N)
y <- matrix(rbinom(N * J, 1, 0.4), N, J)
d <- data.frame(g = g, x = x, z = z)

dump_fit <- function(tag, f) {
  fit <- tobs(formula = f, data = d, y = y, detection = ~ 1,
              family = occu(), engine = "laplace", control = list(verbose = FALSE))
  cat("\n==== ", tag, " ====\n", sep = "")
  cat("names(means):\n"); print(names(fit$means))
  cat("re structure n_coefs:\n")
  if (!is.null(fit$re)) for (r in fit$re) cat("  type=", r$type,
    " intercept=", isTRUE(r$intercept), "\n", sep = "")
  rf <- tryCatch(ranef(fit), error = function(e) paste("ranef ERROR:", conditionMessage(e)))
  cat("ranef():\n"); print(rf)
}
dump_fit("(1+x+z|g) correlated", ~ (1 + x + z | g))
dump_fit("(0+x|g) slope-only",   ~ (0 + x | g))
cat("\n=== done ===\n")
