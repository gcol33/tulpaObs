# Single-process NUTS check, immune to the concurrent tulpa-reinstall race:
# install tulpa (ABI 24) from the bisect worktree, recompile+load tulpaObs
# against it, and fit a NUTS RE model -- all before any other process can swap
# the on-disk tulpa. Distinguishes a genuine NUTS segfault at matched ABI 24
# from the between-process ABI-mismatch artifact.
suppressMessages({
  try(remove.packages("tulpa"), silent = TRUE)
  devtools::install("C:/Users/Gilles Colling/Documents/dev/_tulpa_bisect",
                    quick = TRUE, build = FALSE, upgrade = FALSE, quiet = TRUE)
})
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", recompile = TRUE, quiet = TRUE))
cat("in-process tulpa ABI:",
    trimws(grep("ABI_VERSION =", readLines(system.file(
      "include/tulpa/model_data.h", package = "tulpa")), value = TRUE)), "\n")
set.seed(1)
N <- 40L; J <- 3L
g <- factor(rep(1:8, length.out = N)); h <- factor(rep(1:4, length.out = N))
x <- rnorm(N)
z <- rbinom(N, 1, plogis(0.2 + 0.5 * x))
y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * 0.6)
d <- data.frame(g = g, h = h, x = x)
cat("fitting NUTS, two re() terms...\n"); flush.console()
fit <- tryCatch(
  tobs(~ x + (1 | g) + (1 | h), data = d, y = y, detection = ~ 1,
       family = occu(), method = "nuts",
       control = list(n.iter = 50, n.warmup = 25, seed = 1, verbose = FALSE)),
  error = function(e) { cat("R-ERROR:", conditionMessage(e), "\n"); NULL })
cat("RESULT:", if (is.null(fit)) "NULL" else "tobs_fit OK", "\n")
