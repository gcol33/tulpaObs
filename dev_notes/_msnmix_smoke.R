suppressMessages(library(spAbundance))
cat("sim fns:", paste(grep("sim", ls("package:spAbundance"), value = TRUE, ignore.case = TRUE), collapse = ", "), "\n")

set.seed(1)
S <- 4; J <- 30; nrep <- 4
# y: species x sites x visits
lam <- exp(0.8); p <- 0.5
y <- array(NA_integer_, c(S, J, nrep))
for (s in 1:S) { N <- rpois(J, lam); for (i in 1:J) y[s, i, ] <- rbinom(nrep, N[i], p) }
xa <- rnorm(J); xd <- rnorm(J)
data <- list(y = y, abund.covs = data.frame(xa = xa), det.covs = list(xd = xd))

fit <- msNMix(abund.formula = ~ xa, det.formula = ~ xd, data = data,
              n.batch = 20, batch.length = 25, family = "Poisson",
              n.burn = 200, n.thin = 1, n.chains = 1, verbose = FALSE)
cat("\n=== class ===\n"); print(class(fit))
cat("\n=== names ===\n"); print(names(fit))
cat("\n=== str of *.samples ===\n")
for (nm in grep("samples$", names(fit), value = TRUE)) {
  obj <- fit[[nm]]
  cat(sprintf("%-26s dim=%s  colnames=%s\n", nm,
              paste(dim(as.matrix(obj)), collapse = "x"),
              paste(head(colnames(as.matrix(obj)), 8), collapse = ",")))
}
