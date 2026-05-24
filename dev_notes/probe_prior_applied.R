# Is the penalised default actually breaking the psi-p ridge, or is the prior a
# no-op? Fit the occu-prior test's hard case (N=200, J=4, p_x_det truth 0.8)
# with the prior ON (default) vs OFF (priors=FALSE) and print the detection
# slope + SE for both. If pen == unpen, the prior is not reaching the fit.
WT <- "C:/Users/Gilles Colling/Documents/dev/_tulpa_bisect"
suppressMessages(devtools::install(WT, quick = TRUE, build = FALSE,
                                   upgrade = FALSE, quiet = TRUE))
h <- system.file("include/tulpa/model_data.h", package = "tulpa")
cat("== tulpa:", trimws(grep("TULPA_ABI_VERSION =", readLines(h), value = TRUE)), "==\n")

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))

cat("\n-- default occu_priors() --\n"); print(occu_priors())

truth <- 0.8; N <- 200L; J <- 4L
ep <- eu <- sp <- su <- numeric(0)
for (s in 1:12) {
  set.seed(4000L + s)
  x_occ <- rnorm(N); x_det <- rnorm(N)
  z <- rbinom(N, 1L, plogis(0.5 + 0.5 * x_occ))
  p_tru <- plogis(0.0 + truth * x_det)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
  d <- data.frame(x_occ = x_occ, x_det = x_det)
  fp <- tobs(~ x_occ, data = d, family = occu(), detection = ~ x_det, y = y,
             method = "laplace", control = list(verbose = FALSE))
  fu <- tobs(~ x_occ, data = d, family = occu(), detection = ~ x_det, y = y,
             method = "laplace", priors = FALSE, control = list(verbose = FALSE))
  ep <- c(ep, fp$means[["p_x_det"]]); sp <- c(sp, fp$sds[["p_x_det"]])
  eu <- c(eu, fu$means[["p_x_det"]]); su <- c(su, fu$sds[["p_x_det"]])
}
cat(sprintf("\npenalised : mean est=%.3f (bias %.3f)  mean se=%.3f\n",
            mean(ep), mean(ep) - truth, mean(sp)))
cat(sprintf("unpenalised: mean est=%.3f (bias %.3f)  mean se=%.3f\n",
            mean(eu), mean(eu) - truth, mean(su)))
cat(sprintf("pen==unpen identical? %s  (max|diff|=%.4g)\n",
            isTRUE(all.equal(ep, eu)), max(abs(ep - eu))))
