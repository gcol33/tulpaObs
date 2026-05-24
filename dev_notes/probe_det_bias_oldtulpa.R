# Decisive regression test: detection-slope attenuation at OLD tulpa vs current.
# priors=FALSE works on any tulpa (R-side gate), and at current tulpa
# pen ~= unpen, so the unpenalised bias is a clean cross-version comparison.
# If old tulpa recovers ~0.8 and current gives ~0.58, a tulpa commit caused the
# attenuation. If old tulpa also attenuates, it is inherent / pre-existing.
WT <- "C:/Users/Gilles Colling/Documents/dev/_tulpa_bisect"
suppressMessages(devtools::install(WT, quick = TRUE, build = FALSE,
                                   upgrade = FALSE, quiet = TRUE))
h <- system.file("include/tulpa/model_data.h", package = "tulpa")
cat("== tulpa:", trimws(grep("TULPA_ABI_VERSION =", readLines(h), value = TRUE)),
    "|", as.character(utils::packageDescription("tulpa")$Version), "==\n")

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))
truth <- 0.8; N <- 200L; J <- 4L; e <- numeric(0)
for (s in 1:12) {
  set.seed(4000L + s)
  x_occ <- rnorm(N); x_det <- rnorm(N)
  z <- rbinom(N, 1L, plogis(0.5 + 0.5 * x_occ))
  p_tru <- plogis(0.0 + truth * x_det)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
  d <- data.frame(x_occ = x_occ, x_det = x_det)
  f <- tryCatch(tobs(~ x_occ, data = d, family = occu(), detection = ~ x_det,
                     y = y, method = "laplace", priors = FALSE,
                     control = list(verbose = FALSE)),
                error = function(e) NULL)
  if (!is.null(f)) e <- c(e, f$means[["p_x_det"]])
}
cat(sprintf("n=%d  mean p_x_det=%.3f  bias=%.3f\n", length(e), mean(e), mean(e) - truth))
