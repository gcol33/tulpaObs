setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))
N <- 600L; J <- 6L; e <- s <- numeric(0)
for (k in 1:20) {
  sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = c(0.5, 1.2), beta_det = c(0, 0.8), seed = 7000L + k)
  f <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = occu(),
                     detection = ~ det_cov1, y = sim$y, method = "laplace",
                     control = list(verbose = FALSE)), error = function(e) NULL)
  if (!is.null(f)) { e <- c(e, f$means[["p_det_cov1"]]); s <- c(s, f$sds[["p_det_cov1"]]) }
}
cat("p_det_cov1 estimates (truth 0.8):\n"); print(round(sort(e), 3))
cat(sprintf("\nmean=%.3f  sd=%.3f  mad=%.3f  IQR/1.349=%.3f\n",
            mean(e), sd(e), mad(e), IQR(e) / 1.349))
cat(sprintf("median SE=%.3f\n", median(s)))
cat(sprintf("ratio (medSE/sd)=%.2f   ratio (medSE/mad)=%.2f   ratio(medSE/robustsd)=%.2f\n",
            median(s) / sd(e), median(s) / mad(e), median(s) / (IQR(e) / 1.349)))
