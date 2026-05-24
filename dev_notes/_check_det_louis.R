setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))
N <- 600L; J <- 6L; truth <- c(0.5, 1.2, 0.0, 0.8)
nms <- c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)", "p_det_cov1")
bh <- sh <- matrix(NA_real_, 20L, 4L, dimnames = list(NULL, nms))
for (s in 1:20) {
  sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                       beta_occ = truth[1:2], beta_det = truth[3:4], seed = 7000L + s)
  f <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = occu(),
                     detection = ~ det_cov1, y = sim$y, method = "laplace",
                     control = list(verbose = FALSE)), error = function(e) NULL)
  if (!is.null(f)) { bh[s, ] <- f$means[nms]; sh[s, ] <- f$sds[nms] }
}
keep <- complete.cases(bh) & complete.cases(sh)
ratio <- apply(sh[keep, ], 2, median) / apply(bh[keep, ], 2, sd)
cat("n =", sum(keep), "\n")
for (nm in nms) cat(sprintf("  %-18s ratio(se/sd) = %.2f\n", nm, ratio[nm]))
