setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))
site <- function(N, J, truth = 0.8, ns = 15) {
  e <- s <- numeric(0)
  for (k in seq_len(ns)) {
    sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                         beta_occ = c(0.5, 1.2), beta_det = c(0, truth),
                         seed = 9000L + k)
    f <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = occu(),
                       detection = ~ det_cov1, y = sim$y, method = "laplace",
                       control = list(verbose = FALSE)), error = function(e) NULL)
    if (!is.null(f)) { e <- c(e, f$means[["p_det_cov1"]]); s <- c(s, f$sds[["p_det_cov1"]]) }
  }
  cat(sprintf("SITE N=%d J=%2d  bias=%+.3f  ratio(se/sd)=%.2f  n=%d\n",
              N, J, mean(e) - truth, median(s) / sd(e), length(e))); flush.console()
}
for (J in c(4L, 8L, 16L)) site(400L, J)
site(2000L, 4L)
