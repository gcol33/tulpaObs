# Is the detection-arm attenuation + SE underestimate small-sample (estimator
# consistent, tests under-powered) or structural (estimator bug)? Sweep J (the
# binding axis for detection identifiability) at fixed-ish N and report, for the
# detection slope: mean estimate vs truth (bias) and median(SE)/sd(est) (ratio).
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))

site_run <- function(N, J, truth = 0.8, nseed = 20) {
  e <- s <- numeric(0)
  for (k in seq_len(nseed)) {
    sim <- simulate_occu(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
                         beta_occ = c(0.5, 1.2), beta_det = c(0.0, truth),
                         seed = 9000L + k)
    f <- tryCatch(tobs(~ occ_cov1, data = sim$data, family = occu(),
                       detection = ~ det_cov1, y = sim$y, method = "laplace",
                       control = list(verbose = FALSE)), error = function(e) NULL)
    if (!is.null(f)) { e <- c(e, f$means[["p_det_cov1"]]); s <- c(s, f$sds[["p_det_cov1"]]) }
  }
  cat(sprintf("  site  N=%-5d J=%-3d  bias=%+.3f  ratio(se/sd)=%.2f  (n=%d)\n",
              N, J, mean(e) - truth, median(s) / sd(e), length(e)))
}
visit_run <- function(N, J, p1 = 1.2, nseed = 20) {
  e <- s <- numeric(0)
  for (k in seq_len(nseed)) {
    set.seed(5000L + k)
    zocc <- rbinom(N, 1, plogis(0.4)); eff <- matrix(rnorm(N * J), N, J)
    y <- matrix(0L, N, J)
    for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, zocc[i] * plogis(0.2 + p1 * eff[i, ]))
    df <- data.frame(site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
                     occur = as.vector(t(y)), effort = as.vector(t(eff)))
    od <- tobs_data(df, y = "occur", site = "site_id", visit = "visit", det.covs = "effort")
    f <- tryCatch(tobs(~ 1, data = data.frame(site_id = unique(df$site_id)), y = od$y,
                       detection = ~ effort, visits = od$det.covs, family = occu(),
                       method = "laplace", control = list(verbose = FALSE)),
                  error = function(e) NULL)
    if (!is.null(f)) { e <- c(e, f$means[["p_visit_effort"]]); s <- c(s, f$sds[["p_visit_effort"]]) }
  }
  cat(sprintf("  visit N=%-5d J=%-3d  bias=%+.3f  ratio(se/sd)=%.2f  (n=%d)\n",
              N, J, mean(e) - p1, median(s) / sd(e), length(e)))
}
cat("SITE-LEVEL detection slope (truth 0.8):\n")
for (J in c(6L, 12L, 24L)) site_run(600L, J)
site_run(2000L, 6L)
cat("VISIT-LEVEL detection slope (truth 1.2):\n")
for (J in c(6L, 12L, 24L)) visit_run(300L, J)
visit_run(1500L, 6L)
