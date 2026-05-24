setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(".", quiet = TRUE))

# occu-prior config: site-level detection, N=200, J=4, p slope truth 0.8.
e <- s <- numeric(0)
for (k in 1:20) {
  set.seed(4000L + k)
  N <- 200L; J <- 4L
  xo <- rnorm(N); xd <- rnorm(N)
  z <- rbinom(N, 1, plogis(0.5 + 0.5 * xo)); ptru <- plogis(0.0 + 0.8 * xd)
  y <- matrix(0L, N, J); for (i in 1:N) if (z[i]) y[i, ] <- rbinom(J, 1, ptru[i])
  d <- data.frame(x_occ = xo, x_det = xd)
  f <- tryCatch(tobs(~ x_occ, data = d, family = occu(), detection = ~ x_det,
                     y = y, method = "laplace", control = list(verbose = FALSE)),
                error = function(e) NULL)
  if (!is.null(f)) { e <- c(e, f$means[["p_x_det"]]); s <- c(s, f$sds[["p_x_det"]]) }
}
cov <- mean(abs(e - 0.8) < 1.96 * s)
cat(sprintf("SITE  N=200 J=4: n=%d bias=%+.3f covSE=%.2f coverage=%.2f\n",
            length(e), mean(e) - 0.8, median(s) / sd(e), cov))

# issue8 config: visit-level detection, N=300, J=6, slope truth 1.2.
e <- s <- numeric(0)
for (k in 1:20) {
  set.seed(5000L + k)
  N <- 300L; J <- 6L
  z <- rbinom(N, 1, plogis(0.4)); eff <- matrix(rnorm(N * J), N, J)
  y <- matrix(0L, N, J); for (i in 1:N) y[i, ] <- rbinom(J, 1, z[i] * plogis(0.2 + 1.2 * eff[i, ]))
  df <- data.frame(site_id = rep(1:N, each = J), visit = rep(1:J, times = N),
                   occur = as.vector(t(y)), effort = as.vector(t(eff)))
  od <- tobs_data(df, y = "occur", site = "site_id", visit = "visit", det.covs = "effort")
  f <- tryCatch(tobs(~ 1, data = data.frame(site_id = unique(df$site_id)), y = od$y,
                     detection = ~ effort, visits = od$det.covs, family = occu(),
                     method = "laplace", control = list(verbose = FALSE)),
                error = function(e) NULL)
  if (!is.null(f)) { e <- c(e, f$means[["p_visit_effort"]]); s <- c(s, f$sds[["p_visit_effort"]]) }
}
cov <- mean(abs(e - 1.2) < 1.96 * s)
cat(sprintf("VISIT N=300 J=6: n=%d bias=%+.3f covSE=%.2f coverage=%.2f\n",
            length(e), mean(e) - 1.2, median(s) / sd(e), cov))
