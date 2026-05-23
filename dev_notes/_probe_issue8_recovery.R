suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs"))
set.seed(7)
N <- 300; J <- 6
psi <- plogis(0.4); p0 <- 0.2; p1 <- 1.2
zocc <- rbinom(N, 1, psi)
eff <- matrix(rnorm(N * J), N, J)
y <- matrix(0L, N, J)
for (i in 1:N) for (j in 1:J) y[i, j] <- rbinom(1, 1, zocc[i] * plogis(p0 + p1 * eff[i, j]))

df <- data.frame(
  site_id = rep(seq_len(N), each = J),
  visit   = rep(seq_len(J), times = N),
  occur   = as.vector(t(y)),
  effort  = as.vector(t(eff))
)
od <- tobs_data(df, y = "occur", site = "site_id", visit = "visit",
                det.covs = "effort")
fit <- tobs(~ 1, data = data.frame(site_id = unique(df$site_id)),
            y = od$y, detection = ~ effort, visit_data = od$det.covs,
            family = occu(), engine = "laplace", control = list(verbose = FALSE))

cat("names(means):\n"); print(names(fit$means))
cat("means:\n"); print(round(fit$means, 3))
cat("sds:\n"); print(round(fit$sds, 3))
cat("\np1 truth =", p1, " p0 truth =", p0, "\n")
cat("=== done ===\n")
