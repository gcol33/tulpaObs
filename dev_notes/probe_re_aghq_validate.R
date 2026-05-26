# Validation of the AGHQ debias.
#  (1) n.quad = 1 must reproduce the EM (Laplace) fit -- the quadrature reduces
#      to the Laplace approximation at the mode.
#  (2) Multi-seed mean sigma: EM vs AGHQ vs NUTS at a small cluster size, to
#      check AGHQ reduces the attenuation toward the NUTS (calibrated) value.
devtools::load_all(".", quiet = TRUE)
suppressMessages(library(stats))

sim_re_int <- function(seed, ng, per, J = 6L, b0 = 0.3, b1 = -0.6,
                       sigma = 0.9, p = 0.45) {
  set.seed(seed)
  N <- ng * per; g <- rep(seq_len(ng), each = per); x <- rnorm(N)
  b <- rnorm(ng, 0, sigma)
  z <- rbinom(N, 1, plogis(b0 + b1 * x + b[g]))
  y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x))
}
sig1 <- function(f) f$means[[grep("^sigma_", names(f$means), value = TRUE)[1]]]

fit_one <- function(s, aghq, nq = 9L, method = "laplace") {
  ctrl <- list(verbose = FALSE)
  if (method == "laplace") { ctrl$re.aghq <- aghq; ctrl$n.quad <- nq }
  if (method == "nuts") ctrl <- list(n.iter = 600, n.warmup = 300, seed = 1, verbose = FALSE)
  tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
       family = occu(), method = method, control = ctrl)
}

cat("== (1) n.quad=1 AGHQ vs EM (should match) ==\n")
s <- sim_re_int(seed = 11, ng = 30L, per = 12L)
f_em <- fit_one(s, aghq = FALSE)
f_q1 <- fit_one(s, aghq = TRUE, nq = 1L)
cat(sprintf("  EM sigma=%.4f   AGHQ(Q=1) sigma=%.4f   diff=%.4f\n",
            sig1(f_em), sig1(f_q1), abs(sig1(f_em) - sig1(f_q1))))
cat(sprintf("  EM psi_x=%.4f   AGHQ(Q=1) psi_x=%.4f\n",
            f_em$means[["psi_x"]], f_q1$means[["psi_x"]]))

cat("\n== (2) mean sigma over seeds (truth 0.9), per=8, ng=30 ==\n")
seeds <- 1:15
res <- sapply(seeds, function(sd) {
  s <- sim_re_int(seed = 100 + sd, ng = 30L, per = 8L)
  em <- tryCatch(sig1(fit_one(s, FALSE)), error = function(e) NA)
  aq <- tryCatch(sig1(fit_one(s, TRUE, 9L)), error = function(e) NA)
  nu <- tryCatch({
    fn <- fit_one(s, NA, method = "nuts")
    exp(fn$means[[grep("^log_sigma_", names(fn$means), value = TRUE)]])
  }, error = function(e) NA)
  c(EM = em, AGHQ = aq, NUTS = nu)
})
cat(sprintf("  EM   mean sigma = %.3f  (bias %+.3f)\n", mean(res["EM",], na.rm=TRUE),   mean(res["EM",], na.rm=TRUE) - 0.9))
cat(sprintf("  AGHQ mean sigma = %.3f  (bias %+.3f)\n", mean(res["AGHQ",], na.rm=TRUE), mean(res["AGHQ",], na.rm=TRUE) - 0.9))
cat(sprintf("  NUTS mean sigma = %.3f  (bias %+.3f)\n", mean(res["NUTS",], na.rm=TRUE), mean(res["NUTS",], na.rm=TRUE) - 0.9))
