# Multi-seed correlated-slope recovery: is AGHQ rho unbiased, or does ML push
# it to the boundary? True sigma=(0.775, 0.632), rho=+0.612.
devtools::load_all(".", quiet = TRUE)

Sig <- matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)   # sd=(.775,.632), rho=.612
sim_corr <- function(seed, ng = 40L, per = 12L, J = 6L, p = 0.5) {
  set.seed(seed)
  N <- ng * per; g <- rep(seq_len(ng), each = per); x <- rnorm(N)
  U <- matrix(rnorm(ng * 2), ng, 2) %*% chol(Sig)
  z <- rbinom(N, 1, plogis(0.2 + U[g, 1] + (-0.4 + U[g, 2]) * x))
  y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  data.frame(row.names = NULL, stringsAsFactors = FALSE) -> .
  list(y = y, d = data.frame(g = factor(g), x = x))
}
pull <- function(f) {
  sg <- f$means[grep("^sigma_", names(f$means))]
  rho <- f$means[[grep("^cor_", names(f$means), value = TRUE)]]
  c(s1 = unname(sg[1]), s2 = unname(sg[2]), rho = rho)
}

seeds <- 1:10
runs <- list(
  EM             = list(aghq = FALSE, lkj = 1),
  `AGHQ lkj=1`   = list(aghq = TRUE,  lkj = 1),   # ML, no correlation prior
  `AGHQ lkj=2`   = list(aghq = TRUE,  lkj = 2)    # default LKJ regularization
)
acc <- setNames(vector("list", length(runs)), names(runs))
for (sd in seeds) {
  s <- sim_corr(400 + sd)
  for (nm in names(runs)) {
    cfg <- runs[[nm]]
    f <- tryCatch(tobs(~ x + (1 + x | g), data = s$d, y = s$y, detection = ~ 1,
                       family = occu(), method = "laplace",
                       control = list(re.aghq = cfg$aghq, n.quad = 7L,
                                      re.lkj = cfg$lkj, verbose = FALSE)),
                  error = function(e) NULL)
    if (!is.null(f)) acc[[nm]] <- rbind(acc[[nm]], pull(f))
  }
}
cat(sprintf("truth: s1=%.3f s2=%.3f rho=%.3f\n", sqrt(Sig[1,1]), sqrt(Sig[2,2]),
            Sig[1,2]/sqrt(Sig[1,1]*Sig[2,2])))
for (nm in names(runs)) {
  m <- colMeans(acc[[nm]], na.rm = TRUE)
  cat(sprintf("%-12s mean: s1=%.3f s2=%.3f rho=%.3f  | rho range [%.2f, %.2f]\n",
              nm, m["s1"], m["s2"], m["rho"],
              min(acc[[nm]][,"rho"]), max(acc[[nm]][,"rho"])))
}
