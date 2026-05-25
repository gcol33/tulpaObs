# Probe: areal-spatial Poisson N-mixture (icar/bym2) with calibrated SEs from
# the grid-integrated coefficient covariance.
suppressMessages({
  library(devtools)
  load_all(".", quiet = TRUE)
})

# 6x6 rook-adjacency grid.
side <- 6L; ng <- side * side
coords <- expand.grid(x = seq_len(side), y = seq_len(side))
adj <- matrix(0L, ng, ng)
for (i in seq_len(ng)) for (j in seq_len(ng)) {
  if (i != j && abs(coords$x[i] - coords$x[j]) + abs(coords$y[i] - coords$y[j]) == 1L)
    adj[i, j] <- 1L
}

set.seed(42)
# Spatial field: smooth-ish via averaging random noise over neighbours, demeaned.
phi <- as.numeric(scale(rnorm(ng)))
for (rep in 1:3) {
  phi_new <- phi
  for (i in seq_len(ng)) {
    nb <- which(adj[i, ] == 1L)
    phi_new[i] <- 0.5 * phi[i] + 0.5 * mean(phi[nb])
  }
  phi <- phi_new
}
phi <- 0.6 * as.numeric(scale(phi))      # sd ~ 0.6 spatial offset
phi <- phi - mean(phi)

beta_lambda <- c(log(5), 0.5)            # intercept + 1 abundance slope
beta_p      <- c(0.3, 0.4)               # intercept + 1 detection slope
x_ab <- rnorm(ng)
x_det <- rnorm(ng)
J <- 5L
lambda <- exp(beta_lambda[1] + beta_lambda[2] * x_ab + phi)
p <- plogis(beta_p[1] + beta_p[2] * x_det)
N <- rpois(ng, lambda)
y <- matrix(NA_integer_, ng, J)
for (i in seq_len(ng)) y[i, ] <- rbinom(J, N[i], p[i])
data <- data.frame(abund_cov1 = x_ab, det_cov1 = x_det)

cat("mean count:", round(mean(y), 2), " max:", max(y), "\n\n")

for (model in c("icar", "bym2")) {
  cat("=====", model, "=====\n")
  term <- if (model == "icar") quote(icar(graph = adj)) else quote(bym2(graph = adj))
  f <- as.formula(bquote(~ abund_cov1 + .(term)))
  fit <- tryCatch(
    tobs(formula = f, data = data, family = abun(), detection = ~ det_cov1,
         y = y, method = "nested_laplace", control = list(verbose = FALSE)),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) next
  cat("method:", fit$method, "\n")
  truth <- c(beta_lambda, beta_p)
  names(truth) <- names(fit$means)
  for (nm in names(fit$means)) {
    cat(sprintf("  %-22s est=% .3f  se=%.3f  truth=% .3f  z=% .2f\n",
                nm, fit$means[nm], fit$sds[nm], truth[nm],
                (fit$means[nm] - truth[nm]) / fit$sds[nm]))
  }
  V <- vcov(fit)
  cat("  vcov PD:", isTRUE(all(eigen(V, only.values = TRUE)$values > 0)),
      " finite:", all(is.finite(V)), "\n")
  cat("  hyper:", paste(names(fit$nmix_hyper), collapse = ", "), "\n\n")
}
cat("DONE\n")
