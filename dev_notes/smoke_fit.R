setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all("."))
set.seed(42)
sim <- simulate_occu(N = 100, J = 4, seed = 42)
fit <- tobs(
  formula   = ~ occ_cov1,
  data      = sim$data,
  family    = occu(),
  detection = ~ det_cov1,
  y         = sim$y,
  engine    = "laplace",
  control   = list(verbose = FALSE)
)
print(fit)
cat("class:", paste(class(fit), collapse = "/"), "\n")
cat("n_params:", fit$n_params, "\n")
cat("intercepts: psi =", fit$intercepts$psi, " p =", fit$intercepts$p, "\n")
