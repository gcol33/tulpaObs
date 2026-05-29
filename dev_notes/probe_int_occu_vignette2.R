devtools::load_all(".", quiet = TRUE)

try_one <- function(N_total, J, n_shared, beta_occ, bdet, seed) {
  si <- simulate_int_occu(N_total = N_total, n_data = 2, J = J, n_shared = n_shared,
                          beta_occ = beta_occ, beta_det = bdet, seed = seed)
  fit <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
              y = list(structured = si$y[[1]], casual = si$y[[2]]),
              method = "laplace", control = list(verbose = FALSE))
  ci <- confint(fit)
  cat(sprintf("seed=%d J=(%d,%d) p=(%.2f,%.2f) :: psi_int=%.2f (truth %.2f) psi_x=%.2f (truth %.2f)\n",
              seed, J[1], J[2], plogis(bdet[[1]][1]), plogis(bdet[[2]][1]),
              coef(fit)$psi[1], beta_occ[1], coef(fit)$psi[2], beta_occ[2]))
}

for (s in 1:5) try_one(400, c(6,4), 150, c(0.4,1.0), list(c(1.0,0),c(0.5,0)), s)
cat("---- higher detection, more visits ----\n")
for (s in 1:5) try_one(500, c(8,5), 200, c(0.0,0.8), list(c(1.2,0),c(0.8,0)), s)
