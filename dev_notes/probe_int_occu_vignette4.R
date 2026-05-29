devtools::load_all(".", quiet = TRUE)

widths <- function(seed) {
  si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(6,4), n_shared = 150,
                          beta_occ = c(0.0, 0.9),
                          beta_det = list(c(1.0,0), c(0.6,0)), seed = seed)
  fi <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
             y = list(structured = si$y[[1]], casual = si$y[[2]]),
             method = "laplace", control = list(verbose = FALSE))
  rows1 <- si$site_maps[[1]]
  d1 <- si$data[rows1, , drop = FALSE]
  f1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
             y = si$y[[1]], method = "laplace", control = list(verbose = FALSE))
  ci_i <- confint(fi); ci_1 <- confint(f1)
  wi <- ci_i["psi_x","97.5%"] - ci_i["psi_x","2.5%"]
  w1 <- ci_1["psi_x","97.5%"] - ci_1["psi_x","2.5%"]
  c(int_est = unname(coef(fi)$psi[2]), src1_est = unname(coef(f1)$psi[2]),
    int_width = wi, src1_width = w1, n1 = length(rows1))
}
res <- t(sapply(1:6, widths))
print(round(res, 3))
cat("\nmean int width:", mean(res[,"int_width"]), " mean src1 width:", mean(res[,"src1_width"]), "\n")
