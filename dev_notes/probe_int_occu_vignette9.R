devtools::load_all(".", quiet = TRUE)
si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(2, 5), n_shared = 250,
                        beta_occ = c(0.0, 1.0),
                        beta_det = list(c(0.5, 0), c(0.7, 0)), seed = 11)
fi <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
           y = list(survey = si$y[[1]], casual = si$y[[2]]),
           method = "laplace", control = list(verbose = FALSE))
rows1 <- si$site_maps[[1]]
d1 <- si$data[rows1, , drop = FALSE]
f1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
           y = si$y[[1]], method = "laplace", control = list(verbose = FALSE))
ci_i <- confint(fi); ci_1 <- confint(f1)
cat("LAPLACE integrated psi_x:", round(coef(fi)$psi[2],3), "CI", round(ci_i["psi_x",],3),
    "width", round(ci_i["psi_x","97.5%"]-ci_i["psi_x","2.5%"],3), "\n")
cat("LAPLACE survey-only psi_x:", round(coef(f1)$psi[2],3), "CI", round(ci_1["psi_x",],3),
    "width", round(ci_1["psi_x","97.5%"]-ci_1["psi_x","2.5%"],3), "\n")
