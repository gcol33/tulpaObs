devtools::load_all(".", quiet = TRUE)
si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(2, 5), n_shared = 250,
                        beta_occ = c(0.0, 1.0),
                        beta_det = list(c(0.5, 0), c(0.7, 0)), seed = 11)
fi <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
           y = list(survey = si$y[[1]], casual = si$y[[2]]),
           method = "nuts",
           control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))
rows1 <- si$site_maps[[1]]
d1 <- si$data[rows1, , drop = FALSE]
f1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
           y = si$y[[1]], method = "nuts",
           control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))
si_s <- as.numeric(summary(fi)["psi_x", c("mean","sd","q2.5","q97.5")])
s1_s <- as.numeric(summary(f1)["psi_x", c("mean","sd","q2.5","q97.5")])
cat("integrated psi_x mean/sd/lo/hi:", round(si_s,3), "\n")
cat("survey-only psi_x mean/sd/lo/hi:", round(s1_s,3), "\n")
cat("width int:", round(si_s[4]-si_s[3],3), " width survey:", round(s1_s[4]-s1_s[3],3),"\n")
