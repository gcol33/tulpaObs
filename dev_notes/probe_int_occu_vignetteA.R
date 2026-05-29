devtools::load_all(".", quiet = TRUE)
# clearly different detection between sources: survey p high, casual p low
si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(5, 4), n_shared = 200,
                        beta_occ = c(0.0, 0.9),
                        beta_det = list(c(1.2, 0), c(-0.8, 0)), seed = 21)
cat("p survey =", round(plogis(1.2),2), " p casual =", round(plogis(-0.8),2), "\n")
fi <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
           y = list(survey = si$y[[1]], casual = si$y[[2]]),
           method = "nuts",
           control = list(n.iter = 800, n.warmup = 400, seed = 1, verbose = FALSE))
cat("survey det intercept:", round(coef(fi)$survey,3),
    " -> p =", round(plogis(coef(fi)$survey),2), "\n")
cat("casual det intercept:", round(coef(fi)$casual,3),
    " -> p =", round(plogis(coef(fi)$casual),2), "\n")
print(summary(fi))
