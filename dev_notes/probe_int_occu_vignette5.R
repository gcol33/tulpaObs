devtools::load_all(".", quiet = TRUE)

# Honest payoff with NUTS posterior SD: integrated vs source-1-only.
si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(5,3), n_shared = 150,
                        beta_occ = c(0.0, 0.9),
                        beta_det = list(c(0.8,0), c(0.4,0)), seed = 7)

fi <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
           y = list(structured = si$y[[1]], casual = si$y[[2]]),
           method = "nuts",
           control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))

rows1 <- si$site_maps[[1]]
d1 <- si$data[rows1, , drop = FALSE]
f1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
           y = si$y[[1]], method = "nuts",
           control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))

si_sum <- summary(fi); s1_sum <- summary(f1)
cat("integrated psi_x:\n"); print(si_sum["psi_x", , drop=FALSE])
cat("source1 psi_x:\n"); print(s1_sum["psi_x", , drop=FALSE])
cat("\nn source1 sites:", length(rows1), " n integrated sites:", nrow(si$data), "\n")
cat("integrated CI width:", confint(fi)["psi_x","97.5%"] - confint(fi)["psi_x","2.5%"], "\n")
cat("source1   CI width:", confint(f1)["psi_x","97.5%"] - confint(f1)["psi_x","2.5%"], "\n")
print(colnames(si_sum))
