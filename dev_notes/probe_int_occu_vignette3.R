devtools::load_all(".", quiet = TRUE)

# Does NUTS show the same intercept bias? And does the simulator's z get
# observed correctly? Check: a held-out site observed by both sources.
si <- simulate_int_occu(N_total = 500, n_data = 2, J = c(8,5), n_shared = 200,
                        beta_occ = c(0.0, 0.8),
                        beta_det = list(c(1.2,0), c(0.8,0)), seed = 1)

# naive: fraction of sites with any detection in EITHER source vs true z
# build per-site observed-any
any_det <- rep(NA, 500)
for (s in 1:2) {
  rows <- si$site_maps[[s]]
  ad <- rowSums(si$y[[s]]) > 0
  for (k in seq_along(rows)) {
    r <- rows[k]
    any_det[r] <- isTRUE(any_det[r]) || ad[k]
  }
}
cat("mean true z:", mean(si$truth$z), " mean any-det:", mean(any_det, na.rm=TRUE), "\n")
cat("mean psi truth:", mean(si$truth$psi), "\n")

fitn <- tobs(~ x, data = si$data, family = int_occu(), detection = ~ 1,
             y = list(structured = si$y[[1]], casual = si$y[[2]]),
             method = "nuts",
             control = list(n.iter = 800, n.warmup = 400, seed = 1, verbose = FALSE))
cat("NUTS psi:", coef(fitn)$psi, " (truth 0.0, 0.8)\n")
cat("NUTS det:", coef(fitn)$structured, coef(fitn)$casual, "(truth 1.2, 0.8)\n")
