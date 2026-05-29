devtools::load_all(".", quiet = TRUE)

# Source 1 = small structured survey (few sites, decent detection).
# Source 2 = large opportunistic set (many sites, lower detection).
# Integration should sharpen the SLOPE relative to source-1 alone via more sites.
run <- function(seed) {
  si <- simulate_int_occu(N_total = 600, n_data = 2, J = c(4, 2), n_shared = 40,
                          beta_occ = c(0.0, 1.0),
                          beta_det = list(c(0.8, 0), c(0.3, 0)), seed = seed)
  # source 1 is the small one: cut it down by taking only its own sites
  rows1 <- si$site_maps[[1]]
  # make source 1 genuinely small: keep first 120 of its sites
  keep <- rows1[seq_len(min(120, length(rows1)))]
  # rebuild source-1 y aligned
  idx1 <- match(keep, rows1)
  y1small <- si$y[[1]][idx1, , drop = FALSE]
  rownames(y1small) <- as.character(keep)
  d_all <- si$data; rownames(d_all) <- as.character(seq_len(nrow(d_all)))

  # integrated: small source 1 + full source 2
  y2 <- si$y[[2]]; rownames(y2) <- as.character(si$site_maps[[2]])
  fi <- tobs(~ x, data = d_all, family = int_occu(), detection = ~ 1,
             y = list(structured = y1small, casual = y2),
             method = "nuts",
             control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))
  d1 <- d_all[keep, , drop = FALSE]
  f1 <- tobs(~ x, data = d1, family = occu(), detection = ~ 1,
             y = y1small, method = "nuts",
             control = list(n.iter = 1000, n.warmup = 500, seed = 1, verbose = FALSE))
  si_s <- summary(fi)["psi_x", ]; s1_s <- summary(f1)["psi_x", ]
  c(int_est = si_s["mean"], int_sd = si_s["sd"],
    s1_est = s1_s["mean"], s1_sd = s1_s["sd"],
    n1 = length(keep))
}
print(round(run(3), 3))
